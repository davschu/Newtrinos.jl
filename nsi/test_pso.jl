# PSO settings performance test: injects several random known-truth (optionally
# nuisance-fluctuated) Asimov data realizations for a chosen NSI hypothesis, then runs
# several CLI-specified PSO seed-search settings ("combos") against every truth point and
# checks how well PSO ALONE (no LBFGS polish) recovers the injected truth.
#
# Phase 2's PSO search is seeded BLIND, at the nominal parameter values (not at the
# injected truth) -- passing the truth itself as global_seed_search's seed-anchor would
# warm-start one particle exactly at the answer, trivially "recovering" it regardless of
# PSO settings and defeating the point of testing PSO's own search ability in isolation.
#
# The reference "actual global optimum" for each truth point is a single LBFGS fit seeded
# exactly at that truth point (no PSO search -- truth_param already sits near the
# likelihood optimum by construction, so this reference fit is meant to converge fast, not
# to be blind). This is the ONLY LBFGS call in the whole script -- every PSO combo's own
# result is evaluated directly at the best point PSO itself found, with no local polish
# step, since the point of this script is to test PSO in isolation.
#
# --pso-particles / --global-iterations form a FULL FACTORIAL cross product of settings
# combos (every particle count tested against every generation budget), e.g.
#   --pso-particles 10 20 --global-iterations 10 20
# gives 4 combos: (10,10), (10,20), (20,10), (20,20) -- not just the 2 "matching" pairs.
#
# Diagnostics computed per (combo, truth) task:
#   1. Reference-relative convergence trace: Δχ²_vs_reference(gen) = 2*(pso_trace[gen] -
#      ref_neglogpost), compared against the truth-seeded LBFGS REFERENCE optimum
#      (ref_log_post from Phase 1), NOT the injected truth's own posterior value -- the
#      truth is only one point in parameter space and, especially under
#      --fluctuate-nuisances, is not generally the true local MAP; the LBFGS-from-truth
#      reference fit is a much better estimate of the actual optimum PSO is trying to
#      find. PLOTTED (the only plotted diagnostic; see below).
#   2. Δχ²_optimizer_residual = 2*(log_post - ref_log_post): final PSO result vs. the
#      truth-seeded LBFGS reference optimum (same reference as #1, evaluated only at
#      PSO's final point instead of at every generation) -- saved only.
#   3. nsi_param_diff (fit - truth per scanned NSI param) -- saved only.
#   4. plateau_gen: earliest PSO generation within --plateau-frac-thresh of the trace's
#      total start->end drop -- saved only.
#   5. seed_elapsed: PSO wall time -- saved only.
# Diagnostics 2-5 are only saved to the .jld2 (not plotted) -- inspect them there directly.
#
# Usage: julia --threads=8 test_pso.jl --hypothesis ττ_μμ --n-points 10 \
#            --pso-particles 20 40 60 --global-iterations 20 50 100
#
# Output: nsi/pso_test_results/pso_test_<hypothesis>_<fluctuated|nominal>_<n>points[_<suffix>].{jld2,pdf}

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Distributions
using DensityInterface
using BAT
using DataStructures
using MeasureBase
using ADTypes
using Newtrinos
using FileIO
using Accessors
using ArgParse
using Random
using ValueShapes
using CairoMakie
using Statistics
import Optim

include(joinpath(@__DIR__, "..", "src", "analysis", "optimizer_tests_common.jl"))

const RESULTS_DIR = joinpath(@__DIR__, "pso_test_results")

# LBFGS settings for the Phase-1 truth-seeded reference fit -- the ONLY LBFGS call in this
# script, so these are hardcoded rather than CLI flags (the PSO-settings sweep in Phase 2
# never calls LBFGS at all).
const REF_LBFGS_ITERATIONS = 1000
const REF_LBFGS_GTOL       = 1e-6
const REF_LBFGS_STEP_SIZE  = 1.0
const REF_LBFGS_FTOL       = 0.0
const REF_LBFGS_XTOL       = 0.0

### CLI ###

function parse_command_line()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--hypothesis"
        help = "Which hypothesis to test: eμ, eτ, μτ, ee_μμ, ττ_μμ, all"
        arg_type = String
        default = "ττ_μμ"

        "--n-points"
        help = "Number of independent randomly-generated truth points"
        arg_type = Int
        default = 10

        "--pso-particles"
        help = "Swarm size (n_particles) values to test -- crossed factorially with --global-iterations (every particle count x every generation budget), not zipped elementwise"
        arg_type = Int
        nargs = '+'
        default = [10, 20]

        "--global-iterations"
        help = "PSO generation budget values to test -- crossed factorially with --pso-particles"
        arg_type = Int
        nargs = '+'
        default = [10, 20]

        "--fluctuate-nuisances"
        help = "Additionally draw every nuisance parameter from its own prior per truth point (the \"(fluctuated) asimov\" mode) -- NSI truth is unaffected"
        action = :store_true

        "--vary-physics-truth"
        help = "Also draw θ₂₃'s and Δm²₃₁'s truth randomly from their priors (passed to controlled_truth's vary_physics)"
        action = :store_true

        "--plateau-frac-thresh"
        help = "Plateau-detection threshold: earliest PSO generation within this fraction of the trace's total start->end drop"
        arg_type = Float64
        default = 0.01

        "--seed"
        help = "Base RNG seed for truth generation and PSO starts"
        arg_type = Int
        default = 1234

        "--suffix"
        help = "Extra freeform tag appended to output filenames, e.g. 'run1'"
        arg_type = String
        default = ""
    end
    return parse_args(s)
end

args = parse_command_line()
hypothesis_name          = args["hypothesis"]
n_points                 = args["n-points"]
pso_particles_list       = args["pso-particles"]
global_iterations_list   = args["global-iterations"]
fluctuate_nuisances_flag = args["fluctuate-nuisances"]
vary_physics_truth_flag  = args["vary-physics-truth"]
plateau_frac_thresh      = args["plateau-frac-thresh"]
suffix                   = args["suffix"]

truth_mode_tag = fluctuate_nuisances_flag ? "fluctuated" : "nominal"
output_tag = isempty(suffix) ? "$(hypothesis_name)_$(truth_mode_tag)_$(n_points)points" :
                                "$(hypothesis_name)_$(truth_mode_tag)_$(n_points)points_$(suffix)"

### LOCAL HELPERS ###

"""Full factorial cross product of --pso-particles x --global-iterations settings combos."""
function build_combos(pso_particles::Vector{Int}, global_iterations::Vector{Int})
    [(n_particles=p, global_iterations=g) for p in pso_particles for g in global_iterations]
end

"""Short filename/label tag for one settings combo, e.g. \"p20_g20\"."""
combo_label(n_particles, global_iterations) = "p$(n_particles)_g$(global_iterations)"

"""
    plateau_generation(f_trace; frac_thresh=0.01)

Earliest generation g (1-based) at which f_trace[g] has closed to within `frac_thresh` of
the total f_trace[1]->f_trace[end] drop -- NOT a per-step relative-delta threshold, which
degenerates to reporting generation 1 when the total run-wide improvement is tiny. Returns
`length(f_trace)` if there's no net improvement (total_drop <= 0), or 0 if `f_trace` is
empty.
"""
function plateau_generation(f_trace::Vector{Float64}; frac_thresh::Float64=0.01)
    isempty(f_trace) && return 0
    total_drop = f_trace[1] - f_trace[end]
    total_drop <= 0 && return length(f_trace)
    for g in 1:length(f_trace)
        (f_trace[g] - f_trace[end]) <= frac_thresh * total_drop && return g
    end
    length(f_trace)
end

### PHYSICS CONFIG (mirrors roundtrips.jl / optimizer_tests.jl) ###

function configure_physics()
    osc = Newtrinos.osc.configure(
        Newtrinos.osc.OscillationConfig(
            flavour     = Newtrinos.osc.ThreeFlavour(),
            interaction = Newtrinos.osc.NSI_Standard()
        )
    )
    atm_flux     = Newtrinos.atm_flux.configure()
    earth_layers = Newtrinos.earth_layers.configure(Newtrinos.earth_layers.PREM12())
    xsec         = Newtrinos.xsec.configure()
    physics      = (; osc, atm_flux, earth_layers, xsec)
    (; deepcore = Newtrinos.deepcore.configure(physics))
end

experiments = configure_physics()

ftype = Float64
p = Newtrinos.get_params(experiments)
priors = Newtrinos.get_priors(experiments)
conditional_vars = Dict(:δCP=>0.0, :θ₁₂=>ftype(33.62 * π/180), :θ₁₃=>ftype(8.54  * π/180), :Δm²₂₁=>ftype(7.40e-5),
                         :atm_flux_updown_sigma=>ftype(0.0), :atm_flux_nuenuebar_sigma=>ftype(0.0))
priors = Newtrinos.condition(priors, conditional_vars, p)
@reset priors.θ₂₃   = Uniform(ftype(30*π/180), ftype(60*π/180))
@reset priors.Δm²₃₁ = Uniform(ftype(0.93e-3 + 7.40e-5), ftype(3.93e-3 + 7.40e-5))

@reset priors.deepcore_lifetime        = Uniform(ftype(0), ftype(3.8))
@reset priors.deepcore_ice_scattering  = Truncated(Normal(ftype(1), ftype(0.1)), ftype(0.9), ftype(1.1))
@reset priors.deepcore_ice_absorption  = Truncated(Normal(ftype(1), ftype(0.1)), ftype(0.9), ftype(1.1))
@reset priors.deepcore_opt_eff_lateral = Truncated(Normal(ftype(0), ftype(1.)), ftype(-2), ftype(2.5))

@reset priors.Δ_eμ     = Uniform(-ftype(5),   ftype(5))
@reset priors.Δ_τμ     = Uniform(-ftype(0.1), ftype(0.1))
@reset priors.ε_eμ_abs = Uniform(ftype(0),    ftype(0.3))
@reset priors.ε_eτ_abs = Uniform(ftype(0),    ftype(0.35))
@reset priors.ε_μτ_abs = Uniform(ftype(0),    ftype(0.07))
@reset priors.δ_eμ     = Uniform(ftype(0),    ftype(2π))
@reset priors.δ_eτ     = Uniform(ftype(0),    ftype(2π))
@reset priors.δ_μτ     = Uniform(ftype(0),    ftype(2π))

cp_eμ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eτ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eτ => 0.0, :δ_μτ => 0.0), p)
cp_eτ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eμ => 0.0, :δ_μτ => 0.0), p)
cp_μτ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0), p)
cp_ee_μμ = Newtrinos.condition(priors, Dict(:Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0, :δ_μτ => 0.0), p)
cp_ττ_μμ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0, :δ_μτ => 0.0), p)
cp_all = deepcopy(priors)

all_cps = Dict("eμ" => cp_eμ, "eτ" => cp_eτ, "μτ" => cp_μτ, "ee_μμ" => cp_ee_μμ, "ττ_μμ" => cp_ττ_μμ, "all" => cp_all)
hypothesis_name in keys(all_cps) || error("--hypothesis must be one of $(join(sort(collect(keys(all_cps))), ", ")), got '$hypothesis_name'")
cp_hyp = all_cps[hypothesis_name]

nsi_param_names = (:Δ_eμ, :Δ_τμ, :ε_eμ_abs, :ε_eτ_abs, :ε_μτ_abs, :δ_eμ, :δ_eτ, :δ_μτ)
nuisance_names = (
    :atm_flux_delta_spectral_index, :atm_flux_nuenuebar_sigma, :atm_flux_nuenumu_sigma,
    :atm_flux_numunumubar_sigma, :atm_flux_updown_sigma, :atm_flux_uphorizonzal_sigma,
    :deepcore_atm_muon_scale, :deepcore_ice_absorption, :deepcore_ice_scattering,
    :deepcore_lifetime, :deepcore_opt_eff_headon, :deepcore_opt_eff_lateral,
    :deepcore_opt_eff_overall, :nc_norm, :nutau_cc_norm,
)
physics_names = (:θ₂₃, :Δm²₃₁)

scan_keys = [k for k in nsi_param_names if !(cp_hyp[k] isa ValueShapes.ConstValueDist) && !(cp_hyp[k] isa Number)]

prior_dist = distprod(;cp_hyp...)
algorithm  = make_algorithm("lbfgs", REF_LBFGS_STEP_SIZE)
combos     = build_combos(pso_particles_list, global_iterations_list)
n_combos   = length(combos)

println("Hypothesis '$hypothesis_name': scan_keys=$(Tuple(scan_keys)), n_points=$n_points, combos=$(join([combo_label(c.n_particles,c.global_iterations) for c in combos], ", "))")

### PHASE 1: one truth-seeded LBFGS reference fit per truth point (independent of combo) ###

phase1 = Vector{NamedTuple}(undef, n_points)

println("="^80)
println("Phase 1: building $n_points truth points + reference fits...")

Threads.@threads for t in 1:n_points
    truth_rng = Random.Xoshiro(args["seed"] + t)
    truth_param = controlled_truth(truth_rng, p, cp_hyp, nsi_param_names; vary_physics=vary_physics_truth_flag)
    truth_param = fluctuate_nuisances_flag ? fluctuate_nuisances(truth_rng, truth_param, cp_hyp, nuisance_names) : truth_param

    injected_data    = Newtrinos.generate_asimov_data(experiments, truth_param)
    truth_likelihood = Newtrinos.generate_likelihood(experiments, injected_data)
    truth_log_post   = logdensityof(PosteriorMeasure(truth_likelihood, prior_dist), truth_param)

    t0 = time()
    ref_llh, ref_log_post, ref_fit_param, ref_converged, ref_n_iters = local_find_mle(
        truth_likelihood, prior_dist, truth_param;
        fit_method="optim", algorithm=algorithm, iterations=REF_LBFGS_ITERATIONS,
        g_tol=REF_LBFGS_GTOL, f_tol=REF_LBFGS_FTOL, x_tol=REF_LBFGS_XTOL, objective="posterior",
        ad_backend=ADTypes.AutoForwardDiff())
    ref_elapsed = time() - t0

    phase1[t] = (truth_param=truth_param, truth_likelihood=truth_likelihood, truth_log_post=truth_log_post,
                 ref_llh=ref_llh, ref_log_post=ref_log_post, ref_fit_param=ref_fit_param,
                 ref_converged=ref_converged, ref_n_iters=ref_n_iters, ref_elapsed=ref_elapsed)
    println("  [truth $t/$n_points] truth_log_post=$(round(truth_log_post,digits=4)) " *
            "ref_log_post=$(round(ref_log_post,digits=4)) ref_converged=$ref_converged elapsed=$(round(ref_elapsed,digits=1))s")
end

n_ref_converged = count(r -> r.ref_converged, phase1)
println("Phase 1 done: $n_ref_converged/$n_points reference fits converged")

### PHASE 2: PSO-only sweep, one task per (combo, truth) ###

phase2_tasks = [(combo_idx=c, truth_idx=t) for c in 1:n_combos for t in 1:n_points]
n_tasks = length(phase2_tasks)
phase2 = Vector{NamedTuple}(undef, n_tasks)

println("="^80)
println("Phase 2: running $n_tasks PSO-only tasks ($n_combos combos x $n_points truths)...")

Threads.@threads for k in 1:n_tasks
    combo_idx, truth_idx = phase2_tasks[k].combo_idx, phase2_tasks[k].truth_idx
    combo = combos[combo_idx]
    tp = phase1[truth_idx]
    task_rng = Random.Xoshiro(args["seed"] + 1_000_000 + k)

    t0 = time()
    # Seeded BLIND at nominal params `p`, NOT at tp.truth_param -- global_seed_search uses
    # its `params` argument as the free-dims x0 anchor for the PSO search (extract_bounds),
    # so passing the truth here would warm-start (part of) the swarm exactly at the answer,
    # trivially "recovering" it regardless of PSO settings and making this diagnostic
    # measure almost nothing (see file header).
    fit_param, pso_trace = global_seed_search(tp.truth_likelihood, cp_hyp, p, "pso",
                                               combo.global_iterations, task_rng;
                                               n_particles=combo.n_particles, return_trace=true, objective="posterior")
    seed_elapsed = time() - t0

    llh      = logdensityof(tp.truth_likelihood, fit_param)
    log_post = logdensityof(PosteriorMeasure(tp.truth_likelihood, prior_dist), fit_param)

    Δχ²_optimizer_residual = 2 * (log_post - tp.ref_log_post)
    Δχ²_recovery            = 2 * (log_post - tp.truth_log_post)
    nsi_param_diff = Dict{Symbol,Float64}(sk => Float64(fit_param[sk]) - Float64(tp.truth_param[sk]) for sk in scan_keys)
    plateau_gen = plateau_generation(pso_trace; frac_thresh=plateau_frac_thresh)

    phase2[k] = (combo_idx=combo_idx, truth_idx=truth_idx,
                 combo_label=combo_label(combo.n_particles, combo.global_iterations),
                 fit_param=fit_param, llh=llh, log_post=log_post,
                 Δχ²_optimizer_residual=Δχ²_optimizer_residual, Δχ²_recovery=Δχ²_recovery,
                 nsi_param_diff=nsi_param_diff, plateau_gen=plateau_gen,
                 pso_trace=pso_trace, seed_elapsed=seed_elapsed)
    println("  [task $k/$n_tasks] combo=$(phase2[k].combo_label) truth=$truth_idx " *
            "log_post=$(round(log_post,digits=4)) Δχ²_optimizer_residual=$(round(Δχ²_optimizer_residual,digits=4)) " *
            "plateau_gen=$plateau_gen elapsed=$(round(seed_elapsed,digits=1))s")
end

println("="^80)
println("Phase 2 done. Per-combo summary (mean over $n_points truth points):")
for (c, combo) in enumerate(combos)
    rows = filter(r -> r.combo_idx == c, phase2)
    mean_dchi2    = mean(r.Δχ²_optimizer_residual for r in rows)
    mean_plateau  = mean(r.plateau_gen for r in rows)
    mean_elapsed  = mean(r.seed_elapsed for r in rows)
    println("  $(combo_label(combo.n_particles,combo.global_iterations)): " *
            "mean Δχ²_optimizer_residual=$(round(mean_dchi2,digits=4)) " *
            "mean plateau_gen=$(round(mean_plateau,digits=1)) mean seed_elapsed=$(round(mean_elapsed,digits=1))s")
end

### SAVE RESULTS (.jld2, before plotting) ###
# One self-contained .jld2 PER SETTINGS COMBO (not one combined file across combos) --
# each file carries its own copy of the (combo-independent) Phase-1 truth/reference data
# alongside just that combo's Phase-2 results, so a given combo's results never need to be
# picked back out of interleaved cross-combo arrays.

mkpath(RESULTS_DIR)

for (c, combo) in enumerate(combos)
    rows = filter(r -> r.combo_idx == c, phase2)  # this combo's n_points tasks, already truth_idx-ordered
    combo_tag = combo_label(combo.n_particles, combo.global_iterations)

    nsi_param_diff_flat = Dict{String, Vector{Float64}}(
        "nsi_param_diff_$(k)" => [r.nsi_param_diff[k] for r in rows] for k in scan_keys
    )

    FileIO.save(joinpath(RESULTS_DIR, "pso_test_$(output_tag)_$(combo_tag).jld2"), Dict(
        "hypothesis"  => hypothesis_name,
        "scan_keys"   => collect(scan_keys),
        "output_tag"  => output_tag,
        "n_points"    => n_points,

        "combo_label"        => combo_tag,
        "combo_n_particles"  => combo.n_particles,
        "combo_global_iterations" => combo.global_iterations,

        "truth_params"   => [phase1[r.truth_idx].truth_param for r in rows],
        "truth_log_post" => [phase1[r.truth_idx].truth_log_post for r in rows],
        "ref_llh"        => [phase1[r.truth_idx].ref_llh for r in rows],
        "ref_log_post"   => [phase1[r.truth_idx].ref_log_post for r in rows],
        "ref_fit_params" => [phase1[r.truth_idx].ref_fit_param for r in rows],
        "ref_converged"  => [phase1[r.truth_idx].ref_converged for r in rows],
        "ref_n_iters"    => [phase1[r.truth_idx].ref_n_iters for r in rows],
        "ref_elapsed"    => [phase1[r.truth_idx].ref_elapsed for r in rows],

        "truth_idx"        => [r.truth_idx for r in rows],
        "fit_params"       => [r.fit_param for r in rows],
        "llh"              => [r.llh for r in rows],
        "log_post"         => [r.log_post for r in rows],
        "Δχ²_optimizer_residual" => [r.Δχ²_optimizer_residual for r in rows],
        "Δχ²_recovery"           => [r.Δχ²_recovery for r in rows],
        "nsi_param_diff"   => [r.nsi_param_diff for r in rows],
        nsi_param_diff_flat...,
        "plateau_gen"      => [r.plateau_gen for r in rows],
        "pso_trace"        => [r.pso_trace for r in rows],
        "seed_elapsed"     => [r.seed_elapsed for r in rows],

        "settings" => Dict(
            "hypothesis"          => hypothesis_name,
            "n_points"             => n_points,
            "combo_label"          => combo_tag,
            "combo_n_particles"    => combo.n_particles,
            "combo_global_iterations" => combo.global_iterations,
            "fluctuate_nuisances"  => fluctuate_nuisances_flag,
            "vary_physics_truth"   => vary_physics_truth_flag,
            "plateau_frac_thresh"  => plateau_frac_thresh,
            "seed"                 => args["seed"],
            "ref_lbfgs_iterations" => REF_LBFGS_ITERATIONS,
            "ref_lbfgs_gtol"       => REF_LBFGS_GTOL,
            "ref_lbfgs_step_size"  => REF_LBFGS_STEP_SIZE,
        ),
        "settings_hypothesis"          => hypothesis_name,
        "settings_n_points"            => n_points,
        "settings_combo_label"         => combo_tag,
        "settings_combo_n_particles"   => combo.n_particles,
        "settings_combo_global_iterations" => combo.global_iterations,
        "settings_fluctuate_nuisances" => fluctuate_nuisances_flag,
        "settings_vary_physics_truth"  => vary_physics_truth_flag,
        "settings_plateau_frac_thresh" => plateau_frac_thresh,
        "settings_seed"                => args["seed"],
    ))
    println("Saved results for combo $combo_tag to $(joinpath(RESULTS_DIR, "pso_test_$(output_tag)_$(combo_tag).jld2"))")
end

### PLOT: trace diagnostics only (diagnostic #1) ###
# Diagnostics #2-#5 (final value vs. actual optimum, NSI param recovery, plateau
# generations, runtime) are saved above but not plotted -- inspect them from the .jld2.

n_cols = min(n_combos, 3)
n_rows = cld(n_combos, n_cols)
fig = Figure(size=(500 * n_cols, 400 * n_rows))

for (c, combo) in enumerate(combos)
    row = (c - 1) ÷ n_cols + 1
    col = (c - 1) % n_cols + 1
    # pseudolog10 (sign(x)*log10(1+|x|)) instead of a clipped log10: Δχ²_vs_reference is
    # SIGNED (positive = worse than the reference optimum, negative = PSO's trajectory
    # actually beat the reference -- possible since PSO explores randomly and the
    # reference is only a single local LBFGS fit, not a guaranteed global optimum).
    # Clipping negative values to eps() before log10 (an earlier approach) silently
    # floors every "beat reference" excursion to the same invisible point, which can make
    # an actually-informative trace look like a flat line. pseudolog10 keeps both
    # directions visible on one log-like axis with no data thrown away.
    ax = Axis(fig[row, col], xlabel="PSO generation", ylabel="Δχ²_vs_reference = 2×(trace - ref_neglogpost)",
              yscale=Makie.pseudolog10, title=combo_label(combo.n_particles, combo.global_iterations))

    rows = filter(r -> r.combo_idx == c, phase2)
    for r in rows
        tp = phase1[r.truth_idx]
        # Compared against tp.ref_log_post (the truth-seeded LBFGS reference optimum from
        # Phase 1), NOT tp.truth_log_post (the injected truth's own posterior value) --
        # the reference is the better estimate of the actual optimum PSO is searching for.
        ref_neglogpost = -tp.ref_log_post
        Δχ²_vs_reference = 2 .* (r.pso_trace .- ref_neglogpost)
        color = Makie.wong_colors()[mod1(r.truth_idx, 7)]
        lines!(ax, 0:length(Δχ²_vs_reference)-1, Δχ²_vs_reference; color=(color, 0.8), linewidth=2)
    end
    hlines!(ax, [0.0]; color=:black, linestyle=:dash, linewidth=2)
end

title_suffix = fluctuate_nuisances_flag ? ", fluctuated nuisances" : ""
combos_str = join([combo_label(c.n_particles, c.global_iterations) for c in combos], ", ")
Label(fig[0, 1:n_cols], "PSO trace diagnostics -- hypothesis=$hypothesis_name, n_points=$n_points$title_suffix, combos=$combos_str",
      fontsize=16, font=:bold)

save(joinpath(RESULTS_DIR, "pso_test_$(output_tag).pdf"), fig)
println("Saved plot to $(joinpath(RESULTS_DIR, "pso_test_$(output_tag).pdf"))")
