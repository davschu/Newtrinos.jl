# Thesis-style roundtrip scan: replicates the dissertation's "recovering injected NSI
# hypotheses" test (section 7.6.2, nuisances/theta23/deltam31 pinned at nominal) and,
# with --fluctuate-nuisances, its "recovering fluctuated nuisance parameters" test
# (section 7.6.3, nuisances + theta23/deltam31 redrawn from their own priors per grid
# point). All truth-vs-fit panels (NSI scan panel(s), nuisance/physics panels, and the
# Δχ²_mod histogram) are combined into ONE figure/ONE output file per run -- no separate
# "_nuisances.pdf" file.
#
# Per-hypothesis NSI panel(s):
#   - 1D hypotheses (single real-valued NSI param, e.g. Δ_τμ): one truth-vs-fit diagonal
#     panel.
#   - 2D hypotheses (magnitude+phase pair, e.g. ε_eμ_abs+δ_eμ): THREE panels -- the
#     (magnitude, phase) plane with truth and fit as separate point sets (matching the
#     dissertation's Fig. 7.14b layout), plus two diagonal truth-vs-fit panels (one for
#     magnitude, one for phase), same visual language as the nuisance panels.
#   - "all" (8 free NSI params, Sobol sample): 8 truth-vs-fit diagonal panels, one per
#     free NSI param.
#
# The GRID controls only how truths are injected (via truth_grid, optimizer_tests_common.jl)
# -- it does NOT restrict the fit. Every fit is still free to vary every non-fixed
# parameter in the hypothesis's prior (all physics params + all 15 nuisances), exactly
# like any other roundtrip; only the injected truth is simplified to isolate the scan to
# the NSI parameter(s) being tested (plus, under --fluctuate-nuisances, the nuisances and
# θ₂₃/Δm²₃₁).
#
# --hypothesis: one of eμ, eτ, μτ, ee_μμ, ττ_μμ, all.
#
# Δχ²_mod statistic -- pick ONE via mutually-exclusive flags (--chi2-recovery is default):
#   --chi2-recovery: Δχ²_recovery = 2*(log_post_fit - log_post_truth), i.e. how far the
#     PSO+LBFGS pipeline's best-fit point is (in -2logL units) from the injected truth
#     point. A pipeline/optimizer diagnostic -- "did the fit find the truth" -- NOT a
#     hypothesis-test statistic, and not χ²-distributed with a small fixed dof (its
#     effective dof tracks however many directions are practically constrained by the
#     data at that point).
#   --chi2-mis-modeling: Δχ²_mis-mod = χ²(NSI param(s) fixed at truth, refit) -
#     χ²(NSI param(s) free, best fit) -- the dissertation's actual profile-likelihood
#     Δχ²_mod (section 7.2.1 family): under Wilks' theorem, asymptotically χ²-distributed
#     with dof = number of NSI params fixed (1 for standalone real params, 2 for
#     magnitude+phase pairs). Requires a SECOND LBFGS fit per grid point (NSI param(s)
#     fixed at the injected truth value, everything else re-optimized), roughly doubling
#     per-point cost -- only run when you actually want this statistic.
#
# --prior-effect: additionally decomposes Δχ²_recovery into two pieces, storing both
#   alongside the chosen Δχ²_mod/Δχ²_recovery/Δχ²_mis_modeling statistic(s) (no effect on
#   plots). Under --objective posterior, the injected truth is only guaranteed to maximize
#   the LIKELIHOOD (Asimov data is generated to peak exactly there) -- it is generally NOT
#   the posterior maximum, since the prior can pull the true MAP away from truth (e.g.
#   under --fluctuate-nuisances, nuisances drawn away from their prior mode). This means a
#   positive Δχ²_recovery under --objective posterior conflates two different effects:
#   expected prior pull, and genuine optimizer underperformance. --prior-effect separates
#   them by running one extra CHEAP local-only fit per grid point, seeded exactly at
#   truth_param (no PSO search -- truth_param already sits near the likelihood optimum, so
#   this should converge fast) with objective="posterior" (prior pull is a posterior-only
#   phenomenon; irrelevant under objective="likelihood", so this reference fit always
#   targets posterior regardless of --objective), giving a cheap estimate log_post_ref of
#   the true local posterior maximum near truth:
#     Δχ²_prior_pull        = 2*(log_post_ref - truth_log_post)   -- expected, not a bug
#     Δχ²_optimizer_residual = 2*(log_post      - log_post_ref)    -- should be ≈0 (or
#       negative) for a well-converged main PSO+LBFGS fit; meaningfully negative means the
#       main pipeline underperformed a simple local fit started right at the truth.
#   Both NaN when --prior-effect is not set.
#
# --objective posterior|likelihood: what the local fit maximizes (default posterior, i.e.
#   MAP). likelihood fits against the raw likelihood via a flat prior over the same bounds
#   (i.e. MLE) -- see local_find_mle's docstring in optimizer_tests_common.jl. Δχ²_recovery
#   is computed from llh (not log_post) in likelihood mode, since log_post there is only the
#   flat-prior pseudo-posterior, not the true MAP posterior truth_log_post is computed from.
#
# --global-iterations/--pso-particles/--iterations/--step-size/--g-tol/--f-tol/--x-tol:
#   tunable PSO seed-search and LBFGS polish settings (same flag names as
#   optimizer_tests.jl); all recorded in the saved .jld2's "settings" block.
#
# --toydata: inject Poisson-fluctuated toy data (Newtrinos.generate_toy_data) instead of
#   the default noiseless Asimov data. Also widens the plot to show all nuisance/physics
#   panels (same as --fluctuate-nuisances), since with toy data every free parameter's
#   fit is worth inspecting, not just the scanned NSI param(s).
#
# All grid points are independent (same likelihood-generation machinery, different
# truth) -- runs via Threads.@threads, same index-owned-output pattern as
# test_global_seed.jl. Launch with multiple threads to parallelize, e.g.:
#
# Usage: julia --threads=8 truth_grid_scan.jl --hypothesis ττ_μμ
#
# Output directory: Optimization_tests/roundtrip_results/ (created if missing).
# Output filenames: grid_roundtrip_<asimov|toydata>_<fluctuated|nominal>_<hypothesis>_<posteriorfit|likelihoodfit>[_<suffix>].{pdf,jld2}
# -- "asimov"/"toydata" reflects whether --toydata was set, "fluctuated"/"nominal" reflects
# whether --fluctuate-nuisances was set, "posteriorfit"/"likelihoodfit" reflects --objective;
# --suffix is freeform run-tagging (e.g. "run1").
# Results (.jld2) are saved BEFORE plotting, so they survive even if plotting fails.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "Newtrinos.jl"))

using Distributions
using DensityInterface
using BAT
using DataStructures
using MeasureBase
using ADTypes
using Newtrinos
using FileIO
using Accessors
using Random
using ArgParse
import Optim
import ForwardDiff
using ValueShapes
using CairoMakie
using Sobol
using Statistics

include(joinpath(@__DIR__, "optimizer_tests_common.jl"))

const RESULTS_DIR = joinpath(@__DIR__, "roundtrip_results")

### CLI ###

function parse_command_line()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--hypothesis"
        help = "Which hypothesis to scan: eμ, eτ, μτ, ee_μμ, ττ_μμ, all (all uses Sobol quasi-random sampling over all 8 free NSI params instead of a factorial grid -- use a larger --n-points, e.g. 100-200)"
        arg_type = String
        default = "ττ_μμ"

        "--n-points"
        help = "Total sample points (1D: evenly-spaced values; 2D: round(sqrt(n)) per dimension; all: Sobol quasi-random samples -- recommend 100-200)"
        arg_type = Int
        default = 25

        "--fluctuate-nuisances"
        help = "Replicate the dissertation's section 7.6.3 test: also randomly draw every nuisance parameter, plus θ₂₃/Δm²₃₁, from their own priors (NSI truth from --hypothesis's grid/Sobol point is unaffected). Adds per-parameter truth-vs-fit panels and a Δχ²_mod histogram to the output figure. Also selects the 'fluctuated' vs 'nominal' output filename tag."
        action = :store_true

        "--chi2-recovery"
        help = "Δχ²_mod = pipeline-recovery statistic: 2*(log_post_fit - log_post_truth) (default). Mutually exclusive with --chi2-mis-modeling."
        action = :store_true

        "--chi2-mis-modeling"
        help = "Δχ²_mod = profile-likelihood mis-modeling statistic: χ²(NSI param(s) fixed at truth) - χ²(NSI param(s) free) (dissertation's actual Δχ²_mod). Requires a second LBFGS fit per grid point -- roughly doubles runtime. Mutually exclusive with --chi2-recovery."
        action = :store_true

        "--prior-effect"
        help = "Additionally compute and save Δχ²_prior_pull and Δχ²_optimizer_residual per grid point (decomposes Δχ²_recovery into expected prior-pull vs. genuine optimizer underperformance -- see header comment). Adds one cheap local-only reference fit per grid point, seeded at the truth. Does not affect plots."
        action = :store_true

        "--toydata"
        help = "Generate toy data with bin poisson fluctuation. activates when option --toydata is appended"
        action = :store_true

        "--suffix"
        help = "Extra freeform tag appended to every output filename (grid_roundtrip_<asimov|toydata>_<fluctuated|nominal>_<hypothesis>_<posteriorfit|likelihoodfit>_<suffix>.pdf/.jld2), e.g. 'run1', to keep repeated runs from overwriting each other."
        arg_type = String
        default = ""

        "--objective"
        help = "What bat_findmode maximizes for the local fit: posterior (likelihood x prior, i.e. MAP -- default) or likelihood (likelihood only, via a flat prior over the same bounds -- i.e. MLE)"
        arg_type = String
        default = "posterior"

        "--global-iterations"
        help = "Iteration budget for the PSO global seed-search stage"
        arg_type = Int
        default = 50

        "--pso-particles"
        help = "Swarm size (n_particles) for Optim.ParticleSwarm during the global seed-search stage"
        arg_type = Int
        default = 20

        "--iterations"
        help = "Max iterations for the local LBFGS optimizer (outer_iterations for Fminbox)"
        arg_type = Int
        default = 300

        "--step-size"
        help = "Initial step size (alphaguess) for LBFGS"
        arg_type = Float64
        default = 1.0

        "--g-tol"
        help = "Gradient-norm convergence tolerance"
        arg_type = Float64
        default = 0.0

        "--f-tol"
        help = "Function-value relative convergence tolerance"
        arg_type = Float64
        default = 0.0

        "--x-tol"
        help = "Parameter-value relative convergence tolerance"
        arg_type = Float64
        default = 0.0

    end
    return parse_args(s)
end

args = parse_command_line()
hypothesis_name = args["hypothesis"]
fluctuate_nuisances_flag = args["fluctuate-nuisances"]
n_points = args["n-points"]
objective = lowercase(args["objective"])
global_iterations = args["global-iterations"]
pso_particles = args["pso-particles"]
lbfgs_iterations = args["iterations"]
step_size = args["step-size"]
g_tol = args["g-tol"]
f_tol = args["f-tol"]
x_tol = args["x-tol"]
prior_effect_flag = args["prior-effect"]
toydata_flag = args["toydata"]

args["chi2-recovery"] && args["chi2-mis-modeling"] &&
    error("--chi2-recovery and --chi2-mis-modeling are mutually exclusive -- pick one.")
chi2_mode = args["chi2-mis-modeling"] ? "mis-modeling" : "recovery"  # recovery is default
objective in ("posterior", "likelihood") || error("--objective must be 'posterior' or 'likelihood', got '$objective'")

data_mode_tag = toydata_flag ? "toydata" : "asimov"
truth_mode_tag = fluctuate_nuisances_flag ? "fluctuated" : "nominal"
objective_tag = objective == "likelihood" ? "likelihoodfit" : "posteriorfit"
output_tag = isempty(args["suffix"]) ? "$(data_mode_tag)_$(truth_mode_tag)_$(hypothesis_name)_$(objective_tag)" : "$(data_mode_tag)_$(truth_mode_tag)_$(hypothesis_name)_$(objective_tag)_$(args["suffix"])"

Δχ²_label = chi2_mode == "mis-modeling" ? "Δχ² (mis-mod)" : "Δχ² (recovery)"

### PHYSICS CONFIG (mirrors optimizer_tests.jl / test_global_seed.jl) ###

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

#conditioned priors for 1-by-1 hypotheses (mirrors optimizer_tests.jl exactly)
cp_eμ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eτ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eτ => 0.0, :δ_μτ => 0.0), p)
cp_eτ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eμ => 0.0, :δ_μτ => 0.0), p)
cp_μτ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0), p)
cp_ee_μμ = Newtrinos.condition(priors, Dict(:Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0, :δ_μτ => 0.0), p)
cp_ττ_μμ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0, :δ_μτ => 0.0), p)
cp_all = deepcopy(priors)

all_cps = Dict("eμ" => cp_eμ, "eτ" => cp_eτ, "μτ" => cp_μτ, "ee_μμ" => cp_ee_μμ, "ττ_μμ" => cp_ττ_μμ, "all" => cp_all)
hypothesis_name in keys(all_cps) || error("--hypothesis must be one of $(join(sort(collect(keys(all_cps))), ", ")), got '$hypothesis_name'")
cp_hyp = all_cps[hypothesis_name]
is_all = hypothesis_name == "all"

nsi_param_names = (:Δ_eμ, :Δ_τμ, :ε_eμ_abs, :ε_eτ_abs, :ε_μτ_abs, :δ_eμ, :δ_eτ, :δ_μτ)

# All 15 nuisance parameters (everything free in the prior that isn't an NSI or
# physics-of-interest param). Confirmed directly against Newtrinos.get_priors: 12 have a
# Truncated(Normal(μ,σ)) prior (get a 1σ shaded band in the panel), 3 have a Uniform
# prior (no natural 1σ concept, no band).
nuisance_names = (
    :atm_flux_delta_spectral_index, :atm_flux_nuenuebar_sigma, :atm_flux_nuenumu_sigma,
    :atm_flux_numunumubar_sigma, :atm_flux_updown_sigma, :atm_flux_uphorizonzal_sigma,
    :deepcore_atm_muon_scale, :deepcore_ice_absorption, :deepcore_ice_scattering,
    :deepcore_lifetime, :deepcore_opt_eff_headon, :deepcore_opt_eff_lateral,
    :deepcore_opt_eff_overall, :nc_norm, :nutau_cc_norm,
)

# Physics-of-interest params (kept separate from `nuisance_names` conceptually -- these
# are the oscillation parameters `truth_grid`/`controlled_truth` normally pin at nominal
# -- but fluctuated and plotted alongside the nuisances under --fluctuate-nuisances,
# since the dissertation's 7.6.3 test treats them the same way: both are "not the NSI
# parameter(s) under test" but are free to fluctuate/be refit).
physics_names = (:θ₂₃, :Δm²₃₁)
fluctuated_names = (nuisance_names..., physics_names...)

### BUILD THE TRUTH GRID ###

truths, scan_keys, grid_coords = truth_grid(cp_hyp, p, nsi_param_names; n_points=n_points)
n_grid = length(truths)
is_2d = length(scan_keys) == 2 && !is_all

println("Truth grid built for hypothesis '$hypothesis_name': $(n_grid) points, scanning $(scan_keys)")

algorithm = make_algorithm("lbfgs", step_size)

### ROUNDTRIP SCAN (PSO seed + LBFGS polish at every grid point, parallel) ###

results = Vector{NamedTuple}(undef, n_grid)

println("="^80)
println("Running $(n_grid) grid-point roundtrips on $(Threads.nthreads()) threads (Δχ²_mod mode: $chi2_mode)...")

Threads.@threads for k in 1:n_grid
    seed_rng = Random.Xoshiro(1000 + k)
    truth_param = fluctuate_nuisances_flag ?
        fluctuate_nuisances(seed_rng, truths[k], cp_hyp, fluctuated_names) : truths[k]

    t0 = time()
    # Newtrinos.generate_toy_data takes no RNG argument -- it always draws via plain
    # rand(dist_obj) against Julia's global/task-local RNG, so its Poisson fluctuation is
    # NOT reproducible from seed_rng on its own. Seed the global RNG with an explicit
    # integer (derived from this grid point's own seed, 1000+k, same value seed_rng was
    # built from) immediately before the call, so the toy-data draw becomes deterministic
    # and tied to this grid point. NOTE: Random.seed!(rng::AbstractRNG) (no seed value) is
    # NOT reproducible -- it re-randomizes rng from system entropy (confirmed against
    # Random.seed!'s own docstring example: `rand(Random.seed!(rng), Bool) # not
    # reproducible`) -- only the explicit-integer form `Random.seed!(::Integer)` seeds
    # deterministically, which is why this reseeds the GLOBAL rng by integer rather than
    # trying to copy seed_rng's state into it.
    toydata_flag && Random.seed!(1000 + k)
    injected_data = toydata_flag ?
        Newtrinos.generate_toy_data(experiments, truth_param) : Newtrinos.generate_asimov_data(experiments, truth_param)
    truth_likelihood = Newtrinos.generate_likelihood(experiments, injected_data)
    truth_posterior = PosteriorMeasure(truth_likelihood, prior_dist)
    truth_log_post = logdensityof(truth_posterior, truth_param)

    start_param = global_seed_search(truth_likelihood, cp_hyp, truth_param, "pso", global_iterations, seed_rng; n_particles=pso_particles, objective=objective)
    seed_elapsed = time() - t0

    t1 = time()
    llh, log_post, fit_param, converged, n_iters = local_find_mle(truth_likelihood, prior_dist, start_param;
        fit_method="optim", algorithm=algorithm, iterations=lbfgs_iterations, g_tol=g_tol, f_tol=f_tol, x_tol=x_tol,
        ad_backend=ADTypes.AutoForwardDiff(), objective=objective)
    lbfgs_elapsed = time() - t1

    # χ²_recovery: pipeline-recovery diagnostic (see header comment). Always computed
    # (cheap, already have both values / one extra logdensityof call). Branches on
    # `objective` since with objective="likelihood" the fit target is the flat-prior
    # pseudo-posterior (≈ llh + const), not the true MAP posterior -- log_post/truth_log_post
    # would compare like-for-unlike in that mode (see local_find_mle's docstring).
    Δχ²_recovery = if objective == "likelihood"
        truth_llh = logdensityof(truth_likelihood, truth_param)
        2 * (llh - truth_llh)
    else
        2 * (log_post - truth_log_post)
    end

    prior_pull_elapsed = 0.0
    Δχ²_prior_pull = NaN
    Δχ²_optimizer_residual = NaN
    if prior_effect_flag
        # Cheap local-only reference fit seeded exactly at truth_param (no PSO search --
        # truth_param already sits near the likelihood optimum by construction, so this
        # should converge fast from a good starting point). Always targets objective=
        # "posterior": prior pull is a posterior-only phenomenon (see header comment),
        # so under --objective likelihood this reference fit still estimates the local
        # posterior maximum near truth, independent of what the main fit optimized.
        t3 = time()
        _, log_post_ref, _, converged_ref, _ = local_find_mle(truth_likelihood, prior_dist, truth_param;
            fit_method="optim", algorithm=algorithm, iterations=lbfgs_iterations, g_tol=g_tol, f_tol=f_tol, x_tol=x_tol,
            ad_backend=ADTypes.AutoForwardDiff(), objective="posterior")
        prior_pull_elapsed = time() - t3
        Δχ²_prior_pull = 2 * (log_post_ref - truth_log_post)
        Δχ²_optimizer_residual = 2 * (log_post - log_post_ref)
    end

    profile_elapsed = 0.0
    Δχ²_mis_modeling = NaN
    if chi2_mode == "mis-modeling"
        # Profile fit: NSI param(s) fixed at their injected truth value, everything else
        # (nuisances + θ₂₃/Δm²₃₁) re-optimized from the same PSO-seeded starting point.
        # Δχ²_mis-mod = χ²(fixed) - χ²(free) = -2*log_post_fixed - (-2*log_post_free)
        #             = 2*(log_post_free - log_post_fixed).
        t2 = time()
        fixed_vals = Dict(k => Float64(truth_param[k]) for k in scan_keys)
        cp_fixed = Newtrinos.condition(cp_hyp, fixed_vals, truth_param)
        prior_fixed = distprod(;cp_fixed...)
        # start_param's scan-key value(s) come from the free/unconstrained PSO seed search
        # and generally do NOT equal truth_param's -- but prior_fixed now holds scan_keys
        # constant at the truth value, and bat_findmode's ExplicitInit strictly validates
        # the init point against the prior's constants (same failure mode documented on
        # controlled_truth/_fix_conditioned in optimizer_tests_common.jl). Overwrite the
        # scan key(s) in the seed to match before fitting, or this throws ArgumentError
        # (silently caught by local_find_mle, returning NaN log_post).
        profile_start = merge(start_param, NamedTuple{Tuple(scan_keys)}(Tuple(fixed_vals[k] for k in scan_keys)))
        _, log_post_fixed, _, converged_fixed, _ = local_find_mle(truth_likelihood, prior_fixed, profile_start;
            fit_method="optim", algorithm=algorithm, iterations=lbfgs_iterations, g_tol=g_tol, f_tol=f_tol, x_tol=x_tol,
            ad_backend=ADTypes.AutoForwardDiff(), objective=objective)
        profile_elapsed = time() - t2
        Δχ²_mis_modeling = 2 * (log_post - log_post_fixed)
    end

    total_elapsed = seed_elapsed + lbfgs_elapsed + profile_elapsed + prior_pull_elapsed
    Δχ²_mod = chi2_mode == "mis-modeling" ? Δχ²_mis_modeling : Δχ²_recovery

    # Fit-minus-truth diff for the hypothesis's actually-scanned NSI parameter(s) only
    # (scan_keys -- 1 for 1D hypotheses, 2 for magnitude/phase pairs, up to 8 for "all"),
    # matching what the truth-vs-fit diagonal panels already plot.
    nsi_param_diff = Dict{Symbol,Float64}(k => Float64(fit_param[k]) - Float64(truth_param[k]) for k in scan_keys)

    results[k] = (
        truth_param=truth_param, fit_param=fit_param, converged=converged,
        log_post=log_post, truth_log_post=truth_log_post,
        Δχ²_mod=Δχ²_mod, Δχ²_recovery=Δχ²_recovery, Δχ²_mis_modeling=Δχ²_mis_modeling,
        Δχ²_prior_pull=Δχ²_prior_pull, Δχ²_optimizer_residual=Δχ²_optimizer_residual,
        n_iters=n_iters, grid_coord=grid_coords[k],
        seed_elapsed=seed_elapsed, lbfgs_elapsed=lbfgs_elapsed, profile_elapsed=profile_elapsed,
        prior_pull_elapsed=prior_pull_elapsed,
        total_elapsed=total_elapsed, nsi_param_diff=nsi_param_diff,
    )
    coord_str = if is_all
        "[" * join(round.(grid_coords[k], digits=4), ", ") * "]"
    elseif is_2d
        "($(round(grid_coords[k][1],digits=4)), $(round(grid_coords[k][2],digits=4)))"
    else
        "$(round(grid_coords[k],digits=5))"
    end
    println("  [$k/$n_grid] coord=$coord_str converged=$converged n_iters=$n_iters log_post=$(round(log_post,digits=4)) " *
            "Δχ²_mod=$(round(Δχ²_mod,digits=4)) elapsed=$(round(total_elapsed,digits=1))s")
end

n_converged = count(r -> r.converged, results)
total_wall_s = sum(r.total_elapsed for r in results)
mean_elapsed_s = total_wall_s / n_grid
min_elapsed_s, max_elapsed_s = extrema(r.total_elapsed for r in results)
println("="^80)
println("Converged: $n_converged / $n_grid")
println("Timing: total=$(round(total_wall_s,digits=1))s (summed per-point; wall clock is shorter under threading) " *
        "mean=$(round(mean_elapsed_s,digits=1))s min=$(round(min_elapsed_s,digits=1))s max=$(round(max_elapsed_s,digits=1))s")

### SAVE RESULTS ###
# Saved before plotting so the numeric results survive even if CairoMakie plotting fails.

mkpath(RESULTS_DIR)

# Flat per-key arrays (fit - truth), one Vector{Float64} per scanned NSI param, aligned
# 1:1 with truth_params/fit_params -- alongside the Vector{Dict} form above, since generic
# JLD2/HDF5 viewers can't render a Vector{Dict} nicely (it shows as an opaque reference
# table), but flat numeric vectors display normally.
nsi_param_diff_flat = Dict{String, Vector{Float64}}(
    "nsi_param_diff_$(k)" => [r.nsi_param_diff[k] for r in results] for k in scan_keys
)

FileIO.save(joinpath(RESULTS_DIR, "grid_roundtrip_$(output_tag).jld2"), Dict(
    "hypothesis"      => hypothesis_name,
    "output_tag"      => output_tag,
    "scan_keys"       => collect(scan_keys),
    "grid_coords"     => grid_coords,
    "truth_params"    => [r.truth_param for r in results],
    "fit_params"      => [r.fit_param for r in results],
    "converged"       => [r.converged for r in results],
    "log_post"        => [r.log_post for r in results],
    "truth_log_post"  => [r.truth_log_post for r in results],
    "chi2_mode"       => chi2_mode,
    "Δχ²_mod"         => [r.Δχ²_mod for r in results],
    "Δχ²_recovery"    => [r.Δχ²_recovery for r in results],
    "Δχ²_mis_modeling" => [r.Δχ²_mis_modeling for r in results],
    "Δχ²_prior_pull"        => [r.Δχ²_prior_pull for r in results],
    "Δχ²_optimizer_residual" => [r.Δχ²_optimizer_residual for r in results],
    "nsi_param_diff"  => [r.nsi_param_diff for r in results],
    nsi_param_diff_flat...,
    "n_iters"         => [r.n_iters for r in results],
    "seed_elapsed"    => [r.seed_elapsed for r in results],
    "lbfgs_elapsed"   => [r.lbfgs_elapsed for r in results],
    "profile_elapsed" => [r.profile_elapsed for r in results],
    "prior_pull_elapsed" => [r.prior_pull_elapsed for r in results],
    "total_elapsed"   => [r.total_elapsed for r in results],
    "n_converged"  => n_converged,
    "n_grid"       => n_grid,
    "fluctuate_nuisances" => fluctuate_nuisances_flag,
    "settings" => Dict(
        "objective"         => objective,
        "global_iterations" => global_iterations,
        "pso_particles"     => pso_particles,
        "iterations"        => lbfgs_iterations,
        "step_size"         => step_size,
        "g_tol"             => g_tol,
        "f_tol"             => f_tol,
        "x_tol"             => x_tol,
        "hypothesis"        => hypothesis_name,
        "n_points"          => n_points,
        "fluctuate_nuisances" => fluctuate_nuisances_flag,
        "chi2_mode"         => chi2_mode,
        "prior_effect"      => prior_effect_flag,
        "toydata"           => toydata_flag,
    ),
    # Same settings again, flattened to top-level scalar entries -- a nested Dict (like
    # "settings" above) shows as an opaque reference table in generic JLD2/HDF5 viewers
    # (same issue nsi_param_diff had); flat top-level scalars display normally.
    "settings_objective"           => objective,
    "settings_n_points"            => n_points,
    "settings_global_iterations"   => global_iterations,
    "settings_pso_particles"       => pso_particles,
    "settings_iterations"          => lbfgs_iterations,
    "settings_step_size"           => step_size,
    "settings_g_tol"               => g_tol,
    "settings_f_tol"               => f_tol,
    "settings_x_tol"               => x_tol,
    "settings_prior_effect"        => prior_effect_flag,
    "settings_toydata"             => toydata_flag,
    "settings_fluctuate_nuisances" => fluctuate_nuisances_flag,
))
println("Saved results to $(joinpath(RESULTS_DIR, "grid_roundtrip_$(output_tag).jld2"))")

### PLOT ###
# Single combined figure: NSI panel(s) first, then (if --fluctuate-nuisances) nuisance +
# physics panels, then a Δχ²_mod histogram. All panels share the same small-multiples
# layout (row, 2*col axis / 2*col+1 colorbar), one shared "converged"/"not converged"
# legend in column 1 spanning every row, Δχ²_mod colors the FIT points (converged only --
# a non-converged fit's log_post may be unreliable, so those keep a fixed marker instead).

conv_flags = [r.converged for r in results]
Δχ²_all = [r.Δχ²_mod for r in results]
# colorrange is only ever applied to CONVERGED points' colors (non-converged points always
# render as fixed red X's, never colored by Δχ²_all -- see the per-panel drawing code
# below) -- computing it from the full Δχ²_all (including non-converged points, whose
# Δχ²_mod can be NaN when local_find_mle fails) previously let extrema return (NaN, NaN)
# whenever ANY point failed to converge, even if every OTHER point converged fine. A NaN
# colorrange silently corrupts every Colorbar's internal image and crashes CairoMakie's
# PDF renderer (ArgumentError: start and stop must be finite) -- filtering to converged,
# finite values avoids that regardless of how many points failed to converge.
#
# Separately, a DEGENERATE (zero-width) range -- cmin == cmax, e.g. only 1 converged point,
# or every converged point sharing the exact same Δχ²_mod -- is finite but still crashes
# CairoMakie's colormap interpolation ("Can't interpolate in a range where cmin == cmax").
# Widen a degenerate range by a small symmetric epsilon so the colorbar always has a valid
# (nonzero) span; falls back to the same (-1.0, 1.0) default when there's nothing to plot.
Δχ²_converged_finite = filter(isfinite, Δχ²_all[conv_flags])
colorrange = if isempty(Δχ²_converged_finite)
    (-1.0, 1.0)
else
    lo, hi = extrema(Δχ²_converged_finite)
    lo == hi ? (lo - max(abs(lo), 1.0) * 1e-3 - 1e-6, hi + max(abs(hi), 1.0) * 1e-3 + 1e-6) : (lo, hi)
end

# Build the ordered list of NSI panel specs first (hypothesis-dependent), then append
# nuisance/physics panels (only if fluctuated or toydata), then the histogram panel.
# Each spec is a NamedTuple describing how to draw one panel; `kind` dispatches the draw
# logic below. `plane` panels (2D mag/phase) are handled specially since they overlay
# truth+fit as two point sets rather than a truth-vs-fit diagonal.
show_all_params = fluctuate_nuisances_flag || toydata_flag # show all parameter plots also for toydata runs
panel_specs = NamedTuple[]

if is_all
    for k in scan_keys
        push!(panel_specs, (kind=:diagonal, key=k, xlabel="truth $(k)", ylabel="fit $(k)", title="$(k)"))
    end
elseif is_2d
    mag_key, phase_key = scan_keys
    push!(panel_specs, (kind=:plane, mag_key=mag_key, phase_key=phase_key))
    push!(panel_specs, (kind=:diagonal, key=mag_key, xlabel="truth $(mag_key)", ylabel="fit $(mag_key)", title="$(mag_key) (magnitude)"))
    push!(panel_specs, (kind=:diagonal, key=phase_key, xlabel="truth $(phase_key)", ylabel="fit $(phase_key)", title="$(phase_key) (phase)"))
else
    k1 = scan_keys[1]
    push!(panel_specs, (kind=:diagonal, key=k1, xlabel="truth $(k1)", ylabel="fit $(k1)", title="$(k1)"))
end

if show_all_params
    for k in fluctuated_names
        push!(panel_specs, (kind=:diagonal, key=k, xlabel="truth $(k)", ylabel="fit $(k)", title="$(k)"))
    end
end
push!(panel_specs, (kind=:histogram,))

n_panels = length(panel_specs)
n_cols = is_all || show_all_params || is_2d ? 3 : 1
n_rows = cld(n_panels, n_cols)
fig = Figure(size=(n_cols == 1 ? 950 : 600 * n_cols, max(700, 380 * n_rows)))

plt_conv_for_legend = nothing
for (idx, spec) in enumerate(panel_specs)
    local row = (idx - 1) ÷ n_cols + 1
    local col = (idx - 1) % n_cols + 1

    if spec.kind == :diagonal
        local k = spec.key
        local ax = Axis(fig[row, 2*col], xlabel=spec.xlabel, ylabel=spec.ylabel, title=spec.title)
        local truth_vals = [r.truth_param[k] for r in results]
        local fit_vals   = [r.fit_param[k] for r in results]

        if k in nuisance_names || k in physics_names
            local d = cp_hyp[k]
            if d isa Truncated{<:Normal}
                local μ, σ = d.untruncated.μ, d.untruncated.σ
                hspan!(ax, μ - σ, μ + σ; color=(:grey, 0.25))
            end
        end

        local plt_conv = scatter!(ax, truth_vals[conv_flags], fit_vals[conv_flags];
                             color=Δχ²_all[conv_flags], colormap=:viridis, colorrange=colorrange,
                             markersize=8, label="converged")
        if any(!, conv_flags)
            scatter!(ax, truth_vals[.!conv_flags], fit_vals[.!conv_flags]; color=:red, marker=:xcross, markersize=10, label="not converged")
        end
        ablines!(ax, 0, 1, color=:red, linestyle=:dash)
        Colorbar(fig[row, 2*col + 1], plt_conv, label=Δχ²_label, ticklabelsize=9, labelsize=10)
        if plt_conv_for_legend === nothing
            global plt_conv_for_legend = plt_conv
        end

    elseif spec.kind == :plane
        local mag_key, phase_key = spec.mag_key, spec.phase_key
        local ax = Axis(fig[row, 2*col], xlabel="$(mag_key) (magnitude)", ylabel="$(phase_key) (phase)",
                         title="Truth vs fit (mag/phase plane)")
        local truth_mags   = [r.truth_param[mag_key] for r in results]
        local truth_phases = [r.truth_param[phase_key] for r in results]
        local fit_mags     = [r.fit_param[mag_key] for r in results]
        local fit_phases   = [r.fit_param[phase_key] for r in results]

        scatter!(ax, truth_mags, truth_phases; color=:steelblue, marker=:circle, markersize=10, label="truth")
        local plt_fit = scatter!(ax, fit_mags[conv_flags], fit_phases[conv_flags];
                            color=Δχ²_all[conv_flags], colormap=:viridis, colorrange=colorrange,
                            marker=:diamond, markersize=12, label="fit (converged)")
        if any(!, conv_flags)
            scatter!(ax, fit_mags[.!conv_flags], fit_phases[.!conv_flags]; color=:red, marker=:xcross, markersize=14, label="fit (not converged)")
        end
        Colorbar(fig[row, 2*col + 1], plt_fit, label=Δχ²_label, ticklabelsize=9, labelsize=10)

    elseif spec.kind == :histogram
        local ax_hist = Axis(fig[row, 2*col], xlabel=Δχ²_label, ylabel="count",
                        title="$(Δχ²_label) distribution (converged only)")
        local Δχ²_converged = filter(isfinite, Δχ²_all[conv_flags])
        if !isempty(Δχ²_converged)
            hist!(ax_hist, Δχ²_converged; bins=20, color=(:steelblue, 0.7), strokewidth=1, strokecolor=:black)
            vlines!(ax_hist, [mean(Δχ²_converged)]; color=:red, linestyle=:dash, label="mean")
            axislegend(ax_hist)
        end
    end
end

legend_elems = [
    MarkerElement(color=:steelblue, marker=:circle, markersize=10),
    MarkerElement(color=:red, marker=:xcross, markersize=10),
]
Legend(fig[1:n_rows, 1], legend_elems, ["converged", "not converged"]; tellwidth=true, tellheight=false)

title_suffix = fluctuate_nuisances_flag ? ", fluctuated nuisances" : ""
sample_kind = is_all ? "Sobol sample, " : ""
fit_kind = objective == "likelihood" ? "llh-fit" : "scan"
Label(fig[0, 1:(2*n_cols + 1)], "Truth-vs-fit $fit_kind (hypothesis: $hypothesis_name, $(sample_kind)n=$n_grid$title_suffix)", fontsize=18, font=:bold)

save(joinpath(RESULTS_DIR, "grid_roundtrip_$(output_tag).pdf"), fig)
println("Saved plot to $(joinpath(RESULTS_DIR, "grid_roundtrip_$(output_tag).pdf"))")
