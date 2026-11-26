# Roundtrip scan: inject truths (grid / Sobol sample of the tested parameters of a hypothesis of a setup file), refit with a strategy,
# and plot truth vs fit (Δχ²_recovery or Δχ²_mis-mod). All panels go into ONE figure; the results .jld2 is saved before plotting.
# Usage: julia src/analysis/roundtrips.jl --setup <file> --hypothesis <name> --strategy random lbfgsb [--workers N]  (--help)
# Options, statistics, output files: docs/src/README_optimizer_tests.md and the docs page "Optimization tests".

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

# Driver-only packages: CLI parsing happens before we know --workers/--threads (needed to
# decide whether/how to addprocs), and saving/plotting only ever run on the driver after
# the scan completes -- neither needs to be loaded on worker processes.
using ArgParse
using FileIO
using CairoMakie
using Sobol
using Statistics
using Dates
using Distributed
using SlurmClusterManager

### CLI ###

function parse_command_line()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--setup"
        help = "Setup file defining the physics, experiments, priors and hypotheses (it defines setup(; ftype) -- see nsi/setups/nsi_deepcore_3y.jl for the format)"
        arg_type = String
        required = true

        "--hypothesis"
        help = "Name of the hypothesis (from the setup file) to scan. A hypothesis with more than 2 tested parameters uses Sobol quasi-random sampling instead of a factorial grid -- use a larger --n-points, e.g. 100-200"
        arg_type = String
        default = ""

        "--results-dir"
        help = "Output folder for the saved results and plots (default: analysis_results/results/<setup name>/roundtrip, relative to the working directory)"
        arg_type = String
        default = ""

        "--n-points"
        help = "Total sample points (1 tested parameter: evenly-spaced values; 2: round(sqrt(n)) per dimension; more: Sobol quasi-random samples -- recommend 100-200)"
        arg_type = Int
        default = 25

        "--fluctuate-nuisances"
        help = "Randomly draw every other free parameter of the hypothesis (nuisances, and e.g. θ₂₃/Δm²₃₁) from its own prior (the tested parameters' truth from --hypothesis's grid/Sobol point is unaffected). Adds per-parameter truth-vs-fit panels and a Δχ²_mod histogram to the output figure. Also selects the 'fluctuated' vs 'nominal' output filename tag."
        action = :store_true

        "--chi2-recovery"
        help = "Δχ²_mod = pipeline-recovery statistic: 2*(log_post_fit - log_post_truth) (default). Mutually exclusive with --chi2-mis-modeling."
        action = :store_true

        "--chi2-mis-modeling"
        help = "Δχ²_mod = profile-likelihood mis-modeling statistic: χ²(NSI param(s) fixed at truth) - χ²(NSI param(s) free) (dissertation's actual Δχ²_mod). Requires a second LBFGS fit per grid point -- roughly doubles runtime. Mutually exclusive with --chi2-recovery."
        action = :store_true

        "--prior-effect"
        help = "Additionally compute and save Δχ²_prior_pull and Δχ²_optimizer_residual per grid point (decomposes Δχ²_recovery into expected prior-pull vs. genuine optimizer underperformance -- see docs/src/README_optimizer_tests.md). Adds one cheap local-only reference fit per grid point, seeded at the truth. Does not affect plots."
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

        "--strategy"
        help = "How the fit at every grid point is done: a chain of stages, each starting from the result of the previous one, e.g. --strategy smseed lbfgsb  |  --strategy optim-pso:particles=20,iters=50 lbfgsb:iterations=300  |  --strategy random lbfgsb  |  --strategy mle lbfgsb  |  --strategy manual:ε_eμ_abs=0.1,δ_eμ=1.0 lbfgsb  |  --strategy optim-pso:particles=50,iters=100 (a single global search, not polished). Stage kinds: seed (random, manual), alias (smseed = local fit with the tested parameters pinned at their null values; mle = random + local fit), global (optim-pso, optim-sa, mh-pso, mh-de, mh-abc, bbo, bbo-xnes, bbo-dxnes; stop at `evals` objective evaluations) and local (lbfgsb, optim-lbfgs, ...). Seed, alias and global stages are only allowed as the FIRST stage. Settings of a stage are given as name:key=val,key=val (see optimizer_tests.jl --list-optimizers); the last local stage's settings are also used for the extra fits of --prior-effect/--chi2-mis-modeling."
        arg_type = String
        nargs = '+'
        required = true

        "--start"
        help = "Start point handed to the first stage: nominal (the setup's nominal parameters, blind -- default) or truth (the injected truth, e.g. to anchor a PSO at the truth)"
        arg_type = String
        default = "nominal"

        "--nseeds"
        help = "Number of independent attempts of the whole --strategy per grid point (random stages and global searches use a fresh RNG draw per attempt, so the attempts explore differently); keeps the attempt with the highest log_post. Default 1 is a single attempt. A deterministic strategy (e.g. smseed lbfgsb) gives the same fit in every attempt."
        arg_type = Int
        default = 1

        "--workers"
        help = "Number of distributed worker processes (default: 1, no distributed -- grid points run via Threads.@threads in a single process, same as before). When >1, grid points are dispatched across worker processes via pmap, each worker internally running --threads threads across its own grid point's seed-attempts. Confirmed locally that a MODERATE number of threads per worker (e.g. 4) spread across more workers beats both a single many-threaded process (severe contention -- Julia's per-process stop-the-world GC stalls every thread together under heavy concurrent allocation) and very few many-threaded workers (reintroduces the same per-worker GC contention at smaller scale)."
        arg_type = Int
        default = 1

        "--threads"
        help = "Threads per worker process (only used when --workers > 1; passed to addprocs' exeflags). Ignored (and irrelevant) when --workers is 1 -- that single-process fallback path's own concurrency is fixed by whatever --threads=N the whole `julia` invocation was started with."
        arg_type = Int
        default = 4

        "--multi-node"
        help = "Spawn workers across multiple SLURM-allocated nodes (e.g. LRZ cm4_std) via SlurmClusterManager.jl instead of all on this single machine (plain addprocs(n_workers), the cm4_tiny-compatible default). Only valid when --workers > 1 AND running inside an active SLURM allocation (SLURM_JOB_ID/SLURM_JOBID and SLURM_NTASKS must be set -- this is NOT auto-detected; the flag must be passed explicitly, and the SBATCH script's own --ntasks must equal --workers exactly or this errors out before spawning anything, rather than silently using the wrong worker count. Worker count on multi-node is controlled by the SBATCH script's --ntasks, not by --workers directly -- --workers is only used here as a consistency check."
        action = :store_true

    end
    return parse_args(s)
end

args = parse_command_line()
hypothesis_name = args["hypothesis"]
fluctuate_nuisances_flag = args["fluctuate-nuisances"]
n_points = args["n-points"]
objective = lowercase(args["objective"])
strategy = String[t for s in args["strategy"] for t in split(s)]   # the stages, as separate arguments or in one quoted string (as in sensitivity.jl)
start_mode = lowercase(args["start"])
prior_effect_flag = args["prior-effect"]
toydata_flag = args["toydata"]
nseeds = args["nseeds"]
n_workers = args["workers"]
n_threads = args["threads"]
multi_node_flag = args["multi-node"]
use_distributed = n_workers > 1

multi_node_flag && !use_distributed &&
    error("--multi-node requires --workers > 1 (got --workers $n_workers).")

args["chi2-recovery"] && args["chi2-mis-modeling"] &&
    error("--chi2-recovery and --chi2-mis-modeling are mutually exclusive -- pick one.")
chi2_mode = args["chi2-mis-modeling"] ? "mis-modeling" : "recovery"  # recovery is default
objective in ("posterior", "likelihood") || error("--objective must be 'posterior' or 'likelihood', got '$objective'")
start_mode in ("nominal", "truth") || error("--start must be 'nominal' or 'truth', got '$start_mode'")

Δχ²_label = rich("Δχ", superscript("2"), subscript(chi2_mode == "mis-modeling" ? "mis-mod" : "recovery"))

### DISTRIBUTED WORKERS (optional) ###
# Workers are spawned BEFORE loading packages, which are then loaded with @everywhere. Each worker gets --threads=n_threads: a few multi-threaded
# processes avoid the contention of one process with many threads (Julia's GC stops all threads of a process).
# --multi-node: SlurmManager() reads SLURM_NTASKS itself, so --workers is only checked against it (fails before spawning on a mismatch).
if use_distributed
    if multi_node_flag
        slurm_ntasks_str = get(ENV, "SLURM_NTASKS", nothing)
        slurm_ntasks_str === nothing &&
            error("--multi-node requires an active SLURM allocation with SLURM_NTASKS set " *
                  "(not found in ENV) -- SlurmManager() cannot be constructed outside sbatch/salloc.")
        slurm_ntasks = parse(Int, slurm_ntasks_str)
        slurm_ntasks == n_workers ||
            error("--multi-node: --workers ($n_workers) must equal the SBATCH script's " *
                  "--ntasks (SLURM_NTASKS=$slurm_ntasks) -- fix one or the other, don't let them drift.")
        addprocs(SlurmManager(); exeflags="--threads=$n_threads")
    else
        addprocs(n_workers; exeflags="--threads=$n_threads")
    end
end

@everywhere begin
    using Distributions
    using DensityInterface
    using BAT
    using DataStructures
    using MeasureBase
    using ADTypes
    using Newtrinos
    using Random
    import Optim
    import ForwardDiff
    using ValueShapes
end

# @reset must be expanded in a separate @everywhere block from `using Accessors` (world age).
@everywhere begin
    using Accessors
    include(joinpath(@__DIR__, "optimizer_tests_common.jl"))
end

### SETUP (physics, experiments, priors, hypotheses: defined by the setup file, see --setup) ###

include(abspath(args["setup"]))
S = finalize_setup(setup(; ftype = Float64))
hypothesis_name = isempty(hypothesis_name) ? first(keys(S.hypotheses)) : hypothesis_name   # default: the first hypothesis of the setup
hyp = hypothesis_info(S, hypothesis_name)

data_mode_tag = toydata_flag ? "toydata" : "asimov"
truth_mode_tag = fluctuate_nuisances_flag ? "fluctuated" : "nominal"
objective_tag = objective == "likelihood" ? "likelihoodfit" : "posteriorfit"
output_tag = isempty(args["suffix"]) ? "$(data_mode_tag)_$(truth_mode_tag)_$(hypothesis_name)_$(objective_tag)" : "$(data_mode_tag)_$(truth_mode_tag)_$(hypothesis_name)_$(objective_tag)_$(args["suffix"])"
output_tag = nseeds == 1 ? output_tag : output_tag * "_nseeds$(nseeds)"
output_tag = output_tag * "_" * join((first(split(s, ':')) for s in strategy), "+")   # stage names of the strategy
output_tag = output_tag * "_" * Dates.format(Dates.today(), "ddmmyyyy")             # execution date ends every file name

experiments = S.experiments
p = S.params
cp_hyp = hyp.prior                  # prior of the hypothesis under test
cp_sm = hyp.null_prior              # the same prior with the tested parameters at their null values
tested_names = hyp.tested           # tested parameters = keys of the hypothesis' `null`
fluctuated_names = hyp.fluctuated   # all other free parameters: redrawn by --fluctuate-nuisances, shown in the extra panels
results_dir = isempty(args["results-dir"]) ? joinpath("analysis_results/results", S.name, "roundtrip") : args["results-dir"]

### BUILD THE TRUTH GRID ###

truths, scan_keys, grid_coords = truth_grid(cp_hyp, p, tested_names; n_points=n_points)
n_grid = length(truths)
is_many = length(scan_keys) > 2   # Sobol sample instead of a grid
is_2d = length(scan_keys) == 2

println("Truth grid built for hypothesis '$hypothesis_name': $(n_grid) points, scanning $(scan_keys)")

prior_dist = distprod(;cp_hyp...)

# The strategy is parsed and checked once here (before any expensive work); `final_local_spec` is the local fit that the
# extra fits of --prior-effect/--chi2-mis-modeling use (the last local stage of the strategy, with ad = :forward as in the
# main fits, since the grid points / attempts already run concurrently).
normalize_strategy(strategy)
local_defaults = (; ad = :forward)
extra_spec = with_local_defaults(final_local_spec(strategy), local_defaults)

# `run_grid_point` needs these on whichever process runs it: a plain @everywhere function does not see driver-side variables, so they are
# bundled into one NamedTuple `ctx` and broadcast once (not re-sent per call).
ctx = (; experiments, prior_dist, cp_hyp, cp_sm, p,
       nseeds, strategy, start_mode, local_defaults, extra_spec, objective,
       prior_effect_flag, chi2_mode, scan_keys,
       fluctuated_names, truths, grid_coords, toydata_flag, n_grid, is_many, is_2d,
       fluctuate_nuisances_flag)
use_distributed && @everywhere ctx = $ctx

### ROUNDTRIP SCAN (PSO seed + LBFGS polish at every grid point, parallel) ###

println("="^80)
println("Running $(n_grid) grid-point roundtrips (strategy: $(join(strategy, " ")), $nseeds attempt(s) each) on " *
        (use_distributed ? "$n_workers workers x $n_threads threads each" : "$(Threads.nthreads()) threads") *
        " (Δχ²_mod mode: $chi2_mode)...")

# One grid point's roundtrip: truth-injected likelihood, `nseeds` attempts of the strategy, best attempt, then the optional
# --prior-effect and --chi2-mis-modeling fits. Called from pmap (--workers > 1) or a Threads.@threads loop; defined with @everywhere in both cases.
@everywhere function run_grid_point(k)
    experiments, prior_dist, cp_hyp, cp_sm, p,
        nseeds, strategy, start_mode, local_defaults, extra_spec, objective,
        prior_effect_flag, chi2_mode, scan_keys,
        fluctuated_names, truths, grid_coords, toydata_flag, n_grid, is_many, is_2d,
        fluctuate_nuisances_flag = ctx
    seed_rng = Random.Xoshiro(1000 + k)
    truth_param = fluctuate_nuisances_flag ?
        fluctuate_nuisances(seed_rng, truths[k], cp_hyp, fluctuated_names) : truths[k]

    t0 = time()
    # generate_toy_data draws from the global RNG: seed it with an integer (Random.seed!(rng) without a value is not reproducible).
    toydata_flag && Random.seed!(1000 + k)
    injected_data = toydata_flag ?
        Newtrinos.generate_toy_data(experiments, truth_param) : Newtrinos.generate_asimov_data(experiments, truth_param)
    truth_likelihood = Newtrinos.generate_likelihood(experiments, injected_data)
    truth_posterior = PosteriorMeasure(truth_likelihood, prior_dist)
    truth_log_post = logdensityof(truth_posterior, truth_param)

    # `nseeds` independent attempts of the strategy, concurrently (Threads.@threads); the best by log_post is kept (NaN counts as -Inf).
    # Attempt s has its own RNG Xoshiro(1000k + s): reproducible per (k, s), no shared RNG state between threads.
    start_param = start_mode == "truth" ? truth_param : p
    attempt_results = Vector{Any}(undef, nseeds)
    Threads.@threads for s in 1:nseeds
        attempt_results[s] = find_mle_advanced(truth_likelihood, cp_hyp, start_param; strategy, nseeds = 1,
                                               rng = Random.Xoshiro(1000 * k + s), objective, null_prior = cp_sm, local_defaults)
    end

    # Keep the best attempt by log_post and its own seed point (needed for --chi2-mis-modeling's profile_start); if every attempt failed,
    # fall back to a NaN failure record so that this grid point is recorded as non-convergent.
    best_rec = nothing
    best_log_post = -Inf
    for s in 1:nseeds
        rec = attempt_results[s]
        log_post_s = isnan(rec.log_post) ? -Inf : rec.log_post
        if log_post_s > best_log_post
            best_rec, best_log_post = rec, log_post_s
        end
    end
    if best_rec === nothing
        nan_params = NamedTuple{keys(prior_dist)}(ntuple(_ -> NaN, length(keys(prior_dist))))
        best_result = (NaN, NaN, nan_params, false, missing)
        best_start_param = nan_params
        seed_elapsed, lbfgs_elapsed = 0.0, 0.0
    else
        best_result = (best_rec.llh, best_rec.log_post, best_rec.params, best_rec.converged, best_rec.n_iters)
        best_start_param = best_rec.seed_params
        seed_elapsed = sum((r.time for r in best_rec.stages[1:end-1]); init = 0.0)   # all stages but the last
        lbfgs_elapsed = best_rec.stages[end].time                                      # the last stage
    end
    llh, log_post, fit_param, converged, n_iters = best_result

    # Δχ²_recovery: always computed; with objective="likelihood" the fit target is a flat-prior pseudo-posterior, so llh is compared instead.
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
        # Cheap local fit started at the truth, always with objective="posterior" (prior pull is a posterior-only effect): estimates the local posterior maximum near the truth, whatever the main fit optimized.
        t3 = time()
        log_post_ref = run_optimizer(extra_spec, truth_likelihood, prior_dist, truth_param; objective="posterior").log_post
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
        # The seed's scan key(s) differ from truth_param's, but bat_findmode's ExplicitInit rejects a start that differs from the prior's
        # constants (see controlled_truth): overwrite them with the truth values (using the WINNING attempt's seed).
        profile_start = merge(best_start_param, NamedTuple{Tuple(scan_keys)}(Tuple(fixed_vals[k] for k in scan_keys)))
        log_post_fixed = run_optimizer(extra_spec, truth_likelihood, prior_fixed, profile_start; objective=objective).log_post
        profile_elapsed = time() - t2
        Δχ²_mis_modeling = 2 * (log_post - log_post_fixed)
    end

    total_elapsed = seed_elapsed + lbfgs_elapsed + profile_elapsed + prior_pull_elapsed
    Δχ²_mod = chi2_mode == "mis-modeling" ? Δχ²_mis_modeling : Δχ²_recovery

    # Fit-minus-truth diff for the hypothesis's actually-scanned (tested) parameter(s) only
    # (scan_keys -- 1 for 1D hypotheses, 2 for magnitude/phase pairs, up to 8 for "all"),
    # matching what the truth-vs-fit diagonal panels already plot.
    tested_param_diff = Dict{Symbol,Float64}(k => Float64(fit_param[k]) - Float64(truth_param[k]) for k in scan_keys)

    result = (
        truth_param=truth_param, fit_param=fit_param, converged=converged,
        log_post=log_post, truth_log_post=truth_log_post,
        Δχ²_mod=Δχ²_mod, Δχ²_recovery=Δχ²_recovery, Δχ²_mis_modeling=Δχ²_mis_modeling,
        Δχ²_prior_pull=Δχ²_prior_pull, Δχ²_optimizer_residual=Δχ²_optimizer_residual,
        n_iters=n_iters, grid_coord=grid_coords[k],
        seed_elapsed=seed_elapsed, lbfgs_elapsed=lbfgs_elapsed, profile_elapsed=profile_elapsed,
        prior_pull_elapsed=prior_pull_elapsed,
        total_elapsed=total_elapsed, tested_param_diff=tested_param_diff,
    )
    coord_str = if is_many
        "[" * join(round.(grid_coords[k], digits=4), ", ") * "]"
    elseif is_2d
        "($(round(grid_coords[k][1],digits=4)), $(round(grid_coords[k][2],digits=4)))"
    else
        "$(round(grid_coords[k],digits=5))"
    end
    println("  [$k/$n_grid] coord=$coord_str converged=$converged n_iters=$n_iters log_post=$(round(log_post,digits=4)) " *
            "Δχ²_mod=$(round(Δχ²_mod,digits=4)) elapsed=$(round(total_elapsed,digits=1))s")
    return result
end

# --workers > 1: grid points over the workers with pmap (dynamic load balancing; the per-point cost varies a lot), each worker threads over
# the attempts of its point. --workers 1: Threads.@threads over the grid points.
results = if use_distributed
    pmap(run_grid_point, 1:n_grid)
else
    results_local = Vector{NamedTuple}(undef, n_grid)
    Threads.@threads for k in 1:n_grid
        results_local[k] = run_grid_point(k)
    end
    results_local
end

n_converged = count(r -> r.converged === true, results)
n_noflag = count(r -> ismissing(r.converged), results)   # fits without a convergence flag (global-only strategies)
total_wall_s = sum(r.total_elapsed for r in results)
mean_elapsed_s = total_wall_s / n_grid
min_elapsed_s, max_elapsed_s = extrema(r.total_elapsed for r in results)
println("="^80)
println("Converged: $n_converged / $n_grid" * (n_noflag > 0 ? " ($n_noflag without a convergence flag)" : ""))
println("Timing: total=$(round(total_wall_s,digits=1))s (summed per-point; wall clock is shorter under threading) " *
        "mean=$(round(mean_elapsed_s,digits=1))s min=$(round(min_elapsed_s,digits=1))s max=$(round(max_elapsed_s,digits=1))s")

### SAVE RESULTS ###
# Saved before plotting so the numeric results survive even if CairoMakie plotting fails.

mkpath(results_dir)

# Flat per-key arrays (fit - truth), one Vector{Float64} per scanned NSI param, aligned
# 1:1 with truth_params/fit_params -- alongside the Vector{Dict} form above, since generic
# JLD2/HDF5 viewers can't render a Vector{Dict} nicely (it shows as an opaque reference
# table), but flat numeric vectors display normally.
tested_param_diff_flat = Dict{String, Vector{Float64}}(
    "tested_param_diff_$(k)" => [r.tested_param_diff[k] for r in results] for k in scan_keys
)

FileIO.save(joinpath(results_dir, "grid_roundtrip_$(output_tag).jld2"), Dict(
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
    "tested_param_diff"  => [r.tested_param_diff for r in results],
    "setup"           => S.name,
    tested_param_diff_flat...,
    "n_iters"         => [r.n_iters for r in results],
    "seed_elapsed"    => [r.seed_elapsed for r in results],
    "lbfgs_elapsed"   => [r.lbfgs_elapsed for r in results],
    "profile_elapsed" => [r.profile_elapsed for r in results],
    "prior_pull_elapsed" => [r.prior_pull_elapsed for r in results],
    "total_elapsed"   => [r.total_elapsed for r in results],
    "n_converged"  => n_converged,
    "n_no_convergence_flag" => n_noflag,
    "n_grid"       => n_grid,
    "fluctuate_nuisances" => fluctuate_nuisances_flag,
    "settings" => Dict(
        "objective"         => objective,
        "strategy"          => strategy,
        "start"             => start_mode,
        "hypothesis"        => hypothesis_name,
        "n_points"          => n_points,
        "fluctuate_nuisances" => fluctuate_nuisances_flag,
        "chi2_mode"         => chi2_mode,
        "prior_effect"      => prior_effect_flag,
        "toydata"           => toydata_flag,
        "nseeds"            => nseeds,
        "setup"             => S.name,
    ),
    # Same settings again, flattened to top-level scalar entries -- a nested Dict (like
    # "settings" above) shows as an opaque reference table in generic JLD2/HDF5 viewers
    # (same issue tested_param_diff had); flat top-level scalars display normally.
    "settings_objective"           => objective,
    "settings_n_points"            => n_points,
    "settings_strategy"            => join(strategy, " "),
    "settings_start"               => start_mode,
    "settings_prior_effect"        => prior_effect_flag,
    "settings_toydata"             => toydata_flag,
    "settings_fluctuate_nuisances" => fluctuate_nuisances_flag,
    "settings_nseeds"              => nseeds,
    "settings_setup"               => S.name,
))
println("Saved results to $(joinpath(results_dir, "grid_roundtrip_$(output_tag).jld2"))")

### PLOT ###
# One combined figure: tested-parameter panel(s), then (with --fluctuate-nuisances) nuisance/physics panels, then the Δχ²_mod histogram. Marker:
# circle/diamond = converged, cross = not converged, square = no convergence flag (global-only strategy); failed fits (not finite) are not drawn.

conv_flags = [r.converged === true for r in results]
noflag_flags = [ismissing(r.converged) for r in results]   # no convergence criterion (global-only strategy)
Δχ²_all = [r.Δχ²_mod for r in results]
# colorrange from the FINITE Δχ²_mod values only, widened if degenerate (checked in Float32, as CairoMakie casts down to Float32): a NaN or
# zero-width range crashes CairoMakie's PDF renderer (see "Implementation notes" in docs/src/README_optimizer_tests.md).
fin_flags = isfinite.(Δχ²_all)
sel_conv = conv_flags .& fin_flags        # converged, drawn as circles/diamonds
sel_nonc = (.!conv_flags .& .!noflag_flags) .& fin_flags   # not converged, drawn as crosses in the same color scale
sel_noflag = noflag_flags .& fin_flags    # no convergence flag, drawn as squares in the same color scale
sel_fail = .!fin_flags                    # fit failed (Δχ²_mod not finite), grey crosses
Δχ²_finite = Δχ²_all[fin_flags]
colorrange = if isempty(Δχ²_finite)
    (-1.0, 1.0)
else
    lo, hi = extrema(Δχ²_finite)
    Float32(lo) == Float32(hi) ? (lo - max(abs(lo), 1.0) * 1e-3 - 1e-6, hi + max(abs(hi), 1.0) * 1e-3 + 1e-6) : (lo, hi)
end

# Panel specs: first the tested-parameter panels, then nuisance/physics panels (only if fluctuated or toydata), then the histogram; `plane` panels
# (2D) overlay truth and fit as two point sets.
show_all_params = fluctuate_nuisances_flag || toydata_flag # show all parameter plots also for toydata runs
panel_specs = NamedTuple[]

if is_many
    for k in scan_keys
        push!(panel_specs, (kind=:diagonal, key=k, xlabel="truth $(k)", ylabel="fit $(k)", title="$(k)"))
    end
elseif is_2d
    mag_key, phase_key = scan_keys
    push!(panel_specs, (kind=:plane, mag_key=mag_key, phase_key=phase_key))
    push!(panel_specs, (kind=:diagonal, key=mag_key, xlabel="truth $(mag_key)", ylabel="fit $(mag_key)", title="$(mag_key)"))
    push!(panel_specs, (kind=:diagonal, key=phase_key, xlabel="truth $(phase_key)", ylabel="fit $(phase_key)", title="$(phase_key)"))
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
n_cols = is_many || show_all_params || is_2d ? 3 : 1
n_rows = cld(n_panels, n_cols)
fig = Figure(size=(n_cols == 1 ? 950 : 600 * n_cols, max(700, 380 * n_rows)))

plt_conv_for_legend = nothing
fail_drawn = false   # does any panel show a failed fit?
for (idx, spec) in enumerate(panel_specs)
    local row = (idx - 1) ÷ n_cols + 1
    local col = (idx - 1) % n_cols + 1

    if spec.kind == :diagonal
        local k = spec.key
        local ax = Axis(fig[row, 2*col], xlabel=spec.xlabel, ylabel=spec.ylabel, title=spec.title)
        local truth_vals = [r.truth_param[k] for r in results]
        local fit_vals   = [r.fit_param[k] for r in results]

        if k in fluctuated_names
            local d = cp_hyp[k]
            if d isa Truncated{<:Normal}
                local μ, σ = d.untruncated.μ, d.untruncated.σ
                hspan!(ax, μ - σ, μ + σ; color=(:grey, 0.25))
            end
        end

        local plt_conv = scatter!(ax, truth_vals[sel_conv], fit_vals[sel_conv];
                             color=Δχ²_all[sel_conv], colormap=Reverse(:viridis), colorrange=colorrange,
                             markersize=8, label="converged")
        any(sel_nonc) && scatter!(ax, truth_vals[sel_nonc], fit_vals[sel_nonc]; color=Δχ²_all[sel_nonc], colormap=Reverse(:viridis),
                                  colorrange=colorrange, marker=:xcross, markersize=10, label="not converged")
        any(sel_noflag) && scatter!(ax, truth_vals[sel_noflag], fit_vals[sel_noflag]; color=Δχ²_all[sel_noflag], colormap=Reverse(:viridis),
                                    colorrange=colorrange, marker=:rect, markersize=8, label="no convergence flag")
        # failed fits (Δχ² not finite) usually have NaN fit parameters, which cannot be drawn: only those with finite coordinates are shown
        local failed_here = sel_fail .& isfinite.(truth_vals) .& isfinite.(fit_vals)
        if any(failed_here)
            scatter!(ax, truth_vals[failed_here], fit_vals[failed_here]; color=:red, marker=:xcross, markersize=10)
            global fail_drawn = true
        end
        ablines!(ax, 0, 1, color=:red, linestyle=:dash)
        Colorbar(fig[row, 2*col + 1], plt_conv, label=Δχ²_label, ticklabelsize=9, labelsize=10)
        if plt_conv_for_legend === nothing
            global plt_conv_for_legend = plt_conv
        end

    elseif spec.kind == :plane
        local mag_key, phase_key = spec.mag_key, spec.phase_key
        local ax = Axis(fig[row, 2*col], xlabel="$(mag_key)", ylabel="$(phase_key)",
                         title="Truth vs fit ($(mag_key) - $(phase_key) plane)")
        local truth_mags   = [r.truth_param[mag_key] for r in results]
        local truth_phases = [r.truth_param[phase_key] for r in results]
        local fit_mags     = [r.fit_param[mag_key] for r in results]
        local fit_phases   = [r.fit_param[phase_key] for r in results]

        scatter!(ax, truth_mags, truth_phases; color=:grey60, marker=:circle, markersize=10, label="truth")
        local plt_fit = scatter!(ax, fit_mags[sel_conv], fit_phases[sel_conv];
                            color=Δχ²_all[sel_conv], colormap=Reverse(:viridis), colorrange=colorrange,
                            marker=:diamond, markersize=12, label="fit (converged)")
        any(sel_nonc) && scatter!(ax, fit_mags[sel_nonc], fit_phases[sel_nonc]; color=Δχ²_all[sel_nonc], colormap=Reverse(:viridis),
                                  colorrange=colorrange, marker=:xcross, markersize=14, label="fit (not converged)")
        any(sel_noflag) && scatter!(ax, fit_mags[sel_noflag], fit_phases[sel_noflag]; color=Δχ²_all[sel_noflag], colormap=Reverse(:viridis),
                                    colorrange=colorrange, marker=:rect, markersize=11, label="fit (no convergence flag)")
        local failed_plane = sel_fail .& isfinite.(fit_mags) .& isfinite.(fit_phases)
        if any(failed_plane)
            scatter!(ax, fit_mags[failed_plane], fit_phases[failed_plane]; color=:red, marker=:xcross, markersize=14)
            global fail_drawn = true
        end
        Colorbar(fig[row, 2*col + 1], plt_fit, label=Δχ²_label, ticklabelsize=9, labelsize=10)

    elseif spec.kind == :histogram
        local ax_hist = Axis(fig[row, 2*col], xlabel=Δχ²_label, ylabel="count",
                        title=rich(Δχ²_label, " distribution (all fits with a finite value)"))
        local Δχ²_converged = filter(isfinite, Δχ²_all)   # converged, not converged and no-flag fits; only failed fits are left out
        if !isempty(Δχ²_converged)
            hist!(ax_hist, Δχ²_converged; bins=20, color=(:steelblue, 0.7), strokewidth=1, strokecolor=:black)
            vlines!(ax_hist, [mean(Δχ²_converged)]; color=:red, linestyle=:dash, label="mean")
            axislegend(ax_hist)
        end
    end
end

legend_elems = [
    MarkerElement(color=:steelblue, marker=:circle, markersize=10),
    MarkerElement(color=:steelblue, marker=:xcross, markersize=10),
    MarkerElement(color=:steelblue, marker=:rect, markersize=9),
]
legend_labels = ["converged", "not converged", "no convergence flag"]
if fail_drawn   # the legend entry only appears if a failed fit is actually drawn
    push!(legend_elems, MarkerElement(color=:red, marker=:xcross, markersize=10))
    push!(legend_labels, "fit failed (no Δχ²)")
end
Legend(fig[1:n_rows, 1], legend_elems, legend_labels; tellwidth=true, tellheight=false)

title_suffix = fluctuate_nuisances_flag ? ", fluctuated nuisances" : ""
sample_kind = is_many ? "Sobol sample, " : ""
fit_kind = objective == "likelihood" ? "llh-fit" : "scan"
Label(fig[0, 1:(2*n_cols + 1)], "Truth-vs-fit $fit_kind (hypothesis: $hypothesis_name, $(sample_kind)n=$n_grid$title_suffix)", fontsize=18, font=:bold)

save(joinpath(results_dir, "grid_roundtrip_$(output_tag).pdf"), fig)
println("Saved plot to $(joinpath(results_dir, "grid_roundtrip_$(output_tag).pdf"))")
