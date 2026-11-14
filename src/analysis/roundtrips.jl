# Generic roundtrip scan. Physics, experiments, priors and the hypotheses to test come from a setup
# file (--setup, e.g. nsi/setups/nsi_deepcore_3y.jl; see there for how a hypothesis is specified); this
# script contains no model-specific code. --hypothesis selects one of the setup's hypotheses.
#
# Thesis-style roundtrip scan: replicates the dissertation's "recovering injected NSI
# hypotheses" test (section 7.6.2, nuisances/theta23/deltam31 pinned at nominal) and,
# with --fluctuate-nuisances, its "recovering fluctuated nuisance parameters" test
# (section 7.6.3, nuisances + theta23/deltam31 redrawn from their own priors per grid
# point). All truth-vs-fit panels (NSI scan panel(s), nuisance/physics panels, and the
# Δχ²_mod histogram) are combined into ONE figure/ONE output file per run -- no separate
# "_nuisances.pdf" file.
#
# Per-hypothesis panel(s), by the number of tested parameters of the hypothesis (keys of its `null`):
#   - 1 tested parameter: one truth-vs-fit diagonal panel.
#   - 2 tested parameters (e.g. a magnitude+phase pair): THREE panels -- the plane of the two
#     parameters with truth and fit as separate point sets (matching the dissertation's Fig.
#     7.14b layout), plus two diagonal truth-vs-fit panels (one per parameter), same visual
#     language as the nuisance panels.
#   - more than 2 (Sobol sample): one truth-vs-fit diagonal panel per tested parameter.
#
# The GRID controls only how truths are injected (via truth_grid, optimizer_tests_common.jl)
# -- it does NOT restrict the fit. Every fit is still free to vary every non-fixed
# parameter in the hypothesis's prior (all physics params + all nuisances), exactly
# like any other roundtrip; only the injected truth is simplified to isolate the scan to
# the NSI parameter(s) being tested (plus, under --fluctuate-nuisances, the nuisances and
# θ₂₃/Δm²₃₁).
#
# --setup: setup file defining physics, experiments, priors and hypotheses (required).
# --hypothesis: name of one of the setup's hypotheses.
# --results-dir: output folder (default results/<setup name>/roundtrip, relative to the working directory).
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
# --strategy: the fit at every grid point is a chain of stages (find_mle_advanced in optimizer_tests_common.jl), e.g.
#   --strategy smseed lbfgsb | --strategy optim-pso:particles=20,iters=50 lbfgsb | --strategy random lbfgsb. Every stage starts from
#   the result of the previous one; seed, alias and global stages (random, manual, smseed, mle, optim-pso, ...) are only allowed first.
#   The stage settings (iterations, tolerances, particles, ...) are part of the stage spec and recorded in the saved .jld2's
#   "settings" block. A strategy that is a single global search is not polished and has no convergence flag
#   (`converged = missing`, drawn with its own marker). --start chooses the start point of the first stage (nominal or truth).
#   Local fits should normally use lbfgsb (L-BFGS-B): plain optim-lbfgs jumps onto the prior bounds from a start far from the
#   optimum and then stops with a vanishing gradient.
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
# Usage: julia --threads=8 src/analysis/roundtrips.jl --setup nsi/setups/nsi_deepcore_3y.jl --hypothesis ττ_μμ
#
# Output directory: --results-dir (created if missing).
# Output filenames: grid_roundtrip_<asimov|toydata>_<fluctuated|nominal>_<hypothesis>_<posteriorfit|likelihoodfit>[_<suffix>].{pdf,jld2}
# -- "asimov"/"toydata" reflects whether --toydata was set, "fluctuated"/"nominal" reflects
# whether --fluctuate-nuisances was set, "posteriorfit"/"likelihoodfit" reflects --objective;
# --suffix is freeform run-tagging (e.g. "run1").
# Results (.jld2) are saved BEFORE plotting, so they survive even if plotting fails.

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
        default = "ττ_μμ"

        "--results-dir"
        help = "Output folder for the saved results and plots (default: results/<setup name>/roundtrip, relative to the working directory)"
        arg_type = String
        default = ""

        "--n-points"
        help = "Total sample points (1 tested parameter: evenly-spaced values; 2: round(sqrt(n)) per dimension; more: Sobol quasi-random samples -- recommend 100-200)"
        arg_type = Int
        default = 25

        "--fluctuate-nuisances"
        help = "Replicate the dissertation's section 7.6.3 test: also randomly draw every other free parameter of the hypothesis (nuisances, and e.g. θ₂₃/Δm²₃₁) from its own prior (the tested parameters' truth from --hypothesis's grid/Sobol point is unaffected). Adds per-parameter truth-vs-fit panels and a Δχ²_mod histogram to the output figure. Also selects the 'fluctuated' vs 'nominal' output filename tag."
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
strategy = args["strategy"]
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

data_mode_tag = toydata_flag ? "toydata" : "asimov"
truth_mode_tag = fluctuate_nuisances_flag ? "fluctuated" : "nominal"
objective_tag = objective == "likelihood" ? "likelihoodfit" : "posteriorfit"
output_tag = isempty(args["suffix"]) ? "$(data_mode_tag)_$(truth_mode_tag)_$(hypothesis_name)_$(objective_tag)" : "$(data_mode_tag)_$(truth_mode_tag)_$(hypothesis_name)_$(objective_tag)_$(args["suffix"])"
output_tag = nseeds == 1 ? output_tag : output_tag * "_nseeds$(nseeds)"
output_tag = output_tag * "_" * join((first(split(s, ':')) for s in strategy), "+")   # stage names of the strategy
output_tag = output_tag * "_" * Dates.format(Dates.today(), "ddmmyyyy")             # execution date ends every file name

Δχ²_label = rich("Δχ", superscript("2"), subscript(chi2_mode == "mis-modeling" ? "mis-mod" : "recovery"))

### DISTRIBUTED WORKERS (optional) ###
# Spawn worker processes BEFORE loading any package they'll need, then load those
# packages via @everywhere so both the driver and every worker have them. Each worker
# gets its own --threads=n_threads (passed through addprocs' exeflags), so grid points
# dispatched to it (via pmap below) can still run their own seed-attempts concurrently
# via Threads.@threads internally -- this hybrid shape (a moderate number of separate
# processes, each multi-threaded) was confirmed locally to avoid the severe contention
# that a single process with MANY threads hits under heavy concurrent PSO/LBFGS
# allocation (Julia's current GC is stop-the-world PER PROCESS -- every thread in that
# process stalls together on every GC pause; splitting into a few processes confines
# each pause to that process's own threads only).
#
# --multi-node switches from plain addprocs(n_workers) (all workers on THIS single
# machine -- correct for cm4_tiny) to SlurmClusterManager.jl's SlurmManager(), which
# launches workers across every node SLURM actually allocated to this job (e.g.
# cm4_std's multi-node jobs) via `srun ... --worker` internally. SlurmManager() takes
# NO worker-count argument -- it reads ENV["SLURM_NTASKS"] itself (requires an active
# SLURM allocation; SLURM_JOB_ID/SLURM_JOBID and SLURM_NTASKS must be set, or its own
# constructor errors). Worker count is therefore fixed by the SBATCH script's own
# --ntasks, NOT by --workers directly here -- --workers is only used as an explicit
# consistency check below, so a mismatched SBATCH script (--ntasks != --workers) fails
# loudly before spawning anything, rather than silently running with the wrong worker
# count.
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

# `@reset` (Accessors' macro, used below for prior overrides) must be expanded in a
# SEPARATE top-level @everywhere block from `using Accessors` -- Julia errors with
# "UndefVarError: @reset not defined ... binding may be too new" if both are evaluated
# together as one atomic block (world age hasn't advanced yet to see the macro). In
# practice `@reset` itself is only ever called on the DRIVER (the priors/`cp_*`
# construction below never runs on a worker), but `using Accessors` still needs to be
# loaded everywhere `optimizer_tests_common.jl`'s own code might reference it internally,
# so keep this split defensively rather than relying on that happening to be true today.
@everywhere begin
    using Accessors
    include(joinpath(@__DIR__, "optimizer_tests_common.jl"))
end

### SETUP (physics, experiments, priors, hypotheses: defined by the setup file, see --setup) ###

include(abspath(args["setup"]))
S = finalize_setup(setup(; ftype = Float64))
hyp = hypothesis_info(S, hypothesis_name)

experiments = S.experiments
p = S.params
cp_hyp = hyp.prior                  # prior of the hypothesis under test
cp_sm = hyp.null_prior              # the same prior with the tested parameters at their null values
tested_names = hyp.tested           # tested parameters = keys of the hypothesis' `null`
fluctuated_names = hyp.fluctuated   # all other free parameters: redrawn by --fluctuate-nuisances, shown in the extra panels
results_dir = isempty(args["results-dir"]) ? joinpath("results", S.name, "roundtrip") : args["results-dir"]

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

# `run_grid_point` below needs every one of these on whichever process actually executes
# it. A plain `@everywhere function run_grid_point(k) ... end` does NOT automatically see
# driver-side free variables the way a pmap do-block closure would -- confirmed directly:
# referencing an un-broadcast driver variable from inside an @everywhere function throws
# UndefVarError on the worker. Bundling everything into one NamedTuple and broadcasting it
# ONCE (rather than ~25 individual `@everywhere x = $x` statements, or re-serializing
# large objects like `experiments` on every one of the n_grid calls) is the simplest fix;
# `run_grid_point` below references the global `ctx` directly (resolved to each worker's
# own already-broadcast copy -- confirmed this does NOT re-send ctx per call). When
# `!use_distributed`, this is just a plain global in the single process, no different from
# today's top-level variables.
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

# One grid point's full roundtrip: build its truth-injected likelihood, run n_attempts
# independent PSO-seed + LBFGS-polish attempts CONCURRENTLY via the inner Threads.@threads
# (one Julia thread per attempt -- this is the actual expensive work), reduce to the best
# attempt, then (optionally) the --prior-effect and --chi2-mis-modeling extra fits.
# Called either from pmap (--workers > 1, dispatched across worker processes, each with
# its own --threads threads so this inner Threads.@threads still has something to use) or
# from a plain Threads.@threads loop over k (--workers 1, today's single-process
# behavior) -- identical body either way; only who calls it differs below. Defined via
# @everywhere unconditionally (harmless/equivalent to a plain `function` when there are no
# extra worker processes) so it works identically in both modes without a branch here.
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

    # Run `nseeds` independent attempts of the strategy (find_mle_advanced), keeping whichever reaches
    # the highest log_post (NaN treated as -Inf). Attempts are independent of each other (each runs the
    # strategy from its own RNG against the SAME read-only truth_likelihood/prior for this grid point) so they run
    # concurrently via Threads.@threads, one Julia thread per attempt -- this is the actual expensive work.
    #
    # Each attempt s gets its own Random.Xoshiro(1000*k + s) RNG rather than sharing one seed_rng advanced
    # sequentially across s (which would be inherently serial and would race under concurrent access) --
    # deterministic and reproducible per (k, s). Deterministic strategies (e.g. smseed lbfgsb) give the same
    # fit for every attempt, so nseeds > 1 only repeats it there.
    start_param = start_mode == "truth" ? truth_param : p
    attempt_results = Vector{Any}(undef, nseeds)
    Threads.@threads for s in 1:nseeds
        attempt_results[s] = find_mle_advanced(truth_likelihood, cp_hyp, start_param; strategy, nseeds = 1,
                                               rng = Random.Xoshiro(1000 * k + s), objective, null_prior = cp_sm, local_defaults)
    end

    # Reduce the attempts serially: keep the best by log_post, along with THAT attempt's own seed point (the point
    # its last stage started from; needed below for --chi2-mis-modeling's profile_start). If every attempt failed
    # (log_post stays NaN for all of them), fall back to a NaN-filled failure record, so this grid point is recorded
    # as a non-convergent failure instead of crashing on a `nothing` destructure below.
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
        # best_start_param's scan-key value(s) come from the free/unconstrained PSO seed
        # search and generally do NOT equal truth_param's -- but prior_fixed now holds
        # scan_keys constant at the truth value, and bat_findmode's ExplicitInit strictly
        # validates the init point against the prior's constants (same failure mode
        # documented on controlled_truth/_fix_conditioned in optimizer_tests_common.jl).
        # Overwrite the scan key(s) in the seed to match before fitting, or this throws
        # ArgumentError (silently caught by local_find_mle, returning NaN log_post).
        # Uses best_start_param (the WINNING attempt's seed, tracked above), not just
        # whichever attempt happened to run/finish last.
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

# --workers > 1: dispatch grid points across worker processes via pmap (dynamic
# load-balancing -- handles the significant per-point cost variance seen throughout
# testing better than Threads.@threads's static partition). Each worker still internally
# parallelizes that point's own n_attempts seed-attempts via the inner Threads.@threads
# inside run_grid_point, using the --threads threads it was spawned with.
# --workers 1 (default): today's exact single-process behavior, Threads.@threads over k,
# with results[k] assigned directly (no pmap/serialization involved).
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
# Single combined figure: NSI panel(s) first, then (if --fluctuate-nuisances) nuisance +
# physics panels, then a Δχ²_mod histogram. All panels share the same small-multiples
# layout (row, 2*col axis / 2*col+1 colorbar), one shared "converged"/"not converged"
# legend in column 1 spanning every row, Δχ²_mod colors the FIT points. Converged fits are drawn as
# circles/diamonds, non-converged fits as crosses and fits WITHOUT a convergence flag (`converged = missing`,
# the result of a global-only strategy) as squares, all in the same color scale; only fits whose Δχ²_mod is
# not finite (the fit itself failed) are drawn as grey crosses.

conv_flags = [r.converged === true for r in results]
noflag_flags = [ismissing(r.converged) for r in results]   # no convergence criterion (global-only strategy)
Δχ²_all = [r.Δχ²_mod for r in results]
# colorrange is computed from the FINITE Δχ²_mod values of all points, converged or not (a failed
# fit has Δχ²_mod = NaN when local_find_mle fails; those points are drawn grey). A NaN colorrange
# silently corrupts every Colorbar's internal image and crashes CairoMakie's PDF renderer
# (ArgumentError: start and stop must be finite), so only finite values enter it.
#
# Separately, a DEGENERATE (zero-width) range -- cmin == cmax, e.g. only 1 converged point,
# or every converged point sharing the exact same Δχ²_mod -- is finite but still crashes
# CairoMakie's colormap interpolation ("Can't interpolate in a range where cmin == cmax").
# Widen a degenerate range by a small symmetric epsilon so the colorbar always has a valid
# (nonzero) span; falls back to the same (-1.0, 1.0) default when there's nothing to plot.
#
# The degeneracy check itself must be done in Float32, not Float64: CairoMakie casts
# colorrange (and the plotted colors) down to Float32 internally before rendering, so two
# Float64 values that are merely CLOSE (not bit-equal) -- e.g. 0.0028239744427... vs.
# 0.0028239744568..., differing only ~1e-11 apart -- can still round to the exact same
# Float32 value and silently reproduce the degenerate-range crash even though the Float64
# `lo == hi` check above passed as "not degenerate" (confirmed: this exact pair of
# near-but-not-equal Float64 Δχ²_mod values from a real run crashed CairoMakie's PDF
# renderer with "Can't interpolate in a range where cmin == cmax" despite lo != hi here).
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

# Build the ordered list of NSI panel specs first (hypothesis-dependent), then append
# nuisance/physics panels (only if fluctuated or toydata), then the histogram panel.
# Each spec is a NamedTuple describing how to draw one panel; `kind` dispatches the draw
# logic below. `plane` panels (2D mag/phase) are handled specially since they overlay
# truth+fit as two point sets rather than a truth-vs-fit diagonal.
show_all_params = fluctuate_nuisances_flag || toydata_flag # show all parameter plots also for toydata runs
panel_specs = NamedTuple[]

if is_many
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
                             color=Δχ²_all[sel_conv], colormap=:viridis, colorrange=colorrange,
                             markersize=8, label="converged")
        any(sel_nonc) && scatter!(ax, truth_vals[sel_nonc], fit_vals[sel_nonc]; color=Δχ²_all[sel_nonc], colormap=:viridis,
                                  colorrange=colorrange, marker=:xcross, markersize=10, label="not converged")
        any(sel_noflag) && scatter!(ax, truth_vals[sel_noflag], fit_vals[sel_noflag]; color=Δχ²_all[sel_noflag], colormap=:viridis,
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
        local ax = Axis(fig[row, 2*col], xlabel="$(mag_key) (magnitude)", ylabel="$(phase_key) (phase)",
                         title="Truth vs fit (mag/phase plane)")
        local truth_mags   = [r.truth_param[mag_key] for r in results]
        local truth_phases = [r.truth_param[phase_key] for r in results]
        local fit_mags     = [r.fit_param[mag_key] for r in results]
        local fit_phases   = [r.fit_param[phase_key] for r in results]

        scatter!(ax, truth_mags, truth_phases; color=:steelblue, marker=:circle, markersize=10, label="truth")
        local plt_fit = scatter!(ax, fit_mags[sel_conv], fit_phases[sel_conv];
                            color=Δχ²_all[sel_conv], colormap=:viridis, colorrange=colorrange,
                            marker=:diamond, markersize=12, label="fit (converged)")
        any(sel_nonc) && scatter!(ax, fit_mags[sel_nonc], fit_phases[sel_nonc]; color=Δχ²_all[sel_nonc], colormap=:viridis,
                                  colorrange=colorrange, marker=:xcross, markersize=14, label="fit (not converged)")
        any(sel_noflag) && scatter!(ax, fit_mags[sel_noflag], fit_phases[sel_noflag]; color=Δχ²_all[sel_noflag], colormap=:viridis,
                                    colorrange=colorrange, marker=:rect, markersize=11, label="fit (no convergence flag)")
        local failed_plane = sel_fail .& isfinite.(fit_mags) .& isfinite.(fit_phases)
        if any(failed_plane)
            scatter!(ax, fit_mags[failed_plane], fit_phases[failed_plane]; color=:red, marker=:xcross, markersize=14)
            global fail_drawn = true
        end
        Colorbar(fig[row, 2*col + 1], plt_fit, label=Δχ²_label, ticklabelsize=9, labelsize=10)

    elseif spec.kind == :histogram
        local ax_hist = Axis(fig[row, 2*col], xlabel=Δχ²_label, ylabel="count",
                        title=rich(Δχ²_label, " distribution (converged or no flag)"))
        local Δχ²_converged = filter(isfinite, Δχ²_all[conv_flags .| noflag_flags])
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
