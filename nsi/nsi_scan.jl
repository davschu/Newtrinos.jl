using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Distributions
using Distributed
using SlurmClusterManager
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
import Sobol
using ProgressMeter
import Optim

include(joinpath(@__DIR__, "..", "src", "analysis", "optimizer_tests_common.jl"))

function parse_command_line()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--setup"
        help = "Setup file defining the physics, experiments, priors and hypotheses (it defines setup(; ftype) -- see nsi/setups/nsi_deepcore_3y.jl for the format)"
        arg_type = String
        default = joinpath(@__DIR__, "setups", "nsi_deepcore_3y.jl")

        "--results-dir"
        help = "Output folder for the saved results (default: results/<setup name>/scan, relative to the working directory)"
        arg_type = String
        default = ""

        "--name"
        help = "Name for outputs"
        arg_type = String
        required = true

        "--task"
        help = "Task to perform: Choice of Profile, Scan"
        arg_type = String
        required = true

        "--hypothesis" # maybe change to profile-vars - would be more general - do when modifying for generic usage
        help = "Name of the hypothesis (from the setup file) to test"
        arg_type = String
        required = true

        "--objective"
        help = "What bat_findmode maximizes for the local fit: posterior (likelihood x prior, i.e. MAP -- default) or likelihood (likelihood only, via a flat prior over the same bounds -- i.e. MLE)"
        arg_type = String
        required = true

        "--grid-points"
        help = "Number of grid points for the scan/profile"
        arg_type = Int
        default = 25

        "--seed-strategy"
        help = "seed for the optimizer - available: random, mle, pso, smseed"
        arg_type = String
        default = "random"

        "--nseeds"
        help = "Number of independent local-optimizer starts+fits per gridpoint; keeps the best point"
        arg_type = Int
        default = 1

        "--seed"
        help = "RNG seed for seed-point generation"
        arg_type = Int
        default = 1234

        "--global-iterations"
        help = "Iteration budget for the PSO global seed-search stage (only used with --seed-strategy pso)"
        arg_type = Int
        default = 30

        "--pso-particles"
        help = "Swarm size (n_particles) for Optim.ParticleSwarm during the global seed-search stage (only used with --seed-strategy pso)"
        arg_type = Int
        default = 50

        "--iterations"
        help = "Max iterations for the local LBFGS optimizer"
        arg_type = Int
        default = 2000

        "--step-size"
        help = "Initial step size (alphaguess) for LBFGS"
        arg_type = Float64
        default = 1.0

        "--g-tol"
        help = "Gradient-norm convergence tolerance"
        arg_type = Float64
        default = 1e-6

        "--workers"
        help = "Number of distributed workers (default: 1, no distributed)"
        arg_type = Int
        default = 1

        "--threads"
        help = "Number of threads per worker (only used when --workers > 1)"
        arg_type = Int
        default = 1

        "--multi-node"
        help = "Spawn workers across multiple SLURM-allocated nodes (e.g. LRZ cm4_std) via SlurmClusterManager.jl instead of all on this single machine (plain addprocs(n_workers), the cm4_tiny-compatible default). Only valid when --workers > 1 AND running inside an active SLURM allocation (SLURM_JOB_ID/SLURM_JOBID and SLURM_NTASKS must be set). Worker count on multi-node is controlled by the SBATCH script's --ntasks, not by --workers directly -- --workers is only used here as a consistency check (must equal SLURM_NTASKS exactly, or this errors out before spawning anything)."
        action = :store_true
    end

    return parse_args(s)
end

args = parse_command_line()

name = args["name"]
n_workers = args["workers"]
n_threads = args["threads"]
multi_node_flag = args["multi-node"]
hypothesis_name = args["hypothesis"]
objective = args["objective"]
objective in ("posterior", "likelihood") || error("--objective must be 'posterior' or 'likelihood', got '$objective'")
global_iterations = args["global-iterations"]
pso_particles = args["pso-particles"]
lbfgs_iterations = args["iterations"]
step_size = args["step-size"]
g_tol = args["g-tol"]

use_distributed = n_workers > 1
multi_node_flag && !use_distributed &&
    error("--multi-node requires --workers > 1 (got --workers $n_workers).")
ftype = Float64

map_func = nothing

if use_distributed
    # --multi-node switches from plain addprocs(n_workers) (all workers on THIS single
    # machine -- correct for cm4_tiny) to SlurmClusterManager.jl's SlurmManager(), which
    # launches workers across every node SLURM allocated to this job (e.g. cm4_std's
    # multi-node jobs). SlurmManager() takes NO worker-count argument -- it reads
    # ENV["SLURM_NTASKS"] itself (requires an active SLURM allocation). Worker count is
    # therefore fixed by the SBATCH script's own --ntasks, NOT by --workers directly --
    # --workers is only used as an explicit consistency check here, so a mismatched
    # SBATCH script (--ntasks != --workers) fails loudly before spawning anything.
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

    @everywhere args = $args

    @everywhere begin
        using Distributions
        using DensityInterface
        using BAT
        using DataStructures
        using MeasureBase
        using ADTypes
        using Newtrinos
        using Accessors
        using ValueShapes
        import Sobol
        import Optim

        include(joinpath($(@__DIR__), "..", "src", "analysis", "optimizer_tests_common.jl"))

    end

    # every process (workers and driver) builds its own experiments and likelihood from the setup file;
    # separate @everywhere calls, so that the setup function is defined before it is called
    @everywhere include($(abspath(args["setup"])))
    @everywhere begin
        S = finalize_setup(setup(; ftype = Float64))
        experiments = S.experiments
        likelihood = Newtrinos.generate_likelihood(experiments)
    end

    map_func = pmap
else
    adsel = AutoForwardDiff()
    context = set_batcontext(ad = adsel)

    ##### CONFIGURATION (from the setup file) #####
    include(abspath(args["setup"]))
    S = finalize_setup(setup(; ftype = Float64))
    experiments = S.experiments
    likelihood = Newtrinos.generate_likelihood(experiments)
end

##### PARAMETERS, PRIORS AND HYPOTHESIS (from the setup file) #####

p = S.params
priors = S.priors
hyp = hypothesis_info(S, hypothesis_name)
cp_hyp = hyp.prior       # prior of the hypothesis under test
cp_sm = hyp.null_prior   # the same prior with the tested parameters at their null values (used by --seed-strategy smseed)
results_dir = isempty(args["results-dir"]) ? joinpath("results", S.name, "scan") : args["results-dir"]
mkpath(results_dir)

#ToDo: implement profile scan
#ToDo: implement saving of results to file


##### seed strategy #####

nseeds = args["nseeds"]
seed_strategy = args["seed-strategy"]
rng = Random.Xoshiro(args["seed"])
prior_dist = distprod(;cp_hyp...)
algorithm = make_algorithm("lbfgs", step_size)

seed_params = Vector{Any}(undef, nseeds)
if seed_strategy == "random" # Plain random draw from the (conditioned) prior -- fixed/conditioned entries just resample their own constant
    for s in 1:nseeds
        seed_params[s] = rand(rng, prior_dist)
    end
elseif seed_strategy == "pso"
    # Short bounded PSO global search per seed (global_seed_search, from
    # optimizer_tests_common.jl) -- unlike "random"/"mle" above, PSO seeds are NOT built
    # here: a seed built now, anchored at the nominal params `p`, would search for a good
    # point completely independent of which grid point it's later used for (the scanned
    # key is only overwritten AFTER the PSO search already ran, via `merge(seed_params[s],
    # fixed_here)` in the profile task below), so the PSO search itself never actually
    # explores near its target grid point. Instead, the profile task below runs a fresh
    # PSO search per grid point, anchored at that point's own fixed scan-key value(s) --
    # see its loop body for details.
    #
    # seed_params itself is never read for "pso" (see above), but it's still referenced --
    # in a dead branch -- inside the profile task's `pmap(...) do i ... end` closure (the
    # `else` arm of `start_param = if seed_strategy == "pso" ... else
    # merge(seed_params[s], fixed_here) end`). Julia captures closures based on what the
    # code mentions, not what actually runs, so under --workers>1
    # Distributed.ClusterSerializer still has to serialize/hash `seed_params` as part of
    # shipping that closure to workers -- and hashing a Vector{Any} with undefined
    # (never-assigned) slots throws UndefRefError. Filling every slot with a throwaway
    # value (never read) avoids that without changing "pso"'s actual behavior.
    fill!(seed_params, NamedTuple())
elseif seed_strategy == "mle" # Random start + a full local MLE fit per seed -- the "seed" is the fit result itself.
    # Uses local_find_mle (not Newtrinos.find_mle, which is always MAP/posterior and has no
    # objective knob) so this respects --objective the same way the "pso" branch already
    # does -- global_seed_search's own docstring warns that mixing a posterior-biased seed
    # with a likelihood-only local polish (or vice versa) is exactly the wrong thing to do.
    for s in 1:nseeds
        start_param = rand(rng, prior_dist)
        res = local_find_mle(likelihood, prior_dist, start_param;
                              fit_method="optim", algorithm=algorithm, iterations=lbfgs_iterations,
                              g_tol=g_tol, f_tol=0.0, x_tol=0.0, objective=objective)
        seed_params[s] = res[3]
    end
elseif seed_strategy == "smseed"
    # Null-model seeding: fit the nuisance/physics parameters with the tested parameters of this
    # hypothesis pinned at their null values (e.g. the SM point), then hand that fit's result --
    # nuisances tuned, tested parameters still sitting at their null values -- on as the seed. The idea
    # is that the systematics land at a reasonable point from the null fit alone, so the subsequent fit
    # mostly has to tune the tested parameters themselves rather than re-discovering good nuisance
    # values from scratch. The SM fit doesn't depend on `s`, so it only needs to run once; it's
    # reused for every one of the `nseeds` seeds (matching "mle"/"random"'s per-seed slot
    # convention, just with identical contents across seeds here).
    # `ExplicitInit` strictly validates the init point against every constant entry of
    # `cp_sm`'s prior (its conditioned parameters + the tested parameters pinned at their null values)
    # -- `p` itself (Newtrinos.get_params's raw nominal values) doesn't necessarily agree
    # with those constants, so every key `cp_sm` holds fixed must be forced to its exact
    # conditioned value first, or `bat_findmode` throws `ArgumentError("Cannot set
    # constant value to a different value")` (same issue `_fix_conditioned` solves in
    # optimizer_tests_common.jl).
    sm_overrides = Dict{Symbol, Float64}()
    for k in keys(cp_sm)
        d = cp_sm[k]
        if d isa ValueShapes.ConstValueDist
            sm_overrides[k] = Float64(d.value)
        elseif d isa Number
            sm_overrides[k] = Float64(d)
        end
    end
    sm_init = merge(p, NamedTuple(sm_overrides))
    sm_res = local_find_mle(likelihood, distprod(;cp_sm...), sm_init;
                          fit_method="optim", algorithm=algorithm, iterations=lbfgs_iterations,
                          g_tol=g_tol, f_tol=0.0, x_tol=0.0, objective=objective)
    fill!(seed_params, sm_res[3])
else
    error("Unknown --seed-strategy '$seed_strategy'. Available options: random, mle, pso, smseed")
end


##### specify gridpoints for scan/profile #####

# grid over #grid-points total points in the space of the tested parameters; for more than 2, use a sobol sequence to sample the space quasi-randomly.
scan_keys = collect(hyp.tested)
grid_points = args["grid-points"]

if length(scan_keys) in (1, 2)
    n_per_dim = length(scan_keys) == 1 ? grid_points : round(Int, sqrt(grid_points)) #take sqrt for 2D scan to keep total points ~ grid_points
    vars_to_scan = OrderedDict{Symbol,Int}(k => n_per_dim for k in scan_keys)
    init_values, scanpoints = Newtrinos.generate_scanpoints(vars_to_scan, cp_hyp)
    mesh_flat = vec(Newtrinos._generate_grid(vars_to_scan, cp_hyp)[3])  # grid points in the order of `scanpoints`
else
    # more than 2 tested parameters -- grid-points quasi-random samples via a Sobol low-discrepancy
    # sequence over all tested dimensions jointly, built the same way generate_scanpoints's own
    # make_prior closure builds each scan point (fix the scanned keys, keep everything else free).
    los = Float64[Float64(minimum(cp_hyp[k])) for k in scan_keys]
    his = Float64[Float64(maximum(cp_hyp[k])) for k in scan_keys]
    seq = Sobol.SobolSeq(los, his)
    mesh_flat = [Sobol.next!(seq) for _ in 1:grid_points]
    scanpoints = map(mesh_flat) do v
        p = deepcopy(cp_hyp)
        for (i, k) in enumerate(scan_keys)
            @reset p[k] = v[i]
        end
        distprod(;p...)
    end
    # NOTE: unlike the factorial branch, these are NOT a shared per-axis grid -- values[d][i] and
    # values[d'][i] are paired (same Sobol sample i), but values[d] alone isn't a common axis that
    # combines factorially with other axes. Don't reuse profile()'s `axes=NamedTuple(values)` /
    # reshape-to-grid pattern for this branch.
    init_values = [[pt[i] for pt in mesh_flat] for i in eachindex(scan_keys)]
end

n_grid = length(scanpoints)
println("Grid built for hypothesis '$hypothesis_name': $n_grid points, scanning $(Tuple(scan_keys))")


##### run scan/profile #####

if lowercase(args["task"]) == "profile"
    # Run the local optimizer at every grid point, trying each of the `nseeds` candidate
    # seeds and keeping whichever converges to the highest log_posterior. Uses
    # local_find_mle (optimizer_tests_common.jl, fit_method="optim") directly instead of
    # Newtrinos.profile/multistart_profile: their own nseeds mechanism ignores
    # --seed-strategy, and they're built on find_mle_ext, which never reports an
    # iteration count. Same per-point body inlined in both branches below (no helper
    # function) -- bat_findmode's ExplicitInit strictly validates the init point against
    # the prior's constants, so the seed's scanned key(s) must be overwritten to match
    # this grid point's fixed value first, or it throws ArgumentError.
    # `algorithm` is already built above in ##### seed strategy #####.
    #
    # Under seed_strategy=="pso", each of the `nseeds` seeds is a FRESH PSO search
    # anchored at THIS grid point's own fixed scan-key value(s) (merge(p, fixed_here)),
    # run right here rather than read from the pre-built `seed_params` (see the ##### seed
    # strategy ##### block above) -- this is what makes each PSO search actually explore
    # near its target grid point instead of the unconditioned nominal point. Each grid
    # point gets its own RNG stream (seed_rng, seeded from args["seed"] + this point's
    # linear index) rather than reusing the single top-level `rng`: under
    # Threads.@threads, mutating one shared RNG object from multiple threads concurrently
    # is unsafe, and per-point streams also keep results reproducible independent of
    # pmap/thread scheduling order. global_seed_search consumes seed_rng progressively
    # across the `for s in 1:nseeds` loop, so within one grid point each of the `nseeds`
    # PSO sub-runs still draws a different random start, matching the "several
    # independent PSO runs, keep the best" pattern used elsewhere (e.g. roundtrips.jl).
    # Under seed_strategy=="random"/"mle", seeds are unaffected and come from the
    # pre-built seed_params exactly as before.

    if use_distributed
        opt_results = @showprogress pmap(eachindex(scanpoints)) do i
            fixed_here = NamedTuple{Tuple(scan_keys)}(Tuple(Float64.(mesh_flat[i])))
            seed_rng = Random.Xoshiro(args["seed"] + i)
            best, best_log_post = nothing, -Inf
            t0 = time()
            for s in 1:nseeds
                start_param = if seed_strategy == "pso"
                    global_seed_search(likelihood, scanpoints[i], merge(p, fixed_here), "pso", global_iterations, seed_rng; n_particles=pso_particles, objective=objective)
                else
                    merge(seed_params[s], fixed_here)
                end
                res = local_find_mle(likelihood, scanpoints[i], start_param;
                                      fit_method="optim", algorithm=algorithm, iterations=lbfgs_iterations,
                                      g_tol=g_tol, f_tol=0.0, x_tol=0.0, objective=objective)
                log_post = isnan(res[2]) ? -Inf : res[2]
                if best === nothing || log_post > best_log_post
                    best, best_log_post = res, log_post
                end
            end
            (best..., time() - t0)
        end
    else
        opt_results = Vector{Any}(undef, length(scanpoints))
        @showprogress Threads.@threads for i in eachindex(scanpoints)
            fixed_here = NamedTuple{Tuple(scan_keys)}(Tuple(Float64.(mesh_flat[i])))
            seed_rng = Random.Xoshiro(args["seed"] + i)
            best, best_log_post = nothing, -Inf
            t0 = time()
            for s in 1:nseeds
                start_param = if seed_strategy == "pso"
                    global_seed_search(likelihood, scanpoints[i], merge(p, fixed_here), "pso", global_iterations, seed_rng; n_particles=pso_particles, objective=objective)
                else
                    merge(seed_params[s], fixed_here)
                end
                res = local_find_mle(likelihood, scanpoints[i], start_param;
                                      fit_method="optim", algorithm=algorithm, iterations=lbfgs_iterations,
                                      g_tol=g_tol, f_tol=0.0, x_tol=0.0, objective=objective)
                log_post = isnan(res[2]) ? -Inf : res[2]
                if best === nothing || log_post > best_log_post
                    best, best_log_post = res, log_post
                end
            end
            opt_results[i] = (best..., time() - t0)
        end
    end

    result = Dict(
        "hypothesis"    => hypothesis_name,
        "scan_keys"     => Tuple(scan_keys),
        "mesh_flat"     => mesh_flat,
        "llh"           => [r[1] for r in opt_results],
        "log_posterior" => [r[2] for r in opt_results],
        "fit_params"    => [r[3] for r in opt_results],
        "converged"     => [r[4] for r in opt_results],
        "n_iterations"  => [r[5] for r in opt_results],
        "lbfgs_time"    => [r[6] for r in opt_results],
        "n_grid"        => n_grid,
        "nseeds"        => nseeds,
        "seed_strategy" => seed_strategy,
        "objective"     => objective,
        "setup"         => S.name,
    )
    FileIO.save(joinpath(results_dir, name * "_profile_$(hypothesis_name)_$(grid_points)points_$(seed_strategy)seed$(nseeds)" * ".jld2"), result) 

elseif lowercase(args["task"]) == "scan"
    # No optimization -- straight likelihood evaluation with every non-scanned parameter
    # frozen at its nominal value `p` (matches Newtrinos.scan's own semantics).
    if length(scan_keys) in (1, 2)
        result = Newtrinos.scan(likelihood, cp_hyp, vars_to_scan, p)
    else
        llhs = Vector{Float64}(undef, n_grid)
        @showprogress Threads.@threads for i in 1:n_grid
            point = merge(p, NamedTuple{Tuple(scan_keys)}(Tuple(Float64.(mesh_flat[i]))))
            llhs[i] = logdensityof(likelihood, point)
        end
        result = Dict("hypothesis"=>hypothesis_name, "scan_keys"=>Tuple(scan_keys), "mesh_flat"=>mesh_flat, "llh"=>llhs, "n_grid"=>n_grid, "setup"=>S.name)
    end
    FileIO.save(joinpath(results_dir, name * "_scan_$(hypothesis_name)_$(grid_points)points" * ".jld2"), result)
else
    error("--task must be 'Profile' or 'Scan' (case-insensitive), got '$(args["task"])'")
end