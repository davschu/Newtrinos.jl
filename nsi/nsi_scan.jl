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

        "--strategy"
        help = "How the fit at every grid point is done: a chain of stages, each starting from the result of the previous one, e.g. --strategy random lbfgsb  |  --strategy smseed lbfgsb  |  --strategy mle lbfgsb  |  --strategy optim-pso:particles=50,iters=30 lbfgsb:iterations=2000  |  --strategy manual:ε_eμ_abs=0.1 lbfgsb. Stage kinds: seed (random, manual), alias (smseed = local fit with the tested parameters pinned at their null values; mle = random + local fit), global (optim-pso, optim-sa, mh-pso, mh-de, mh-abc, bbo, bbo-xnes, bbo-dxnes; stop at `evals` objective evaluations) and local (lbfgsb, optim-lbfgs, ...). Seed, alias and global stages are only allowed as the FIRST stage. A seed or alias first stage does not depend on the grid point and is computed ONCE per scan (the seed's scanned parameters are overwritten with each grid point's values); a global first stage is run at every grid point. Settings of a stage: name:key=val,key=val (see optimizer_tests.jl --list-optimizers)."
        arg_type = String
        nargs = '+'
        required = true

        "--nseeds"
        help = "Number of independent attempts of the strategy per gridpoint (for a seed first stage: number of different seeds); keeps the best point"
        arg_type = Int
        default = 1

        "--seed"
        help = "RNG seed for seed-point generation"
        arg_type = Int
        default = 1234

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
strategy = args["strategy"]
strategy_tag = join((first(split(s, ':')) for s in strategy), "+")   # stage names, for the output file name

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
cp_sm = hyp.null_prior   # the same prior with the tested parameters at their null values (used by the smseed stage of --strategy)
results_dir = isempty(args["results-dir"]) ? joinpath("results", S.name, "scan") : args["results-dir"]
mkpath(results_dir)

#ToDo: implement profile scan
#ToDo: implement saving of results to file


##### strategy and seeds #####

nseeds = args["nseeds"]
rng = Random.Xoshiro(args["seed"])
normalize_strategy(strategy)   # checks the stage names/settings/order before any expensive work

# A seed or alias FIRST stage (random, manual, smseed, mle) does not depend on the grid point: it is computed once here
# (nseeds seeds; smseed is deterministic, so one fit is shared) and every grid point starts its remaining stages from
# the seed with the scanned parameters overwritten by the point's own values (bat_findmode's ExplicitInit strictly
# validates the start against the prior's constants). A GLOBAL first stage (optim-pso, mh-de, ...) is run at every
# grid point instead, anchored at that point's fixed scan values -- a seed computed now would search independently of
# the point it is used for. With nseeds attempts per point, the best fit (by log_post) is kept.
first_kind = OPTIMIZERS[parse_optimizer_spec(strategy[1]).name].kind
hoisted = first_kind in (:seed, :alias)
rest_strategy = strategy[2:end]
hoisted && isempty(rest_strategy) && error("--strategy: '$(strategy[1])' only creates a seed; add a fitting stage after it, e.g. --strategy $(strategy[1]) lbfgsb")

# seed_params is also referenced (in a dead branch) inside the pmap closure below, so it must have every slot defined
# for Distributed's serialization (hashing a Vector{Any} with undefined slots throws UndefRefError) -- hence the fill.
seed_params = Vector{Any}(undef, nseeds)
fill!(seed_params, NamedTuple())
if hoisted
    for s in 1:nseeds
        seed_params[s] = if parse_optimizer_spec(strategy[1]).name == "smseed" && s > 1
            seed_params[1]   # the null-model fit does not depend on s
        else
            find_mle_advanced(likelihood, cp_hyp, p; strategy = [strategy[1]], rng, objective, null_prior = cp_sm).params
        end
    end
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
    # Run the strategy at every grid point (find_mle_advanced), keeping the best of `nseeds` attempts by log_posterior.
    # Uses the optimizer registry (optimizer_tests_common.jl) directly instead of Newtrinos.profile/multistart_profile:
    # their own nseeds mechanism has no strategy and they are built on find_mle_ext, which never reports an iteration
    # count. For a global first stage (see above) each point gets its own RNG stream (seed_rng, seeded from args["seed"]
    # + this point's linear index) rather than reusing the single top-level `rng`: under Threads.@threads, mutating one
    # shared RNG object from multiple threads concurrently is unsafe, and per-point streams also keep results
    # reproducible independent of pmap/thread scheduling order.

    if use_distributed
        opt_results = @showprogress pmap(eachindex(scanpoints)) do i
            fixed_here = NamedTuple{Tuple(scan_keys)}(Tuple(Float64.(mesh_flat[i])))
            seed_rng = Random.Xoshiro(args["seed"] + i)
            t0 = time()
            rec = if hoisted
                find_mle_advanced(likelihood, scanpoints[i], p; strategy = rest_strategy, starts = [merge(seed_params[s], fixed_here) for s in 1:nseeds],
                                  rng = seed_rng, objective = objective)
            else
                find_mle_advanced(likelihood, scanpoints[i], merge(p, fixed_here); strategy = strategy, nseeds = nseeds,
                                  rng = seed_rng, objective = objective)
            end
            best = (rec.llh, rec.log_post, rec.params, rec.converged, rec.n_iters)
            (best..., time() - t0)
        end
    else
        opt_results = Vector{Any}(undef, length(scanpoints))
        @showprogress Threads.@threads for i in eachindex(scanpoints)
            fixed_here = NamedTuple{Tuple(scan_keys)}(Tuple(Float64.(mesh_flat[i])))
            seed_rng = Random.Xoshiro(args["seed"] + i)
            t0 = time()
            rec = if hoisted
                find_mle_advanced(likelihood, scanpoints[i], p; strategy = rest_strategy, starts = [merge(seed_params[s], fixed_here) for s in 1:nseeds],
                                  rng = seed_rng, objective = objective)
            else
                find_mle_advanced(likelihood, scanpoints[i], merge(p, fixed_here); strategy = strategy, nseeds = nseeds,
                                  rng = seed_rng, objective = objective)
            end
            best = (rec.llh, rec.log_post, rec.params, rec.converged, rec.n_iters)
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
        "strategy"      => strategy,
        "objective"     => objective,
        "setup"         => S.name,
    )
    FileIO.save(joinpath(results_dir, name * "_profile_$(hypothesis_name)_$(grid_points)points_$(strategy_tag)_nseeds$(nseeds)" * ".jld2"), result) 

elseif lowercase(args["task"]) == "scan"
    # No optimization -- straight likelihood evaluation with every non-scanned parameter
    # frozen at its nominal value `p` (matches Newtrinos.scan's own semantics).
    if length(scan_keys) in (1, 2)
        r = Newtrinos.scan(likelihood, cp_hyp, vars_to_scan, p)
        # FileIO/JLD2 cannot save the NewtrinosResult struct directly -> store its fields in a Dict
        result = Dict("hypothesis"=>hypothesis_name, "scan_keys"=>Tuple(scan_keys), "axes"=>r.axes,
                      "values"=>r.values, "meta"=>r.meta, "mesh_flat"=>mesh_flat, "n_grid"=>n_grid, "setup"=>S.name)
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