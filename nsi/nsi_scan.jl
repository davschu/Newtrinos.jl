# Profile (or plain likelihood scan) of one hypothesis of a setup file, fitted with a strategy at every grid point.
# Usage: julia nsi/nsi_scan.jl --setup <file> --name <n> --task profile --hypothesis <name> --objective posterior --strategy random lbfgsb [--workers N]
# Saves a Newtrinos.NewtrinosResult (key "result"). Options and output: nsi/README.md.

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
using Dates
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
    # --multi-node: SlurmManager() reads SLURM_NTASKS itself, so --workers is only a consistency check (a mismatch fails before spawning).
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

##### strategy #####

nseeds = args["nseeds"]
normalize_strategy(strategy)   # checks the stage names/settings/order before any expensive work
run_date = Dates.format(Dates.today(), "ddmmyyyy")   # execution date, ends every file name

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

# The result is a Newtrinos.NewtrinosResult (like Newtrinos.profile / scan). A grid (1 or 2 tested parameters) has one axis per tested
# parameter; for more tested parameters (Sobol sample, no common axes) the only axis is the point index and the tested parameters are
# ordinary entries of `values`.
sobol = !(length(scan_keys) in (1, 2))
result_axes = sobol ? (; point = collect(1:n_grid)) : NamedTuple{Tuple(scan_keys)}(init_values)
result_meta = Dict{String, Any}("hypothesis" => hypothesis_name, "setup" => S.name, "scan_keys" => Tuple(scan_keys), "n_grid" => n_grid,
                                "mesh_flat" => mesh_flat, "sobol" => sobol)
sobol || (result_meta["vars_to_scan"] = vars_to_scan)

if lowercase(args["task"]) == "profile"
    # Seed/alias first stages are computed ONCE (plan), a global first stage runs per point; every point has its own RNG (seed + index), so results do not
    # depend on scheduling. The pieces are used instead of profile_advanced so that `likelihood` stays a worker-global inside pmap (a closure would ship it).
    plan = _profile_plan(likelihood, cp_hyp, p, scanpoints, mesh_flat, scan_keys; strategy, nseeds, seed = args["seed"], objective, null_prior = cp_sm)
    t_start = time()
    if use_distributed
        @everywhere plan = $plan
        recs = @showprogress pmap(i -> _profile_point_fit(likelihood, plan, i), eachindex(scanpoints))
    else
        recs = Vector{Any}(undef, n_grid)
        @showprogress Threads.@threads for i in eachindex(scanpoints)
            recs[i] = _profile_point_fit(likelihood, plan, i)
        end
    end
    result = _profile_result(plan, recs, time() - t_start; axes = result_axes, meta = result_meta)
    FileIO.save(joinpath(results_dir, name * "_profile_$(hypothesis_name)_$(grid_points)points_$(strategy_tag)_nseeds$(nseeds)_$(run_date).jld2"), "result", result)

elseif lowercase(args["task"]) == "scan"
    # No optimization -- straight likelihood evaluation with every non-scanned parameter
    # frozen at its nominal value `p` (matches Newtrinos.scan's own semantics).
    if !sobol
        result = Newtrinos.scan(likelihood, cp_hyp, vars_to_scan, p)
        merge!(result.meta, result_meta)
    else
        t_start = time()
        llhs = Vector{Float64}(undef, n_grid)
        @showprogress Threads.@threads for i in 1:n_grid
            point = merge(p, NamedTuple{Tuple(scan_keys)}(Tuple(Float64.(mesh_flat[i]))))
            llhs[i] = logdensityof(likelihood, point)
        end
        vals = OrderedDict{Symbol, Any}(k => fill(p[k], n_grid) for k in setdiff(keys(p), scan_keys))
        for (j, k) in enumerate(scan_keys)
            vals[k] = Float64[Float64(mesh_flat[i][j]) for i in 1:n_grid]
        end
        vals[:llh] = llhs
        vals[:log_posterior] = llhs
        meta = Dict{String, Any}("task" => "scan", "priors" => cp_hyp, "params" => p, "exec_time" => time() - t_start)
        merge!(meta, result_meta)
        Newtrinos.add_meta!(meta)
        result = NewtrinosResult(axes = result_axes, values = NamedTuple(vals), meta = meta)
    end
    FileIO.save(joinpath(results_dir, name * "_scan_$(hypothesis_name)_$(grid_points)points_$(run_date)" * ".jld2"), "result", result)
else
    error("--task must be 'Profile' or 'Scan' (case-insensitive), got '$(args["task"])'")
end
