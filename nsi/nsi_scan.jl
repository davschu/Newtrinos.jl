using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Distributions
using Distributed
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

include(joinpath(@__DIR__, "..", "src", "analysis", "cli_common.jl"))
include(joinpath(@__DIR__, "..", "src", "analysis", "optimizer_tests_common.jl"))

function parse_command_line()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--experiments"
        help = "List of experiments to run"
        nargs = '+'
        required = true

        "--name"
        help = "Name for outputs"
        arg_type = String
        required = true

        "--task"
        help = "Task to perform: Choice of Profile, Scan"
        arg_type = String
        required = true

        "--hypothesis" # maybe change to profile-vars - would be more general - do when modifying for generic usage
        help = "Hypothesis to test: eμ, eτ, μτ, ee_μμ, ττ_μμ, all"
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
        help = "seed for the optimizer - available: random, mle, pso"
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

        "--workers"
        help = "Number of distributed workers (default: 1, no distributed)"
        arg_type = Int
        default = 1

        "--threads"
        help = "Number of threads per worker (only used when --workers > 1)"
        arg_type = Int
        default = 1
    end

    return parse_args(s)
end

args = parse_command_line()

name = args["name"]
n_workers = args["workers"]
n_threads = args["threads"]
hypothesis_name = args["hypothesis"]
objective = args["objective"]
objective in ("posterior", "likelihood") || error("--objective must be 'posterior' or 'likelihood', got '$objective'")

use_distributed = n_workers > 1
ftype = Float64
mkpath("scan_results")

map_func = nothing

if use_distributed
    addprocs(n_workers; exeflags="--threads=$n_threads")

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

        include(joinpath($(@__DIR__), "..", "src", "analysis", "cli_common.jl"))
        include(joinpath($(@__DIR__), "..", "src", "analysis", "optimizer_tests_common.jl"))

        ##### CONFIGURATION (built independently on every worker) #####
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
        experiments  = Newtrinos.configure_experiments(args["experiments"], physics)

        likelihood = Newtrinos.generate_likelihood(experiments)
    end

    map_func = pmap
else
    adsel = AutoForwardDiff()
    context = set_batcontext(ad = adsel)

    ##### CONFIGURATION #####
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
    experiments  = Newtrinos.configure_experiments(args["experiments"], physics)

    likelihood = Newtrinos.generate_likelihood(experiments)
end

##### SPECIFY PARAMETER SPACE, PRIORS #####

p = Newtrinos.get_params(experiments)
priors = Newtrinos.get_priors(experiments)

#specifying priors

# Standard oscillation priors fixed to paper values
conditional_vars = Dict(:δCP=>0.0, :θ₁₂=>ftype(33.62 * π/180), :θ₁₃=>ftype(8.54  * π/180), :Δm²₂₁=>ftype(7.40e-5))
priors = Newtrinos.condition(priors, conditional_vars, p)
@reset priors.θ₂₃   = Uniform(ftype(30*π/180), ftype(60*π/180))
@reset priors.Δm²₃₁ = Uniform(ftype(0.93e-3 + 7.40e-5), ftype(3.93e-3 + 7.40e-5))#paper sets range for deltam^2_23 so just use flat prior from sum in NO: Δm²₃₁ = Δm²₂₃ + Δm²₂₁ with Δm²₂₁ fixed

# Detector-systematic prior ranges corrected to match paper Table III 
@reset priors.deepcore_lifetime        = Uniform(ftype(0), ftype(3.8)) #lifetime: paper allows [0, 3.8] years
@reset priors.deepcore_ice_scattering  = Truncated(Normal(ftype(1), ftype(0.1)), ftype(0.9), ftype(1.1)) # ice scattering/absorption: paper truncates at ±1σ (Normal(1,0.1) -> [0.9,1.1])
@reset priors.deepcore_ice_absorption  = Truncated(Normal(ftype(1), ftype(0.1)), ftype(0.9), ftype(1.1))
@reset priors.deepcore_opt_eff_lateral = Truncated(Normal(ftype(0), ftype(1.)), ftype(-2), ftype(2.5)) #paper's truncation is asymmetric, +2.5σ/-2σ

# NSI parameter bounds — Standard parameterization
@reset priors.Δ_eμ     = Uniform(-ftype(5),   ftype(5))
@reset priors.Δ_τμ     = Uniform(-ftype(0.1), ftype(0.1))
@reset priors.ε_eμ_abs = Uniform(ftype(0),    ftype(0.3))
@reset priors.ε_eτ_abs = Uniform(ftype(0),    ftype(0.35))
@reset priors.ε_μτ_abs = Uniform(ftype(0),    ftype(0.07))
@reset priors.δ_eμ     = Uniform(ftype(0),    ftype(2π)) #used [0,π] for 0,180 degree with two gridpoints, but use [0,2π] for general grid scan
@reset priors.δ_eτ     = Uniform(ftype(0),    ftype(2π))
@reset priors.δ_μτ     = Uniform(ftype(0),    ftype(2π))

##### SPECIFY HYPOTHESIS AND GRIDPOINTS #####

cp_eμ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eτ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eτ => 0.0, :δ_μτ => 0.0), p)
cp_eτ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eμ => 0.0, :δ_μτ => 0.0), p)
cp_μτ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0), p)
cp_ee_μμ = Newtrinos.condition(priors, Dict(:Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0, :δ_μτ => 0.0), p)
cp_ττ_μμ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0, :δ_μτ => 0.0), p)
cp_all = deepcopy(priors)

all_cps = Dict("eμ" => cp_eμ, "eτ" => cp_eτ, "μτ" => cp_μτ, "ee_μμ" => cp_ee_μμ, "ττ_μμ" => cp_ττ_μμ, "all" => cp_all)
hypothesis_name in keys(all_cps) || error("--hypothesis must be one of $(join(sort(collect(keys(all_cps))), ", ")), got '$hypothesis_name'")
cp_hyp = all_cps[hypothesis_name]
is_all = hypothesis_name == "all" # bool for special case of scanning all nsi params

#ToDo: implement profile scan
#ToDo: implement saving of results to file


##### seed strategy #####

nseeds = args["nseeds"]
seed_strategy = args["seed-strategy"]
rng = Random.Xoshiro(args["seed"])
prior_dist = distprod(;cp_hyp...)
algorithm = make_algorithm("lbfgs", 1.0)

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
elseif seed_strategy == "mle" # Random start + a full local MLE fit per seed -- the "seed" is the fit result itself.
    # Uses local_find_mle (not Newtrinos.find_mle, which is always MAP/posterior and has no
    # objective knob) so this respects --objective the same way the "pso" branch already
    # does -- global_seed_search's own docstring warns that mixing a posterior-biased seed
    # with a likelihood-only local polish (or vice versa) is exactly the wrong thing to do.
    for s in 1:nseeds
        start_param = rand(rng, prior_dist)
        res = local_find_mle(likelihood, prior_dist, start_param;
                              fit_method="optim", algorithm=algorithm, iterations=2000,
                              g_tol=1e-6, f_tol=0.0, x_tol=0.0, objective=objective)
        seed_params[s] = res[3]
    end
else
    error("Unknown --seed-strategy '$seed_strategy'. Available options: random, mle, pso")
end


##### specify gridpoints for scan/profile #####

# grid over #grid-points total points in the free NSI parameter space, for >2 free nsi params, use a sobol sequence to sample the space quasi-randomly.
nsi_param_names = (:Δ_eμ, :Δ_τμ, :ε_eμ_abs, :ε_eτ_abs, :ε_μτ_abs, :δ_eμ, :δ_eτ, :δ_μτ)
scan_keys = [k for k in nsi_param_names if !(cp_hyp[k] isa ValueShapes.ConstValueDist) && !(cp_hyp[k] isa Number)]
grid_points = args["grid-points"]

if length(scan_keys) in (1, 2)
    n_per_dim = length(scan_keys) == 1 ? grid_points : round(Int, sqrt(grid_points)) #take sqrt for 2D scan to keep total points ~ grid_points
    vars_to_scan = OrderedDict{Symbol,Int}(k => n_per_dim for k in scan_keys)
    init_values, scanpoints, mesh_flat = Newtrinos.generate_scanpoints(vars_to_scan, cp_hyp)
else
    # "all": 8 free NSI params -- grid-points quasi-random samples via a Sobol low-discrepancy
    # sequence over all 8 dimensions jointly, built the same way generate_scanpoints's own
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
                    global_seed_search(likelihood, scanpoints[i], merge(p, fixed_here), "pso", 30, seed_rng; n_particles=50, objective=objective)
                else
                    merge(seed_params[s], fixed_here)
                end
                res = local_find_mle(likelihood, scanpoints[i], start_param;
                                      fit_method="optim", algorithm=algorithm, iterations=2000,
                                      g_tol=1e-6, f_tol=0.0, x_tol=0.0, objective=objective)
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
                    global_seed_search(likelihood, scanpoints[i], merge(p, fixed_here), "pso", 30, seed_rng; n_particles=50, objective=objective)
                else
                    merge(seed_params[s], fixed_here)
                end
                res = local_find_mle(likelihood, scanpoints[i], start_param;
                                      fit_method="optim", algorithm=algorithm, iterations=2000,
                                      g_tol=1e-6, f_tol=0.0, x_tol=0.0, objective=objective)
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
    )
    FileIO.save(joinpath("scan_results/", name * "_profile_$(hypothesis_name)_$(grid_points)points_$(seed_strategy)seed$(nseeds)" * ".jld2"), result) 

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
        result = Dict("hypothesis"=>hypothesis_name, "scan_keys"=>Tuple(scan_keys), "mesh_flat"=>mesh_flat, "llh"=>llhs, "n_grid"=>n_grid)
    end
    FileIO.save(joinpath("scan_results/", name * "_scan_$(hypothesis_name)_$(grid_points)points" * ".jld2"), result)
else
    error("--task must be 'Profile' or 'Scan' (case-insensitive), got '$(args["task"])'")
end