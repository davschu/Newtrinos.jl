using Distributed
using DataStructures
using FileIO
using Accessors
using ArgParse

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

        "--plot"
        help = "Enable plotting"
        action = :store_true

        "--workers"
        help = "Number of workers to use"
        arg_type = Int
        default = 4
    end

    return parse_args(s)
end

args = parse_command_line()

addprocs(args["workers"])
name = args["name"]

@everywhere args = $args

@everywhere begin
    using Distributions
    using DensityInterface
    using BAT
    using MeasureBase
    using ADTypes
    using Newtrinos

    include(joinpath(@__DIR__, "cli_common.jl"))

    ##### PHYSICS CONFIG #####
    # To use defaults:
    experiments = configure_experiments(args["experiments"])

    # To override physics, uncomment and modify:
    # osc = Newtrinos.osc.configure(Newtrinos.osc.OscillationConfig(
    #     flavour=Newtrinos.osc.ThreeFlavour(ordering=:IO),
    #     interaction=Newtrinos.osc.SI(),
    # ))
    # atm_flux = Newtrinos.atm_flux.configure()
    # earth_layers = Newtrinos.earth_layers.configure()
    # xsec = Newtrinos.xsec.configure()
    # physics = (; osc, atm_flux, earth_layers, xsec)
    # experiments = configure_experiments(args["experiments"], physics)

    likelihood = Newtrinos.generate_likelihood(experiments)
end

params = Newtrinos.get_params(experiments)
priors = Newtrinos.get_priors(experiments)
conditional_vars = Dict(:θ₁₂=>params.θ₁₂, :δCP=>-1.89, :Δm²₂₁=>params.Δm²₂₁)
priors = Newtrinos.condition(priors, conditional_vars, params)

@reset priors.Δm²₃₁ = Uniform(0.0018, 0.0028)
@reset priors.θ₂₃ = Uniform(pi/4-0.2, pi/4+0.2)
vars_to_scan = OrderedDict()
vars_to_scan[:θ₂₃] = 7
vars_to_scan[:Δm²₃₁] = 7
values, scanpoints = Newtrinos.generate_scanpoints(vars_to_scan, priors)

cache_dir = joinpath(@__DIR__, name)

if !isdir(cache_dir)
    mkdir(cache_dir)
end

t1 = time()
opt_results = pmap(x -> Newtrinos.find_mle_cached(likelihood, x, deepcopy(params), cache_dir), scanpoints)
res = Newtrinos.assemble_profile_results(opt_results, size(scanpoints))
t2 = time()

meta = Dict("task"=> "profile", "priors"=>priors, "vars_to_scan"=>vars_to_scan, "params"=>params, "exec_time"=>t2-t1, "cache_dir"=>cache_dir)
Newtrinos.add_meta!(meta)
axes = NamedTuple{tuple(keys(vars_to_scan)...)}(values)
result = NewtrinosResult(axes=axes, values=res, meta=meta)

FileIO.save(name * ".jld2", Dict("result" => result))

if args["plot"]
    using CairoMakie
    fig = Figure()
    ax = Axis(fig[1,1])
    plot!(ax, result)
    ax.xlabel = String(collect(keys(vars_to_scan))[1])

    if length(vars_to_scan) == 1
        ax.ylabel = "-2ΔLLH"
    else
        ax.ylabel = String(collect(keys(vars_to_scan))[2])
    end
    save(name * ".png", fig)
end
