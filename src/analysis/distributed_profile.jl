using Distributed
using DataStructures
using FileIO
using Accessors
using ArgParse

function parse_command_line()
    s = ArgParseSettings()

    @add_arg_table s begin
        # "--experiments"
        # help = "List of experiments to run"
        # nargs = '+'
        # required = true

        # "--ordering"
        # help = "NMO: either NO or IO"
        # arg_type = String
        # default = "NO"      

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

@everywhere begin
    using LinearAlgebra
    using Distributions
    using DensityInterface
    using Base
    using ForwardDiff
    using BAT
    using IterTools
    using MeasureBase
    using ADTypes
    using Newtrinos

    osc_cfg = Newtrinos.osc.OscillationConfig(
        flavour=Newtrinos.osc.ThreeFlavour(ordering=Symbol(:NO)),
        propagation=Newtrinos.osc.Basic(),
        states=Newtrinos.osc.All(),
        interaction=Newtrinos.osc.SI()
        )
    osc = Newtrinos.osc.configure(osc_cfg)

    atm_flux = Newtrinos.atm_flux.configure(
            Newtrinos.atm_flux.AtmFluxConfig(nominal_model=Newtrinos.atm_flux.HKKM("kam-ally-20-01-mtn-solmin.d")
            )
        )
    earth_layers = Newtrinos.earth_layers.configure()
    xsec=Newtrinos.xsec.configure(Newtrinos.xsec.Differential_H2O())

    physics = (; osc, atm_flux, earth_layers, xsec);

    # dynamically construct named tuple from experiment names
    function configure_experiments(experiment_list, physics)
        pairs = (Symbol(lowercase(exp)) => getproperty(getproperty(Newtrinos, Symbol(lowercase(exp))), :configure)(physics) for exp in experiment_list)
        return (; pairs...)
    end

    experiments = configure_experiments(["super_k",], physics)
    likelihood = Newtrinos.generate_likelihood(experiments)
end

    
params = Newtrinos.get_params(experiments)
priors = Newtrinos.get_priors(experiments)
conditional_vars = Dict(:θ₁₂=>params.θ₁₂, :δCP=>-1.89, :Δm²₂₁=>params.Δm²₂₁)
priors = Newtrinos.condition(priors, conditional_vars, params)

@reset priors.Δm²₃₁ = Uniform(0.0018, 0.0028)
@reset priors.θ₂₃ = Uniform(pi/4-0.2, pi/4+0.2)
vars_to_scan = OrderedDict()
vars_to_scan[:θ₂₃] = 11
vars_to_scan[:Δm²₃₁] = 11
values, scanpoints = Newtrinos.generate_scanpoints(vars_to_scan, priors)

cache_dir = joinpath(@__DIR__, name)

if isnothing(cache_dir)
    mkdir(cache_dir)
end

t1 = time()
results = Array{Any}(undef, size(scanpoints))
llhs = Array{Any}(undef, size(scanpoints))
log_posteriors = Array{Any}(undef, size(scanpoints))

opt_results = pmap(x -> Newtrinos.find_mle_cached(likelihood, x, deepcopy(params), cache_dir), scanpoints)

for (i, opt_result) in enumerate(opt_results)
    llhs[i] = opt_result[1]
    log_posteriors[i] = opt_result[2]
    results[i] = opt_result[3]
end

s = OrderedDict(key=>[x[key] for x in results] for key in keys(first(results)))
s[:llh] = llhs
s[:log_posterior] = log_posteriors
res = NamedTuple(s)

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