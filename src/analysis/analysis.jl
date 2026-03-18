using LinearAlgebra
using Distributions
using DensityInterface
using ForwardDiff
using BAT
using IterTools
using DataStructures
using MeasureBase
using ADTypes
using Newtrinos
using FileIO
using Accessors
using ArgParse

include("cli_common.jl")

function parse_command_line()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--experiments"
        help = "List of experiments to run"
        nargs = '+'
        required = true

        "--ordering"
        help = "NMO: either NO or IO"
        arg_type = String
        default = "NO"

        "--name"
        help = "Name for outputs"
        arg_type = String
        required = true

        "--task"
        help = "Task to perform: Choice of NestedSampling, ImportanceSampling, Profile, Scan"
        arg_type = String
        required = true

        "--plot"
        help = "Enable plotting"
        action = :store_true
    end

    return parse_args(s)
end

args = parse_command_line()

adsel = AutoForwardDiff()
context = set_batcontext(ad = adsel)

name = args["name"]

physics = configure_physics(args["ordering"])
experiments = configure_experiments(args["experiments"], physics)
p = Newtrinos.get_params(experiments)
priors = Newtrinos.get_priors(experiments)

# Variables to condition on (=fix)
conditional_vars = Dict(:θ₁₂=>p.θ₁₂, :δCP=>-1.89, :Δm²₂₁=>p.Δm²₂₁)

# For profile / scan task only: choose scan grid
vars_to_scan = OrderedDict()
vars_to_scan[:θ₂₃] = 11
vars_to_scan[:Δm²₃₁] = 11

###### END CONFIG ######

likelihood = Newtrinos.generate_likelihood(experiments);

priors = Newtrinos.condition(priors, conditional_vars, p)

@reset priors.Δm²₃₁ = Uniform(0.0018, 0.0028)
@reset priors.θ₂₃ = Uniform(pi/4-0.2, pi/4+0.2)

if lowercase(args["task"]) == "nestedsampling"
    import UltraNest
    prior = distprod(;priors...)
    posterior = PosteriorMeasure(likelihood, prior)
    samples = bat_sample(posterior, ReactiveNestedSampling()).result
    FileIO.save(name * ".jld2", Dict("samples" => samples))

elseif lowercase(args["task"]) == "importancesampling"
    prior = distprod(;priors...)
    posterior = PosteriorMeasure(likelihood, prior)

    seed_points = load("darkdim_seeds.jld2")["df"]
    seed_points = seed_points[seed_points.ca3 .< 0, :]
    init_samples = make_init_samples(posterior, seed_points[1:10, :], 10_000)

    FileIO.save(name * "_init_samples.jld2", Dict(String(a)=>init_samples[a] for a in keys(init_samples)))
    whack_samples = whack_many_moles(posterior, init_samples, target_samplesize=100_000, cache_dir=name)
    FileIO.save(name * ".jld2", Dict(String(a)=>whack_samples[a] for a in keys(whack_samples)))
else
    if lowercase(args["task"]) == "profile"
        result = Newtrinos.profile(likelihood, priors, vars_to_scan, p, cache_dir=name)
    elseif lowercase(args["task"]) == "scan"
        result = Newtrinos.scan(likelihood, priors, vars_to_scan, p)
    end
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
        ax.title = args["task"] * ": " * join(args["experiments"], " + ")
        save(name * ".png", fig)
    end
end
