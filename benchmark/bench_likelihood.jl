using ArgParse
using BenchmarkTools
using LinearAlgebra
using Distributions
using DensityInterface
using BAT
using MeasureBase
using ADTypes
using DataStructures
using Accessors
using Newtrinos
import ForwardDiff

include(joinpath(@__DIR__, "..", "src", "analysis", "cli_common.jl"))

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
    end

    return parse_args(s)
end

args = parse_command_line()

println("Configuring physics and experiments...")
physics = configure_physics(args["ordering"])
experiments = configure_experiments(args["experiments"], physics)

params = Newtrinos.get_params(experiments)
priors = Newtrinos.get_priors(experiments)

likelihood = Newtrinos.generate_likelihood(experiments)

println()
println("=" ^60)
println("LIKELIHOOD BENCHMARK")
println("=" ^60)
println("  Experiments: ", join(args["experiments"], ", "))
println("  Ordering:    ", args["ordering"])
println("  Parameters:  ", length(keys(params)))
println()

# Warmup
logdensityof(likelihood, params)

println("=== Likelihood evaluation ===")
b_llh = @benchmark logdensityof($likelihood, $params)
display(b_llh)
println()

# Gradient via ForwardDiff
println("=== Gradient (ForwardDiff) ===")
grad_f(p) = logdensityof(likelihood, p)
ForwardDiff.gradient(grad_f, params)  # warmup

b_grad = @benchmark ForwardDiff.gradient($grad_f, $params)
display(b_grad)
println()

# Summary
t_llh = median(b_llh).time / 1e6
t_grad = median(b_grad).time / 1e6
n_params = length(keys(params))

println("=== Summary ===")
println("  Likelihood:  $(round(t_llh, digits=2)) ms (median)")
println("  Gradient:    $(round(t_grad, digits=2)) ms (median)")
println("  Ratio grad/llh: $(round(t_grad / t_llh, digits=1))×")
println("  Parameters:  $n_params")
println("  Grad cost per param: $(round(t_grad / n_params, digits=2)) ms")
