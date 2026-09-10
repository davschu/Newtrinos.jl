#this file is for testing different optimizer settings for the nsi analysis via roundtrips with injected truths

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "Newtrinos.jl"))

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
using CairoMakie
using Random
import Optim
import ForwardDiff
using ProgressMeter
using ValueShapes

include(joinpath(@__DIR__, "optimizer_tests_common.jl"))

### CLI ###

function parse_command_line()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--experiments"
        help = "List of experiments to run"
        nargs = '+'
        default = ["deepcore"]

        "--hypotheses"
        help = "Subset of hypotheses to run, by name: eμ, eτ, μτ, ee_μμ, ττ_μμ, all. Space-separated; default runs all 6. Note 'all' here selects the single cp_all hypothesis (every NSI param free), not \"every hypothesis\" -- omit --hypotheses entirely to run all 6."
        nargs = '+'
        default = ["eμ", "eτ", "μτ", "ee_μμ", "ττ_μμ", "all"]

        "--name"
        help = "Name for outputs"
        arg_type = String
        default = "optimizer_tests"

        "--n-truths"
        help = "Number of injected truths per hypothesis"
        arg_type = Int
        default = 20

        "--seed"
        help = "RNG seed for truth injection and random starts"
        arg_type = Int
        default = 1234

        "--truth-mode"
        help = "How injected truths are generated: controlled (default -- only the hypothesis's free NSI parameter(s) [+ theta23/deltam31 if --vary-physics-truth] are drawn randomly, everything else pinned at nominal values) or random (every free parameter, including all nuisances, drawn randomly from its prior)"
        arg_type = String
        default = "controlled"

        "--vary-physics-truth"
        help = "In --truth-mode controlled, also draw theta23's and deltam31's truths randomly from their priors instead of pinning them at their nominal values. No effect in --truth-mode random (they already vary there)."
        action = :store_true

        "--data-mode"
        help = "Injected-data mode: asimov (expected values) or toy (Poisson-fluctuated)"
        arg_type = String
        default = "asimov"

        "--fit-method"
        help = "Fit method: optim (tunable BAT.bat_findmode + OptimAlg, see --optimizer etc.), newtrinos (standard Newtrinos.find_mle), or ext (Newtrinos.find_mle_ext) -- newtrinos/ext are unmodified, for comparison"
        arg_type = String
        default = "optim"

        "--optimizer"
        help = "Optim.jl algorithm: lbfgs, gradientdescent, conjugategradient, neldermead (only used when --fit-method optim)"
        arg_type = String
        default = "lbfgs"

        "--objective"
        help = "What bat_findmode maximizes for --fit-method optim: posterior (likelihood x prior, i.e. MAP -- default) or likelihood (likelihood only, via a flat prior over the same bounds -- i.e. MLE). Ignored for --fit-method newtrinos/ext."
        arg_type = String
        default = "posterior"

        "--iterations"
        help = "Max iterations for the local optimizer (outer_iterations for Fminbox)"
        arg_type = Int
        default = 2000

        "--step-size"
        help = "Initial step size (alphaguess) for first-order algorithms"
        arg_type = Float64
        default = 1.0

        "--g-tol"
        help = "Gradient-norm convergence tolerance"
        arg_type = Float64
        default = 1e-6

        "--f-tol"
        help = "Function-value relative convergence tolerance"
        arg_type = Float64
        default = 0.0

        "--x-tol"
        help = "Parameter-value relative convergence tolerance"
        arg_type = Float64
        default = 0.0

        "--seed-strategy"
        help = "Starting-point strategy: random, pso, or sa"
        arg_type = String
        default = "random"

        "--global-iterations"
        help = "Iteration budget for the PSO/SA global seed-search stage"
        arg_type = Int
        default = 100

        "--nseeds"
        help = "Number of independent local-optimizer starts+fits per roundtrip; keeps the best log_posterior (mirrors Newtrinos.profile's nseeds)"
        arg_type = Int
        default = 1

        "--pso-particles"
        help = "Swarm size (n_particles) for Optim.ParticleSwarm during the global seed-search stage (only used with --seed-strategy pso)"
        arg_type = Int
        default = 10

        "--plot"
        help = "Enable plotting"
        action = :store_true

        "--trace-plots"
        help = "Save PSO/LBFGS convergence-trace plots per hypothesis (adds trace-collection overhead to do_roundtrip; requires --plot, otherwise has no effect)"
        action = :store_true

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
run_start_time = time()

name           = args["name"]
n_truths       = args["n-truths"]
data_mode      = lowercase(args["data-mode"])
truth_mode     = lowercase(args["truth-mode"])
vary_physics   = args["vary-physics-truth"]
fit_method     = lowercase(args["fit-method"])
optimizer_name = lowercase(args["optimizer"])
objective      = lowercase(args["objective"])
iterations     = args["iterations"]
step_size      = args["step-size"]
g_tol          = args["g-tol"]
f_tol          = args["f-tol"]
x_tol          = args["x-tol"]
seed_strategy  = lowercase(args["seed-strategy"])
global_iters   = args["global-iterations"]
nseeds         = args["nseeds"]
pso_particles  = args["pso-particles"]
trace_plots    = args["trace-plots"]
n_workers      = args["workers"]
n_threads      = args["threads"]
use_distributed = n_workers > 1

data_mode in ("asimov", "toy") || error("--data-mode must be 'asimov' or 'toy', got '$data_mode'")
fit_method in ("optim", "newtrinos", "ext") || error("--fit-method must be 'optim', 'newtrinos' or 'ext', got '$fit_method'")
seed_strategy in ("random", "pso", "sa") || error("--seed-strategy must be 'random', 'pso', or 'sa', got '$seed_strategy'")
objective in ("posterior", "likelihood") || error("--objective must be 'posterior' or 'likelihood', got '$objective'")
truth_mode in ("controlled", "random") || error("--truth-mode must be 'controlled' or 'random', got '$truth_mode'")

algorithm = make_algorithm(optimizer_name, step_size)

### PHYSICS CONFIG ###

function configure_physics()
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
    (; deepcore = Newtrinos.deepcore.configure(physics))
end

# Set up distributed workers if requested
if use_distributed
    addprocs(n_workers; exeflags="--threads=$n_threads")

    @everywhere begin
        using Distributions
        using DensityInterface
        using BAT
        using DataStructures
        using MeasureBase
        using ADTypes
        using Newtrinos
        using ValueShapes
        using Accessors
        import Optim
        import ForwardDiff
        include(joinpath($(@__DIR__), "optimizer_tests_common.jl"))
    end
end

experiments = configure_physics()

ftype=Float64
p= Newtrinos.get_params(experiments)
priors = Newtrinos.get_priors(experiments)
conditional_vars = Dict(:δCP=>0.0, :θ₁₂=>ftype(33.62 * π/180), :θ₁₃=>ftype(8.54  * π/180), :Δm²₂₁=>ftype(7.40e-5))
priors = Newtrinos.condition(priors, conditional_vars, p)
@reset priors.θ₂₃   = Uniform(ftype(30*π/180), ftype(60*π/180))
@reset priors.Δm²₃₁ = Uniform(ftype(0.93e-3 + 7.40e-5), ftype(3.93e-3 + 7.40e-5))

@reset priors.deepcore_lifetime        = Uniform(ftype(0), ftype(3.8))
@reset priors.deepcore_ice_scattering  = Truncated(Normal(ftype(1), ftype(0.1)), ftype(0.9), ftype(1.1))
@reset priors.deepcore_ice_absorption  = Truncated(Normal(ftype(1), ftype(0.1)), ftype(0.9), ftype(1.1))
@reset priors.deepcore_opt_eff_lateral = Truncated(Normal(ftype(0), ftype(1.)), ftype(-2), ftype(2.5))

# NSI parameter bounds — Standard parameterization
@reset priors.Δ_eμ     = Uniform(-ftype(5),   ftype(5))
@reset priors.Δ_τμ     = Uniform(-ftype(0.1), ftype(0.1))
@reset priors.ε_eμ_abs = Uniform(ftype(0),    ftype(0.3))
@reset priors.ε_eτ_abs = Uniform(ftype(0),    ftype(0.35))
@reset priors.ε_μτ_abs = Uniform(ftype(0),    ftype(0.07))
@reset priors.δ_eμ     = Uniform(ftype(0),    ftype(2π)) #used [0,π] for 0,180 degree with two gridpoints, but use [0,2π] for general grid scan
@reset priors.δ_eτ     = Uniform(ftype(0),    ftype(2π))
@reset priors.δ_μτ     = Uniform(ftype(0),    ftype(2π))

#conditioned priors for 1-by-1 hypotheses
cp_eμ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eτ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eτ => 0.0, :δ_μτ => 0.0), p)
cp_eτ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eμ => 0.0, :δ_μτ => 0.0), p)
cp_μτ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0), p)
cp_ee_μμ = Newtrinos.condition(priors, Dict(:Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0, :δ_μτ => 0.0), p)
cp_ττ_μμ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :ε_μτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0, :δ_μτ => 0.0), p)

#all params non-zero hypothesis
cp_all = deepcopy(priors) 

all_cps   = [cp_eμ, cp_eτ, cp_μτ, cp_ee_μμ, cp_ττ_μμ, cp_all]
all_names = ["eμ", "eτ", "μτ", "ee_μμ", "ττ_μμ", "all"] # 6 hypotheses
interesting_nuisances = [:Δ_eμ, :Δ_τμ, :Δm²₃₁, :δ_eμ, :δ_eτ, :δ_μτ, :ε_eμ_abs, :ε_eτ_abs, :ε_μτ_abs, :θ₂₃]
nsi_param_names = (:Δ_eμ, :Δ_τμ, :ε_eμ_abs, :ε_eτ_abs, :ε_μτ_abs, :δ_eμ, :δ_eτ, :δ_μτ)

# Filter to the CLI-selected subset (default: all 6), always in the canonical
# eμ/eτ/μτ/ee_μμ/ττ_μμ/all order regardless of the order flags were passed in.
selected_hypotheses = args["hypotheses"]
for h in selected_hypotheses
    h in all_names || error("--hypotheses: unknown hypothesis '$h'. Choose from: $(join(all_names, ", "))")
end
keep_idx = [idx for (idx, nm) in enumerate(all_names) if nm in selected_hypotheses]
cps  = all_cps[keep_idx]
scan = all_names[keep_idx]
n_hyp = length(cps)


### ROUNDTRIPS WITH INJECTED TRUTH ###

rng = Random.Xoshiro(args["seed"])

# Pre-generate all truths and per-task RNGs from the single seeded `rng` before
# parallelizing, so results are reproducible regardless of thread/worker scheduling order.
work = Tuple{Int,Int}[(j, i) for j in 1:n_hyp for i in 1:n_truths]
n_work = length(work)

truth_seeds = Vector{Any}(undef, n_work)
start_rngs  = Vector{Random.Xoshiro}(undef, n_work)
for k in 1:n_work
    j, _ = work[k]
    prior_dist = distprod(;cps[j]...)
    truth_seeds[k] = truth_mode == "controlled" ?
        controlled_truth(rng, p, cps[j], nsi_param_names; vary_physics=vary_physics) :
        rand(rng, prior_dist)
    start_rngs[k]  = Random.Xoshiro(rand(rng, UInt))
end

random_truth_seeds                 = Matrix{Any}(undef, n_hyp, n_truths)
likelihoods_for_random_truth_seeds = Matrix{Any}(undef, n_hyp, n_truths)
scan_results                       = Matrix{Any}(undef, n_hyp, n_truths)
random_truth_posteriors            = Matrix{Any}(undef, n_hyp, n_truths)
fit_quality_A                      = Matrix{Float64}(undef, n_hyp, n_truths)
fit_quality_per_param              = Matrix{Any}(undef, n_hyp, n_truths)
all_traces                         = Matrix{Any}(undef, n_hyp, n_truths)

if use_distributed
    @everywhere experiments = $experiments
    results_flat = @showprogress pmap(1:n_work) do k
        do_roundtrip(k, work, truth_seeds, data_mode, experiments, cps, seed_strategy,
                     global_iters, start_rngs, fit_method, algorithm, iterations, g_tol, f_tol, x_tol;
                     nseeds=nseeds, pso_particles=pso_particles, trace_plots=trace_plots, objective=objective)
    end
else
    results_flat = Vector{Any}(undef, n_work)
    @showprogress Threads.@threads for k in 1:n_work
        results_flat[k] = do_roundtrip(k, work, truth_seeds, data_mode, experiments, cps, seed_strategy,
                                        global_iters, start_rngs, fit_method, algorithm, iterations, g_tol, f_tol, x_tol;
                                        nseeds=nseeds, pso_particles=pso_particles, trace_plots=trace_plots, objective=objective)
    end
end

for (j, i, truth_param, truth_likelihood, fit_result, posterior, A, per_param, all_seed_traces) in results_flat
    random_truth_seeds[j,i]                 = truth_param
    likelihoods_for_random_truth_seeds[j,i]  = truth_likelihood
    scan_results[j,i]                        = fit_result
    random_truth_posteriors[j,i]             = posterior
    fit_quality_A[j,i]                       = A
    fit_quality_per_param[j,i]               = per_param
    all_traces[j,i]                          = all_seed_traces
end

for j in 1:n_hyp
    println("scans for $(scan[j]) completed")
end


### SUMMARY OUTPUT ###

# quantify fit quality with A = 1/n_params * Σ[|param_i_fit - param_i_truth|/(|param_i_truth| + |param_i_fit| + ϵ)]^2     for each injected truth
# build mean for combined information of several roundtrips with same settings A_bar = 1/n_truths * Σ[A]
# (computed per-roundtrip in do_roundtrip/fit_quality; here we aggregate across all n_hyp*n_truths roundtrips)

n_total = n_hyp * n_truths
converged_flags = [scan_results[j,i][4] for j in 1:n_hyp, i in 1:n_truths]
n_converged     = count(identity, converged_flags)
n_not_converged = n_total - n_converged

iters_all = [scan_results[j,i][5] for j in 1:n_hyp, i in 1:n_truths]
iters_valid = filter(x -> !ismissing(x), vec(iters_all))
mean_iterations = isempty(iters_valid) ? missing : sum(iters_valid) / length(iters_valid)

A_bar = sum(fit_quality_A) / length(fit_quality_A)

# Per-parameter quality: mean squared relative error for each parameter, averaged over
# every roundtrip where that parameter was free (which hypotheses vary by parameter).
per_param_sums   = Dict{Symbol, Float64}()
per_param_counts = Dict{Symbol, Int}()
for j in 1:n_hyp, i in 1:n_truths
    for (k, v) in fit_quality_per_param[j,i]
        per_param_sums[k]   = get(per_param_sums, k, 0.0) + v
        per_param_counts[k] = get(per_param_counts, k, 0) + 1
    end
end
per_param_quality = Dict(k => per_param_sums[k] / per_param_counts[k] for k in keys(per_param_sums))

# A_bar restricted to physics/NSI parameters (excludes nuisances). With controlled
# truths, nuisance truths are pinned at nominal (often 0.0), where fit_quality's
# relative-error formula is degenerate (any nonzero fit reads as ~100% error regardless
# of absolute closeness) -- A_bar_physical avoids that artifact by only averaging over
# the parameters actually varied by --truth-mode controlled.
physical_param_names = (nsi_param_names..., :θ₂₃, :Δm²₃₁)
physical_keys = [k for k in keys(per_param_quality) if k in physical_param_names]
A_bar_physical = isempty(physical_keys) ? NaN :
    sum(per_param_quality[k] for k in physical_keys) / length(physical_keys)

total_elapsed_s = time() - run_start_time

summary_lines = String[]
push!(summary_lines, "=== optimizer_tests summary ===")
push!(summary_lines, "name: $name")
push!(summary_lines, "total roundtrips: $n_total ($n_hyp hypotheses [$(join(scan, ", "))] x $n_truths truths)")
push!(summary_lines, "total wall time: $(round(total_elapsed_s, digits=1))s ($(round(total_elapsed_s/60, digits=1)) min)")
push!(summary_lines, "converged: $n_converged / $n_total")
push!(summary_lines, "not converged: $n_not_converged / $n_total")
push!(summary_lines, "tolerances: g_tol=$g_tol f_tol=$f_tol x_tol=$x_tol")
push!(summary_lines, "mean iterations: $(mean_iterations)" * (fit_method == "optim" ? "" : " (n/a for --fit-method $fit_method)"))
push!(summary_lines, "fit quality A_bar (mean over all roundtrips): $A_bar")
push!(summary_lines, "fit quality A_bar_physical (mean over NSI/theta23/deltam31 params only): $A_bar_physical")
push!(summary_lines, "fit quality per parameter (mean squared relative error):")
for k in sort(collect(keys(per_param_quality)))
    push!(summary_lines, "  $k: $(per_param_quality[k])")
end
summary_text = join(summary_lines, "\n")

println(summary_text)


### PLOTTING ###

if args["plot"]
    fig = Figure(size=(2000, 400 * (2 + length(interesting_nuisances))))

    for j in 1:n_hyp
        #log_posterior -- when --objective likelihood, scan_results[j,i][2] is the
        #flat-prior pseudo-posterior (≈ llh + const, see local_find_mle's docstring),
        #NOT the true MAP posterior log_p_truth is computed from -- label this clearly
        #so the mismatch is visible rather than silently misread as MAP recovery.
        fit_ylabel = objective == "likelihood" ? "best fit flat-prior log_post (∝ llh)" : "best fit value log_p"
        ax = Axis(fig[1, j], xlabel = "truth value log_p", ylabel = fit_ylabel, title = "$(scan[j])")
        log_p_truth = [logdensityof(random_truth_posteriors[j,i], random_truth_seeds[j,i]) for i in 1:n_truths]
        log_p_scan  = [scan_results[j,i][2] for i in 1:n_truths]
        scatter!(ax, log_p_truth, log_p_scan)
        ablines!(ax, 0, 1, color=:red) # 0 = x-intercept, 1 = slope

        #log likelihood
        ax2 = Axis(fig[2, j], xlabel = "truth value llh", ylabel = "best fit value llh")
        llh_truth = [logdensityof(likelihoods_for_random_truth_seeds[j,i], random_truth_seeds[j,i]) for i in 1:n_truths]
        llh_scan  = [scan_results[j,i][1] for i in 1:n_truths]
        scatter!(ax2, llh_truth, llh_scan)
        ablines!(ax2, 0, 1, color=:red)

        for (nui, nuisance) in enumerate(interesting_nuisances)
            ax_nui = Axis(fig[2 + nui, j], xlabel="truth value $(nuisance)", ylabel = "best fit value $(nuisance)")
            idxs = [i for i in 1:n_truths if haskey(random_truth_seeds[j,i], nuisance) && haskey(scan_results[j,i][3], nuisance)]
            if !isempty(idxs)
                nui_truth = [random_truth_seeds[j,i][nuisance] for i in idxs]
                nui_scan  = [scan_results[j,i][3][nuisance] for i in idxs]
                scatter!(ax_nui, nui_truth, nui_scan)
                ablines!(ax_nui, 0, 1, color=:red)
            end
        end
    end

    save(name * "_roundtrips.pdf", fig)
end


### OPTIMIZER TRACE PLOTTING ###

# Per hypothesis: one combined PSO panel and one combined LBFGS panel, overlaying every
# truth (one color per truth) and every seed's individual trajectory (nseeds>1) in that
# truth's color, with a dashed horizontal line at that truth's -logposterior.
if args["plot"] && trace_plots
    fig2 = Figure(size=(1200, 400 * n_hyp))
    truth_colors = Makie.wong_colors()

    for j in 1:n_hyp
        ax_pso   = Axis(fig2[j, 1], xlabel="PSO generation", ylabel="-logposterior (obj)",
                         yscale=log10, title="$(scan[j]) -- PSO seed search")
        ax_lbfgs = Axis(fig2[j, 2], xlabel="LBFGS iteration", ylabel="-logposterior",
                         yscale=log10, title="$(scan[j]) -- LBFGS polish")

        for i in 1:n_truths
            color = truth_colors[mod1(i, length(truth_colors))]
            truth_neglogpost = -logdensityof(random_truth_posteriors[j,i], random_truth_seeds[j,i])
            for (pso_trace, lbfgs_trace) in all_traces[j,i]
                !isempty(pso_trace)   && lines!(ax_pso,   0:length(pso_trace)-1,   max.(pso_trace, eps());   color=color, label="truth $i")
                !isempty(lbfgs_trace) && lines!(ax_lbfgs, 0:length(lbfgs_trace)-1, max.(lbfgs_trace, eps()); color=color, label="truth $i")
            end
            hlines!(ax_pso,   [max(truth_neglogpost, eps())]; color=color, linestyle=:dash)
            hlines!(ax_lbfgs, [max(truth_neglogpost, eps())]; color=color, linestyle=:dash)
        end
        axislegend(ax_lbfgs, merge=true, unique=true, position=:rt)
    end

    save(name * "_optimizer_traces.pdf", fig2)
end


### SAVE RESULTS ###

results_dict = Dict{String, Any}(
    "settings" => Dict(
        "data_mode"          => data_mode,
        "truth_mode"         => truth_mode,
        "vary_physics"       => vary_physics,
        "fit_method"         => fit_method,
        "optimizer"          => optimizer_name,
        "objective"          => objective,
        "iterations"         => iterations,
        "step_size"          => step_size,
        "g_tol"              => g_tol,
        "f_tol"              => f_tol,
        "x_tol"              => x_tol,
        "seed_strategy"      => seed_strategy,
        "global_iterations"  => global_iters,
        "nseeds"             => nseeds,
        "pso_particles"      => pso_particles,
        "trace_plots"        => trace_plots,
        "n_truths"           => n_truths,
        "seed"               => args["seed"],
        "experiments"        => args["experiments"],
    ),
    "hypotheses"           => scan,
    "truth_params"         => random_truth_seeds,
    "fitted_llh"           => [scan_results[j,i][1] for j in 1:n_hyp, i in 1:n_truths],
    "fitted_log_posterior" => [scan_results[j,i][2] for j in 1:n_hyp, i in 1:n_truths],
    "fitted_params"        => [scan_results[j,i][3] for j in 1:n_hyp, i in 1:n_truths],
    "converged"            => [scan_results[j,i][4] for j in 1:n_hyp, i in 1:n_truths],
    "n_iterations"         => iters_all,
    "fit_quality_A"        => fit_quality_A,
    "fit_quality_per_param_per_roundtrip" => fit_quality_per_param,
    "optimizer_traces"     => all_traces,
    "summary" => Dict(
        "n_total"                  => n_total,
        "n_converged"              => n_converged,
        "n_not_converged"          => n_not_converged,
        "mean_iterations"          => mean_iterations,
        "fit_quality_A_bar"        => A_bar,
        "fit_quality_A_bar_physical" => A_bar_physical,
        "fit_quality_per_param"    => per_param_quality,
        "total_elapsed_s"          => total_elapsed_s,
        "text"                     => summary_text,
    ),
)

FileIO.save(name * "_roundtrip_results.jld2", results_dict)
