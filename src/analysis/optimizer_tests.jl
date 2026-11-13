# Compares optimizers on the same problems. For each of several injected truth points of one hypothesis of a setup
# file (--setup, e.g. nsi/setups/nsi_deepcore_3y.jl) every optimizer given in --optimizers is run (--nrepeats times,
# each with its own random seed) with the same objective, bounds and start point and the same evaluation budget, and
# compared with
#   * the log-posterior of the injected truth (a successful optimizer ends at or above it), and
#   * a truth-seeded reference fit (--reference, by default a long LBFGS fit started at the truth).
# Timing, number of objective evaluations, per-parameter errors and the convergence trace (best objective value so far
# against evaluations) are recorded for every run.
#
# Optimizers are named entries of the registry in optimizer_tests_common.jl (`OPTIMIZERS`); a spec is "name" or
# "name:key=val,key=val", e.g.
#   --optimizers optim-lbfgs "optim-pso:particles=50" mh-de "mh-de:N=30,F=0.5" bbo-xnes
# (run with --list-optimizers to see the names and their settings). Every run is one task; with --workers N the tasks
# are distributed over N single-threaded worker processes, so timings are comparable.
#
# Start point of every optimizer: the setup's nominal parameters (blind, NOT the truth). Gradient-free global methods
# search the box of the free parameters' prior bounds and start from nothing (optim-pso / optim-sa from that start).
#
# Output (in --results-dir, default results/<setup name>/optimizer_tests):
#   <name>_<hypothesis>_summary.txt      table per optimizer spec
#   <name>_<hypothesis>_results.jld2     every run (parameters, scores, calls, time, trace), truths, reference fits
#   <name>_<hypothesis>_traces.pdf       best objective value so far vs evaluations, one panel per truth
#   <name>_<hypothesis>_summary.pdf      distance to the reference per optimizer, and time vs distance

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

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
using Statistics
import Optim
import ForwardDiff
using ProgressMeter
using ValueShapes

include(joinpath(@__DIR__, "optimizer_tests_common.jl"))

### CLI ###

function parse_command_line()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--setup"
        help = "Setup file defining the physics, experiments, priors and hypotheses (it defines setup(; ftype) -- see nsi/setups/nsi_deepcore_3y.jl for the format)"
        arg_type = String
        default = ""

        "--list-optimizers"
        help = "Print the registered optimizers with their settings and exit"
        action = :store_true

        "--optimizers"
        help = "Optimizer specs to compare, space-separated, each \"name\" or \"name:key=val,key=val\" (see --list-optimizers); the same optimizer may appear several times with different settings"
        arg_type = String
        nargs = '+'
        default = ["optim-lbfgs", "optim-pso:particles=20", "mh-de"]

        "--hypothesis"
        help = "Name of the hypothesis (from the setup file) to test; default: the first one of the setup"
        arg_type = String
        default = ""

        "--budget"
        help = "Default maximum number of objective evaluations per run of a global optimizer that gives no `evals` of its own (e.g. mh-de:evals=3000 overrides it); local optimizers are limited by their own `iterations`"
        arg_type = Int
        default = 2000

        "--n-truths"
        help = "Number of injected truth points"
        arg_type = Int
        default = 3

        "--nrepeats"
        help = "Independent runs (different random seeds) of every optimizer on every truth point"
        arg_type = Int
        default = 1

        "--reference"
        help = "Optimizer spec of the truth-seeded reference fit (started AT the truth)"
        arg_type = String
        default = "optim-lbfgs:iterations=1000"

        "--tol"
        help = "A run counts as 'reached the reference' if its score is within this value of the reference's score (log-posterior units)"
        arg_type = Float64
        default = 0.1

        "--fluctuate-nuisances"
        help = "Also draw every other free parameter of the hypothesis from its prior for each truth point (the tested parameters' truth is unaffected)"
        action = :store_true

        "--vary-truth-params"
        help = "Names of additional parameters whose truth is drawn randomly from their priors, besides the tested parameters of the hypothesis, e.g. --vary-truth-params θ₂₃ Δm²₃₁"
        arg_type = String
        nargs = '*'
        default = String[]

        "--data-mode"
        help = "Injected-data mode: asimov (expected values) or toy (Poisson-fluctuated)"
        arg_type = String
        default = "asimov"

        "--objective"
        help = "What the optimizers maximize: posterior (likelihood x prior, i.e. MAP -- default) or likelihood (likelihood only). Scores are log_posterior resp. llh."
        arg_type = String
        default = "posterior"

        "--seed"
        help = "Base RNG seed for the truth points and the optimizer runs"
        arg_type = Int
        default = 1234

        "--name"
        help = "Name prefix for the output files"
        arg_type = String
        default = "optimizer_tests"

        "--results-dir"
        help = "Output folder (default: results/<setup name>/optimizer_tests, relative to the working directory)"
        arg_type = String
        default = ""

        "--no-warmup"
        help = "Skip the warm-up run of every optimizer on every worker (without it the first run of each optimizer includes compilation time in its timing)"
        action = :store_true

        "--workers"
        help = "Number of distributed single-threaded worker processes (default: 1, no distributed: tasks run one after the other in this process)"
        arg_type = Int
        default = 1
    end

    return parse_args(s)
end

args = parse_command_line()

if args["list-optimizers"]
    for nm in sort(collect(keys(OPTIMIZERS)))
        e = OPTIMIZERS[nm]
        println(rpad(nm, 26), rpad(String(e.kind), 8), isempty(e.defaults) ? "(no settings)" : join(("$k=$(repr(v))" for (k, v) in pairs(e.defaults)), ", "))
    end
    exit(0)
end
isempty(args["setup"]) && error("--setup is required (e.g. --setup nsi/setups/nsi_deepcore_3y.jl)")

run_start_time = time()
name       = args["name"]
n_truths   = args["n-truths"]
nrepeats   = args["nrepeats"]
budget     = args["budget"]
tol        = args["tol"]
objective  = lowercase(args["objective"])
data_mode  = lowercase(args["data-mode"])
vary_truth_params = Tuple(Symbol.(args["vary-truth-params"]))
n_workers  = args["workers"]
use_distributed = n_workers > 1

data_mode in ("asimov", "toy") || error("--data-mode must be 'asimov' or 'toy', got '$data_mode'")
objective in ("posterior", "likelihood") || error("--objective must be 'posterior' or 'likelihood', got '$objective'")

specs    = [parse_optimizer_spec(s) for s in args["optimizers"]]
ref_spec = parse_optimizer_spec(args["reference"])
OPTIMIZERS[ref_spec.name].kind == :local || @warn "the reference is normally a local fit started at the truth; got a global optimizer ($(ref_spec.name))"
n_specs  = length(specs)
labels   = [s.label for s in specs]
length(unique(labels)) == n_specs || error("--optimizers: every spec must be different (got duplicates: $labels)")

# Distributed workers: single-threaded, so that wall times of different optimizers are comparable
if use_distributed
    addprocs(n_workers; exeflags = "--threads=1")

    @everywhere begin
        using Distributions
        using DensityInterface
        using BAT
        using DataStructures
        using MeasureBase
        using ADTypes
        using Newtrinos
        using ValueShapes
        import Optim
        import ForwardDiff
    end
    # SEPARATE block: `using Accessors` must be loaded before the helper file defines methods using it
    @everywhere begin
        using Accessors
        include(joinpath($(@__DIR__), "optimizer_tests_common.jl"))
    end
end

### SETUP (physics, experiments, priors, hypotheses: defined by the setup file, see --setup) ###

include(abspath(args["setup"]))
S = finalize_setup(setup(; ftype = Float64))

hyp_name = isempty(args["hypothesis"]) ? first(keys(S.hypotheses)) : args["hypothesis"]
hyp = hypothesis_info(S, hyp_name)
experiments = S.experiments
p = S.params
prior = hyp.prior
prior_dist = distprod(; prior...)
tested_names = hyp.tested
fluctuated_names = hyp.fluctuated
for k in vary_truth_params
    k in keys(prior) || error("--vary-truth-params: unknown parameter $k")
end
plotted_params = unique([collect(tested_names); collect(vary_truth_params)])   # parameters whose recovery is reported
results_dir = isempty(args["results-dir"]) ? joinpath("results", S.name, "optimizer_tests") : args["results-dir"]
mkpath(results_dir)
out_prefix = joinpath(results_dir, "$(name)_$(hyp_name)")

println("setup $(S.name), hypothesis '$hyp_name' (tested: $(join(tested_names, ", "))), $n_truths truths x $n_specs optimizers x $nrepeats repeats, budget $budget evaluations")

### TRUTH POINTS ###

rng = Random.Xoshiro(args["seed"])
truths = Vector{Any}(undef, n_truths)
for i in 1:n_truths
    t = controlled_truth(rng, p, prior, tested_names; extra_names = vary_truth_params)
    truths[i] = args["fluctuate-nuisances"] ? fluctuate_nuisances(rng, t, prior, fluctuated_names) : t
end

### TASKS ###

start = _fix_conditioned(p, prior)   # nominal parameters with the prior's pinned parameters forced to their pinned values
ctx = (; experiments, prior, prior_dist, start, truths, data_mode, objective, budget, specs, ref_spec, seed = args["seed"])
use_distributed && @everywhere ctx = $ctx

# One task = one optimizer run on one truth point: t = (; i, s, r) -- truth i, spec s (0: the reference fit), repeat r.
# `ctx` is a per-process global (see roundtrips.jl for why it is broadcast once instead of captured).
@everywhere function opt_task(t)
    c = ctx
    truth = c.truths[t.i]
    data = c.data_mode == "toy" ? Newtrinos.generate_toy_data(c.experiments, truth) : Newtrinos.generate_asimov_data(c.experiments, truth)
    L = Newtrinos.generate_likelihood(c.experiments, data)
    score(r) = c.objective == "likelihood" ? r.llh : r.log_post
    seed = c.seed + 1_000_003 * t.i + 10_007 * t.s + t.r
    Random.seed!(seed)                       # the (task-local) global RNG: used by Optim's PSO and BlackBoxOptim
    if t.s == 0                              # reference fit, started at the truth; also the truth's own score
        r = run_optimizer(c.ref_spec, L, c.prior_dist, truth; rng = Random.Xoshiro(seed), objective = c.objective)
        truth_score = c.objective == "likelihood" ? logdensityof(L, truth) : logdensityof(PosteriorMeasure(L, c.prior_dist), truth)
        return merge(r, (; i = t.i, s = 0, r = t.r, score = score(r), truth_score))
    end
    r = run_optimizer(c.specs[t.s], L, c.prior, c.start; budget = c.budget, rng = Random.Xoshiro(seed), objective = c.objective)
    merge(r, (; i = t.i, s = t.s, r = t.r, score = score(r), truth_score = NaN))
end

# Warm-up: run every optimizer once, tiny, on every worker (and here) so that compilation does not enter the timings
@everywhere function warmup_task(w)
    c = ctx
    truth = c.truths[1]
    L = Newtrinos.generate_likelihood(c.experiments, Newtrinos.generate_asimov_data(c.experiments, truth))
    for sp in vcat(c.specs, [c.ref_spec])
        d = OPTIMIZERS[sp.name].defaults
        st = sp.settings
        :iterations in keys(d) && (st = merge(st, (; iterations = 2)))
        :iters in keys(d) && (st = merge(st, (; iters = 1)))
        run_optimizer((; sp..., settings = st), L, c.prior, c.start; budget = 40, rng = Random.Xoshiro(1), objective = c.objective)
    end
    nothing
end

if !args["no-warmup"]
    println("warm-up (compilation) ...")
    if use_distributed
        @sync for w in workers()
            @async remotecall_wait(warmup_task, w, 1)
        end
    else
        warmup_task(1)
    end
end

ref_tasks = [(; i, s = 0, r = 1) for i in 1:n_truths]
opt_tasks = [(; i, s, r) for i in 1:n_truths for s in 1:n_specs for r in 1:nrepeats]
run_tasks(tasks) = use_distributed ? @showprogress(pmap(opt_task, tasks)) : @showprogress(map(opt_task, tasks))
println("running $(length(ref_tasks)) reference fits + $(length(opt_tasks)) optimizer runs ...")
ref_records = run_tasks(ref_tasks)
records     = run_tasks(opt_tasks)

### EVALUATION ###

b = extract_bounds(prior, p)
width = Dict(k => b.hi_vec[b.free_x_index[findfirst(==(k), b.all_keys)]] - b.lo_vec[b.free_x_index[findfirst(==(k), b.all_keys)]]
             for k in plotted_params if k in b.all_keys && b.is_free_mask[findfirst(==(k), b.all_keys)])
ref_by_truth = Dict(r.i => r for r in ref_records)

# excess_*: how much worse than the truth / the reference a run is (>= 0: worse or equal; < 0: better)
excess_truth(r) = ref_by_truth[r.i].truth_score - r.score
excess_ref(r)   = ref_by_truth[r.i].score - r.score
param_err(r, k) = haskey(width, k) ? abs(Float64(r.params[k]) - Float64(truths[r.i][k])) / width[k] : NaN
fin(x) = filter(isfinite, x)
med(x) = isempty(fin(x)) ? NaN : median(fin(x))

summary_lines = String[]
push!(summary_lines, "=== optimizer_tests summary ===")
push!(summary_lines, "setup: $(S.name)   hypothesis: $hyp_name (tested: $(join(tested_names, ", ")))")
push!(summary_lines, "$n_truths truths x $nrepeats repeats, budget $budget evaluations, objective $objective, data $data_mode" *
                     (args["fluctuate-nuisances"] ? ", nuisances fluctuated" : "") *
                     (isempty(vary_truth_params) ? "" : ", varied truth params: $(join(vary_truth_params, ", "))"))
push!(summary_lines, "reference: $(ref_spec.label) at the truth; score = " * (objective == "likelihood" ? "llh" : "log_posterior"))
push!(summary_lines, "excess = reference score - score of the run (0: as good as the reference; > 0: worse; < 0: better); the plots show Δχ²_ref = 2 × excess")
push!(summary_lines, "below truth = run ended with a score lower than the injected truth's (a failure of the optimizer)")
push!(summary_lines, "reached = excess <= $tol")
push!(summary_lines, "param err = |fit - truth| / prior width, median over runs")
push!(summary_lines, "")
header = rpad("optimizer", 34) * rpad("runs", 6) * rpad("med excess", 12) * rpad("max excess", 12) * rpad("reached", 9) *
         rpad("below truth", 12) * rpad("med calls", 11) * rpad("med time/s", 11) * rpad("failed", 7) * "param err (" * join(plotted_params, ", ") * ")"
push!(summary_lines, header)
spec_summaries = Dict{String, Any}()
for (s, sp) in enumerate(specs)
    rs = [r for r in records if r.s == s]
    ex = [excess_ref(r) for r in rs]
    below = count(r -> isfinite(r.score) && r.score < ref_by_truth[r.i].truth_score - 1e-6, rs)
    reached = count(e -> isfinite(e) && e <= tol, ex)
    nfail = count(r -> !isfinite(r.score), rs)
    perr = [med([param_err(r, k) for r in rs]) for k in plotted_params]
    spec_summaries[sp.label] = (; n = length(rs), median_excess = med(ex), max_excess = isempty(fin(ex)) ? NaN : maximum(fin(ex)), reached, below_truth = below,
                                median_calls = med([Float64(coalesce(r.n_calls, NaN)) for r in rs]), median_time = med([r.time for r in rs]), failed = nfail, param_err = perr)
    sm = spec_summaries[sp.label]
    push!(summary_lines, rpad(sp.label, 34) * rpad(string(sm.n), 6) * rpad(string(round(sm.median_excess, sigdigits = 3)), 12) *
        rpad(string(round(sm.max_excess, sigdigits = 3)), 12) * rpad("$(reached)/$(sm.n)", 9) * rpad(string(below), 12) *
        rpad(string(round(sm.median_calls, digits = 0)), 11) * rpad(string(round(sm.median_time, digits = 1)), 11) * rpad(string(nfail), 7) *
        join(round.(perr, sigdigits = 2), "  "))
end
push!(summary_lines, "")
push!(summary_lines, "reference fits: median time $(round(med([r.time for r in ref_records]), digits = 1)) s, converged $(count(r -> r.converged === true, ref_records))/$(length(ref_records))" *
                     ", reference score minus truth score: " * join(round.([r.score - r.truth_score for r in ref_records], sigdigits = 3), ", ") *
                     " (negative: the reference did not reach the truth's own value)")
push!(summary_lines, "total wall time: $(round(time() - run_start_time, digits = 1)) s")
summary_text = join(summary_lines, "\n")
println(summary_text)
write(out_prefix * "_summary.txt", summary_text * "\n")

### PLOTS ###

colors = Makie.to_colormap(:tab10)
speccolor(s) = colors[mod1(s, length(colors))]
Δ_label = rich("Δχ", superscript("2"), subscript("ref"))   # Δχ²_ref = 2 (reference score - score of the run)
ysc = Makie.Symlog10(0.1)   # linear within ±0.1, logarithmic outside: 0 (the reference) and negative values (better than the reference) fit on the axis

# Trace of a run in Δχ²_ref = 2 × (excess above the reference) against evaluations. Global optimizers: the trace of the best
# objective value so far per evaluation. Local fits: Optim's per-iteration trace (in BAT's transformed space, shifted to
# end at the run's final objective value), with the iterations spread evenly over the run's likelihood calls.
function trace_xy(r)
    ref_obj = -ref_by_truth[r.i].score
    if r.trace_evals === missing
        tr = collect(Float64, r.trace_best)
        isempty(tr) && return Float64[], Float64[]
        tr = tr .- tr[end] .- r.score
        n = length(tr)
        calls = r.n_calls === missing ? n : r.n_calls
        x = n == 1 ? [Float64(calls)] : [1 + (calls - 1) * (j - 1) / (n - 1) for j in 1:n]
        return x, 2 .* (tr .- ref_obj)
    end
    Float64.(r.trace_evals), 2 .* (r.trace_best .- ref_obj)
end

ncols = min(n_truths, 3)
nrows = cld(n_truths, ncols)
fig = Figure(size = (550 * ncols, 420 * nrows + 60))
for i in 1:n_truths
    ax = Axis(fig[(i - 1) ÷ ncols + 1, (i - 1) % ncols + 1], xlabel = "objective evaluations", ylabel = Δ_label, title = "truth $i",
              xscale = log10, yscale = ysc)
    for r in records
        r.i == i || continue
        x, y = trace_xy(r)
        isempty(x) && continue
        lines!(ax, x, y; color = speccolor(r.s), label = labels[r.s])
    end
    hlines!(ax, [0.0]; color = :black, linestyle = :dash)                                                        # the reference fit
    hlines!(ax, [2 * (ref_by_truth[i].score - ref_by_truth[i].truth_score)]; color = :gray, linestyle = :dot)   # the injected truth's own value
end
Legend(fig[nrows + 1, 1:ncols], [LineElement(color = speccolor(s)) for s in 1:n_specs], labels, orientation = :horizontal, tellwidth = false)
save(out_prefix * "_traces.pdf", fig)

fig2 = Figure(size = (1300, 480))
ax1 = Axis(fig2[1, 1], ylabel = Δ_label, yscale = ysc,
           xticks = (1:n_specs, labels), xticklabelrotation = π / 6, title = "distance to the reference")
for (s, sp) in enumerate(specs)
    rs = [r for r in records if r.s == s]
    ex = [2 * excess_ref(r) for r in rs]
    ok = isfinite.(ex)
    any(ok) && scatter!(ax1, fill(s, count(ok)) .+ 0.12 .* (rand(count(ok)) .- 0.5), ex[ok]; color = speccolor(s))
end
hlines!(ax1, [0.0]; color = :black, linestyle = :dash)   # the reference fit
hlines!(ax1, [2 * tol]; color = :gray, linestyle = :dot)   # 'reached' threshold
ax2 = Axis(fig2[1, 2], xlabel = "wall time / s", ylabel = Δ_label, xscale = log10, yscale = ysc, title = "time vs distance")
for (s, sp) in enumerate(specs)
    rs = [r for r in records if r.s == s && isfinite(excess_ref(r))]
    isempty(rs) || scatter!(ax2, [r.time for r in rs], [2 * excess_ref(r) for r in rs]; color = speccolor(s), label = labels[s])
end
hlines!(ax2, [0.0]; color = :black, linestyle = :dash)
Legend(fig2[1, 3], ax2, framevisible = false)
save(out_prefix * "_summary.pdf", fig2)

### SAVE RESULTS ###

results_dict = Dict{String, Any}(
    "settings" => Dict("setup" => S.name, "hypothesis" => hyp_name, "tested" => collect(string.(tested_names)), "optimizers" => labels,
                       "reference" => ref_spec.label, "budget" => budget, "n_truths" => n_truths, "nrepeats" => nrepeats,
                       "objective" => objective, "data_mode" => data_mode, "fluctuate_nuisances" => args["fluctuate-nuisances"],
                       "vary_truth_params" => collect(string.(vary_truth_params)), "seed" => args["seed"], "tol" => tol,
                       "workers" => n_workers),
    "truths"          => truths,
    "reference_fits"  => ref_records,
    "records"         => records,
    "excess_ref"      => [excess_ref(r) for r in records],
    "excess_truth"    => [excess_truth(r) for r in records],
    "summary"         => spec_summaries,
    "summary_text"    => summary_text,
)
FileIO.save(out_prefix * "_results.jld2", results_dict)
println("saved results to $results_dir")
