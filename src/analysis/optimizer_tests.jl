# Compares optimizer STRATEGIES (chains of stages, e.g. "optim-pso lbfgsb @3") on injected truths of one hypothesis of a setup file:
# accuracy against the injected truth and a truth-seeded reference fit, cost (evaluations, time), convergence traces.
# Usage: julia src/analysis/optimizer_tests.jl --setup <file> --strategies lbfgsb "mh-de lbfgsb" ... [--workers N]  (--list-optimizers, --help)
# Options, output files and plots: docs/src/README_optimizer_tests.md and the docs page "Optimization tests".

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
using Dates
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

        "--strategies", "--optimizers"
        dest_name = "strategies"
        help = "Strategies to compare, one QUOTED argument each: stage specs separated by spaces, each \"name\" or \"name:key=val,key=val\" (see --list-optimizers; stage rules as in roundtrips.jl --strategy: seed/alias/global stages only first), optionally ending in @N = N independent attempts of that strategy, the best being kept, e.g. --strategies lbfgsb \"mh-de:evals=3000\" \"optim-pso:particles=20 lbfgsb @3\" \"smseed lbfgsb\". A single optimizer is a one-stage strategy. Time and calls are those of ALL attempts and stages"
        arg_type = String
        nargs = '+'
        default = ["optim-lbfgs", "optim-pso:particles=20", "mh-de"]

        "--nseeds"
        help = "Default number of attempts of a strategy that gives no @N of its own"
        arg_type = Int
        default = 1

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

specs    = [parse_strategy_arg(s; nseeds = args["nseeds"]) for s in args["strategies"]]   # (; stages, nseeds, label)
ref_spec = parse_optimizer_spec(args["reference"])
OPTIMIZERS[ref_spec.name].kind == :local || @warn "the reference is normally a local fit started at the truth; got a global optimizer ($(ref_spec.name))"
n_specs  = length(specs)
labels   = [s.label for s in specs]
length(unique(labels)) == n_specs || error("--strategies: every strategy must be different (got duplicates: $labels)")

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
run_date = Dates.format(Dates.today(), "ddmmyyyy")   # execution date, ends every file name
out_prefix = joinpath(results_dir, "$(name)_$(hyp_name)")

println("setup $(S.name), hypothesis '$hyp_name' (tested: $(join(tested_names, ", "))), $n_truths truths x $n_specs strategies x $nrepeats repeats, default budget $budget evaluations")

### TRUTH POINTS ###

rng = Random.Xoshiro(args["seed"])
truths = Vector{Any}(undef, n_truths)
for i in 1:n_truths
    t = controlled_truth(rng, p, prior, tested_names; extra_names = vary_truth_params)
    truths[i] = args["fluctuate-nuisances"] ? fluctuate_nuisances(rng, t, prior, fluctuated_names) : t
end

### TASKS ###

start = _fix_conditioned(p, prior)   # nominal parameters with the prior's pinned parameters forced to their pinned values
ctx = (; experiments, prior, prior_dist, null_prior = hyp.null_prior, start, truths, data_mode, objective, budget, specs, ref_spec, seed = args["seed"])
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
    st = c.specs[t.s]
    r = find_mle_advanced(L, c.prior, c.start; strategy = st.stages, nseeds = st.nseeds, budget = c.budget, rng = Random.Xoshiro(seed),
                          objective = c.objective, null_prior = c.null_prior)
    # time and calls of the record are those of ALL attempts and stages (the summary and the plots use these fields)
    merge(r, (; time = r.total_time, n_calls = r.total_calls, i = t.i, s = t.s, r = t.r, score = score(r), truth_score = NaN))
end

# a stage spec with tiny settings (for the warm-up; defined on the workers too)
@everywhere function tiny_stage(str)
    sp = parse_optimizer_spec(str)
    d = OPTIMIZERS[sp.name].defaults
    st = sp.settings
    :iterations in keys(d) && (st = merge(st, (; iterations = 2)))
    :iters in keys(d) && (st = merge(st, (; iters = 1)))
    :iters in keys(d) || :evals in keys(d) && (st = merge(st, (; evals = 40)))
    (; sp..., settings = st)
end

# Warm-up: run every optimizer once, tiny, on every worker (and here) so that compilation does not enter the timings
@everywhere function warmup_task(w)
    c = ctx
    truth = c.truths[1]
    L = Newtrinos.generate_likelihood(c.experiments, Newtrinos.generate_asimov_data(c.experiments, truth))
    for sp in c.specs
        find_mle_advanced(L, c.prior, c.start; strategy = [tiny_stage(s) for s in sp.stages], budget = 40, rng = Random.Xoshiro(1),
                          objective = c.objective, null_prior = c.null_prior)
    end
    run_optimizer(tiny_stage(c.ref_spec.label), L, c.prior, c.start; budget = 40, rng = Random.Xoshiro(1), objective = c.objective)
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
push!(summary_lines, "calls and time are summed over all stages and attempts (@N) of a strategy; the score is that of the best attempt")
push!(summary_lines, "below truth = run ended with a score lower than the injected truth's (a failure of the optimizer)")
push!(summary_lines, "reached = excess <= $tol")
push!(summary_lines, "param err = |fit - truth| / prior width, median over runs")
push!(summary_lines, "")
label_width = max(44, maximum(length, labels) + 2)   # the strategy column is as wide as the longest label
header = rpad("strategy", label_width) * rpad("runs", 6) * rpad("med excess", 12) * rpad("max excess", 12) * rpad("reached", 9) *
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
    push!(summary_lines, rpad(sp.label, label_width) * rpad(string(sm.n), 6) * rpad(string(round(sm.median_excess, sigdigits = 3)), 12) *
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
write(out_prefix * "_summary_$(run_date).txt", summary_text * "\n")

### PLOTS ###

colors = Makie.to_colormap(:tab10)
speccolor(s) = colors[mod1(s, length(colors))]
Δ_label = rich("Δχ", superscript("2"), subscript("ref"))   # Δχ²_ref = 2 (reference score - score of the run)
# Symlog scale by hand (linear within ±noise_band, logarithmic outside; 0 -- the reference -- and negative values, i.e. better than the
# reference, fit on the axis). The data are transformed before plotting and the axes are linear: with Makie's own nonlinear axis scales
# CairoMakie resamples lines into tiny segments and dashed/dotted line styles come out ragged.
noise_band = 0.1   # |Δχ²_ref| below this: only optimizer noise is left
ysym(y) = sign(y) * (abs(y) <= noise_band ? abs(y) / noise_band : 1 + log10(abs(y) / noise_band))

# Trace of a run in Δχ²_ref = 2 × (excess above the reference) against evaluations: the best objective value so far through
# all stages and attempts of the strategy (`trace_x`/`trace_y` of find_mle_advanced; the evaluation counter continues from stage to
# stage and from attempt to attempt). Global stages record every evaluation; local fits record per iteration (in BAT's transformed
# space, shifted to end at the stage's final objective value), spread evenly over the stage's likelihood calls (over its
# iteration count for `newtrinos`, which does not count calls).
# Returned per record: one (x, Δχ²_ref, stage) per attempt, every attempt drawn from its own start (same colour in the plot).
function trace_xy(r)
    ref_obj = -ref_by_truth[r.i].score
    map(r.attempt_traces) do t
        ok = isfinite.(t.y)
        (t.x[ok], 2 .* (t.y[ok] .- ref_obj), t.stage[ok])
    end
end

# One line per run of a stage of the strategy (the line style changes with the stage); a segment starts at the last point of the previous
# run so that the line is continuous.
const STAGE_STYLES = [:solid, :dashdot, :dot, :dash]
stagestyle(k) = STAGE_STYLES[mod1(k, length(STAGE_STYLES))]
function stage_segments(x, y, st)
    segs = Tuple{Int, Vector{Float64}, Vector{Float64}}[]
    a = 1
    for b in 1:length(st)
        if b == length(st) || st[b + 1] != st[a]   # end of a run of equal stage index (attempts after the first start a new run)
            idx = (a > 1 ? a - 1 : a):b
            push!(segs, (st[a], x[idx], y[idx]))
            a = b + 1
        end
    end
    segs
end

# y ticks of the symlog axis: 0, ±0.1 and the decades above, mirrored below zero when there are negative values
function symlog_ticks(ymin, ymax)
    kmax = max(0, ceil(Int, log10(max(ymax, 1.0))))
    pos = [10.0^k for k in -1:kmax]
    neg = ymin < -0.1 ? [10.0^k for k in -1:max(0, ceil(Int, log10(-ymin)))] : Float64[]
    vals = sort(vcat(-neg, 0.0, pos))
    lab(v) = v == 0 ? rich("0") : rich(v < 0 ? "-10" : "10", superscript(string(round(Int, log10(abs(v))))))
    (ysym.(vals), [lab(v) for v in vals])
end
# x ticks of a log10 axis (the data are transformed with log10 before plotting)
function log_ticks(xmin, xmax)
    ks = floor(Int, log10(xmin)):ceil(Int, log10(xmax))
    (Float64.(ks), [rich("10", superscript(string(k))) for k in ks])
end

ncols = min(n_truths, 3)
nrows = cld(n_truths, ncols)
trace_data = [trace_xy(r) for r in records]   # per record: (x, Δχ²_ref, stage) per attempt
flat_traces = reduce(vcat, trace_data; init = Tuple{Vector{Float64}, Vector{Float64}, Vector{Int}}[])
all_y = reduce(vcat, (d[2] for d in flat_traces); init = Float64[])
all_x = reduce(vcat, (d[1] for d in flat_traces); init = Float64[])
yticks_sym = symlog_ticks(isempty(all_y) ? 0.0 : minimum(all_y), isempty(all_y) ? 1.0 : maximum(all_y))
xlo, xhi = isempty(all_x) ? (1.0, 10.0) : (max(minimum(all_x), 1.0), max(maximum(all_x), 10.0))
xticks_log = log_ticks(xlo, xhi)
xspan = [log10(xlo), log10(xhi)]
fig = Figure(size = (550 * ncols, 420 * nrows + 40 + 22 * cld(n_specs, clamp(floor(Int, 550 * ncols / (6.5 * maximum(length, labels; init = 10) + 50)), 1, n_specs))))
for i in 1:n_truths
    ax = Axis(fig[(i - 1) ÷ ncols + 1, (i - 1) % ncols + 1], xlabel = "objective evaluations", ylabel = Δ_label, title = "truth $i",
              xticks = xticks_log, yticks = yticks_sym)
    hspan!(ax, -1, 1; color = (:grey, 0.3))                                                                      # |Δχ²| <= noise_band: only optimizer noise left in here
    for (ri, r) in enumerate(records)
        r.i == i || continue
        for (x, y, st) in trace_data[ri]   # one trace per attempt, all in the colour of the strategy
            isempty(x) && continue
            for (k, xs, ys) in stage_segments(log10.(max.(x, 1.0)), ysym.(y), st)
                lines!(ax, xs, ys; color = speccolor(r.s), linestyle = stagestyle(k))
            end
        end
    end
    lines!(ax, xspan, [0.0, 0.0]; color = :black, linestyle = :dash)                                             # the reference fit
    truth_level = ysym(2 * (ref_by_truth[i].score - ref_by_truth[i].truth_score))
    lines!(ax, xspan, [truth_level, truth_level]; color = :gray, linestyle = :dot)                               # the injected truth's own value
end
# legend: only the strategies (colours)
leg_elems = Any[LineElement(color = speccolor(s)) for s in 1:n_specs]
leg_cols = clamp(floor(Int, 550 * ncols / (6.5 * maximum(length, labels; init = 10) + 50)), 1, n_specs)   # as many columns as fit the figure width
Legend(fig[nrows + 1, 1:ncols], leg_elems, String[labels...], orientation = :horizontal, nbanks = cld(n_specs, leg_cols), tellwidth = false)
save(out_prefix * "_traces_$(run_date).pdf", fig)

excesses = [2 * excess_ref(r) for r in records]   # Δχ²_ref of every run
ex_fin = filter(isfinite, excesses)
yticks_sum = symlog_ticks(isempty(ex_fin) ? 0.0 : minimum(ex_fin), isempty(ex_fin) ? 1.0 : maximum(ex_fin))
fig2 = Figure(size = (1300, 480))
ax1 = Axis(fig2[1, 1], ylabel = Δ_label, xlabel = "strategy (see legend)", yticks = yticks_sum,
           xticks = (1:n_specs, string.(1:n_specs)), title = "distance to the reference")
for (s, sp) in enumerate(specs)
    rs = [r for r in records if r.s == s]
    ex = [2 * excess_ref(r) for r in rs]
    ok = isfinite.(ex)
    any(ok) && scatter!(ax1, fill(s, count(ok)) .+ 0.12 .* (rand(count(ok)) .- 0.5), ysym.(ex[ok]); color = speccolor(s))
end
hspan!(ax1, -1, 1; color = (:grey, 0.3))
hlines!(ax1, [0.0]; color = :black, linestyle = :dash)   # the reference fit
hlines!(ax1, [ysym(2 * tol)]; color = :gray, linestyle = :dot)   # 'reached' threshold
ax2 = Axis(fig2[1, 2], xlabel = "wall time / s", ylabel = Δ_label, yticks = yticks_sum, title = "time vs distance")
for (s, sp) in enumerate(specs)
    rs = [r for r in records if r.s == s && isfinite(excess_ref(r))]
    isempty(rs) || scatter!(ax2, [r.time for r in rs], [ysym(2 * excess_ref(r)) for r in rs]; color = speccolor(s), label = "$s: $(labels[s])")
end
hspan!(ax2, -1, 1; color = (:grey, 0.3))
hlines!(ax2, [0.0]; color = :black, linestyle = :dash)
Legend(fig2[1, 3], ax2, framevisible = false)
save(out_prefix * "_summary_$(run_date).pdf", fig2)

### SAVE RESULTS ###

results_dict = Dict{String, Any}(
    "settings" => Dict("setup" => S.name, "hypothesis" => hyp_name, "tested" => collect(string.(tested_names)), "strategies" => labels, "nseeds" => [s.nseeds for s in specs],
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
FileIO.save(out_prefix * "_results_$(run_date).jld2", results_dict)
println("saved results to $results_dir")
