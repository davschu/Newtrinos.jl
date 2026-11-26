# Sensitivity studies on one hypothesis of a setup file, with data injected at the setup's parameters (Asimov, or toys for --method toys):
#   --method fisher  expected uncertainties and the a priori correlation matrix from the Hessian at the truth (no fits); --truth selects the truth(s):
#                    nominal | grid | range:<intervals> | point:<name=value,...> | bf:<result.jld2>; one overview figure (sensitivity_plots.jl)
#   --method scan    one-at-a-time scans over the full prior range of every free parameter: others fixed and (unless --no-profile) others re-fitted
#   --method toys    distribution of the test statistic q = 2 (log post. of the full fit - log post. of the null fit) for toys injected at the null
# Usage: julia src/analysis/sensitivity.jl --setup <file> --hypothesis <name> --method fisher|scan|toys [--workers N]  (--help)
# Options, output files and plots: docs/src/README_optimizer_tests.md (section "sensitivity.jl").

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
using LinearAlgebra
using Dates
using ProgressMeter
using ValueShapes
using Sobol
import Optim
import ForwardDiff

include(joinpath(@__DIR__, "optimizer_tests_common.jl"))
include(joinpath(@__DIR__, "sensitivity_plots.jl"))

### CLI ###

function parse_command_line()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--setup"
        help = "Setup file defining the physics, experiments, priors and hypotheses (it defines setup(; ftype) -- see nsi/setups/nsi_deepcore_3y.jl for the format). The truth is the setup's parameters (with the pinned parameters of the hypothesis forced)"
        arg_type = String
        required = true

        "--hypothesis"
        help = "Name of the hypothesis (from the setup file); its free parameters are the ones studied, its tested parameters are the parameters of interest. Default: the first one of the setup"
        arg_type = String
        default = ""

        "--method"
        help = "fisher (Hessian at the truth, no fits), scan (one-at-a-time scans over the full prior range) or toys (test-statistic distribution under the null)"
        arg_type = String
        default = "fisher"

        "--params"
        help = "scan: names of the free parameters to scan. One name: a 1D scan of it (the truth is a grid point). Two or more names: ONE joint scan of all of them (2: a 2D grid, more: a Sobol sample; the others are fixed or profiled at every point; the nominal truth is not added as a point). Default (no --params): all free scalar parameters of the hypothesis, one at a time"
        arg_type = String
        nargs = '*'
        default = String[]

        "--n-points"
        help = "scan: grid points over the prior range (infinite prior ranges: quantiles 0.001 to 0.999). With one parameter (or none: each one separately): this many points per parameter. With a joint scan (--params with 2 or more names): the TOTAL number of points: round(sqrt(n)) per axis for 2 parameters (so 64 gives 8 x 8), n Sobol points for more"
        arg_type = Int
        default = 11

        "--range"
        help = "scan with 2 or more --params: restrict the joint scan to these intervals, name=lo,hi;name2=lo,hi (parameters not listed keep their whole prior range; needed for a parameter with an infinite prior range)"
        arg_type = String
        default = ""

        "--no-profile"
        help = "scan: only the scans with all other parameters fixed (likelihood evaluations, no fits); default: also re-fit all other free parameters at every grid point"
        action = :store_true

        "--strategy"
        help = "scan/toys: how every fit is done: ONE strategy, its stages (optimizer specs as in roundtrips.jl, e.g. lbfgsb:iterations=2000) given as separate arguments or in one quoted string, e.g. \"mle optim-lbfgs:iterations=2000\". Unlike roundtrips.jl several arguments are the stages of one strategy, not several strategies to compare. A last token @N (scan only) makes N independent attempts per fit and keeps the best, which only matters for a random first stage (random, mle, a global search). Local optimizers only make sense here"
        arg_type = String
        nargs = '+'
        default = ["lbfgsb"]

        "--n-truths"
        help = "fisher with --truth grid or range:...: number of truths, spread over the tested parameters (1 or 2 tested: grid, more: Sobol)"
        arg_type = Int
        default = 1

        "--truth"
        help = "fisher: where the data are injected. nominal (the setup's parameters, default); grid (truth_grid over the whole prior ranges of the tested parameters, --n-truths points; use it if the Fisher matrix is singular at the nominal point, e.g. a phase whose magnitude is zero); range:lo,hi or range:name=lo,hi;name2=lo,hi (the same grid restricted to these intervals, e.g. a suspicious region; with one tested parameter lo,hi is enough); point:name=x,name2=y (ONE truth: these parameters, tested or nuisances, set on top of the nominal ones); bf:file.jld2 (ONE truth: the best fit of a NewtrinosResult stored under the key \"result\", as written by save_result; parameters it does not contain stay nominal)"
        arg_type = String
        default = "nominal"

        "--drop-boundaries"
        help = "fisher with --truth grid or range:...: keep the truths off the interval bounds (cell-centred grid values lo + (i-1/2)/n * width instead of values from lo to hi), so no tested parameter sits at a bound or has no effect (e.g. a phase at zero magnitude). All parameters stay in the matrix"
        action = :store_true

        "--top-n"
        help = "fisher: number of nuisances drawn in the overview figure (the ones with the largest |impact| on a tested parameter at any truth); all of them are in the summary text and the results file"
        arg_type = Int
        default = 25

        "--information"
        help = "fisher: expected (Fisher information of the Poisson model at the truth, no data, always positive semi-definite; default) or observed (Hessian of the log likelihood of the rounded Asimov counts, which can be indefinite where the data are insensitive)"
        arg_type = String
        default = "expected"

        "--n-toys"
        help = "toys: number of toy data sets injected at the null of the hypothesis"
        arg_type = Int
        default = 100

        "--seed"
        help = "Base RNG seed (toy data, optimizer runs)"
        arg_type = Int
        default = 1234

        "--name"
        help = "Name prefix for the output files"
        arg_type = String
        default = "sensitivity"

        "--results-dir"
        help = "Output folder (default: results/<setup name>/sensitivity, relative to the working directory)"
        arg_type = String
        default = ""

        "--workers"
        help = "Number of distributed single-threaded worker processes (default: 1, no distributed: tasks run one after the other in this process)"
        arg_type = Int
        default = 1
    end

    return parse_args(s)
end

args = parse_command_line()
run_start_time = time()
method = lowercase(args["method"])
method in ("fisher", "scan", "toys") || error("--method must be fisher, scan or toys, got '$method'")
first(split(args["truth"], ':'; limit = 2)) in ("nominal", "grid", "range", "point", "bf") ||
    error("--truth must be nominal, grid, range:<spec>, point:<spec> or bf:<file.jld2>, got '$(args["truth"])'")
args["top-n"] >= 1 || error("--top-n must be at least 1")
args["information"] in ("expected", "observed") || error("--information must be expected or observed, got '$(args["information"])'")
n_workers = args["workers"]
use_distributed = n_workers > 1

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
        using LinearAlgebra
        using Sobol
        import Optim
        import ForwardDiff
    end
    # SEPARATE block: `using Accessors` must be loaded before the helper file defines methods using it
    @everywhere begin
        using Accessors
        include(joinpath($(@__DIR__), "optimizer_tests_common.jl"))
    end
end

### SETUP ###

include(abspath(args["setup"]))
S = finalize_setup(setup(; ftype = Float64))
hyp_name = isempty(args["hypothesis"]) ? first(keys(S.hypotheses)) : args["hypothesis"]
hyp = hypothesis_info(S, hyp_name)
experiments = S.experiments
prior = hyp.prior
tested = hyp.tested
free_keys = [k for k in keys(prior) if !_is_fixed(prior[k])]
nominal = _fix_conditioned(S.params, prior)   # nominal parameters with the prior's pinned parameters forced to their pinned values
strategy = String[t for s in args["strategy"] for t in split(s)]   # one strategy: its stages, as separate arguments or in one quoted string
nseeds = 1   # a last token @N: N independent attempts of the strategy per fit, the best is kept (scan only)
if !isempty(strategy) && (m = match(r"^@(\d+)$", strategy[end])) !== nothing
    nseeds = parse(Int, m.captures[1])
    pop!(strategy)
end
nseeds >= 1 || error("--strategy: @N needs N >= 1")
(nseeds == 1 || method == "scan") || error("--strategy: @N (several attempts) is only supported for --method scan")
normalize_strategy(strategy)   # check the strategy before anything expensive
results_dir = isempty(args["results-dir"]) ? joinpath("results", S.name, "sensitivity") : args["results-dir"]
mkpath(results_dir)
run_date = Dates.format(Dates.today(), "ddmmyyyy")
out_prefix = joinpath(results_dir, "$(args["name"])_$(hyp_name)_$(method)")
# a figure must never cost the results of a long run: its failure (e.g. the results folder was moved while the run was going) is reported, the summary and the results file are still written
function save_figure(path, fig)
    try
        mkpath(dirname(path))
        save(path, fig)
    catch e
        @error "the figure $path could not be saved (the summary and the results file are still written)" exception = (e, catch_backtrace())
    end
end
println("setup $(S.name), hypothesis '$hyp_name' (tested: $(join(tested, ", "))), method $method, $(length(free_keys)) free parameters")

# width of a prior (finite range, else 6 standard deviations), per scalar parameter or vector component
prior_width(d::Distributions.UnivariateDistribution) = isfinite(maximum(d) - minimum(d)) ? maximum(d) - minimum(d) : 6 * std(d)
prior_width(d, i) = d isa Distributions.UnivariateDistribution ? prior_width(d) : 6 * sqrt(var(d)[i])
function name_width(n)   # width of the prior of a name of `correlation` (`k` or `k[i]`)
    m = match(r"^(.*)\[(\d+)\]$", string(n))
    isnothing(m) ? prior_width(prior[Symbol(n)]) : prior_width(prior[Symbol(m.captures[1])], parse(Int, m.captures[2]))
end

summary_lines = String[]
say(s = "") = (println(s); push!(summary_lines, s))
results = Dict{String, Any}("settings" => Dict("setup" => S.name, "hypothesis" => hyp_name, "method" => method, "tested" => collect(string.(tested)),
                                                "seed" => args["seed"], "strategy" => strategy, "nseeds" => nseeds, "information" => args["information"], "drop_boundaries" => args["drop-boundaries"]))

### FISHER ###

if method == "fisher"
    truth_sel = truth_from_spec(args["truth"], prior, S.params, tested; n_points = args["n-truths"], interior = args["drop-boundaries"])
    truths = truth_sel.truths
    results["settings"]["truth"] = args["truth"]
    results["settings"]["truth_description"] = truth_sel.description
    ctx = (; experiments, prior, truths, bound_tol = 0.0, information = args["information"])
    use_distributed && @everywhere ctx = $ctx

    # inverse of the Hessian of the log posterior at the truth = expected covariance: of the Fisher information of the model (expected),
    # or of the Hessian of the log likelihood of the (rounded) Asimov counts (observed)
    @everywhere function fisher_task(i)
        c = ctx
        t0 = time()
        out = try
            if c.information == "expected"
                Newtrinos.correlation(c.experiments, c.prior, c.truths[i]; bound_tol = c.bound_tol)
            else
                L = Newtrinos.generate_likelihood(c.experiments, Newtrinos.generate_asimov_data(c.experiments, c.truths[i]))
                Newtrinos.correlation(L, c.prior, c.truths[i]; bound_tol = c.bound_tol)
            end
        catch e
            e isa Union{LinearAlgebra.SingularException, LinearAlgebra.PosDefException} || rethrow(e)
            nothing
        end
        println("truth $i of $(length(c.truths)): ", isnothing(out) ? "no result (not positive definite)" : "done", " in ", round(time() - t0, digits = 1), " s")
        out
    end
    cs = use_distributed ? pmap(fisher_task, 1:length(truths)) : map(fisher_task, 1:length(truths))

    say("=== sensitivity (fisher) summary ===")
    say("setup: $(S.name)   hypothesis: $hyp_name (tested: $(join(tested, ", ")))   information: $(args["information"])")
    say("truth: $(truth_sel.description)")
    say(args["information"] == "expected" ?
        "Fisher information of the Poisson model at the truth (J' diag(1/mu) J from the Jacobian of the expected counts, no data) plus the curvature of the priors; sigma_marg = expected uncertainty with all other parameters free, sigma_cond = with all others fixed;" :
        "Hessian of the log posterior at the truth on (rounded) Asimov data; sigma_marg = expected uncertainty with all other parameters free, sigma_cond = with all others fixed;")
    say("width = prior width (6 sigma for unbounded priors); the impact of nuisance j on tested parameter i = shift of i for a +1 sigma_marg move of j (first order)")
    fisher_results = []
    for (i, c) in enumerate(cs)
        say()
        say("--- truth $i: " * join(("$k = $(round.(truths[i][k], sigdigits = 4))" for k in tested), ", "))
        if isnothing(c)
            say("the information matrix is not positive definite at this truth (see the warning above for the smallest eigenvalues and the parameters in them: a parameter without effect or prior curvature, e.g. a phase whose magnitude is zero)" *
                (args["information"] == "observed" ? ", or the rounded Asimov counts make the Hessian indefinite where the data are insensitive (try --information expected)" : "") * ". No result for this truth")
            push!(fisher_results, nothing)
            continue
        end
        names = c.names
        F = inv(c.cov)
        σ_marg = sqrt.(diag(c.cov))
        σ_cond = 1 ./ sqrt.(diag(F))
        width = [name_width(n) for n in names]
        say(args["information"] == "expected" ?
            "offset of the truth from the stationary point of the log posterior: $(round(c.offset, sigdigits = 3)) sigma (the pull of the priors only); dropped at a prior bound: $(isempty(c.dropped) ? "none" : join(c.dropped, ", "))" :
            "offset of the truth from the stationary point of the log posterior: $(round(c.offset, sigdigits = 3)) sigma (should be << 1; Asimov counts are rounded); dropped at a prior bound: $(isempty(c.dropped) ? "none" : join(c.dropped, ", "))")
        w = maximum(length ∘ string, names) + 2
        say(rpad("parameter", w) * rpad("sigma_marg", 13) * rpad("sigma_cond", 13) * rpad("width", 13) * "sigma_marg / width")
        for (j, n) in enumerate(names)
            say(rpad(string(n), w) * rpad(string(round(σ_marg[j], sigdigits = 3)), 13) * rpad(string(round(σ_cond[j], sigdigits = 3)), 13) *
                rpad(string(round(width[j], sigdigits = 3)), 13) * string(round(σ_marg[j] / width[j], sigdigits = 3)))
        end
        # impact of the nuisance parameters on the tested parameters
        ti = [j for (j, n) in enumerate(names) if n in tested]
        impact = [c.cov[a, j] / σ_marg[j] for a in ti, j in 1:length(names)]
        σ_fixed = [begin   # sigma of the tested parameter a with parameter j held fixed
                       keep = [m for m in 1:length(names) if m != j]
                       a == j ? NaN : sqrt(inv(F[keep, keep])[findfirst(==(a), keep), findfirst(==(a), keep)])
                   end for a in ti, j in 1:length(names)]
        for (r, a) in enumerate(ti)
            say()
            say("impact on $(names[a]) (sigma_marg = $(round(σ_marg[a], sigdigits = 3))), nuisances ranked by |shift|:")
            for j in sortperm(abs.(impact[r, :]); rev = true)
                j == a && continue
                say("  " * rpad(string(names[j]), w) * "shift " * rpad(string(round(impact[r, j], sigdigits = 3)), 12) *
                    "sigma with $(names[j]) fixed: " * string(round(σ_fixed[r, j], sigdigits = 3)))
            end
        end
        push!(fisher_results, (; names, cov = c.cov, corr = c.corr, σ_marg, σ_cond, width, impact, σ_fixed, tested_idx = ti, offset = c.offset, dropped = c.dropped))
    end
    results["truths"] = truths
    results["fisher"] = fisher_results

    # the overview figure (sensitivity_plots.jl); the summary and the results file are written below even if the plot fails
    if all(isnothing, fisher_results)
        say("\nno truth gave a result: no plots")
    else
        try
            save("$(out_prefix)_overview_$(run_date).pdf", static_overview(sensitivity_data(results); top_n = args["top-n"]))
        catch e
            @error "the overview plot failed (the summary and the results file are still written)" exception = (e, catch_backtrace())
        end
    end
end

### SCAN ###

if method == "scan"
    scan_keys = isempty(args["params"]) ? [k for k in free_keys if !(nominal[k] isa AbstractVector)] : Symbol.(args["params"])
    for k in scan_keys
        k in free_keys || error("--params: $k is not a free parameter of hypothesis '$hyp_name'")
        nominal[k] isa AbstractVector && error("--params: $k is a vector parameter, only scalar parameters can be scanned")
    end
    length(unique(scan_keys)) == length(scan_keys) || error("--params: a parameter is given twice")
    joint = length(args["params"]) >= 2   # several named parameters: one joint scan (2D grid or Sobol sample) instead of one scan per parameter
    n_points = args["n-points"]
    if joint
        n_points >= (length(scan_keys) == 2 ? 4 : 2) || error("--n-points is the total number of points of a joint scan: at least 4 for 2 parameters (2 per axis), got $n_points")
    else
        isempty(args["range"]) || error("--range is only used with 2 or more --params (a joint scan)")
        n_points >= 2 || error("--n-points must be at least 2 (a grid from the lower to the upper end of the prior range), got $n_points")
    end
    data = Newtrinos.generate_asimov_data(experiments, nominal)
    function grid(d)
        lo, hi = minimum(d), maximum(d)
        (isfinite(lo) && isfinite(hi)) || ((lo, hi) = (quantile(d, 0.001), quantile(d, 0.999)))
        collect(range(lo, hi, length = n_points))
    end
    if joint
        listed = isempty(args["range"]) ? Symbol[] : [Symbol(strip(first(split(item, '=')))) for item in split(args["range"], ';'; keepempty = false)]
        for k in scan_keys
            k in listed || isfinite(maximum(prior[k]) - minimum(prior[k])) || error("--params: $k has an infinite prior range, give an interval with --range $k=lo,hi")
        end
        jtruths, _, jcoords = if isempty(args["range"])
            truth_grid(prior, nominal, scan_keys; n_points)
        else
            try
                r = truth_from_spec("range:" * args["range"], prior, nominal, scan_keys; n_points)
                (r.truths, r.scan_keys, r.grid_coords)
            catch e
                e isa ErrorException ? error(replace(e.msg, "--truth range" => "--range")) : rethrow()
            end
        end
        tasks = [(; i, pin = NamedTuple{Tuple(scan_keys)}(Tuple(Float64(t[k]) for k in scan_keys))) for (i, t) in enumerate(jtruths)]
    else
        grids = Dict(k => sort(unique(vcat(grid(prior[k]), Float64(nominal[k])))) for k in scan_keys)   # the truth itself is a grid point (Δχ² = 0)
        tasks = [(; i, pin = NamedTuple{(k,)}((v,))) for (i, (k, v)) in enumerate(Iterators.flatten([(k, v) for v in grids[k]] for k in scan_keys))]
    end
    profile = !args["no-profile"]
    ctx = (; experiments, data, prior, nominal, strategy, nseeds, profile, seed = args["seed"])
    use_distributed && @everywhere ctx = $ctx

    @everywhere function scan_task(t)
        c = ctx
        L = Newtrinos.generate_likelihood(c.experiments, c.data)
        p = merge(c.nominal, t.pin)   # t.pin: the scanned parameter(s) of this point
        llh_fixed = logdensityof(L, p)
        c.profile || return (; llh_fixed, llh = NaN, params = p)
        pinned = Newtrinos.condition(c.prior, Dict{Symbol, Float64}(k => Float64(v) for (k, v) in pairs(t.pin)), c.nominal)
        r = find_mle_advanced(L, pinned, p; strategy = c.strategy, nseeds = c.nseeds, rng = Random.Xoshiro(c.seed + t.i), objective = "posterior")
        (; llh_fixed, llh = r.llh, params = r.params)
    end

    L0 = Newtrinos.generate_likelihood(experiments, data)
    llh_truth = logdensityof(L0, nominal)
    println((joint ? "joint scan of $(join(scan_keys, ", ")): $(length(tasks)) points" : "$(length(scan_keys)) parameters x $(n_points) points") * (profile ? ", each profiled" : "") * " (reference llh at the truth: $(round(llh_truth, digits = 3))) ...")
    recs = use_distributed ? @showprogress(pmap(scan_task, tasks)) : @showprogress(map(scan_task, tasks))

    # Δχ² = 2 (llh(truth) - llh(point)); the 1σ crossing is interpolated linearly on the grid, going out from the truth
    function crossing(xs, ys, level)
        for i in 2:length(xs)
            if ys[i] >= level
                return xs[i-1] + (level - ys[i-1]) * (xs[i] - xs[i-1]) / (ys[i] - ys[i-1])
            end
        end
        NaN
    end
    function interval(xs, ys, x0)
        i0 = argmin(abs.(xs .- x0))
        (crossing(xs[i0:-1:1] .- x0, ys[i0:-1:1], 1.0), crossing(xs[i0:end] .- x0, ys[i0:end], 1.0))   # offsets from the truth (lower side negative)
    end

    if joint
    ndof = length(scan_keys)
    coords = [Float64[t.pin[k] for k in scan_keys] for t in tasks]
    dj_fixed = 2 .* (llh_truth .- [r.llh_fixed for r in recs])
    dj_prof = profile ? 2 .* (llh_truth .- [r.llh for r in recs]) : fill(NaN, length(recs))
    levels = [0.68, 0.9, 0.95]
    thr = quantile.(Chisq(ndof), levels)   # Δχ² of the 68 / 90 / 95 % regions for ndof scanned parameters
    dj = profile ? dj_prof : dj_fixed
    ok = findall(isfinite, dj)
    others = [a for a in tested if !(a in scan_keys)]
    jshifts = profile ? [maximum(abs(Float64(r.params[a]) - Float64(nominal[a])) / prior_width(prior[a]) for r in recs) for a in others] : Float64[]

    say("=== sensitivity (scan, joint) summary ===")
    say("setup: $(S.name)   hypothesis: $hyp_name (tested: $(join(tested, ", ")))   joint scan of $(join(scan_keys, ", ")): $(length(tasks)) points ($(ndof == 2 ? "2D grid" : "Sobol sample"), --n-points $n_points is the total), Asimov data at the setup's parameters")
    say("Δχ² = 2 (llh(truth) - llh(point)); 'fixed': all other parameters at the truth, 'profiled': all other free parameters re-fitted ($(join(strategy, " "))$(nseeds > 1 ? " @$nseeds" : ""))")
    say("regions: Δχ² below the Chisq($ndof) quantiles " * join(("$(round(Int, 100l)) % = $(round(t, sigdigits = 3))" for (l, t) in zip(levels, thr)), ", ") * " (the 1σ level of one parameter, Δχ² = 1, does not apply)")
    say()
    say(rpad("", 12) * rpad("max Δχ²", 12) * rpad("min Δχ²", 12) * join((rpad("n < $(round(t, sigdigits = 3))", 12) for t in thr), ""))
    for (name, d) in (profile ? [("fixed", dj_fixed), ("profiled", dj_prof)] : [("fixed", dj_fixed)])
        f = filter(isfinite, d)
        say(rpad(name, 12) * rpad(isempty(f) ? "NaN" : string(round(maximum(f), sigdigits = 3)), 12) * rpad(isempty(f) ? "NaN" : string(round(minimum(f), sigdigits = 3)), 12) *
            join((rpad(string(count(<(t), f)), 12) for t in thr), "") * "  (of $(length(f)) finite points)")
    end
    if !isempty(ok)
        ib = ok[argmin(dj[ok])]
        say("lowest Δχ² ($(profile ? "profiled" : "fixed")): $(round(dj[ib], sigdigits = 3)) at " * join(("$k = $(round(coords[ib][j], sigdigits = 4))" for (j, k) in enumerate(scan_keys)), ", "))
    end
    profile && !isempty(others) && say("shift = largest change of a tested parameter in the profiled fits, in units of its prior width: " * join(("$a: $(round(s, sigdigits = 2))" for (a, s) in zip(others, jshifts)), ", "))
    results["scan"] = Dict{Symbol, Any}(:joint => (; keys = scan_keys, coords, delta_chi2_fixed = dj_fixed, delta_chi2_profiled = dj_prof, params = [r.params for r in recs], ndof, levels, thresholds = thr, shifts = jshifts))
    results["llh_truth"] = llh_truth
    results["settings"]["scan_params"] = collect(string.(scan_keys))
    results["settings"]["range"] = args["range"]

    cmax = max(maximum(filter(isfinite, vcat(dj_fixed, dj_prof)); init = 0.0), 1e-6)
    panels = profile ? [("others fixed", dj_fixed), ("others profiled", dj_prof)] : [("others fixed", dj_fixed)]
    if ndof == 2
        k1, k2 = scan_keys
        gx = sort(unique(c[1] for c in coords)); gy = sort(unique(c[2] for c in coords))
        fig = Figure(size = (480 * length(panels) + 120, 420))
        for (n, (ttl, d)) in enumerate(panels)
            Z = fill(NaN, length(gx), length(gy))
            for (c, v) in zip(coords, d)
                Z[searchsortedfirst(gx, c[1]), searchsortedfirst(gy, c[2])] = v
            end
            ax = Axis(fig[1, n]; title = ttl, xlabel = string(k1), ylabel = string(k2))
            heatmap!(ax, gx, gy, Z; colorrange = (0, cmax), colormap = :viridis)
            try
                contour!(ax, gx, gy, Z; levels = thr, color = :white, linewidth = 1.5)   # 68 / 90 / 95 % regions
            catch e
                @warn "contours not drawn" exception = e
            end
            scatter!(ax, [Float64(nominal[k1])], [Float64(nominal[k2])]; color = :red, marker = :xcross, markersize = 12)
        end
        Colorbar(fig[1, length(panels) + 1]; colorrange = (0, cmax), colormap = :viridis, label = "Δχ²")
    else
        prs = [(a, b) for a in 1:ndof for b in a+1:ndof]
        ncols = min(length(prs), 3); nrows = cld(length(prs), ncols)
        fig = Figure(size = (400 * ncols + 120, 340 * nrows))
        for (n, (a, b)) in enumerate(prs)
            ax = Axis(fig[(n - 1) ÷ ncols + 1, (n - 1) % ncols + 1]; xlabel = string(scan_keys[a]), ylabel = string(scan_keys[b]))
            scatter!(ax, [c[a] for c in coords], [c[b] for c in coords]; color = dj, colorrange = (0, cmax), colormap = :viridis, markersize = 9)
            scatter!(ax, [Float64(nominal[scan_keys[a]])], [Float64(nominal[scan_keys[b]])]; color = :red, marker = :xcross, markersize = 12)
        end
        Colorbar(fig[1:nrows, ncols + 1]; colorrange = (0, cmax), colormap = :viridis, label = "Δχ² ($(profile ? "profiled" : "fixed"))")
    end
    save_figure("$(out_prefix)_$(run_date).pdf", fig)
    else
    say("=== sensitivity (scan) summary ===")
    say("setup: $(S.name)   hypothesis: $hyp_name (tested: $(join(tested, ", ")))   $(n_points) points per parameter over the prior range, Asimov data at the setup's parameters")
    say("Δχ² = 2 (llh(truth) - llh(point)); 'fixed': all other parameters at the truth, 'profiled': all other free parameters re-fitted ($(join(strategy, " "))$(nseeds > 1 ? " @$nseeds" : ""))")
    say("max Δχ² = Δχ² at the worst grid point (how strongly the data exclude the extremes of the prior range); 1σ = distance from the truth to Δχ² = 1 (NaN: not reached within the prior range)")
    say("shift = largest change of a tested parameter in the profiled fits, in units of its prior width")
    say()
    say(rpad("parameter", 26) * rpad("max Δχ² fixed", 15) * (profile ? rpad("max Δχ² profiled", 18) : "") * rpad("1σ fixed (-, +)", 26) * (profile ? rpad("1σ profiled (-, +)", 26) : "") *
        (profile ? "shift of " * join(tested, ", ") : ""))
    scan_results = Dict{Symbol, Any}()
    for k in scan_keys
        sel = [r for (t, r) in zip(tasks, recs) if first(keys(t.pin)) == k]
        xs = grids[k]
        d_fixed = 2 .* (llh_truth .- [r.llh_fixed for r in sel])
        d_prof = profile ? 2 .* (llh_truth .- [r.llh for r in sel]) : fill(NaN, length(xs))
        shifts = profile ? [maximum(abs(Float64(r.params[a]) - Float64(nominal[a])) / prior_width(prior[a]) for r in sel) for a in tested if a != k] : Float64[]
        scan_results[k] = (; x = xs, delta_chi2_fixed = d_fixed, delta_chi2_profiled = d_prof, params = [r.params for r in sel], shifts)
        fmt(t) = "(" * join(round.(t, sigdigits = 3), ", ") * ")"
        say(rpad(string(k), 26) * rpad(string(round(maximum(d_fixed), sigdigits = 3)), 15) * (profile ? rpad(string(round(maximum(d_prof), sigdigits = 3)), 18) : "") *
            rpad(fmt(interval(xs, d_fixed, nominal[k])), 26) * (profile ? rpad(fmt(interval(xs, d_prof, nominal[k])), 26) : "") *
            (profile ? join(round.(shifts, sigdigits = 2), ", ") : ""))
    end
    results["scan"] = scan_results
    results["llh_truth"] = llh_truth

    ncols = min(length(scan_keys), 4)
    nrows = cld(length(scan_keys), ncols)
    fig = Figure(size = (320 * ncols, 260 * nrows))
    for (n, k) in enumerate(scan_keys)
        ax = Axis(fig[(n - 1) ÷ ncols + 1, (n - 1) % ncols + 1], title = string(k), xlabel = string(k), ylabel = "Δχ²")
        lines!(ax, scan_results[k].x, scan_results[k].delta_chi2_fixed; color = :gray, linestyle = :dash, label = "others fixed")
        profile && lines!(ax, scan_results[k].x, scan_results[k].delta_chi2_profiled; color = :steelblue, label = "others profiled")
        hlines!(ax, [1.0]; color = :black, linestyle = :dot)
        vlines!(ax, [Float64(nominal[k])]; color = :red, linestyle = :dot)
    end
    Legend(fig[nrows + 1, 1:ncols], [LineElement(color = :gray, linestyle = :dash), LineElement(color = :steelblue)], ["others fixed", "others profiled"],
           orientation = :horizontal, tellwidth = false)
    save_figure("$(out_prefix)_$(run_date).pdf", fig)
    end   # joint / one at a time
end

### TOYS ###

if method == "toys"
    truth_null = _fix_conditioned(merge(S.params, S.hypotheses[hyp_name].null), prior)   # toys are injected at the null: tested parameters at their no-effect values
    ctx = (; experiments, prior, tested, null_prior = hyp.null_prior, truth = truth_null, strategy, seed = args["seed"])
    use_distributed && @everywhere ctx = $ctx

    @everywhere function toy_task(i)
        c = ctx
        rng = Random.Xoshiro(c.seed + 1_000_003 * i)
        Random.seed!(c.seed + 1_000_003 * i)   # generate_toy_data draws from the global RNG
        L = Newtrinos.generate_likelihood(c.experiments, Newtrinos.generate_toy_data(c.experiments, c.truth))
        null = find_mle_advanced(L, c.null_prior, _fix_conditioned(c.truth, c.null_prior); strategy = c.strategy, rng, objective = "posterior", null_prior = c.null_prior)
        # the full fit starts from the null fit (its tested parameters at the null values), so it cannot end below it
        full = find_mle_advanced(L, c.prior, null.params; strategy = c.strategy, rng, objective = "posterior", null_prior = c.null_prior)
        # the prior density of the tested parameters themselves (e.g. -log(width) of a Uniform) is not part of the statistic: removed from both log posteriors
        tested_prior(p) = sum(logpdf(c.prior[k], p[k]) for k in c.tested)
        q = 2 * ((full.log_post - tested_prior(full.params)) - (null.log_post - tested_prior(null.params)))
        (; i, q, full = full.params, null = null.params, log_post_full = full.log_post, log_post_null = null.log_post,
           converged_full = full.converged, converged_null = null.converged)
    end
    println("$(args["n-toys"]) toys injected at the null of '$hyp_name' ...")
    recs = use_distributed ? @showprogress(pmap(toy_task, 1:args["n-toys"])) : @showprogress(map(toy_task, 1:args["n-toys"]))
    q = [r.q for r in recs]
    qf = filter(isfinite, q)
    k = length(tested)
    levels = [0.68, 0.9, 0.95, 0.99]

    say("=== sensitivity (toys) summary ===")
    say("setup: $(S.name)   hypothesis: $hyp_name (tested: $(join(tested, ", ")))   $(args["n-toys"]) toys at the null, fits: $(join(strategy, " "))")
    say("q = 2 (log posterior of the full fit - of the null fit), without the prior density of the tested parameters themselves; Wilks: q ~ Chisq($k)  (the tested parameters' number)")
    say("finite q: $(length(qf)) / $(length(q)), negative q (failed fit): $(count(<(0), qf)); fits not converged: full $(count(r -> r.converged_full !== true, recs)), null $(count(r -> r.converged_null !== true, recs))")
    say(rpad("level", 8) * rpad("q quantile (toys)", 20) * "Chisq($k) quantile (Wilks)")
    for l in levels
        say(rpad(string(l), 8) * rpad(string(round(quantile(qf, l), sigdigits = 3)), 20) * string(round(quantile(Chisq(k), l), sigdigits = 3)))
    end
    results["toys"] = (; q, records = recs, truth = truth_null, ndof = k)

    fig = Figure(size = (600, 420))
    ax = Axis(fig[1, 1], xlabel = "q", ylabel = "density", title = "test statistic under the null ($(length(qf)) toys)")
    hist!(ax, qf; normalization = :pdf, bins = 20, color = (:steelblue, 0.6), label = "toys")
    xs = range(max(1e-3, minimum(qf)), max(maximum(qf), 1.0), length = 200)
    lines!(ax, xs, pdf.(Chisq(k), xs); color = :black, label = "Chisq($k) (Wilks)")
    axislegend(ax)
    save_figure("$(out_prefix)_$(run_date).pdf", fig)
end

### SAVE ###

say()
say("total wall time: $(round(time() - run_start_time, digits = 1)) s")
mkpath(results_dir)   # again: the folder may have been moved or deleted during a long run
write("$(out_prefix)_summary_$(run_date).txt", join(summary_lines, "\n") * "\n")
results["summary_text"] = join(summary_lines, "\n")
FileIO.save("$(out_prefix)_results_$(run_date).jld2", results)
println("saved results to $results_dir")
