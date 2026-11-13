# Shared optimizer helpers for optimizer_tests.jl. `include`d both in the main process
# and (via @everywhere) on distributed workers, so it must not reference any variable
# that isn't passed in as an argument.

import Random
import Metaheuristics
import BlackBoxOptim
import OptimizationLBFGSB

"""Build an Optim.jl algorithm object from the CLI-selected name and step size."""
function make_algorithm(optimizer_name, step_size)
    alphaguess = Optim.LineSearches.InitialStatic(alpha=step_size)
    if optimizer_name == "lbfgs"
        Optim.LBFGS(alphaguess=alphaguess)
    elseif optimizer_name == "gradientdescent"
        Optim.GradientDescent(alphaguess=alphaguess)
    elseif optimizer_name == "conjugategradient"
        Optim.ConjugateGradient(alphaguess=alphaguess)
    elseif optimizer_name == "neldermead"
        Optim.NelderMead()
    else
        error("Unknown --optimizer '$optimizer_name'. Choose from: lbfgs, gradientdescent, conjugategradient, neldermead")
    end
end

# ---- analysis setups ----------------------------------------------------------------------------
# A setup file (see nsi/setups/nsi_deepcore_3y.jl) defines `setup(; ftype)` returning
# `(; name, experiments, params, priors, hypotheses)`, with `hypotheses` a Dict name => (; prior, null).

"""True for a parameter prior that holds the parameter fixed (`ConstValueDist` / raw `Number`)."""
_is_fixed(d) = d isa ValueShapes.ConstValueDist || d isa Number

"""
    finalize_setup(S) -> S

Check that the NamedTuple `S` returned by a setup file has the required fields
(`name, experiments, params, priors, hypotheses`) and that every hypothesis is a
`(; prior, null)` whose tested parameters (`keys(null)`) are free in its `prior`.
"""
function finalize_setup(S)
    for f in (:name, :experiments, :params, :priors, :hypotheses)
        f in keys(S) || error("setup is missing the field `$f` (required: name, experiments, params, priors, hypotheses)")
    end
    for (hname, h) in S.hypotheses
        (:prior in keys(h) && :null in keys(h)) || error("hypothesis '$hname' must be a NamedTuple (; prior, null)")
        for k in keys(h.null)
            k in keys(h.prior) || error("hypothesis '$hname': tested parameter $k is not a parameter of its prior")
            _is_fixed(h.prior[k]) && error("hypothesis '$hname': tested parameter $k is fixed in its prior, it must be free")
        end
    end
    S
end

"""
    hypothesis_info(S, hypothesis_name) -> NamedTuple

Everything the scripts derive from one hypothesis of a setup `S`:
`(; name, prior, tested, null_prior, fluctuated)` with
- `tested`: the tested parameters, `keys(null)` (their order is the order of the scan axes),
- `null_prior`: `prior` with the tested parameters fixed at their null values (the model without
  the effect under test),
- `fluctuated`: all other free parameters of `prior`; the ones redrawn from their priors by
  `--fluctuate-nuisances`.
"""
function hypothesis_info(S, hypothesis_name)
    haskey(S.hypotheses, hypothesis_name) ||
        error("--hypothesis must be one of $(join(sort(collect(keys(S.hypotheses))), ", ")), got '$hypothesis_name'")
    h = S.hypotheses[hypothesis_name]
    tested = Tuple(keys(h.null))
    null_prior = Newtrinos.condition(h.prior, Dict{Symbol, Float64}(k => Float64(h.null[k]) for k in tested), S.params)
    fluctuated = Tuple(k for k in keys(h.prior) if !_is_fixed(h.prior[k]) && !(k in tested))
    (; name = hypothesis_name, prior = h.prior, tested, null_prior, fluctuated)
end

"""
    extract_bounds(prior, params)

Separate fixed (ConstValueDist / raw Number) from free parameters in `prior`,
returning a NamedTuple `(all_keys, is_free_mask, fixed_float, free_x_index, free_dists, lo_vec, hi_vec, x0_vec)`.
"""
function extract_bounds(prior, params)
    all_keys = tuple(keys(prior)...)
    n = length(all_keys)

    is_free_mask = Vector{Bool}(undef, n)
    fixed_float  = Vector{Float64}(undef, n)
    free_x_index = Vector{Int}(undef, n)
    free_dists   = Any[]
    lo_vec       = Float64[]
    hi_vec       = Float64[]
    x0_vec       = Float64[]
    n_free       = 0

    for (i, k) in enumerate(all_keys)
        d = prior[k]
        if d isa ValueShapes.ConstValueDist
            is_free_mask[i] = false
            fixed_float[i]  = Float64(d.value)
            free_x_index[i] = 0
        elseif d isa Number
            is_free_mask[i] = false
            fixed_float[i]  = Float64(d)
            free_x_index[i] = 0
        else
            n_free          += 1
            is_free_mask[i]  = true
            fixed_float[i]   = 0.0
            free_x_index[i]  = n_free
            push!(free_dists, d)
            lo = Float64(minimum(d))
            hi = Float64(maximum(d))
            push!(lo_vec, lo)
            push!(hi_vec, hi)
            θ0   = Float64(params[k])
            ε    = 1e-8
            lo_s = isfinite(lo) ? lo + ε : θ0 - 1e10
            hi_s = isfinite(hi) ? hi - ε : θ0 + 1e10
            push!(x0_vec, clamp(θ0, lo_s, hi_s))
        end
    end

    (; all_keys, is_free_mask, fixed_float, free_x_index, free_dists, lo_vec, hi_vec, x0_vec)
end

"""
    flat_prior(prior)

Return a `distprod`-wrapped `NamedTupleDist` (same as `prior` itself, e.g. from
`distprod(;cps[j]...)`) with the same keys as `prior`: fixed entries (`ConstValueDist` /
raw `Number`) unchanged, every free distribution replaced by `Uniform(minimum(d),
maximum(d))` over the same bounds. A `Uniform` prior has constant log-density on its
support (confirmed: `logpdf(Uniform(-0.1,0.1), x)` is the same for every `x` in range),
so it contributes nothing to the shape of `logdensityof(PosteriorMeasure(likelihood,
flat_prior(prior)), x)` beyond an additive constant -- maximizing that is equivalent to
maximizing the likelihood alone (MLE), while still giving `bat_findmode`/`OptimAlg`'s
`PriorToNormal` reparametrization the bounds it needs.

`prior` must itself be a `distprod`-wrapped `NamedTupleDist` (as `local_find_mle`
receives it), not a bare `NamedTuple` of distributions -- `PosteriorMeasure` requires an
`AbstractMeasure`/`Distribution`, which a plain `NamedTuple` is not; re-wrapping via
`distprod` is what makes the result usable there.
"""
function flat_prior(prior)
    kvs = [
        (d isa ValueShapes.ConstValueDist || d isa Number) ? (k => d) : (k => Uniform(Float64(minimum(d)), Float64(maximum(d))))
        for (k, d) in pairs(prior)
    ]
    distprod(; NamedTuple(kvs)...)
end

"""
    randomize_free_params(rng, params, prior)

Like `Newtrinos.randomize_params`, but also leaves parameters fixed via a raw `Number`
prior untouched (not just `ValueShapes.ConstValueDist`) — `Newtrinos.condition` fixes
parameters by setting the prior to a plain `Float64`, which `Newtrinos.randomize_params`
doesn't recognize as fixed and would otherwise try (and fail) to `rand` from.
"""
function randomize_free_params(rng, params, prior)
    p = deepcopy(params)
    for key in keys(p)
        d = prior[key]
        if !(d isa ValueShapes.ConstValueDist) && !(d isa Number)
            @reset p[key] = rand(rng, d)
        end
    end
    return p
end

"""
    controlled_truth(rng, nominal_params, prior, tested_names; extra_names=())

Build an injected-truth NamedTuple starting from `nominal_params` (e.g.
`Newtrinos.get_params(experiments)`), overriding: (1) every free (non-fixed) parameter
in `prior` whose key is in `tested_names` (the parameters under test of a hypothesis),
drawn randomly from its prior distribution, (2) every free parameter in `prior` whose key
is in `extra_names` (e.g. `(:θ₂₃, :Δm²₃₁)` to also vary those), drawn randomly from its prior,
and (3) every parameter `prior` holds **fixed** (`ConstValueDist` /
`Number`, e.g. oscillation parameters conditioned in a setup file), forced to that
exact fixed value -- `nominal_params` (`Newtrinos.get_params`'s own defaults) can
disagree with the specific constants a hypothesis's `prior` was conditioned to (verified
directly: `p[:δCP]=1.0` vs. the conditioned `0.0`), and `bat_findmode`'s `ExplicitInit`
strictly validates the init point against the prior's constants, raising
`ArgumentError("Cannot set constant value to a different value")` if they differ.

Every other parameter -- nuisances and the `extra_names` not given -- stays at its nominal
value. This isolates truth variation to the parameter(s) actually being tested by
`prior`'s hypothesis, avoiding prior-pull artifacts from randomly-fluctuated nuisance
parameters landing in low-density prior regions by chance (see `randomize_free_params`
for the same fixed-parameter check used here).
"""
function controlled_truth(rng, nominal_params, prior, tested_names; extra_names=())
    overrides = Dict{Symbol, Float64}()
    for k in (tested_names..., extra_names...)
        d = prior[k]
        _is_fixed(d) && continue
        overrides[k] = rand(rng, d)
    end
    _fix_conditioned(merge(nominal_params, NamedTuple(overrides)), prior)
end

"""
    fluctuate_nuisances(rng, truth_param, prior, nuisance_names)

Redraw every parameter in `nuisance_names` that's free in `prior` from its own prior
distribution, leaving `truth_param`'s other entries (NSI parameters, θ₂₃/Δm²₃₁,
conditioned-fixed params) untouched. Matches the dissertation's "recovering fluctuated
nuisance parameters" test (section 7.6.3): fluctuated nuisances, non-fluctuated
(Asimov) bin counts, NSI truth unchanged -- confirmed directly from the dissertation's
own wording ("pseudo-data are generated with fluctuated nuisance parameter values and
non-fluctuated bin counts").
"""
function fluctuate_nuisances(rng, truth_param, prior, nuisance_names)
    overrides = Dict{Symbol, Float64}()
    for k in nuisance_names
        d = prior[k]
        (d isa ValueShapes.ConstValueDist || d isa Number) && continue
        overrides[k] = rand(rng, d)
    end
    merge(truth_param, NamedTuple(overrides))
end

"""
    _fix_conditioned(params, prior)

Force every parameter `prior` holds fixed (`ConstValueDist` / raw `Number`) to that
exact fixed value in `params`, leaving every other key untouched. Factored out of
`controlled_truth` so `truth_grid` can reuse the same fixed-parameter-forcing logic
(see `controlled_truth`'s docstring for why this matters: `nominal_params` can disagree
with the specific constants a hypothesis's `prior` was conditioned to).
"""
function _fix_conditioned(params, prior)
    overrides = Dict{Symbol, Float64}()
    for k in keys(prior)
        d = prior[k]
        if d isa ValueShapes.ConstValueDist
            overrides[k] = Float64(d.value)
        elseif d isa Number
            overrides[k] = Float64(d)
        end
    end
    merge(params, NamedTuple(overrides))
end


"""
    truth_grid(prior, nominal_params, tested_names; n_points=25)

Generate a grid of `n_points` Asimov-truth NamedTuples scanning the free tested parameters of
a hypothesis (`tested_names`, e.g. `keys(null)` of a setup's hypothesis). All other parameters
stay at `nominal_params`, except the ones `prior` holds fixed, which are forced to their fixed
value (`_fix_conditioned`, same as `controlled_truth`). The grid only controls how the injected
truth is built; the subsequent fit is still free to vary every non-fixed parameter in `prior`.

The sampling strategy follows from the number of FREE tested parameters:
- 1 parameter: `n_points` evenly spaced values over its prior range (1D).
- 2 parameters: a 2D grid with `round(sqrt(n_points))` points per dimension (e.g. 5x5=25 for
  `n_points=25`), in the order of `tested_names` (first name = first coordinate).
- more than 2: a quasi-random Sobol sequence (`SobolSeq`) over all free dimensions
  simultaneously, far better space-filling than a factorial grid in that many dimensions.
- 0 free tested parameters is not a valid hypothesis and errors.

Returns `(truths, scan_keys, grid_coords)`: `truths` is the vector of truth NamedTuples,
`scan_keys` the tuple of scanned parameter names, `grid_coords` the scanned coordinates with
`truths` (`Vector{Float64}` for 1D, `Vector{Tuple{Float64,Float64}}` for 2D,
`Vector{Vector{Float64}}` -- one vector per sample, ordered to match `scan_keys` -- for
the >2D Sobol case).
"""
function truth_grid(prior, nominal_params, tested_names; n_points=25)
    free = [k for k in tested_names if !_is_fixed(prior[k])]

    if length(free) == 1
        k = free[1]
        d = prior[k]
        lo, hi = Float64(minimum(d)), Float64(maximum(d))
        grid_coords = collect(range(lo, hi, length=n_points))
        truths = [_fix_conditioned(merge(nominal_params, NamedTuple{(k,)}((v,))), prior) for v in grid_coords]
        return truths, (k,), grid_coords
    elseif length(free) == 2
        key1, key2 = free[1], free[2]
        d1, d2 = prior[key1], prior[key2]
        n_per_dim = Int(round(sqrt(n_points)))
        vals1 = collect(range(Float64(minimum(d1)), Float64(maximum(d1)), length=n_per_dim))
        vals2 = collect(range(Float64(minimum(d2)), Float64(maximum(d2)), length=n_per_dim))
        grid_coords = [(a, b) for a in vals1 for b in vals2]
        truths = [_fix_conditioned(merge(nominal_params, NamedTuple{(key1,key2)}((a,b))), prior)
                  for (a,b) in grid_coords]
        return truths, (key1, key2), grid_coords
    elseif length(free) > 2
        scan_keys = Tuple(free)
        los = Float64[Float64(minimum(prior[k])) for k in scan_keys]
        his = Float64[Float64(maximum(prior[k])) for k in scan_keys]
        seq = Sobol.SobolSeq(los, his)
        grid_coords = Vector{Vector{Float64}}(undef, n_points)
        truths = Vector{typeof(_fix_conditioned(nominal_params, prior))}(undef, n_points)
        for i in 1:n_points
            v = Sobol.next!(seq)
            grid_coords[i] = v
            truths[i] = _fix_conditioned(merge(nominal_params, NamedTuple{scan_keys}(Tuple(v))), prior)
        end
        return truths, scan_keys, grid_coords
    else
        error("truth_grid: none of the tested parameters $(collect(tested_names)) is free in the prior -- " *
              "a hypothesis needs at least one free tested parameter.")
    end
end

# ---- optimizer registry ----------------------------------------------------------------------------
# One table of named optimizers. A spec string "name:key=val,key=val" (see `parse_optimizer_spec`) selects an
# entry and overrides some of its defaults. Two kinds:
#   :global  gradient-free search in the box [lo, hi] of the free parameters, backend
#            `fn(counter, b, cfg, budget, rng) -> (; converged, n_iters, native)`; run by `global_find_mle`.
#   :local   gradient-based fit through `local_find_mle`, run by `run_optimizer`.
# To add an optimizer write a backend function and call `register_optimizer!` (end of this section).

const OPTIMIZERS = Dict{String, NamedTuple}()

"""
    register_optimizer!(name, kind, fn, defaults)

Add (or replace) an entry of the optimizer registry. `kind` is `:global` or `:local`; `defaults` is the NamedTuple
of settings the entry understands (a spec may override only these keys, the type of a default decides how a
command-line value is parsed; a default of `nothing` means "unset, parse as number or string").
"""
function register_optimizer!(name::AbstractString, kind::Symbol, fn, defaults::NamedTuple)
    kind in (:global, :local) || error("register_optimizer!: kind must be :global or :local, got $kind")
    OPTIMIZERS[String(name)] = (; kind, fn, defaults)
    nothing
end

function _parse_setting(default, v::AbstractString)
    v = strip(v)
    default isa Bool    && return parse(Bool, v)
    default isa Symbol  && return Symbol(v)
    default isa Integer && return something(tryparse(Int, v), Int(round(parse(Float64, v))))
    default isa AbstractFloat && return parse(Float64, v)
    default isa AbstractString && return String(v)
    # default === nothing: number if it looks like one, else string
    n = tryparse(Int, v);  n === nothing || return n
    x = tryparse(Float64, v);  x === nothing || return x
    String(v)
end

"""
    parse_optimizer_spec(str) -> (; name, label, settings)

Parse `"name"` or `"name:key=val,key=val"`. `settings` holds only the keys given in `str`; an unknown optimizer name
or setting key is an error that lists the valid ones. `label` is the spec string itself (used to tell apart several
specs of the same optimizer in tables and plots).
"""
function parse_optimizer_spec(str::AbstractString)
    name, rest = occursin(':', str) ? String.(split(str, ':'; limit = 2)) : (String(str), "")
    haskey(OPTIMIZERS, name) || error("unknown optimizer '$name'. Registered: $(join(sort(collect(keys(OPTIMIZERS))), ", "))")
    defaults = OPTIMIZERS[name].defaults
    kv = Pair{Symbol, Any}[]
    for item in split(rest, ','; keepempty = false)
        occursin('=', item) || error("optimizer spec '$str': expected key=value, got '$item'")
        k, v = split(item, '='; limit = 2)
        sym = Symbol(strip(k))
        sym in keys(defaults) || error("optimizer '$name' has no setting '$sym'. Valid settings: $(join(keys(defaults), ", "))")
        push!(kv, sym => _parse_setting(defaults[sym], v))
    end
    (; name, label = String(str), settings = NamedTuple(kv))
end

function _merge_settings(name, defaults, settings)
    for k in keys(settings)
        k in keys(defaults) || error("optimizer '$name' has no setting '$k'. Valid settings: $(join(keys(defaults), ", "))")
    end
    merge(defaults, settings)
end

# ---- shared objective and evaluation counter ----

"""
    make_objective(likelihood, prior, params; objective="posterior") -> (obj, b)

The function that the gradient-free optimizers minimize, `obj(x)` for the vector `x` of the FREE parameters of `prior`
(in prior order; fixed parameters sit at their fixed value, `b = extract_bounds(prior, params)`):
`-(llh + log prior)` for `objective="posterior"`, `-llh` for `"likelihood"`, `Inf` where not finite. It is the
quantity `local_find_mle` maximizes (up to the sign), so values from both kinds of optimizer are comparable.
"""
function make_objective(likelihood, prior, params; objective = "posterior")
    b = extract_bounds(prior, params)
    function obj(x::AbstractVector)
        all_x = [b.is_free_mask[i] ? x[b.free_x_index[i]] : b.fixed_float[i] for i in 1:length(b.all_keys)]
        full_params = NamedTuple{b.all_keys}(Tuple(all_x))
        llh_val = logdensityof(likelihood, full_params)
        if objective == "likelihood"
            isfinite(llh_val) ? -llh_val : Inf
        else
            prior_val = sum(logpdf(b.free_dists[j], x[j]) for j in 1:length(b.free_dists))
            isfinite(llh_val) && isfinite(prior_val) ? -(llh_val + prior_val) : Inf
        end
    end
    obj, b
end

"""Thrown by a `CountingObjective` whose evaluation budget is used up (caught by `global_find_mle`)."""
struct BudgetExhausted <: Exception end

"""
    CountingObjective(f, budget)

Callable wrapper around the objective `f` that records the cumulative number of evaluations and the best value so far
at every call (`evals`, `bests`) plus the best point (`best_x`), whatever library calls it, and throws
`BudgetExhausted` once `budget` evaluations are used (`budget = nothing`: no limit). One instance per optimizer run
(not thread-safe).
"""
mutable struct CountingObjective{F}
    f::F
    budget::Union{Int, Nothing}
    n::Int
    best::Float64
    best_x::Vector{Float64}
    evals::Vector{Int}
    bests::Vector{Float64}
end
CountingObjective(f, budget) = CountingObjective(f, budget, 0, Inf, Float64[], Int[], Float64[])

function (c::CountingObjective)(x::AbstractVector)
    !isnothing(c.budget) && c.n >= c.budget && throw(BudgetExhausted())
    v = Float64(c.f(x))
    c.n += 1
    if isempty(c.best_x) || v < c.best
        c.best = v
        c.best_x = Vector{Float64}(x)
    end
    push!(c.evals, c.n)
    push!(c.bests, c.best)
    v
end

# ---- global backends: fn(counter, b, cfg, budget, rng) ----

function _backend_optim_pso(counter, b, cfg, budget, rng)
    iters = !isnothing(cfg.iters) ? cfg.iters :
            !isnothing(budget)    ? max(1, budget ÷ cfg.particles) :
            error("optim-pso needs `iters` or a budget")
    x0  = clamp.(b.x0_vec, b.lo_vec, b.hi_vec)
    alg = Optim.ParticleSwarm(lower = b.lo_vec, upper = b.hi_vec, n_particles = cfg.particles)
    res = Optim.optimize(counter, x0, alg, Optim.Options(iterations = iters, store_trace = true))
    (; converged = missing, n_iters = Optim.iterations(res), native = Optim.f_trace(res))
end

function _backend_optim_sa(counter, b, cfg, budget, rng)
    iters = !isnothing(cfg.iters) ? cfg.iters : !isnothing(budget) ? budget : error("optim-sa needs `iters` or a budget")
    x0  = clamp.(b.x0_vec, b.lo_vec, b.hi_vec)
    res = Optim.optimize(counter, x0, Optim.SimulatedAnnealing(), Optim.Options(iterations = iters, store_trace = true))
    (; converged = missing, n_iters = Optim.iterations(res), native = Optim.f_trace(res))
end

function _mh_run(counter, b, budget, alg_fn, rng)
    isnothing(budget) && error("Metaheuristics optimizers need a budget (function evaluations)")
    bounds = permutedims(hcat(b.lo_vec, b.hi_vec))
    # `iterations` is large so that f_calls_limit (not the default iteration cap) is the stopping rule
    opts = Metaheuristics.Options(f_calls_limit = budget, iterations = 10_000_000, seed = rand(rng, 1:10^9))
    state = Metaheuristics.optimize(x -> counter(x), bounds, alg_fn(opts))   # (it wants a plain function, not a callable struct)
    (; converged = missing, n_iters = state.iteration, native = Float64[])
end
_backend_mh_pso(counter, b, cfg, budget, rng) = _mh_run(counter, b, budget, o -> Metaheuristics.PSO(N = cfg.N, C1 = cfg.C1, C2 = cfg.C2, ω = cfg.ω, options = o), rng)
_backend_mh_de(counter, b, cfg, budget, rng)  = _mh_run(counter, b, budget, o -> Metaheuristics.DE(N = cfg.N, F = cfg.F, CR = cfg.CR, options = o), rng)
_backend_mh_abc(counter, b, cfg, budget, rng) = _mh_run(counter, b, budget, o -> Metaheuristics.ABC(N = cfg.N, Ne = cfg.Ne, No = cfg.No, limit = cfg.limit, options = o), rng)

function _backend_bbo(counter, b, cfg, budget, rng)
    isnothing(budget) && error("BlackBoxOptim optimizers need a budget (function evaluations)")
    kw = Pair{Symbol, Any}[:SearchRange => [(b.lo_vec[i], b.hi_vec[i]) for i in eachindex(b.lo_vec)],
                           :Method => cfg.method, :MaxFuncEvals => budget, :TraceMode => :silent]
    isnothing(cfg.population) || push!(kw, :PopulationSize => cfg.population)
    # BlackBoxOptim draws from the (task-local) global RNG: seed it before the call (see `run_optimizer`)
    res = BlackBoxOptim.bboptimize(counter; kw...)
    (; converged = missing, n_iters = BlackBoxOptim.f_calls(res), native = Float64[])
end

"""
    global_find_mle(likelihood, prior, params, name; budget=nothing, settings=(;), rng, objective="posterior",
                    return_trace=false)

Gradient-free global search for the best point of the free parameters of `prior` (bounded box), using the registered
`:global` optimizer `name` (`optim-pso`, `optim-sa`, `mh-pso`, `mh-de`, `mh-abc`, `bbo`, ... see `OPTIMIZERS`) with
`settings` overriding its defaults. `budget` is the maximum number of objective evaluations (`nothing`: the
optimizer's own limit, e.g. `iters` of `optim-pso`). Not seeded from `params` except that `optim-pso`/`optim-sa`
start from it.

Returns `(llh, log_posterior, params, converged, n_iters)` like `local_find_mle`, except that `converged` is always `missing`:
a global search has no convergence criterion (it normally ends on its budget); a failed fit is signalled by `NaN` values. With `return_trace=true` a 6th element
`(; evals, best, native)`: cumulative objective evaluations and best objective value so far (`-log posterior` for the
posterior objective) at every call, and the library's own per-iteration trace if it has one. `log_posterior` is
`llh` for `objective="likelihood"` (flat prior; compare `llh`).
"""
function global_find_mle(likelihood, prior, params, name::AbstractString; budget = nothing, settings = (;),
                         rng = Random.default_rng(), objective = "posterior", return_trace = false)
    haskey(OPTIMIZERS, name) || error("unknown optimizer '$name'. Registered: $(join(sort(collect(keys(OPTIMIZERS))), ", "))")
    entry = OPTIMIZERS[name]
    entry.kind == :global || error("optimizer '$name' is a local optimizer; use `run_optimizer`/`local_find_mle`")
    cfg = _merge_settings(name, entry.defaults, settings)
    obj, b = make_objective(likelihood, prior, params; objective)
    all(isfinite, b.lo_vec) && all(isfinite, b.hi_vec) ||
        error("global optimizer '$name' needs finite bounds for every free parameter")
    counter = CountingObjective(obj, budget)
    info = try
        entry.fn(counter, b, cfg, budget, rng)
    catch e
        e isa BudgetExhausted || rethrow(e)
        (; converged = missing, n_iters = missing, native = Float64[])
    end
    if isempty(counter.best_x) || !isfinite(counter.best)
        nan_params = NamedTuple{b.all_keys}(ntuple(_ -> NaN, length(b.all_keys)))
        result = (NaN, NaN, nan_params, missing, missing)
        return return_trace ? (result..., (; evals = counter.evals, best = counter.bests, native = info.native)) : result
    end
    all_x = [b.is_free_mask[i] ? counter.best_x[b.free_x_index[i]] : b.fixed_float[i] for i in 1:length(b.all_keys)]
    fit = NamedTuple{b.all_keys}(Tuple(all_x))
    llh = logdensityof(likelihood, fit)
    log_post = objective == "likelihood" ? llh : -counter.best
    result = (llh, log_post, fit, info.converged, info.n_iters)
    return_trace ? (result..., (; evals = counter.evals, best = counter.bests, native = info.native)) : result
end

"""
    generate_seed(likelihood, prior, params, strategy, global_iters, rng; n_particles=10, return_trace=false, objective="posterior")

Start point for the local optimizer. `strategy == "random"`: a fresh random draw from the priors
(`randomize_free_params`). `"pso"` / `"sa"`: a short bounded global search (`optim-pso` with `n_particles` and
`global_iters` generations, resp. `optim-sa`) via `global_find_mle`; its best point (free entries; fixed entries at their
fixed values) is returned as a `params` NamedTuple. Any other registered `:global` optimizer name works too, with a
budget of `global_iters * n_particles` evaluations.

`objective` ("posterior" or "likelihood") selects the target of the global search; pass the same `objective` to the
subsequent `local_find_mle`, otherwise a posterior-biased seed is handed to a likelihood-only polish.

`return_trace=true` additionally returns the best objective value so far per generation (Optim's own trace) for
`"pso"`/`"sa"`, or per evaluation for the other optimizers.
"""
function generate_seed(likelihood, prior, params, strategy, global_iters, rng; n_particles = 10, return_trace = false, objective = "posterior")
    if strategy == "random"
        seeded = randomize_free_params(rng, params, prior)
        return return_trace ? (seeded, Float64[]) : seeded
    end
    res = if strategy == "pso"
        global_find_mle(likelihood, prior, params, "optim-pso"; settings = (; particles = n_particles, iters = global_iters),
                        rng, objective, return_trace = true)
    elseif strategy == "sa"
        global_find_mle(likelihood, prior, params, "optim-sa"; settings = (; iters = global_iters),
                        rng, objective, return_trace = true)
    else
        global_find_mle(likelihood, prior, params, strategy; budget = global_iters * n_particles,
                        rng, objective, return_trace = true)
    end
    seeded = res[3]
    if !return_trace
        return seeded
    end
    tr = res[6]
    (seeded, strategy in ("pso", "sa") ? tr.native : tr.best)
end

# ---- run any registered optimizer on one problem ----

"""
    run_optimizer(spec, likelihood, prior, params; budget=nothing, rng, objective="posterior") -> NamedTuple

Run the optimizer described by `spec` (a string, see `parse_optimizer_spec`, or its parsed NamedTuple) on one
problem and return a uniform record
`(; label, llh, log_post, params, converged, n_iters, n_calls, time, trace_evals, trace_best)`:
(`converged` is `missing` for the `:global` optimizers: they end on their budget, there is no convergence criterion),
wall time of the fit, the number of objective evaluations (`n_calls`; for `:local` optimizers the number of
likelihood calls including every ForwardDiff pass, `missing` if it could not be counted), and the trace of the best
objective value so far (`-log posterior` for the posterior objective) against cumulative evaluations (for `:local`
optimizers per iteration, with `trace_evals = missing`). `prior` is a NamedTuple of distributions or its `distprod`.
"""
function run_optimizer(spec, likelihood, prior, params; budget = nothing, rng = Random.default_rng(), objective = "posterior")
    spec isa AbstractString && (spec = parse_optimizer_spec(spec))
    entry = OPTIMIZERS[spec.name]
    cfg = _merge_settings(spec.name, entry.defaults, spec.settings)
    t0 = time()
    if entry.kind == :global
        r = global_find_mle(likelihood, prior, params, spec.name; budget, settings = spec.settings, rng, objective, return_trace = true)
        return (; label = spec.label, llh = r[1], log_post = r[2], params = r[3], converged = r[4], n_iters = r[5],
                n_calls = isempty(r[6].evals) ? 0 : r[6].evals[end], time = time() - t0,
                trace_evals = r[6].evals, trace_best = r[6].best)
    end
    prior_local = prior isa NamedTuple ? distprod(; prior...) : prior   # BAT needs the distprod-wrapped prior
    r = entry.fn(likelihood, prior_local, params, cfg, objective)   # (llh, log_post, params, converged, n_iters, trace, n_calls)
    (; label = spec.label, llh = r[1], log_post = r[2], params = r[3], converged = r[4], n_iters = r[5],
       n_calls = r[7], time = time() - t0, trace_evals = missing, trace_best = r[6])
end

# local backends: fn(likelihood, prior, params, cfg, objective) -> (llh, log_post, params, converged, n_iters, trace, n_calls)

# counts the likelihood calls of a local fit (every ForwardDiff pass counts as one call)
struct CountedLikelihood{L}
    l::L
    n::Base.RefValue{Int}
end
DensityInterface.DensityKind(c::CountedLikelihood) = DensityInterface.DensityKind(c.l)   # BAT requires an `IsDensity` likelihood
DensityInterface.logdensityof(c::CountedLikelihood, x) = (c.n[] += 1; DensityInterface.logdensityof(c.l, x))

function _ad_backend(name::Symbol)
    name == :auto && return nothing
    name == :forward && return ADTypes.AutoForwardDiff()
    name == :polyester && return ADTypes.AutoPolyesterForwardDiff()
    error("unknown ad setting '$name' (auto, forward, polyester)")
end

function _local_optim(alg_name)
    function (likelihood, prior, params, cfg, objective)
        algorithm = make_algorithm(alg_name, cfg.step_size)
        counted = CountedLikelihood(likelihood, Ref(0))
        r = local_find_mle(counted, prior, params; fit_method = "optim", algorithm, iterations = cfg.iterations,
                           g_tol = cfg.g_tol, f_tol = cfg.f_tol, x_tol = cfg.x_tol, return_trace = true, objective,
                           ad_backend = _ad_backend(cfg.ad))
        (r[1], r[2], r[3], r[4], r[5], r[6], counted.n[])
    end
end
_local_newtrinos(likelihood, prior, params, cfg, objective) =
    (r = local_find_mle(likelihood, prior, params; fit_method = "newtrinos", return_trace = true); (r[1], r[2], r[3], r[4], r[5], r[6], missing))
_local_ext(likelihood, prior, params, cfg, objective) =
    (r = local_find_mle(likelihood, prior, params; fit_method = "ext", iterations = cfg.iterations, g_tol = cfg.g_tol, return_trace = true);
     (r[1], r[2], r[3], r[4], r[5], r[6], missing))

# ---- registered optimizers ----
const _LOCAL_DEFAULTS = (; iterations = 2000, g_tol = 1e-6, f_tol = 0.0, x_tol = 0.0, step_size = 1.0, ad = :auto)
register_optimizer!("optim-lbfgs",             :local, _local_optim("lbfgs"),             _LOCAL_DEFAULTS)
register_optimizer!("optim-gradientdescent",   :local, _local_optim("gradientdescent"),   _LOCAL_DEFAULTS)
register_optimizer!("optim-conjugategradient", :local, _local_optim("conjugategradient"), _LOCAL_DEFAULTS)
register_optimizer!("optim-neldermead",        :local, _local_optim("neldermead"),        _LOCAL_DEFAULTS)
function _local_lbfgsb(likelihood, prior, params, cfg, objective)
    counted = CountedLikelihood(likelihood, Ref(0))
    r = local_find_mle(counted, prior, params; fit_method = "lbfgsb", iterations = cfg.iterations, g_tol = cfg.reltol, lbfgsb_m = cfg.m,
                       return_trace = true, objective, ad_backend = _ad_backend(cfg.ad))
    (r[1], r[2], r[3], r[4], r[5], r[6], counted.n[])
end
register_optimizer!("lbfgsb", :local, _local_lbfgsb, (; iterations = 1000, reltol = 1e-7, m = 10, ad = :auto))
register_optimizer!("newtrinos", :local, _local_newtrinos, (;))
register_optimizer!("ext",       :local, _local_ext,       (; iterations = 2000, g_tol = 1e-6))
register_optimizer!("optim-pso", :global, _backend_optim_pso, (; particles = 20, iters = nothing))
register_optimizer!("optim-sa",  :global, _backend_optim_sa,  (; iters = nothing))
register_optimizer!("mh-pso", :global, _backend_mh_pso, (; N = 100, C1 = 2.0, C2 = 2.0, ω = 0.8))
register_optimizer!("mh-de",  :global, _backend_mh_de,  (; N = 100, F = 0.7, CR = 0.5))
register_optimizer!("mh-abc", :global, _backend_mh_abc, (; N = 50, Ne = 25, No = 25, limit = 10))
register_optimizer!("bbo",      :global, _backend_bbo, (; method = :adaptive_de_rand_1_bin_radiuslimited, population = nothing))
register_optimizer!("bbo-xnes",  :global, _backend_bbo, (; method = :xnes, population = nothing))
register_optimizer!("bbo-dxnes", :global, _backend_bbo, (; method = :dxnes, population = nothing))

"""
    _nudge_into_bounds(prior, params; ε=1e-8)

Return `params` with every free parameter moved at least `ε` inside the bounds of its prior
(`[lo+ε, hi-ε]`; an infinite side is left alone), and fixed parameters untouched. A start
exactly on a bound makes BAT's transformation of bounded priors to an unconstrained space
non-finite, so the fit would fail; interior start values are returned unchanged.
"""
function _nudge_into_bounds(prior, params; ε=1e-8)
    nudged = Dict{Symbol, Any}()
    for k in keys(params)
        k in keys(prior) || continue
        d = prior[k]
        _is_fixed(d) && continue
        lo, hi = minimum(d), maximum(d)
        v = params[k]
        nv = v
        isfinite(lo) && nv < lo + ε && (nv = lo + ε)
        isfinite(hi) && nv > hi - ε && (nv = hi - ε)
        nv == v || (nudged[k] = convert(typeof(v), nv))
    end
    isempty(nudged) ? params : merge(params, NamedTuple(nudged))
end

"""
    local_find_mle(likelihood, prior, params; fit_method, algorithm, iterations, g_tol, f_tol, x_tol)

Find the MLE using either:
- `fit_method == "optim"` (default): `BAT.bat_findmode` + `OptimAlg`, i.e. the same
  machinery `Newtrinos.find_mle` uses (BAT's `PriorToNormal` reparametrization of bounded
  priors into unconstrained space before running the chosen Optim.jl algorithm), but with
  the algorithm/iterations/tolerances CLI-tunable instead of hardcoded.

  Earlier versions of this function used `Optim.optimize` + `Optim.Fminbox` directly
  (box constraints in the raw parameter space, matching `Newtrinos.find_mle_ext`'s
  approach). That turned out to be dramatically slower than `Newtrinos.find_mle` on this
  likelihood: `Fminbox`'s barrier method re-solves a full inner optimization at each of
  several shrinking-penalty outer steps, multiplying the number of (expensive) likelihood
  + gradient evaluations, and the barrier term distorts curvature near the many
  zero-lower-bound NSI parameters (`ε_eμ_abs`, etc.), which is exactly where these fits
  tend to sit. Routing through `OptimAlg`'s reparametrization avoids all of that, the same
  way `Newtrinos.find_mle` already does.
- `fit_method == "lbfgsb"`: L-BFGS-B (`OptimizationLBFGSB`, through BAT's `OptimizationAlg`) in the same transformed space, like
  `Newtrinos.find_mle` but tunable: `iterations` -> `maxiters`, `g_tol` -> `reltol` (projected-gradient tolerance), `lbfgsb_m` -> memory.
  Returns the per-iteration `-log posterior` as trace when `return_trace=true`; `converged` is the solver's `Success` return code.
- `fit_method == "newtrinos"`: the standard `Newtrinos.find_mle` unmodified (hardcoded
  LBFGS/maxiters=2000/f_abstol=1e-2), included here so it can be compared against the
  tunable path on the same roundtrips.
- `fit_method == "ext"`: `Newtrinos.find_mle_ext` unmodified.

All three return `(llh, log_posterior, params, converged, n_iterations)`. `n_iterations`
is only available for `"optim"` (`Optim.iterations(res.info)`); it's `missing` for
`"newtrinos"`/`"ext"`, since those don't expose an iteration count in their own return value.

`return_trace=true` additionally returns, as a 6th element, the best-objective-value-per-
iteration trace (`Optim.f_trace`, in `-logposterior` units) for `"optim"`; an empty
`Float64[]` for `"newtrinos"`/`"ext"`, which don't expose Optim's internals. Default
`false` keeps the normal 5-tuple return for all existing callers.

`objective` (only used for `fit_method == "optim"`) selects what `bat_findmode`
maximizes: `"posterior"` (default, likelihood x prior, i.e. MAP -- unchanged behavior)
or `"likelihood"` (likelihood alone, via `flat_prior(prior)` -- see that function's
docstring). `llh` (the first return value) is always `logdensityof(likelihood, ...)`
either way; only the target measure `bat_findmode` searches over changes, so
`log_posterior` (2nd return value) is the *flat-prior* posterior (≈ llh + a constant)
when `objective="likelihood"`, not the true MAP posterior -- callers comparing across
`objective` settings should use `llh`, not `log_posterior`, for a like-for-like
likelihood comparison.
"""
function local_find_mle(likelihood, prior, params; fit_method="optim", algorithm=nothing, iterations=nothing, g_tol=nothing, f_tol=nothing, x_tol=nothing, return_trace=false, objective="posterior", ad_backend=nothing, lbfgsb_m=10)
    if fit_method == "newtrinos"
        llh_, lp_, res_ = Newtrinos.find_mle(likelihood, prior, params)  # no convergence flag returned: finite posterior counts as converged
        result = (llh_, lp_, res_, isfinite(lp_), missing)
        return return_trace ? (result..., Float64[]) : result
    elseif fit_method == "ext"
        result = (Newtrinos.find_mle_ext(likelihood, prior, params; g_tol=g_tol, maxiters=iterations)..., missing)
        return return_trace ? (result..., Float64[]) : result
    end

    # fit_method == "optim" or "lbfgsb"
    fit_method in ("optim", "lbfgsb") || error("local_find_mle: unknown fit_method '$fit_method' (optim, lbfgsb, newtrinos, ext)")
    params = _nudge_into_bounds(prior, params)  # a start on a prior bound would make BAT's transformation non-finite
    target_prior = objective == "likelihood" ? flat_prior(prior) : prior
    posterior = PosteriorMeasure(likelihood, target_prior)
    n_free = count(!(d isa ValueShapes.ConstValueDist) && !(d isa Number) for d in values(prior))

    @info "Running Optimization ($fit_method, $(fit_method == "lbfgsb" ? "L-BFGS-B" : typeof(algorithm).name.name)) for point" n_free=n_free iterations=iterations objective=objective

    try
        t0 = time()
        # Newtrinos.find_mle sets this internally (adsel = Newtrinos.select_ad(length(params)));
        # bat_findmode needs an AD backend on the context too, or gradient-based algorithms
        # (LBFGS etc.) fail immediately with "requires automatic differentiation".
        # set_batcontext(ad = Newtrinos.select_ad(n_free)) this caused thread racing / threading conflicts
	ad = ad_backend === nothing ? Newtrinos.select_ad(n_free) : ad_backend # for choosing ad backend when calling local_find_mle
	context = BATContext(ad = ad) # for setting the bat context locally inside bat_findmode and not globally 
        @info "starting optimization"
        # BAT 5: the Optim backend sits inside TransformedMaxDensity (which also takes the start point); `store_trace` is an OptimAlg field
        # `bat_findmode` only returns the mode, so go through `evalmeasure` to also get the Optim result (convergence, iterations, trace)
        # fit_method "lbfgsb": L-BFGS-B (OptimizationLBFGSB, via BAT's OptimizationAlg) in the same transformed space, as Newtrinos.find_mle does but
        # tunable: `iterations` -> maxiters, `g_tol` -> reltol (L-BFGS-B's projected-gradient tolerance), `lbfgsb_m` -> memory size.
        # Its first step is scaled by the gradient norm, unlike LBFGS with a static first step. `store_trace` records the objective per iteration.
        optalg_obj = fit_method == "lbfgsb" ?
            OptimizationAlg(optalg=OptimizationLBFGSB.LBFGSB(m=lbfgsb_m), maxiters=iterations, reltol=Float64(g_tol), store_trace=return_trace) :
            OptimAlg(optalg=algorithm, maxiters=iterations, store_trace=return_trace, kwargs=(g_tol=g_tol, f_reltol=f_tol, x_reltol=x_tol))
        em = BAT.evalmeasure(posterior, TransformedMaxDensity(optalg=optalg_obj, init=ExplicitInit([params])), context) # set bat context locally
        res = (result = BAT.mode(em), info = em.evalinfo.result.info, trace = fit_method == "lbfgsb" ? em.evalinfo.result.trace_trafo : nothing)
        @info "extracting fit results"
        llh      = logdensityof(likelihood, res.result)
        log_post = logdensityof(posterior, res.result)
        if fit_method == "lbfgsb"
            converged = Symbol(res.info.retcode) == :Success    # (MaxIters: ran into `iterations`)
            n_iters   = res.info.stats.iterations
            f_trace   = isnothing(res.trace) ? Float64[] : -collect(Float64, res.trace.logd)   # -log posterior per iteration, as Optim's f_trace
        else
            converged = Optim.converged(res.info)
            n_iters   = Optim.iterations(res.info)
            f_trace   = return_trace ? Optim.f_trace(res.info) : Float64[]
        end
        @info "Optimization ($fit_method) finished" converged=converged elapsed_s=round(time()-t0, digits=1) iterations=n_iters
        result = (llh, log_post, res.result, converged, n_iters)
        return return_trace ? (result..., f_trace) : result
    catch e
	(e isa ArgumentError || e isa AssertionError) || rethrow(e)
        @warn "local_find_mle caught $(typeof(e)), returning NaN" exception=(e, catch_backtrace())
        nan_params = NamedTuple{keys(prior)}(ntuple(_ -> NaN, length(keys(prior))))
        result = (NaN, NaN, nan_params, false, missing)
        return return_trace ? (result..., Float64[]) : result
    end
end
