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

# One table of named optimizers; a spec "name:key=val,key=val" (`parse_optimizer_spec`) selects an entry and overrides defaults. Kinds: :global (gradient-free
# search in the prior box, `fn(counter, b, cfg, budget, rng)`, run by `global_find_mle`), :local (through `local_find_mle`), and the stages :seed/:alias below.
# To add an optimizer: write a backend and call `register_optimizer!` (end of this section).

const OPTIMIZERS = Dict{String, NamedTuple}()

"""
    register_optimizer!(name, kind, fn, defaults)

Add (or replace) an entry of the optimizer registry. `kind` is `:global`, `:local`, `:seed` or `:alias` (see the strategy section); `defaults` is the NamedTuple
of settings the entry understands (a spec may override only these keys, the type of a default decides how a
command-line value is parsed; a default of `nothing` means "unset, parse as number or string").
"""
function register_optimizer!(name::AbstractString, kind::Symbol, fn, defaults::NamedTuple; open_keys::Bool = false)
    kind in (:global, :local, :seed, :alias) || error("register_optimizer!: kind must be :global, :local, :seed or :alias, got $kind")
    OPTIMIZERS[String(name)] = (; kind, fn, defaults, open_keys)   # open_keys: any setting key is allowed (manual seed: parameter names)
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
    entry = OPTIMIZERS[name]
    defaults = entry.defaults
    kv = Pair{Symbol, Any}[]
    for item in split(rest, ','; keepempty = false)
        occursin('=', item) || error("optimizer spec '$str': expected key=value, got '$item'")
        k, v = split(item, '='; limit = 2)
        sym = Symbol(strip(k))
        entry.open_keys || sym in keys(defaults) || error("optimizer '$name' has no setting '$sym'. Valid settings: $(join(keys(defaults), ", "))")
        push!(kv, sym => _parse_setting(sym in keys(defaults) ? defaults[sym] : nothing, v))
    end
    (; name, label = String(str), settings = NamedTuple(kv))
end

function _merge_settings(name, entry, settings)
    entry.open_keys && return merge(entry.defaults, settings)
    for k in keys(settings)
        k in keys(entry.defaults) || error("optimizer '$name' has no setting '$k'. Valid settings: $(join(keys(entry.defaults), ", "))")
    end
    merge(entry.defaults, settings)
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
_backend_mh_abc(counter, b, cfg, budget, rng) = _mh_run(counter, b, budget, o -> Metaheuristics.ABC(N = cfg.N, Ne = something(cfg.Ne, (cfg.N + 1) ÷ 2), No = something(cfg.No, (cfg.N + 1) ÷ 2),
                                                                                                         limit = cfg.limit, options = o), rng)   # Ne/No default to (N+1)÷2 like the library

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
    cfg = _merge_settings(name, entry, settings)
    isnothing(get(cfg, :evals, nothing)) || (budget = cfg.evals)   # `evals` (a setting of every global optimizer) overrides the default budget
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
    entry.kind in (:global, :local) || error("run_optimizer: '$(spec.name)' is a $(entry.kind) stage, not an optimizer; use `find_mle_advanced`")
    cfg = _merge_settings(spec.name, entry, spec.settings)
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
const _LBFGSB_DEFAULTS = (; iterations = 1000, reltol = 1e-7, m = 10, ad = :auto)
register_optimizer!("lbfgsb", :local, _local_lbfgsb, _LBFGSB_DEFAULTS)
register_optimizer!("newtrinos", :local, _local_newtrinos, (;))
register_optimizer!("ext",       :local, _local_ext,       (; iterations = 2000, g_tol = 1e-6))
register_optimizer!("optim-pso", :global, _backend_optim_pso, (; particles = 20, iters = nothing, evals = nothing))
register_optimizer!("optim-sa",  :global, _backend_optim_sa,  (; iters = nothing, evals = nothing))
register_optimizer!("mh-pso", :global, _backend_mh_pso, (; N = 100, C1 = 2.0, C2 = 2.0, ω = 0.8, evals = nothing))
register_optimizer!("mh-de",  :global, _backend_mh_de,  (; N = 100, F = 0.7, CR = 0.5, evals = nothing))
register_optimizer!("mh-abc", :global, _backend_mh_abc, (; N = 50, Ne = nothing, No = nothing, limit = 10, evals = nothing))
register_optimizer!("bbo",      :global, _backend_bbo, (; method = :adaptive_de_rand_1_bin_radiuslimited, population = nothing, evals = nothing))
register_optimizer!("bbo-xnes",  :global, _backend_bbo, (; method = :xnes, population = nothing, evals = nothing))
register_optimizer!("bbo-dxnes", :global, _backend_bbo, (; method = :dxnes, population = nothing, evals = nothing))

# A strategy is an ordered list of stages (spec strings); every stage starts from the result of the previous one. Kinds: seed (`random`, `manual`),
# alias (`smseed`, `mle`; `lbfgsb` settings), global, local. Seed, alias and global stages may only be FIRST (see docs/src/README_optimizer_tests.md).
_seed_random(cfg, prior, params, rng) = rand(rng, prior isa NamedTuple ? distprod(; prior...) : prior)

function _seed_manual(cfg, prior, params, rng)
    for k in keys(cfg)
        k in keys(params) || error("manual stage: unknown parameter '$k'. Parameters: $(join(keys(params), ", "))")
    end
    merged = merge(params, NamedTuple{keys(cfg)}(Float64.(values(cfg))))
    prior isa NamedTuple ? _fix_conditioned(merged, prior) : merged   # pinned parameters keep their pinned value
end

# alias expansion: spec -> vector of stages (name, label, settings, null)
_expand_smseed(sp) = [(; name = "lbfgsb", label = sp.label, settings = sp.settings, null = true)]
_expand_mle(sp)    = [(; name = "random", label = sp.label * " [random]", settings = (;), null = false),
                      (; name = "lbfgsb", label = sp.label * " [fit]", settings = sp.settings, null = false)]

register_optimizer!("random", :seed, _seed_random, (;))
register_optimizer!("manual", :seed, _seed_manual, (;); open_keys = true)
register_optimizer!("smseed", :alias, _expand_smseed, _LBFGSB_DEFAULTS)
register_optimizer!("mle",    :alias, _expand_mle,    _LBFGSB_DEFAULTS)

"""
    normalize_strategy(strategy) -> Vector{NamedTuple}

Parse a strategy (a vector of spec strings and/or parsed stage NamedTuples `(; name, label, settings)`), expand the
aliases `smseed`/`mle` into their stages and check the order: seed, alias and global stages only as the first stage,
all later stages local. Returns stages `(; name, label, settings, null)`; `null = true` marks a local fit with the
null prior (`smseed`).
"""
function normalize_strategy(strategy)
    isempty(strategy) && error("strategy: needs at least one stage")
    specs = [s isa AbstractString ? parse_optimizer_spec(s) : s for s in strategy]
    stages = NamedTuple[]
    for (i, sp) in enumerate(specs)
        haskey(OPTIMIZERS, sp.name) || error("unknown optimizer '$(sp.name)'. Registered: $(join(sort(collect(keys(OPTIMIZERS))), ", "))")
        kind = OPTIMIZERS[sp.name].kind
        if i > 1 && kind != :local
            error("strategy: '$(sp.label)' is $(kind == :alias ? "an alias" : "a $kind") stage and may only be the first stage; every later stage must be a local optimizer " *
                  "(a global search or a seed ignores the previous result)")
        end
        if kind == :alias
            append!(stages, OPTIMIZERS[sp.name].fn(sp))
        else
            push!(stages, (; name = sp.name, label = sp.label, settings = sp.settings, null = false))
        end
    end
    stages
end

"""
    parse_strategy_arg(str; nseeds=1) -> (; stages, nseeds, label)

A strategy from one command-line argument: stage specs separated by spaces, optionally followed by `@N` = N independent attempts of
this strategy (the best is kept), e.g. `"optim-pso:particles=20 lbfgsb @3"`. Without `@N` the given default `nseeds` applies.
The stages are checked with `normalize_strategy`.
"""
function parse_strategy_arg(str; nseeds = 1)
    toks = split(strip(str))
    if !isempty(toks) && (m = match(r"^@(\d+)$", toks[end])) !== nothing
        nseeds = parse(Int, m.captures[1])
        pop!(toks)
    end
    nseeds >= 1 || error("strategy '$str': the number of attempts must be at least 1")
    stages = String.(toks)
    normalize_strategy(stages)   # validates names, settings and order
    (; stages, nseeds, label = join(stages, " ") * (nseeds > 1 ? " @$nseeds" : ""))
end

"""Stage spec of a manual seed for the Julia API: `manual_stage((; a = 1.0, b = 2.0))` = the CLI stage `manual:a=1.0,b=2.0`."""
manual_stage(params::NamedTuple) = (; name = "manual", label = "manual", settings = params)

"""
    with_local_defaults(spec, defaults) -> spec

Spec of a local optimizer with `defaults` (a NamedTuple, e.g. `(; ad = :forward)`) filled in for the settings it knows and the spec does not
set itself. Other kinds of optimizer are returned unchanged.
"""
function with_local_defaults(spec, defaults)
    entry = OPTIMIZERS[spec.name]
    entry.kind == :local || return spec
    (; spec..., settings = merge(NamedTuple(k => v for (k, v) in pairs(defaults) if k in keys(entry.defaults)), spec.settings))
end

"""
    final_local_spec(strategy) -> spec

The local fit that an extra fit should use (e.g. the reference fits of `roundtrips.jl`): the last stage of `strategy` if it is a local
optimizer on the hypothesis prior, otherwise a default `lbfgsb`.
"""
function final_local_spec(strategy)
    last = normalize_strategy(strategy)[end]
    OPTIMIZERS[last.name].kind == :local && !last.null ? (; name = last.name, label = last.label, settings = last.settings) : parse_optimizer_spec("lbfgsb")
end

function _stage_record(label, params; llh = NaN, log_post = NaN)
    (; label, llh, log_post, params, converged = missing, n_iters = missing, n_calls = 0, time = 0.0, trace_evals = Int[], trace_best = Float64[])
end

function _run_strategy(likelihood, prior, start, stages, strategy_label; rng, objective, null_prior, budget, local_defaults)
    cur = start
    records = NamedTuple[]
    seed_params = start
    for (i, st) in enumerate(stages)
        entry = OPTIMIZERS[st.name]
        i == length(stages) && (seed_params = cur)
        if entry.kind == :seed
            t0 = time()
            cfg = _merge_settings(st.name, entry, st.settings)
            cur = entry.fn(cfg, prior, cur, rng)
            push!(records, merge(_stage_record(st.label, cur), (; time = time() - t0)))
        else
            stage_prior = prior
            start_here = cur
            if st.null
                isnothing(null_prior) && error("strategy: stage '$(st.label)' needs a null prior (`null_prior` = the hypothesis' prior with the tested parameters at their null values)")
                stage_prior = null_prior
                start_here = _fix_conditioned(cur, null_prior)   # the null prior holds the tested parameters at their null values
            end
            rec = run_optimizer(with_local_defaults((; name = st.name, label = st.label, settings = st.settings), local_defaults), likelihood, stage_prior, start_here; budget, rng, objective)
            push!(records, rec)
            cur = rec.params
        end
        if !all(isfinite, (Float64(v) for v in values(cur)))   # a failed stage: do not feed NaN to the next one
            break
        end
    end
    last = records[end]
    xs, ys, stg, span = _strategy_trace(records, objective)
    total_calls = any(r -> r.n_calls === missing, records) ? missing : sum(r.n_calls for r in records)
    merge(last, (; time = sum(r.time for r in records), seed_params, stages = records, strategy = strategy_label,
                 trace_x = xs, trace_y = ys, trace_stage = stg, trace_span = span, total_calls))
end

# Trace of one stage as (evaluations, objective value so far). Global stages record the best objective value per evaluation;
# local fits record per iteration (in BAT's transformed space), shifted here to end at the stage's final objective value, with the
# iterations spread evenly over the stage's likelihood calls (over the iteration count if the calls were not counted).
function _stage_trace(rec, objective)
    score = objective == "likelihood" ? rec.llh : rec.log_post
    if rec.trace_evals === missing
        tr = collect(Float64, rec.trace_best)
        (isempty(tr) || !isfinite(score)) && return Float64[], Float64[]
        tr = tr .- tr[end] .- score
        n = length(tr)
        calls = rec.n_calls === missing ? n : rec.n_calls
        return (n == 1 ? [Float64(calls)] : [1 + (calls - 1) * (j - 1) / (n - 1) for j in 1:n]), tr
    end
    Float64.(rec.trace_evals), Float64.(rec.trace_best)
end

_stage_span(rec) = rec.n_calls === missing ? Float64(length(rec.trace_best)) : Float64(rec.n_calls)

# the stages' traces one after the other (each stage's evaluation counter continues where the previous one ended) as the best
# objective value so far; also returns the total number of evaluations
function _strategy_trace(records, objective)
    xs, ys, stg, off = Float64[], Float64[], Int[], 0.0
    for (k, r) in enumerate(records)
        x, y = _stage_trace(r, objective)
        append!(xs, x .+ off)
        append!(ys, y)
        append!(stg, fill(k, length(x)))   # which stage of the strategy each trace point belongs to
        off += _stage_span(r)
    end
    xs, _running_best(ys), stg, off
end
_running_best(ys) = accumulate(min, replace(ys, NaN => Inf))

"""
    find_mle_advanced(likelihood, prior, params; strategy, nseeds=1, starts=nothing, rng, objective="posterior",
                      null_prior=nothing, budget=nothing, local_defaults=(;)) -> NamedTuple

Run a strategy (see `normalize_strategy`: a chain of stages such as `["smseed", "lbfgsb"]` or
`["optim-pso:particles=20,iters=50", "lbfgsb:iterations=300"]`) from the start point `params` and return the record of
the last stage (the fields of `run_optimizer`) extended by `seed_params` (the point the last stage started from),
`stages` (the records of all stages) and `strategy` (the stage labels); `time` is the total time of all stages.

With `nseeds > 1` the whole strategy is repeated (random stages and global searches draw from `rng` progressively) and the
best fit by `log_post` is kept (`NaN` counts as -Inf); with `starts` (a vector of start points) there is one attempt per start
instead. The record also has `n_attempts`, `total_time` and `total_calls` (summed over ALL attempts and stages; `missing` if a local stage
did not count its calls; `time` stays that of the kept attempt) and the combined trace `trace_x` (evaluations) / `trace_y` (best
objective value so far, `-score`) through all stages and attempts, with `trace_stage` (the stage index of every trace point). `null_prior` is needed by `smseed`; `budget` is the default number of evaluations of global stages without `evals`;
`local_defaults` overrides the defaults of the local stages' settings (only keys the stage knows), e.g. `(; ad = :forward)`.
A stage that returns non-finite parameters ends the chain (the record then holds NaN).
"""
function find_mle_advanced(likelihood, prior, params; strategy, nseeds = 1, starts = nothing, rng = Random.default_rng(),
                           objective = "posterior", null_prior = nothing, budget = nothing, local_defaults = (;))
    stages = normalize_strategy(strategy)
    strategy_label = join((s isa AbstractString ? s : s.label for s in strategy), " ")
    n_attempts = isnothing(starts) ? nseeds : length(starts)
    best, best_score = nothing, -Inf
    attempts = NamedTuple[]
    for s in 1:n_attempts
        start = isnothing(starts) ? params : starts[s]
        rec = _run_strategy(likelihood, prior, start, stages, strategy_label; rng, objective, null_prior, budget, local_defaults)
        push!(attempts, rec)
        score = isnan(rec.log_post) ? -Inf : rec.log_post
        if best === nothing || score > best_score
            best, best_score = rec, score
        end
    end
    # the cost of the strategy is that of ALL attempts; the trace runs through the attempts one after the other
    xs, ys, stg, off = Float64[], Float64[], Int[], 0.0
    for a in attempts
        append!(xs, a.trace_x .+ off)
        append!(ys, a.trace_y)
        append!(stg, a.trace_stage)
        off += a.trace_span
    end
    merge(best, (; n_attempts, total_time = sum(a.time for a in attempts),
                 total_calls = any(a -> a.total_calls === missing, attempts) ? missing : sum(a.total_calls for a in attempts),
                 trace_x = xs, trace_y = _running_best(ys), trace_stage = stg, trace_span = off,
                 # the same trace attempt by attempt, each from its own start (x from 0, best objective so far within the attempt)
                 attempt_traces = [(; x = a.trace_x, y = _running_best(a.trace_y), stage = a.trace_stage) for a in attempts]))
end

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
    if n_free == 0   # nothing to fit (e.g. every free parameter is a scanned one): L-BFGS-B would crash on zero dimensions
        result = (logdensityof(likelihood, params), logdensityof(posterior, params), params, true, 0)
        return return_trace ? (result..., Float64[]) : result
    end

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
        # `bat_findmode` only returns the mode: go through `evalmeasure` to also get the Optim/Optimization result (convergence, iterations, trace).
        # "lbfgsb" = L-BFGS-B (OptimizationLBFGSB via BAT's OptimizationAlg, tunable: `iterations`, `g_tol` -> reltol, `lbfgsb_m`), same transformed space as Newtrinos.find_mle.
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

# `profile_advanced` = `Newtrinos.profile` with a strategy for the fit at every grid point, warm starting and a `NewtrinosResult` with the optimizer
# run data in `meta`; its pieces (`_profile_plan`, `_profile_point_fit`, `_profile_result`) are also used directly by nsi/nsi_scan.jl.

import ProgressMeter

_spec_of(s) = s isa AbstractString ? parse_optimizer_spec(s) : s

"""The stages that run at a warm-started grid point by default: the strategy's local stages (a seed, alias or global first stage
would ignore the warm start), `lbfgsb` if there are none."""
function _default_warm_strategy(strategy)
    strategy = collect(strategy)
    OPTIMIZERS[_spec_of(strategy[1]).name].kind == :local && return strategy
    rest = strategy[2:end]
    isempty(rest) ? ["lbfgsb"] : rest
end

"""
    _profile_plan(likelihood, priors, params, scanpoints, mesh_flat, scan_keys; strategy, nseeds=1, seed=1234, objective="posterior",
                  null_prior=nothing, budget=nothing, local_defaults=(;), starts=[(;)], warm_strategy=nothing) -> NamedTuple

Everything a profile scan needs besides the likelihood: checks the strategy and, for a seed or alias FIRST stage (`random`, `manual`,
`smseed`, `mle`), computes the seeds ONCE (they do not depend on the grid point; `nseeds` seeds from `Xoshiro(seed)`, `smseed` is
deterministic so one fit is shared) from the hypothesis prior `priors` and the start `params`. A global first stage is run at
every grid point instead. `scanpoints` are the per-point priors (scanned parameters fixed), `mesh_flat` the scanned values of every
point and `scan_keys` their names.
"""
function _profile_plan(likelihood, priors, params, scanpoints, mesh_flat, scan_keys; strategy, nseeds = 1, seed = 1234, objective = "posterior",
                       null_prior = nothing, budget = nothing, local_defaults = (;), starts = [(;)], warm_strategy = nothing)
    strategy = collect(strategy)
    normalize_strategy(strategy)   # checks the stage names, settings and order before any expensive work
    nseeds >= 1 || error("nseeds must be at least 1")
    objective in ("posterior", "likelihood") || error("objective must be 'posterior' or 'likelihood', got '$objective'")
    first_spec = _spec_of(strategy[1])
    hoisted = OPTIMIZERS[first_spec.name].kind in (:seed, :alias)
    rest_strategy = strategy[2:end]
    hoisted && isempty(rest_strategy) && error("strategy: '$(first_spec.label)' only creates a seed; add a fitting stage after it, e.g. [\"$(first_spec.label)\", \"lbfgsb\"]")
    warm = isnothing(warm_strategy) ? _default_warm_strategy(strategy) : collect(warm_strategy)
    all(s -> OPTIMIZERS[_spec_of(s).name].kind == :local, warm) || error("warm_strategy: every stage must be a local optimizer (a seed or global stage would ignore the warm start)")
    # every slot is defined (the plan is shipped to distributed workers; undefined slots cannot be serialized)
    seed_params = Vector{Any}(undef, nseeds)
    fill!(seed_params, NamedTuple())
    if hoisted
        rng = Random.Xoshiro(seed)
        for s in 1:nseeds
            seed_params[s] = if first_spec.name == "smseed" && s > 1
                seed_params[1]
            else
                find_mle_advanced(likelihood, priors, params; strategy = [strategy[1]], rng, objective, null_prior, budget, local_defaults).params
            end
        end
    end
    (; scanpoints, mesh_flat, scan_keys = collect(Symbol.(scan_keys)), params, priors, strategy, rest_strategy, hoisted, seed_params, nseeds, seed,
       objective, null_prior, budget, local_defaults, starts = collect(starts), warm_strategy = warm, shape = size(scanpoints))
end

"""
    _profile_point_fit(likelihood, plan, i; warm=nothing) -> record of find_mle_advanced

The fit at grid point `i`. Without `warm` the whole strategy is run (a hoisted seed, if any, replaces the first stage), the best of
`nseeds` attempts by `log_post` is kept; each point has its own RNG stream `Xoshiro(plan.seed + i)`, so results do not depend on
the scheduling order. With `warm` (the fitted parameters of a neighbouring point or an interpolated seed) only
`plan.warm_strategy` is run, started from `warm` with the scanned parameters set to this point's values.
"""
function _profile_point_fit(likelihood, plan, i; warm = nothing)
    sp = plan.scanpoints[i]
    fixed = NamedTuple{Tuple(plan.scan_keys)}(Tuple(Float64.(plan.mesh_flat[i])))
    rng = Random.Xoshiro(plan.seed + i)
    kw = (; rng, objective = plan.objective, null_prior = plan.null_prior, budget = plan.budget, local_defaults = plan.local_defaults)
    if warm !== nothing
        start = _fix_conditioned(merge(warm, fixed), sp)
        return find_mle_advanced(likelihood, sp, start; strategy = plan.warm_strategy, kw...)
    end
    # start-value overrides of multi-start fits (start values of the point's scanned parameters do not matter: dropped)
    free(st) = NamedTuple(k => v for (k, v) in pairs(st) if !(sp[k] isa ValueShapes.ConstValueDist))
    overrides = unique(free.(plan.starts))
    explicit = length(overrides) > 1 || !isempty(first(overrides))
    if plan.hoisted
        starts = [merge(plan.seed_params[s], o, fixed) for s in 1:plan.nseeds for o in overrides]
        find_mle_advanced(likelihood, sp, plan.params; strategy = plan.rest_strategy, starts, kw...)
    elseif explicit
        starts = [merge(plan.params, o, fixed) for _ in 1:plan.nseeds for o in overrides]
        find_mle_advanced(likelihood, sp, plan.params; strategy = plan.strategy, starts, kw...)
    else
        find_mle_advanced(likelihood, sp, merge(plan.params, fixed); strategy = plan.strategy, nseeds = plan.nseeds, kw...)
    end
end

# the optimizer run data of all points, as arrays shaped like the grid (stored in `meta["optimizer_runs"]`)
function _optimizer_runs(recs, shape)
    arr(f) = reshape([f(r) for r in recs], shape)
    Dict{String, Any}(
        "converged"    => arr(r -> r.converged),
        "n_iterations" => arr(r -> r.n_iters),
        "time"         => arr(r -> r.total_time),     # all stages and attempts
        "n_calls"      => arr(r -> r.total_calls),
        "n_attempts"   => arr(r -> r.n_attempts),
        "strategy"     => arr(r -> r.strategy),
        "stages"       => arr(r -> [(; label = s.label, llh = s.llh, log_post = s.log_post, converged = s.converged, n_iters = s.n_iters, n_calls = s.n_calls, time = s.time)
                                    for s in r.stages]),
        "trace_x"      => arr(r -> r.trace_x),         # objective evaluations, best objective so far (-score), stage index
        "trace_y"      => arr(r -> r.trace_y),
        "trace_stage"  => arr(r -> r.trace_stage),
    )
end

"""
    _profile_result(plan, recs, exec_time; axes, meta=Dict()) -> NewtrinosResult

Assemble the records of all grid points (`recs`, in the order of `plan.scanpoints`) like `Newtrinos.profile` does (one array per
parameter plus `llh` and `log_posterior`) and attach the settings and the optimizer run data (`meta["optimizer_runs"]`).
"""
function _profile_result(plan, recs, exec_time; axes, meta = Dict())
    values = Newtrinos.assemble_profile_results([(r.llh, r.log_post, r.params) for r in recs], plan.shape)
    m = Dict{String, Any}("task" => "profile", "priors" => plan.priors, "params" => plan.params, "exec_time" => exec_time,
                          "strategy" => plan.strategy, "nseeds" => plan.nseeds, "objective" => plan.objective, "seed" => plan.seed,
                          "starts" => plan.starts, "optimizer_runs" => _optimizer_runs(recs, plan.shape))
    merge!(m, meta)
    Newtrinos.add_meta!(m)
    NewtrinosResult(axes = axes, values = values, meta = m)
end

# start values of every grid point interpolated from earlier profile results (see `make_seed_interpolants`)
function _warm_params_from_results(vars_to_scan, mesh_flat, params, seed_result, seed_results)
    n_axes = length(keys(vars_to_scan))
    coords(i) = Tuple(Float64(mesh_flat[i][d]) for d in 1:n_axes)
    if !isnothing(seed_result)
        itps = Newtrinos.make_seed_interpolants(seed_result)
        ks = collect(keys(itps))
        @info "Using per-point seeds from smoothed result ($(length(ks)) params interpolated)"
        return [merge(params, NamedTuple{Tuple(ks)}(Tuple(itps[k](coords(i)...) for k in ks))) for i in eachindex(mesh_flat)]
    elseif !isnothing(seed_results)
        itps = Newtrinos.make_merged_seed_interpolants(seed_results, keys(vars_to_scan))
        ks = collect(keys(itps))
        @info "Using merged per-point seeds from $(length(seed_results)) files ($(length(ks)) params interpolated)"
        return [merge(params, NamedTuple{Tuple(ks)}(Tuple(itps[k](coords(i)...) for k in ks))) for i in eachindex(mesh_flat)]
    end
    nothing
end

"""
    profile_advanced(likelihood, priors, vars_to_scan, params; strategy, null_prior=nothing, nseeds=1, starts=[(;)],
                     seed_result=nothing, seed_results=nothing, sequential=false, start_from=nothing, warm_strategy=nothing,
                     map_func=nothing, objective="posterior", seed=1234, budget=nothing, local_defaults=nothing) -> NewtrinosResult

Profile likelihood scan like `Newtrinos.profile` (grid from `vars_to_scan`, all other free parameters optimized at every point), but
the fit at each point is `find_mle_advanced` with `strategy` (a vector of stage specs, e.g. `["smseed", "lbfgsb"]` or
`["optim-pso:particles=20,iters=30", "lbfgsb"]`). The returned `NewtrinosResult` has the usual `axes`/`values` and the settings plus
the optimizer run data of every point in `meta` (`meta["optimizer_runs"]`: converged, iterations, time and calls of all stages and
attempts, the stage records and the combined trace, as arrays shaped like the grid).

- `strategy`, `nseeds`, `null_prior` (needed by `smseed`), `budget`, `objective`, `seed`: as in `find_mle_advanced` / `_profile_plan`. A seed or alias
  first stage (`random`, `smseed`, `mle`, ...) is computed once for the whole scan; a global first stage is run at every point.
- `starts`: start-value overrides, e.g. `[(δCP = π/4,), (δCP = 3π/4,)]`: one attempt per override, the best fit is kept.
- Warm starting: `sequential=true` fits the points in BFS order from `start_from` (a NamedTuple of axis values; default the central
  point), each point started from its already fitted neighbour; `seed_result`/`seed_results` (a `NewtrinosResult` or several of earlier
  profiles) start every point from interpolated best-fit parameters. Warm-started points run only `warm_strategy` (default: the
  strategy's local stages, since a seed or global stage would ignore the warm start); the root point of a sequential scan without
  seed results runs the full strategy.
- `map_func`: custom map for parallelism (e.g. `pmap`); default `Threads.@threads`. `local_defaults` fills settings of the local stages
  (default `(; ad = :forward)` when threaded). Sequential scans run in this process.
"""
function profile_advanced(likelihood, priors, vars_to_scan, params; strategy, null_prior = nothing, nseeds = 1, starts = [(;)],
                          seed_result = nothing, seed_results = nothing, sequential = false, start_from = nothing, warm_strategy = nothing,
                          map_func = nothing, objective = "posterior", seed = 1234, budget = nothing, local_defaults = nothing)
    t1 = time()
    if all([isa(priors[var], Number) for var in setdiff(keys(priors), keys(vars_to_scan))])
        return Newtrinos.scan(likelihood, priors, vars_to_scan, params)
    end
    axis_values, scanpoints = Newtrinos.generate_scanpoints(vars_to_scan, priors)
    mesh_flat = vec(Newtrinos._generate_grid(vars_to_scan, priors)[3])   # grid points in the order of `scanpoints`
    threaded = isnothing(map_func) && Threads.nthreads() > 1
    ld = isnothing(local_defaults) ? (threaded ? (; ad = :forward) : (;)) : local_defaults
    plan = _profile_plan(likelihood, priors, params, scanpoints, mesh_flat, collect(keys(vars_to_scan)); strategy, nseeds, seed, objective,
                         null_prior, budget, local_defaults = ld, starts, warm_strategy)
    warm_params = _warm_params_from_results(vars_to_scan, mesh_flat, params, seed_result, seed_results)
    recs = Vector{Any}(undef, length(scanpoints))

    if sequential
        isnothing(map_func) || @warn "sequential=true ignores map_func; running in this process"
        grid_shape = size(scanpoints)
        center = div(prod(grid_shape), 2) + 1
        if !isnothing(start_from)
            ks = collect(keys(vars_to_scan))
            idx = [argmin(abs.(axis_values[findfirst(==(k), ks)] .- Float64(start_from[k]))) for k in ks]
            center = LinearIndices(grid_shape)[CartesianIndex(Tuple(idx))]
            @info "Sequential profile: starting from grid indices $idx (nearest to start_from)"
        else
            @info "Sequential profile: starting from the central grid point (linear index $center)"
        end
        order, parent = Newtrinos._sequential_order(grid_shape, center)
        prog = ProgressMeter.Progress(length(order))
        for i in order
            w = if parent[i] > 0
                pp = recs[parent[i]].params
                all(v -> isfinite(Float64(v)), Tuple(pp)) ? pp : nothing   # a failed neighbour is no start
            else
                isnothing(warm_params) ? nothing : warm_params[i]
            end
            recs[i] = _profile_point_fit(likelihood, plan, i; warm = w)
            ProgressMeter.next!(prog)
        end
    else
        work(i) = _profile_point_fit(likelihood, plan, i; warm = isnothing(warm_params) ? nothing : warm_params[i])
        if isnothing(map_func)
            prog = ProgressMeter.Progress(length(scanpoints))
            Threads.@threads for i in eachindex(scanpoints)
                recs[i] = work(i)
                ProgressMeter.next!(prog)
            end
        else
            recs .= collect(map_func(work, collect(eachindex(scanpoints))))
        end
    end
    meta = Dict{String, Any}("vars_to_scan" => vars_to_scan, "sequential" => sequential, "warm_strategy" => plan.warm_strategy)
    _profile_result(plan, recs, time() - t1; axes = NamedTuple{Tuple(keys(vars_to_scan))}(axis_values), meta)
end
