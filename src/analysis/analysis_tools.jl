using LinearAlgebra
using TOML
using Distributions
using DensityInterface
using InverseFunctions
import ForwardDiff
import PolyesterForwardDiff
using BAT
using Optimization
using OptimizationLBFGSB
import Optim
using IterTools
using DataStructures
using ADTypes
using AutoDiffOperators
using ContentHashes
using ValueShapes
using FileIO
using FillArrays
import JLD2
using MeasureBase
using FunctionChains
using Accessors
using Logging
using ProgressMeter
using Dates
using LibGit2
import Mooncake
import Random
using Interpolations
using Statistics
using ..Newtrinos

const AD_BACKEND = Ref{Symbol}(:auto)

"""
    set_ad_backend(backend::Symbol)

Set the global AD backend. Valid values: `:auto`, `:forwarddiff`, `:polyester`, `:mooncake`.
"""
function set_ad_backend(backend::Symbol)
    backend in (:auto, :forwarddiff, :polyester, :mooncake) || error("Unknown AD backend: $backend. Choose from: auto, forwarddiff, polyester, mooncake")
    AD_BACKEND[] = backend
end

"""
    select_ad(n_params; threshold=12)

Choose AD backend for gradient-based optimization based on the global setting
(see [`set_ad_backend`](@ref)).

- `:auto` — ForwardDiff for ≤ threshold params, PolyesterForwardDiff otherwise
- `:forwarddiff` — standard chunked ForwardDiff
- `:polyester` — threaded chunked ForwardDiff (PolyesterForwardDiff)
- `:mooncake` — reverse-mode AD via Mooncake (constant overhead, lower memory)
"""
function select_ad(n_params; threshold=12)
    backend = AD_BACKEND[]
    if backend == :auto
        return n_params > threshold ? AutoPolyesterForwardDiff() : AutoForwardDiff()
    elseif backend == :forwarddiff
        return AutoForwardDiff()
    elseif backend == :polyester
        return AutoPolyesterForwardDiff()
    elseif backend == :mooncake
        return AutoMooncake(Mooncake.Config())
    end
end

# AD backend is set via set_ad_backend() + set_batcontext(ad = select_ad(n))
# Do NOT set a hardcoded default here — it overrides user selection

adsel = AutoForwardDiff()
set_batcontext(ad = adsel)

# ── Core Types ──────────────────────────────────────────────────────

"""
    NewtrinosResult(; axes, values, meta=Dict())

Container for the results of a [`scan`](@ref) or [`profile`](@ref) run.

# Fields
- `axes::NamedTuple`: scan grid coordinates, one entry per scanned parameter.
  Each entry is a vector of grid-point values.
- `values::NamedTuple`: per-grid-point outputs — likelihood values (`llh`,
  `log_posterior`) and optimized nuisance parameter arrays.
- `meta::Dict`: execution metadata (hostname, date, git commit, timing).
  Populated by [`add_meta!`](@ref).
"""
@kwdef struct NewtrinosResult
    axes::NamedTuple
    values::NamedTuple
    meta::Dict = Dict()
end

# ── Utility Functions ───────────────────────────────────────────────

"""
    sort_nt(nt::NamedTuple) -> NamedTuple

Return a copy of a named tuple `nt` with keys sorted alphabetically.

Used internally by [`safe_merge`](@ref), [`get_params`](@ref), and
[`get_priors`](@ref) to ensure a canonical key ordering across modules.

# Arguments
- `nt::NamedTuple`: input NamedTuple.

# Returns
A NamedTuple with the same entries as `nt`, sorted alphabetically by key name.
"""
function sort_nt(nt::NamedTuple)
    keys_sorted = sort(collect(keys(nt)))
    values_sorted = getindex.(Ref(nt), keys_sorted)
    return NamedTuple{Tuple(keys_sorted)}(values_sorted)
end

"""
    safe_merge(nt_list::NamedTuple...) -> NamedTuple

Merge NamedTuples, asserting that any key shared by two or more inputs has the
same value in all of them.

Raises an `error` if a key conflict is detected. The result is sorted
alphabetically by key.

# Arguments
- `nt_list::NamedTuple...`: one or more NamedTuples to merge.

# Returns
A sorted NamedTuple containing all key-value pairs from the inputs.

# Examples
```julia
a = (x = 1.0, y = 2.0)
b = (y = 2.0, z = 3.0)
safe_merge(a, b)   # (x = 1.0, y = 2.0, z = 3.0)
```
"""
function safe_merge(nt_list::NamedTuple...)
    merged = NamedTuple()
    for nt in nt_list
        for (k, v) in pairs(nt)
            if haskey(merged, k)
                if merged[k] != v
                    error("Conflict on key '$k': $(merged[k]) ≠ $v")
                end
            end
        end
        merged = merge(merged, nt)
    end
    sort_nt(merged)
end

# ── Wrapper Type ─────────────────────────────────────────────────────
# Defined before accessors since get_params/get_priors dispatch on it

"""
    Wrapper <: Newtrinos.Experiment

Wraps a [`Newtrinos.Experiment`](@ref) with parameter name aliases.

When two experiments share a physics parameter under different names (e.g.
`theta13` vs `sin2theta13`), a `Wrapper` renames parameters seen by the
outer analysis while transparently translating back to the original names
before calling the inner experiment's `forward_model` and `plot` functions.

Construct with `Wrapper(experiment, aliases)` where `aliases` is a
`Dict{Symbol,Symbol}` mapping original parameter names to their aliases.

# Fields
- `x::Newtrinos.Experiment`: the wrapped experiment.
- `aliases::Dict{Symbol,Symbol}`: map from original name → alias.
- `translated_keys::Vector{Symbol}`: all parameter keys after alias substitution.
- `reverse_lookup::Dict{Symbol,Symbol}`: map from alias → original name.

# Examples
```julia
aliases = Dict(:dm21_sq => :Dm21)
wrapped = Wrapper(dayabay_experiment, aliases)
```
"""
struct Wrapper <: Newtrinos.Experiment
    x::Newtrinos.Experiment
    aliases::Dict{Symbol, Symbol}  # actual -> alias
    translated_keys::Vector{Symbol}
    reverse_lookup::Dict{Symbol, Symbol} # alias -> actual
end

"""
    Base.getproperty(wrapper::Wrapper, name::Symbol)

Property accessor for [`Wrapper`](@ref) that intercepts `forward_model` and
`plot` to inject parameter name translation.

When `name` is `:forward_model` or `:plot`, returns a closure that translates
aliased parameter names back to their original names before delegating to the
inner experiment. All other property accesses are forwarded directly to the
wrapped experiment `wrapper.x`.

# Arguments
- `wrapper::Wrapper`: the wrapping experiment.
- `name::Symbol`: the property name to access.

# Returns
The requested property value of `name`, with translation logic applied for
`:forward_model` and `:plot`.
"""
function Base.getproperty(wrapper::Wrapper, name::Symbol)
    if name ∈ (:x, :aliases, :translated_keys, :reverse_lookup)
        return getfield(wrapper, name)
    end
    if name == :forward_model
        function forward_model(params)
            own_params = NamedTuple{Tuple(wrapper.translated_keys)}(params)
            orig_param_names = Tuple([get(wrapper.reverse_lookup, k, k) for k in keys(own_params)])
            orig_params = NamedTuple{orig_param_names}(values(own_params))
            return wrapper.x.forward_model(orig_params)
        end
        return forward_model
    end
    if name == :plot
        function plot(params, data=wrapper.x.assets.observed)
            own_params = NamedTuple{Tuple(wrapper.translated_keys)}(params)
            orig_param_names = Tuple([get(wrapper.reverse_lookup, k, k) for k in keys(own_params)])
            orig_params = NamedTuple{orig_param_names}(values(own_params))
            return wrapper.x.plot(orig_params, data)
        end
        return plot
    end
    return getfield(wrapper.x, name)
end

# ── Parameter & Prior Accessors ─────────────────────────────────────

"""
    get_params(x; start = :nominal) -> NamedTuple

Extract parameter values as a sorted NamedTuple.

Dispatches on the type of `x`:
- `Physics`: returns `x.params` sorted by key.
- `Experiment`: merges `x.params` with the parameters of its `physics` module.
- `Wrapper`: applies the alias translation from [`Wrapper`](@ref) to the keys.
- `NamedTuple` of modules: merges parameters from all modules via
  [`safe_merge`](@ref), raising an error on conflicts.

# Arguments
- `x`: a [`Newtrinos.Physics`](@ref), [`Newtrinos.Experiment`](@ref),
  [`Wrapper`](@ref), or `NamedTuple` of such objects.
- `start`: `:nominal` (default) for the nominal values, or `:best_fit` to replace each experiment's own parameters
  (`x.params`, never those of its physics modules, e.g. the oscillation parameters) by the values stored for the
  experiment's best fit in its default setup (see [`stored_best_fit`](@ref)). A good start point for joint fits;
  parameters without a stored value keep their nominal value.

# Returns
A sorted `NamedTuple` of parameter name → value.
"""
function get_params(x::Newtrinos.Physics; start::Symbol = :nominal)
    sort_nt(x.params)
end

function get_params(x::Newtrinos.Experiment; start::Symbol = :nominal)
    own = x.params
    if start === :best_fit
        bf = stored_best_fit(x)
        own = merge(own, NamedTuple(k => convert_like(own[k], bf[k]) for k in keys(own) if haskey(bf, k)))
    elseif start !== :nominal
        throw(ArgumentError("start must be :nominal or :best_fit, got :$start"))
    end
    safe_merge(own, get_params(x.physics))
end

function get_params(w::Newtrinos.Wrapper; start::Symbol = :nominal)
    NamedTuple{Tuple(w.translated_keys)}(values(get_params(w.x; start)))
end

function get_params(modules::NamedTuple; start::Symbol = :nominal)
    all_params = [get_params(m; start) for m in modules]
    safe_merge(all_params...)
end

convert_like(nominal::Number, v) = convert(typeof(nominal), v)
convert_like(nominal::AbstractArray, v) = convert(typeof(nominal), reshape(collect(v), size(nominal)))

"""
    module_dir(m::Module) -> String

Source directory of the experiment module `m` (the directory of the file defining its `configure`), resolved
relative to the installed Newtrinos package; `nothing` if `m` has no `configure`.
"""
function module_dir(m::Module)
    isdefined(m, :configure) || return nothing
    file = replace(string(first(methods(getproperty(m, :configure))).file), '\\' => '/')
    i = findlast("/src/", file)
    isnothing(i) && return dirname(file)
    dirname(joinpath(pkgdir(Newtrinos), file[first(i)+1:end]))
end

"""
    stored_best_fit(x::Experiment) -> NamedTuple

Best-fit values of the experiment's own parameters in its default setup, read from `bestfit.toml` next to the module
defining the experiment type (see [`module_dir`](@ref)); empty if there is no such file. The files are written by
`src/utils/generate_bestfit.jl`.
"""
function stored_best_fit(x::Newtrinos.Experiment)
    m = parentmodule(typeof(x))
    dir = module_dir(m)
    isnothing(dir) && return NamedTuple()
    f = joinpath(dir, "bestfit.toml")
    isfile(f) || return NamedTuple()
    toml = TOML.parsefile(f)
    name = get(get(toml, "meta", Dict()), "experiment", string(nameof(m)))
    name == string(nameof(m)) || error("$f was written for experiment `$name`, not `$(nameof(m))`")
    d = toml["params"]
    NamedTuple(Symbol(k) => v for (k, v) in d)
end

"""
    get_priors(x) -> NamedTuple

Extract prior distributions as a sorted NamedTuple.

Dispatches on the type of `x`:
- `Physics`: returns `x.priors` sorted by key.
- `Experiment`: merges `x.priors` with priors from its `physics` module.
- `Wrapper`: applies the alias translation from [`Wrapper`](@ref) to the keys.
- `NamedTuple` of modules: merges priors from all modules via
  [`safe_merge`](@ref), raising an error on conflicts.

# Arguments
- `x`: a [`Newtrinos.Physics`](@ref), [`Newtrinos.Experiment`](@ref),
  [`Wrapper`](@ref), or `NamedTuple` of such objects.

# Returns
A sorted `NamedTuple` of parameter name → prior distribution.
"""
function get_priors(x::Newtrinos.Physics)
    sort_nt(x.priors)
end

function get_priors(x::Newtrinos.Experiment)
    safe_merge(x.priors, get_priors(x.physics))
end

function get_priors(w::Newtrinos.Wrapper)
    NamedTuple{Tuple(w.translated_keys)}(values(get_priors(w.x)))
end

function get_priors(modules::NamedTuple)
    all_priors = [get_priors(m) for m in modules]
    safe_merge(all_priors...)
end

"""
    condition(priors, conditional_vars, params) -> NamedTuple

Fix a subset of parameters to specific values by replacing their prior
distributions with constants (`ConstValueDist`).

Two dispatch methods:

**Array form** (conditional_vars is an array) — fixes each symbol in `conditional_vars` to its value in `params`:
```julia
condition(priors, [:theta23, :dm31_sq], params)
```

**Dict form** (conditional_vars is a dict) — fixes each key to the paired value, or to `params[key]` when the
value is `nothing`:
```julia
condition(priors, Dict(:theta23 => pi/4, :dm31_sq => nothing), params)
```

# Arguments
- `priors::NamedTuple`: prior distributions to modify.
- `conditional_vars`: `AbstractArray{Symbol}` or `AbstractDict{Symbol}` of
  parameters to fix.
- `params::NamedTuple`: fallback parameter values used when a dict entry is
  `nothing` or when using the array form.

# Returns
The modified `priors` NamedTuple with the specified parameters frozen.
"""
function condition(priors::NamedTuple, conditional_vars::AbstractArray, p)
    for var in conditional_vars
        @reset priors[var] = p[var]
    end
    priors
end

function condition(priors::NamedTuple, conditional_vars::AbstractDict, p)
    for var in keys(conditional_vars)
        if isnothing(conditional_vars[var])
            @reset priors[var] = p[var]
        else
            @reset priors[var] = conditional_vars[var]
        end
    end
    priors
end

"""
    Wrapper(x::Newtrinos.Experiment, aliases::Dict{Symbol,Symbol}) -> Wrapper

Construct a [`Wrapper`](@ref) around `x` with the given parameter name aliases.

Derives `translated_keys` by applying `aliases` to the full key list returned
by `get_params(x)`, and builds the `reverse_lookup` for inverse translation.

# Arguments
- `x::Newtrinos.Experiment`: experiment to wrap.
- `aliases::Dict{Symbol,Symbol}`: map from original parameter name → alias.

# Returns
A [`Wrapper`](@ref) instance.
"""
function Wrapper(x::Newtrinos.Experiment, aliases::Dict{Symbol, Symbol})
    original_keys = keys(get_params(x))
    translated_keys = [get(aliases, k, k) for k in original_keys]
    reverse_lookup = Dict(value => key for (key, value) in aliases)
    return Wrapper(x, aliases, translated_keys, reverse_lookup)
end

# ── Experiment Utilities ────────────────────────────────────────────

"""
    get_observed(experiments::NamedTuple) -> NamedTuple

Extract the observed data from each experiment in a NamedTuple of experiments.

# Arguments
- `experiments::NamedTuple`: named collection of [`Newtrinos.Experiment`](@ref) objects.

# Returns
A NamedTuple with the same keys as `experiments`, each entry holding the
corresponding `experiment.assets.observed` data.
"""
function get_observed(experiments::NamedTuple)
    NamedTuple{keys(experiments)}(e.assets.observed for e in experiments)
end

"""
    get_fwd_model(experiments::NamedTuple) -> Function

Compose the forward models of all experiments into a single joint model.

Builds a `NamedTuple` of per-experiment forward models and composes them
with `ffanout` and `distprod` so that calling the result with a parameter
NamedTuple returns a joint product distribution over all experiments.

# Arguments
- `experiments::NamedTuple`: named collection of [`Newtrinos.Experiment`](@ref) objects.

# Returns
A callable `params -> joint_distribution`.
"""
function get_fwd_model(experiments::NamedTuple)
    fwd_models = NamedTuple{keys(experiments)}(e.forward_model for e in experiments)
    distprod ∘ ffanout(fwd_models)
end

"""
    generate_likelihood(experiments[, observed]) -> likelihood

Construct a joint likelihood from a NamedTuple of configured experiments.

Combines all experiment forward models into a product distribution and
pairs it with the observed data. The returned object supports
`DensityInterface.logdensityof(likelihood, params)`.

# Arguments
- `experiments::NamedTuple`: named collection of configured
  [`Newtrinos.Experiment`](@ref) objects.
- `observed` (optional): NamedTuple of observed data arrays matching
  `experiments`. Defaults to the observed data stored in each experiment's
  `assets.observed` field.

# Returns
A likelihood object compatible with `DensityInterface.logdensityof`.

# Examples
```julia
experiments = (deepcore = Newtrinos.deepcore.configure(),
               dayabay  = Newtrinos.dayabay.configure())
llh     = generate_likelihood(experiments)
logdensityof(llh, params)
```
"""
function generate_likelihood(experiments::NamedTuple, observed=get_observed(experiments))
    likelihoodof(get_fwd_model(experiments), observed)
end

"""
    correlated_priors_vars(priors, vars, dist) -> (corr_prior, other_prior)

Replace the independent priors for `vars` with a single correlated
multivariate distribution, splitting the prior NamedTuple into two parts.

Useful when a set of parameters has a known covariance structure (e.g. from
an external fit) that cannot be captured by independent marginals.

# Arguments
- `priors::NamedTuple`: the full prior NamedTuple.
- `vars::Union{AbstractArray, Tuple}`: parameter names to correlate.
- `dist::Distributions.Distribution`: a multivariate distribution over `vars`
  (must match the length of `vars`).

# Returns
A 2-tuple `(corr_prior, other_prior)` where:
- `corr_prior`: a callable returning a `ReshapedDist` over `vars`.
- `other_prior`: a `distprod` of the remaining (uncorrelated) priors.
"""
function correlated_priors_vars(priors::NamedTuple, vars::Union{AbstractArray, Tuple}, dist::Distribution)
    named_shapes = NamedTuple(var => ValueShapes.ScalarShape{Real}() for var in vars)
    corr_prior = Returns(ValueShapes.ReshapedDist(dist, ValueShapes.NamedTupleShape(named_shapes)))
    keys_to_keep = Tuple(k for k in keys(priors) if k ∉ vars)
    other_prior = distprod(;NamedTuple{keys_to_keep}(priors)...)
    return corr_prior, other_prior
end

"""
    generate_toy_data(experiment, params) -> data
    generate_toy_data(experiments::NamedTuple, params) -> NamedTuple

Generate random toy data by sampling from the experiment forward model.

Calls `experiment.forward_model(params)` to obtain the predictive
distribution and draws one sample from it.

# Arguments
- `experiment`: a single [`Newtrinos.Experiment`](@ref) or a NamedTuple of them.
- `params::NamedTuple`: parameter values at which to evaluate the forward model.

# Returns
- Single-experiment form: one sample from the forward model distribution.
- Multi-experiment form: a NamedTuple with the same keys as `experiments`,
  each entry holding one sample.
"""
function generate_toy_data(experiment::Newtrinos.Experiment, params::NamedTuple)
    dist_obj = experiment.forward_model(params)
    rand(dist_obj)
end

function generate_toy_data(experiments::NamedTuple, params::NamedTuple)
    map(experiments) do experiment
        dist_obj = experiment.forward_model(params)
        rand(dist_obj)
    end
end

"""
    generate_asimov_data(experiment, params) -> data
    generate_asimov_data(experiments::NamedTuple, params) -> NamedTuple

Generate Asimov (expected-value) data from the experiment forward model.

Evaluates `mean(forward_model(params))`. For Poisson-distributed observables
the mean is rounded to the nearest integer so that the data type matches what
a real run would produce.

# Arguments
- `experiment`: a single [`Newtrinos.Experiment`](@ref) or a NamedTuple of them.
- `params::NamedTuple`: parameter values at which to evaluate the forward model.

# Returns
- Single-experiment form: the expected-value data array (or NamedTuple of arrays).
- Multi-experiment form: a NamedTuple with the same keys as `experiments`.
"""
function generate_asimov_data(experiment::Newtrinos.Experiment, params::NamedTuple)
    dist_obj = experiment.forward_model(params)
    asimov_data_flt = mean(dist_obj)
    check_dist(d) = (d isa Distributions.Poisson) |
        (d isa Distributions.ProductDistribution && !isempty(d.dists) && first(d.dists) isa Distributions.Poisson) |
        (d isa Distributions.Product && !isempty(d.v) && first(d.v) isa Distributions.Poisson)

    if dist_obj isa ValueShapes.NamedTupleDist
        for key in keys(dist_obj)
            if check_dist(dist_obj[key])
                @info "Poisson-based model for $(key). Rounding Asimov data to nearest integer."
                @reset asimov_data_flt[key] = round.(Int, asimov_data_flt[key])
            end
        end
        return asimov_data_flt
    end
    if check_dist(dist_obj)
        @info "Poisson-based model. Rounding Asimov data to nearest integer."
        return round.(Int, asimov_data_flt)
    end

    @info "Not Poisson-based model. Returning std floating-point Asimov data."
    return asimov_data_flt
end

function generate_asimov_data(experiments::NamedTuple, params::NamedTuple)
    map(experiments) do experiment
        generate_asimov_data(experiment, params)
    end
end

# ── Optimization ────────────────────────────────────────────────────

"""
    find_mle(likelihood, prior, params; adsel=AutoPolyesterForwardDiff())
        -> (llh, log_posterior, result)

Find the Maximum Likelihood Estimator (MLE) using L-BFGS optimization via BAT.jl.

Parameters whose prior is a `ConstValueDist` are held fixed during optimization.
If the optimizer raises an `ArgumentError` (e.g. due to a degenerate starting
point), returns `(NaN, NaN, NaN-filled NamedTuple)` instead of propagating the
error.

# Arguments
- `likelihood`: a likelihood object supporting `DensityInterface.logdensityof`.
- `prior`: a `distprod`-style prior NamedTuple or distribution.
- `params::NamedTuple`: starting parameter values for the optimizer.
- `adsel`: AD backend selector (default: `AutoPolyesterForwardDiff()`).

# Returns
A 3-tuple `(llh, log_posterior, result)` where:
- `llh::Float64`: log-likelihood at the optimum.
- `log_posterior::Float64`: log-posterior at the optimum.
- `result::NamedTuple`: optimized parameter values.
"""
function find_mle(likelihood, prior, params; adsel = AutoPolyesterForwardDiff())
    try
        set_batcontext(ad = adsel)
        posterior = PosteriorMeasure(likelihood, prior)

        msg = "Running Optimization for point "
        for key in keys(prior)
            if prior[key] isa ValueShapes.ConstValueDist
                value = prior[key].value
                @reset params[key] = value
                msg *= " $(key): $(value)"
            end
        end

        @info msg
        res = bat_findmode(posterior, TransformedMaxDensity(optalg=OptimizationAlg(optalg=OptimizationLBFGSB.LBFGSB(), kwargs = (reltol=1e-7, maxiters=1000)), init = ExplicitInit([params])))

        return logdensityof(likelihood, res.result), logdensityof(posterior, res.result), res.result
    catch e
        # BAT v5's checked_logdensityof wraps a caught error in BAT.EvalException,
        # whose constructor requires an `AbstractMeasure` target; a plain `Likelihood`
        # (as produced by `likelihoodof`) doesn't qualify, so the wrapping itself
        # throws a MethodError instead of propagating the original ArgumentError
        # (see https://github.com/bat/BAT.jl, BAT.EvalException in v5.1.0).
        is_wrapped_argerror = e isa MethodError && e.f === BAT.EvalException &&
            length(e.args) >= 4 && e.args[4] isa ArgumentError
        if e isa ArgumentError || is_wrapped_argerror
            return NaN, NaN, (; (k => NaN for k in keys(params))... )
        else
            rethrow(e)
        end
    end
end

"""
    find_mle_ext(likelihood, prior, params; g_tol=1e-6, maxiters=2000)

Like [`find_mle`](@ref) but bypasses BAT entirely. Uses `Optim.Fminbox(LBFGS(m=20))`
with bounds derived directly from the prior distributions and ForwardDiff gradients.
More robust for combined fits where one experiment dominates the gradient.

Parameters with `ConstValueDist` (or a raw number) in the prior are treated as fixed
and excluded from the optimization. All other parameters are optimized with bounds
`[minimum(d), maximum(d)]` from their prior distribution (Uniform/Truncated give finite
bounds; Normal gives ±Inf, which Fminbox treats as unconstrained).
"""
function find_mle_ext(likelihood, prior, params; g_tol=1e-6, maxiters=2000)
    all_keys = tuple(keys(prior)...)
    n = length(all_keys)

    # Separate fixed (ConstValueDist or raw Number) from free params
    is_free_mask  = Vector{Bool}(undef, n)
    fixed_float   = Vector{Float64}(undef, n)  # sentinel 0.0 for free slots
    free_x_index  = Vector{Int}(undef, n)       # free_x_index[i] = position in x vector (0 if fixed)
    free_dists    = Any[]
    lo_vec        = Float64[]
    hi_vec        = Float64[]
    x0_vec        = Float64[]
    n_free        = 0

    for (i, k) in enumerate(all_keys)
        d = prior[k]
        if d isa ValueShapes.ConstValueDist
            is_free_mask[i]  = false
            fixed_float[i]   = Float64(d.value)
            free_x_index[i]  = 0
        elseif d isa Number
            is_free_mask[i]  = false
            fixed_float[i]   = Float64(d)
            free_x_index[i]  = 0
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

    function obj(x::AbstractVector{T}) where T
        all_x = [is_free_mask[i] ? x[free_x_index[i]] : T(fixed_float[i]) for i in 1:n]
        full_params = NamedTuple{all_keys}(Tuple(all_x))
        llh_val   = logdensityof(likelihood, full_params)
        prior_val = sum(logpdf(free_dists[j], x[j]) for j in 1:n_free)
        return -(llh_val + prior_val)
    end

    msg = "Running Optimization (ext) for point"
    for (i, k) in enumerate(all_keys)
        is_free_mask[i] || (msg *= " $(k): $(fixed_float[i])")
    end
    @info msg

    try
        result = Optim.optimize(
            obj,
            (G, x) -> ForwardDiff.gradient!(G, obj, x),
            lo_vec, hi_vec, x0_vec,
            Optim.Fminbox(Optim.LBFGS(m=20)),
            Optim.Options(g_tol=g_tol, outer_iterations=maxiters)
        )

        x_opt      = Optim.minimizer(result)
        opt_params = NamedTuple{all_keys}(Tuple(
            is_free_mask[i] ? x_opt[free_x_index[i]] : fixed_float[i] for i in 1:n
        ))

        converged = Optim.converged(result)
        converged || @warn "Optimizer (ext) did not converge (iters=$(Optim.iterations(result)), g_residual=$(Optim.g_residual(result)))"

        llh      = logdensityof(likelihood, opt_params)
        log_post = -(Optim.minimum(result))
        return llh, log_post, opt_params, converged
    catch e
        if e isa ArgumentError
            nan_params = NamedTuple{all_keys}(ntuple(_ -> NaN, n))
            return NaN, NaN, nan_params, false
        else
            rethrow(e)
        end
    end
end

# Cache key of a fit. ContentHashes cannot serialise `Char`, which every `MvNormal` prior contains (the `uplo` field of its
# Cholesky factor, e.g. in the correlated solar-flux prior `solar_norms`): such priors are hashed with their multivariate
# normals replaced by mean and covariance. Priors without them keep their previous keys.
_hash_key(prior) = prior
function _hash_key(prior::ValueShapes.NamedTupleDist)
    any(d -> d isa Distributions.AbstractMvNormal, values(prior)) || return prior
    NamedTuple{keys(prior)}(map(d -> d isa Distributions.AbstractMvNormal ? (:MvNormal, collect(mean(d)), Matrix(cov(d))) : d, values(prior)))
end

"""
    find_mle_cached(likelihood, prior, params, cache_dir)
        -> (llh, log_posterior, result)

Like [`find_mle`](@ref) but caches results to disk using content hashing.

A hash of `(prior, params)` is computed with `ContentHashes.hash`; if a
matching `.jld2` file exists in `cache_dir` the optimization is skipped and
the cached result is returned. Otherwise the result is computed and saved.

# Arguments
- `likelihood`: a likelihood object supporting `DensityInterface.logdensityof`.
- `prior`: prior distribution (same as [`find_mle`](@ref)).
- `params::NamedTuple`: starting parameter values.
- `cache_dir::Union{String,Nothing}`: directory for cached `.jld2` files.
  Pass `nothing` to disable caching.

# Returns
A 3-tuple `(llh, log_posterior, result)` — see [`find_mle`](@ref).
"""
function find_mle_cached(likelihood, prior, params, cache_dir; adsel = AutoPolyesterForwardDiff())
    opt_result = nothing

    h = ContentHashes.hash([_hash_key(prior), params])

    if !isnothing(cache_dir)
        fname = joinpath(cache_dir, "$h.jld2")
        if isfile(fname)
            @info "using cached file $fname"
            cached = FileIO.load(fname)
            opt_result = (cached["llh"], cached["log_posterior"], cached["result"])
        end
    end

    if isnothing(opt_result)
        opt_result = find_mle(likelihood, prior, params; adsel = adsel)
    end

    if !isnothing(cache_dir)
        fname = joinpath(cache_dir, "$h.jld2")
        FileIO.save(fname, OrderedDict("llh"=>opt_result[1], "log_posterior"=>opt_result[2], "result"=>opt_result[3]))
    end

    opt_result
end

# ── Scanning & Profiling ────────────────────────────────────────────

"""
    _generate_grid(vars_to_scan, priors) -> (vars, values, mesh)

Construct the Cartesian grid used by [`generate_scanpoints`](@ref) and
[`scan`](@ref).

Grid points for each variable are placed at evenly-spaced quantiles of its
prior distribution. All combinations are formed via `IterTools.product`.

# Arguments
- `vars_to_scan::OrderedDict{Symbol,Int}`: maps parameter names to grid sizes.
- `priors::NamedTuple`: prior distributions used to compute quantiles.

# Returns
A 3-tuple `(vars, values, mesh)` where:
- `vars`: `Vector{Symbol}` of scanned parameter names.
- `values`: vector of grid-point vectors (one per variable).
- `mesh`: Cartesian product array of grid-point tuples.

# Examples
```julia
    priors = (a=Uniform(0.0, 2.0), b=Uniform(-1.0, 1.0))
    vars_to_scan = OrderedDict(:a => 3, :b => 4)
    vars, values, mesh = Newtrinos._generate_grid(vars_to_scan, priors)
    # => mesh is a 3x4 grid of evenly spaced points from within the prior range of a and b 
```
"""
function _generate_grid(vars_to_scan, priors)
    vars = collect(keys(vars_to_scan))
    values = [quantile(priors[var], collect(range(0,1,vars_to_scan[var]))) for var in vars]
    mesh = collect(IterTools.product(values...))
    vars, values, mesh
end

"""
    generate_scanpoints(vars_to_scan, priors) -> (values, scanpoints)

Build a grid of fixed-parameter priors for a profile likelihood scan.

For each variable in `vars_to_scan`, grid points are placed at evenly-spaced
quantiles of its prior distribution. All combinations are formed via a
Cartesian product. At each grid point, the prior for the scanned variables
is replaced with a `ConstValueDist` fixing them to the grid values.

# Arguments
- `vars_to_scan::OrderedDict{Symbol,Int}`: maps each parameter name to the
  desired number of evenly spaced grid points within the prior range.
- `priors::NamedTuple`: prior distributions for all parameters.

# Returns
A 2-tuple `(values, scanpoints)` where:
- `values`: vector of grid-point vectors (one per scanned parameter).
- `scanpoints`: array of `distprod` priors, one per grid point.
"""
function generate_scanpoints(vars_to_scan, priors)
    vars, values, mesh = _generate_grid(vars_to_scan, priors)
    scanpoints = Array{Any}(undef, size(mesh))

    function make_prior(vals)
        p = deepcopy(priors)
        for i in 1:length(vars)
            @reset p[vars[i]] = vals[i]
        end
        distprod(;p...)
    end

    for i in eachindex(mesh)
        scanpoints[i] = make_prior(mesh[i])
    end

    values, scanpoints
end

"""
    assemble_profile_results(opt_results, result_size) -> NamedTuple

Collect the 'raw' per-grid-point optimization results and structure them into flat arrays within a NamedTuple.

Stacks the `(llh, log_posterior, params)` 3-tuples returned by
[`find_mle_cached`](@ref) into a NamedTuple of arrays with shape
`result_size`.

# Arguments
- `opt_results`: iterable of `(llh, log_posterior, NamedTuple)` 3-tuples.
- `result_size`: shape of the output arrays (matching the scan grid).

# Returns
A `NamedTuple` with one array-entry per parameter (values across the grid) plus
`:llh` and `:log_posterior` arrays.

# Example
```julia
opt_results = [
    (-10.0, -12.0, (a=1.0, b=2.0)), # llh, log_posterior, params output of find_mle
    (-5.0,  -7.0,  (a=3.0, b=4.0))
]
res = Newtrinos.assemble_profile_results(opt_results, (2,))
# returns res = (a = [1.0, 3.0], b = [2.0, 4.0], llh = [-10.0, -5.0], log_posterior = [-12.0, -7.0])
```
"""
function assemble_profile_results(opt_results, result_size)
    results = Array{Any}(undef, result_size)
    llhs = Array{Float64}(undef, result_size)
    log_posteriors = Array{Float64}(undef, result_size)
    for (i, opt_result) in enumerate(opt_results)
        llhs[i] = opt_result[1]
        log_posteriors[i] = opt_result[2]
        results[i] = opt_result[3]
    end
    s = OrderedDict(key=>[x[key] for x in results] for key in keys(first(results)))
    s[:llh] = llhs
    s[:log_posterior] = log_posteriors
    NamedTuple(s)
end

"""
    randomize_params(rng, params, prior)

Sample free parameters from their prior distributions to create randomized starting points.
Parameters fixed via `ConstValueDist` are left unchanged.
"""
function randomize_params(rng, params, prior)
    p = deepcopy(params)
    for key in keys(p)
        d = prior[key]
        if !(d isa ValueShapes.ConstValueDist)
            @reset p[key] = rand(rng, d)
        end
    end
    return p
end

"""
    _sequential_order(grid_shape, center_linear)

Compute a BFS traversal order of grid points (as linear indices) starting from `center_linear`,
walking to Manhattan-adjacent neighbours. Returns `(order, parent)`:
- `order`: vector of linear indices in visit order
- `parent`: for each linear index, the linear index of its BFS parent (-1 for the start)
"""
function _sequential_order(grid_shape, center_linear::Int)
    li = LinearIndices(grid_shape)
    ci = CartesianIndices(grid_shape)
    n  = prod(grid_shape)
    nd = length(grid_shape)

    visited = falses(n)
    order   = Int[]
    sizehint!(order, n)
    parent  = fill(-1, n)

    visited[center_linear] = true
    push!(order, center_linear)
    queue = [center_linear]
    head  = 1

    while head <= length(queue)
        curr    = queue[head]; head += 1
        curr_ci = ci[curr]
        for d in 1:nd, delta in (-1, 1)
            coords = ntuple(i -> i == d ? curr_ci[i] + delta : curr_ci[i], nd)
            all(1 .<= coords .<= grid_shape) || continue
            nb = li[CartesianIndex(coords)]
            visited[nb] && continue
            visited[nb]  = true
            parent[nb]   = curr
            push!(order, nb)
            push!(queue, nb)
        end
    end

    order, parent
end

# AD backend for the per-point fits of a profile scan. With the default threaded loop the threads
# are already busy with scan points, so a threaded AD backend (Polyester) inside each fit only adds
# contention; use plain ForwardDiff there. Otherwise (one thread, or a custom map_func such as pmap
# whose workers have their own threads) keep the find_mle default.
_profile_adsel(map_func) = isnothing(map_func) && Threads.nthreads() > 1 ? AutoForwardDiff() : AutoPolyesterForwardDiff()

"""
    _profile(likelihood, scanpoints, params, cache_dir; map_func=nothing, adsel=...)
        -> NamedTuple

Execute [`find_mle_cached`](@ref) at every point in `scanpoints`, collect
the results, and assemble them into a NamedTuple of result arrays with [`assemble_profile_results`](@ref).

Uses `Threads.@threads` by default. Pass a custom `map_func` (e.g. `pmap`)
to override the parallelism strategy.

# Arguments
- `likelihood`: a likelihood object supporting `DensityInterface.logdensityof`.
- `scanpoints`: array of prior objects, one per grid point.
- `params::NamedTuple`: starting values for nuisance optimization.
- `cache_dir::Union{String,Nothing}`: passed to [`find_mle_cached`](@ref).
- `map_func`: optional custom mapping function; default=`nothing` uses threaded loops.
- `adsel`: AD backend for the per-point fits; plain `AutoForwardDiff()` for the threaded loop,
  `AutoPolyesterForwardDiff()` otherwise.

# Returns
A NamedTuple of result arrays as produced by [`assemble_profile_results`](@ref).
"""
function _profile(likelihood, scanpoints, params, cache_dir; map_func=nothing, adsel=_profile_adsel(map_func), starts=[(;)])
    # one fit per start (start values override `params`), keep the one with the largest log posterior
    function do_work(i)
        # start values of parameters that are fixed at this grid point (scanned) do not matter: drop them
        free(st) = NamedTuple(k => v for (k, v) in pairs(st) if !(scanpoints[i][k] isa ValueShapes.ConstValueDist))
        sts = unique(free.(starts))
        fits = [find_mle_cached(likelihood, scanpoints[i], merge(deepcopy(params), s), cache_dir; adsel = adsel) for s in sts]
        lp = [isnan(f[2]) ? -Inf : f[2] for f in fits]
        fits[argmax(lp)]
    end

    if isnothing(map_func)
        # Default: threaded execution
        opt_results = Array{Any}(undef, size(scanpoints))
        @showprogress Threads.@threads for i in eachindex(scanpoints)
            opt_results[i] = do_work(i)
        end
    else
        # Custom map (e.g., pmap for distributed)
        work = collect(eachindex(scanpoints))
        opt_results_flat = map_func(do_work, work)
        opt_results = reshape(opt_results_flat, size(scanpoints))
    end
    assemble_profile_results(opt_results, size(scanpoints))
end

"""
    profile(likelihood, priors, vars_to_scan, params;
            cache_dir=nothing, map_func=nothing, adsel=...) -> NewtrinosResult

Run a profile likelihood scan over a parameter grid.

Creates grid points defined by `vars_to_scan` with [`generate_scanpoints`](@ref).
At each grid point, all nuisance parameters are optimized with [`_profile`](@ref) , i.e. via [`find_mle_cached`](@ref). 
If all non-scanned parameters have fixed (`Number`) priors, falls back to [`scan`](@ref) automatically.
Also collects meta data of the profile process and attaches it the returned object.

# Arguments
- `likelihood`: a likelihood object supporting `DensityInterface.logdensityof`.
- `priors::NamedTuple`: prior distributions for all parameters.
- `vars_to_scan::OrderedDict{Symbol,Int}`: parameters to scan and grid sizes.
- `params::NamedTuple`: starting values for nuisance parameter optimization.
- `cache_dir::Union{String,Nothing}`: directory for caching MLE results
  (created if absent). Pass `nothing` to disable (this is the default).
- `map_func`: custom mapping function for parallelism (e.g. `pmap` for
  distributed workers). Defaults to `Threads.@threads`.
- `adsel`: AD backend for the per-point fits. Defaults to `AutoForwardDiff()` for the threaded
  loop (the threads are busy with scan points, so a threaded AD backend would only contend), and
  to `AutoPolyesterForwardDiff()` otherwise.
- `starts`: start-value overrides for multi-start fits, e.g.
  `[(δCP = π/4,), (δCP = 3π/4,), (δCP = 5π/4,), (δCP = 7π/4,)]`: each grid point is fitted from every
  start (merged into `params`) and the fit with the largest log posterior is kept. Needed when a
  profiled parameter has several local optima (δCP, the θ₂₃ octant). Default: one start, `params`.

# Returns
A [`NewtrinosResult`](@ref) NewtrinosResult(axes, values, meta) 
with the scan grid axes, per-point profiling results, and meta data.
"""
function profile(likelihood, priors, vars_to_scan, params; cache_dir=nothing, map_func=nothing, adsel=_profile_adsel(map_func), starts=[(;)])
    t1 = time()
    # check if there is actually any variable to be profiled over, or if they are all just Numbers
    if all([isa(priors[var], Number) for var in setdiff(keys(priors), keys(vars_to_scan))])
        return scan(likelihood, priors, vars_to_scan, params)
    end

    values, scanpoints = generate_scanpoints(vars_to_scan, priors)
    if !isnothing(cache_dir)
        if isdir(cache_dir)
            @info "Reusing cache dir `$(cache_dir)`"
        else
            mkdir(cache_dir)
        end
    end
    res = _profile(likelihood, scanpoints, params, cache_dir; map_func=map_func, adsel=adsel, starts=starts)
    t2 = time()
    meta = Dict("task"=> "profile", "priors"=>priors, "vars_to_scan"=>vars_to_scan, "params"=>params, "starts"=>starts, "exec_time"=>t2-t1, "cache_dir"=>cache_dir)
    add_meta!(meta)
    axes = NamedTuple{tuple(keys(vars_to_scan)...)}(values)
    result = NewtrinosResult(axes=axes, values=res, meta=meta)
end

"""
    multistart_profile(likelihood, priors, vars_to_scan, starts; cache_dir=nothing, map_func=nothing)
        -> NewtrinosResult

Run [`profile`](@ref) once per starting point in `starts` over the same grid, then combine
grid-point-by-grid-point, keeping whichever start achieved the higher `log_posterior` at each
point. Guards against the local optimizer used by `find_mle`/`find_mle_cached` (gradient-based,
single fixed starting point) getting stuck in the wrong basin — e.g. the θ₂₃ octant degeneracy
in atmospheric oscillation fits — for part of the grid.

# Arguments
- `likelihood`, `priors`, `vars_to_scan`: same as [`profile`](@ref).
- `starts::AbstractVector{<:NamedTuple}`: candidate starting parameter sets to try at every grid
  point (e.g. one per θ₂₃ octant).
- `cache_dir`, `map_func`: passed through to each underlying `profile` call. Safe to share one
  `cache_dir` across all starts — `find_mle_cached`'s cache key hashes `(prior, params)`, and
  `params` differs per start, so entries never collide.

# Returns
A [`NewtrinosResult`](@ref) with the same `axes` as a single `profile` call, and `values` taken
pointwise from whichever start had the higher `log_posterior` at each grid point.
"""
function multistart_profile(likelihood, priors, vars_to_scan, starts; cache_dir=nothing, map_func=nothing)
    results = [profile(likelihood, priors, vars_to_scan, s; cache_dir=cache_dir, map_func=map_func) for s in starts]

    idx = CartesianIndices(results[1].values.log_posterior)
    best_start = map(i -> argmax([r.values.log_posterior[i] for r in results]), idx)

    combined = NamedTuple(
        key => map(i -> results[best_start[i]].values[key][i], idx)
        for key in keys(results[1].values)
    )
    NewtrinosResult(axes=results[1].axes, values=combined, meta=Dict("task"=>"multistart_profile", "n_starts"=>length(starts)))
end

"""
    scan(likelihood, priors, vars_to_scan, params;
         gradient_map=false) -> NewtrinosResult

Run a simple likelihood scan on a parameter grid (NO nuisance optimization).

Generates grid with [`_generate_grid`](@ref) via `vars_to_scan` and then 
evaluates `logdensityof(likelihood, params)` at each grid point without
optimizing over nuisance parameters. Faster than [`profile`](@ref) but
does not account for nuisance parameter variations. 
Also collects meta data of the scan process and attaches it the returned object.

# Arguments
- `likelihood`: a likelihood object supporting `DensityInterface.logdensityof`.
- `priors::NamedTuple`: prior distributions used to place grid points at
  quantiles.
- `vars_to_scan::OrderedDict{Symbol,Int}`: parameters to scan and grid sizes.
- `params::NamedTuple`: nominal values; non-scanned parameters are held fixed.
- `gradient_map::Bool`: default=`false`; if `true`, also evaluates `ForwardDiff.gradient` at
  each grid point and includes per-parameter gradient arrays in the result.

# Returns
A [`NewtrinosResult`](@ref) NewtrinosResult(axes, values, meta) 
with the scan grid axes, per-point scan results, and meta data.
"""
function scan(likelihood, priors, vars_to_scan, params; gradient_map=false)
    t1 = time()
    vars, values, mesh = _generate_grid(vars_to_scan, priors)
    scanpoints = Array{Any}(undef, size(mesh))

    function make_params(vals)
        p = deepcopy(params)
        for i in 1:length(vars)
            @reset p[vars[i]] = vals[i]
        end
        return p
    end

    for i in eachindex(mesh)
        scanpoints[i] = make_params(mesh[i])
    end

    llhs = Array{Float64}(undef, size(scanpoints))
    if gradient_map
        grads = Array{Any}(undef, size(scanpoints))
    end

    @showprogress Threads.@threads for i in eachindex(scanpoints)
        p = scanpoints[i]
        llhs[i] = logdensityof(likelihood, p)
        if gradient_map
            grads[i] = ForwardDiff.gradient(x -> logdensityof(likelihood, x), p)
        end
    end

    s = OrderedDict{Symbol, Array}(key=>Fill(params[key], size(mesh)) for key in setdiff(keys(params), keys(vars_to_scan)))
    if gradient_map
        g = OrderedDict(Symbol(key, "_grad")=>[x[key] for x in grads] for key in keys(first(grads)))
        s = merge(s, g)
    end
    s[:llh] = llhs
    s[:log_posterior] = llhs
    res = NamedTuple(s)
    t2 = time()
    meta = Dict("task"=>"scan", "priors"=>priors, "vars_to_scan"=>vars_to_scan, "params"=>params, "exec_time"=>t2-t1,)
    add_meta!(meta)
    axes = NamedTuple{tuple(keys(vars_to_scan)...)}(values)
    result = NewtrinosResult(axes=axes, values=res, meta=meta)
end

# ── Results ─────────────────────────────────────────────────────────

"""
    bestfit(result::NewtrinosResult) -> NamedTuple

Extract the best-fit parameter values from a scan or profile result.

Finds the grid point with the highest `log_posterior` and returns both
the scanned parameter values (from `result.axes`) and all optimized
nuisance parameter values at that point.

# Arguments
- `result::NewtrinosResult`: output of [`scan`](@ref) or [`profile`](@ref).

# Returns
A `NamedTuple` containing all parameter values at the best-fit point,
including the scanned axis values, nuisance parameters, `llh`, and
`log_posterior`.
"""
function bestfit(result::NewtrinosResult)
    idx = argmax(result.values.log_posterior)
    bf = OrderedDict(var=>result.values[var][idx] for var in keys(result.values))
    for i in 1:length(result.axes)
        bf[keys(result.axes)[i]] = result.axes[i][idx[i]]
    end
    NamedTuple(bf)
end

"""
    add_meta!(meta::Dict)

Populate a metadata dictionary with execution environment information. When Newtrinos 
is loaded from the julia package index, it marks the git-repo information as missing.

Adds the following keys in-place:
- `"hostname"`: result of `gethostname()`.
- `"username"`: from `ENV["USER"]` or `ENV["USERNAME"]`.
- `"date"`: current date-time formatted as `"yyyy-mm-dd HH:MM:SS"`.
- `"repo"`: path to the Newtrinos.jl root.
- `"package_version"`: loaded version of the Newtrinos.jl package.
- `"commit_hash"`: current HEAD commit hash.
- `"repo_clean"`: `true` if the repository has no uncommitted changes.

Called automatically by [`profile`](@ref) and [`scan`](@ref) before
returning their [`NewtrinosResult`](@ref). 

# Arguments
- `meta::Dict`: dictionary to populate in-place.

# Returns
`nothing` (mutates `meta`).
"""
function add_meta!(meta)
    meta["hostname"] = gethostname()
    meta["username"] = get(ENV, "USER", get(ENV, "USERNAME", "unknown"))
    meta["date"] = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    repo = dirname(dirname(pathof(Newtrinos)))
    meta["repo"] = repo
    meta["package_version"] = pkgversion(Newtrinos)
    try
        meta["commit_hash"] = LibGit2.head(repo)
        meta["repo_clean"] = !LibGit2.isdirty(LibGit2.GitRepo(repo))
    catch e    
        meta["commit_hash"] = missing
        meta["repo_clean"] = missing
    end
end

function _gaussian_blur(arr::Array, sigma::Float64)
    sigma <= 0 && return float.(arr)
    result = float.(arr)
    for d in 1:ndims(arr)
        n = size(result, d)
        radius = min(ceil(Int, 3 * sigma), n - 1)
        k_range = -radius:radius
        kernel = [exp(-k^2 / (2 * sigma^2)) for k in k_range]
        old = copy(result)
        for idx in CartesianIndices(size(old))
            s = 0.0; w = 0.0
            for (ki, delta) in enumerate(k_range)
                ni = idx[d] + delta
                1 <= ni <= n || continue
                nidx = CartesianIndex(ntuple(i -> i == d ? ni : idx[i], ndims(arr)))
                s += kernel[ki] * old[nidx]
                w += kernel[ki]
            end
            result[idx] = s / w
        end
    end
    result
end

"""
    smooth_result(result; outlier_threshold=10.0, gaussian_sigma=1.0)

Prepare a `NewtrinosResult` for use as per-point seeds in `refine_profile`:

1. Detect outlier scan points where `log_posterior` is more than `outlier_threshold` below
   the maximum of its direct grid neighbours (indicating a failed optimisation).
2. Replace each outlier's parameter values with the mean of its valid (non-outlier) neighbours.
3. Apply a Gaussian blur of width `gaussian_sigma` (in grid units) to all parameter arrays.

Returns a new `NewtrinosResult` with smoothed values.
"""
function smooth_result(result::NewtrinosResult; outlier_threshold=10.0, gaussian_sigma=1.0)
    lp = result.values.log_posterior        # already grid-shaped
    grid_shape = size(lp)
    ndim = ndims(lp)
    param_keys = [k for k in keys(result.values) if k ∉ (:llh, :log_posterior)]

    # Helper: direct grid neighbours of idx (in-bounds only)
    function neighbor_indices(idx)
        nbrs = CartesianIndex[]
        for d in 1:ndim, delta in (-1, 1)
            nidx = CartesianIndex(ntuple(i -> i == d ? idx[i] + delta : idx[i], ndim))
            checkbounds(Bool, lp, nidx) && push!(nbrs, nidx)
        end
        nbrs
    end

    # 1. Detect outliers
    is_outlier = falses(grid_shape)
    for idx in CartesianIndices(grid_shape)
        nbrs = neighbor_indices(idx)
        isempty(nbrs) && continue
        lp[idx] < maximum(lp[n] for n in nbrs) - outlier_threshold && (is_outlier[idx] = true)
    end
    @info "smooth_result: $(count(is_outlier)) outlier(s) detected in $(prod(grid_shape)) scan points"

    # 2. Replace outlier param values with mean of valid neighbours, then Gaussian-blur
    new_params = Dict{Symbol, Vector{Float64}}()
    for k in param_keys
        arr = reshape(copy(float.(result.values[k])), grid_shape)
        for idx in CartesianIndices(grid_shape)
            is_outlier[idx] || continue
            valid_vals = [arr[n] for n in neighbor_indices(idx) if !is_outlier[n]]
            isempty(valid_vals) || (arr[idx] = mean(valid_vals))
        end
        new_params[k] = vec(_gaussian_blur(arr, gaussian_sigma))
    end

    # Also repair log_posterior at outlier points (for reference only, not used as seeds)
    new_lp = copy(float.(lp))
    for idx in CartesianIndices(grid_shape)
        is_outlier[idx] || continue
        valid_vals = [lp[n] for n in neighbor_indices(idx) if !is_outlier[n]]
        isempty(valid_vals) || (new_lp[idx] = mean(valid_vals))
    end

    new_values = NamedTuple{keys(result.values)}(
        [k == :log_posterior ? new_lp :
         k == :llh           ? result.values.llh :
         new_params[k]
         for k in keys(result.values)])
    meta = merge(result.meta, Dict("smoothed"=>true,
                                   "outlier_threshold"=>outlier_threshold,
                                   "gaussian_sigma"=>gaussian_sigma))
    NewtrinosResult(axes=result.axes, values=new_values, meta=meta)
end

"""
    make_seed_interpolants(result)

Build a linear interpolant (with linear extrapolation outside the grid) for each
optimised parameter in `result`, indexed by the scan axis values.
Returns a `Dict{Symbol, AbstractExtrapolation}`.
"""
function make_seed_interpolants(result::NewtrinosResult)
    axis_vals = Tuple(collect(Float64, v) for v in values(result.axes))
    grid_shape = size(result.values.log_posterior)
    param_keys = [k for k in keys(result.values) if k ∉ (:llh, :log_posterior)]
    itps = Dict{Symbol, Any}()
    for k in param_keys
        arr = Float64.(reshape(result.values[k], grid_shape))
        itp = interpolate(axis_vals, arr, Gridded(Linear()))
        itps[k] = extrapolate(itp, Linear())
    end
    itps
end

"""
    make_merged_seed_interpolants(seed_results, joint_axis_keys)

Build merged per-parameter seed evaluators from multiple `NewtrinosResult` files.

For each parameter:
- Present in one file: use that file's interpolant.
- Present in multiple files: evaluate each and take the mean.

`joint_axis_keys`: the scan axis symbols of the new joint scan (e.g. `keys(vars_to_scan)`).
Source files scanned over a subset of these axes are handled via coordinate projection.

Returns a `Dict{Symbol, Function}` where each function takes `(coords...)` positional
coordinates aligned to `joint_axis_keys`.
"""
function make_merged_seed_interpolants(
    seed_results::AbstractVector,
    joint_axis_keys
)
    joint_keys = collect(joint_axis_keys)
    source_groups = map(seed_results) do r
        src_keys = collect(keys(r.axes))
        idxs = [findfirst(==(ak), joint_keys) for ak in src_keys]
        bad = findfirst(isnothing, idxs)
        isnothing(bad) || error("Seed file axis $(src_keys[bad]) not in joint axes $joint_keys")
        (idxs, make_seed_interpolants(r))
    end

    all_param_keys = unique(reduce(vcat, [collect(keys(sg[2])) for sg in source_groups]))

    merged = Dict{Symbol, Function}()
    for k in all_param_keys
        sources_k = [(idxs, itps[k]) for (idxs, itps) in source_groups if haskey(itps, k)]
        if length(sources_k) == 1
            (idxs, itp) = sources_k[1]
            merged[k] = (coords...) -> itp(coords[idxs]...)
        else
            merged[k] = (coords...) -> mean(itp(coords[idxs]...) for (idxs, itp) in sources_k)
        end
    end
    merged
end

"""
    correlation(likelihood, prior, params; bound_tol = 1e-6) -> (; cov, corr, names, dropped, offset)
    correlation(result::NewtrinosResult, likelihood; bound_tol = 1e-6)
    correlation(experiments::NamedTuple, prior, params; bound_tol = 1e-6)

Parameter covariance and Pearson correlation from the inverse Hessian of the negative log posterior (log likelihood plus log prior of the
free parameters) at the point `params`, i.e. the Gaussian (Laplace) approximation: `cov = inv(-H)`, `corr = cov_ij / sqrt(cov_ii cov_jj)`.
Each parameter's correlation with the others is the marginal one, all other free parameters being profiled.

The first two methods use the *observed* information: the Hessian of the log likelihood of the data. The result method evaluates at the
best fit of a profile result (`bestfit(result)`, using the priors stored in `result.meta["priors"]`). At a grid point of a dense profile this
is close enough to the true minimum; `offset` (the Newton step to the stationary point of the log posterior, `cov * gradient`, as the
largest component in units of its parameter's uncertainty) tells how far from it the point is. The likelihood method works at any point,
e.g. at the truth on Asimov data (`generate_asimov_data`), but the rounding of the Asimov counts then enters the Hessian.

The `experiments` method uses the *expected* (Fisher) information at `params` of Poisson experiments, without any data: the Hessian of
`sum(mu0 * log(mu(q)) - mu(q))` over all bins, with `mu0` the expected counts at `params` fixed (the Asimov log likelihood without
rounding), equal to `J' diag(1/mu) J` with `J = d mu / d q`. It is positive semi-definite by construction. `offset` then only measures the
pull of the priors. Other than Poisson forward models raise an error.

Only free parameters enter: those with a `Number` or `ConstValueDist` prior are fixed; scalar parameters closer than `bound_tol` (relative
to the prior width) to a prior bound are dropped as well (listed in `dropped`, `bound_tol = 0` keeps all), as the Gaussian approximation
fails there. Vector parameters (e.g. `solar_norms`) enter with all their components, named `solar_norms[1]`, ...
An error (`PosDefException`) is thrown if the negative Hessian is not positive definite (no maximum, or a flat direction).

# Arguments
- `likelihood`: a likelihood object supporting `DensityInterface.logdensityof` (see [`generate_likelihood`](@ref)).
- `experiments`: NamedTuple of experiments, as for [`generate_likelihood`](@ref).
- `prior`: NamedTuple of prior distributions (fixed parameters as numbers or `ConstValueDist`).
- `params::NamedTuple`: the point at which the Hessian is taken, all parameters.

# Returns
A `NamedTuple` `(; cov, corr, names, dropped, offset)`: the matrices of the entering parameters in the order of `names`.
Plot with `Newtrinos.plotting.plot_correlation`.
"""
function correlation(likelihood, prior, params::NamedTuple; bound_tol = 1e-6)
    _correlation(q -> logdensityof(likelihood, merge(params, q)), prior, params, bound_tol)
end

function correlation(experiments::NamedTuple, prior, params::NamedTuple; bound_tol = 1e-6)
    counts(d) = d isa ValueShapes.NamedTupleDist ? reduce(vcat, [counts(d[k]) for k in keys(d)]) : [mean(p) for p in _poisson_leaves(d)]
    ks = map(e -> keys(get_params(e)), values(experiments))
    expected(p) = reduce(vcat, [counts(e.forward_model(NamedTuple{k}(p))) for (e, k) in zip(values(experiments), ks)])
    μ0 = expected(params)
    _correlation(q -> (μ = expected(merge(params, q)); sum(m0 > 0 ? m0 * log(m) - m : -m for (m0, m) in zip(μ0, μ))), prior, params, bound_tol)
end

# shared part: `data_logl(q)` is the data log likelihood at the free parameters `q` (the others fixed at `params`)
function _correlation(data_logl, prior, params::NamedTuple, bound_tol)
    is_free(k) = k in keys(prior) && !(prior[k] isa Number || prior[k] isa ValueShapes.ConstValueDist)
    function at_bound(k)
        d = prior[k]
        d isa Distributions.UnivariateDistribution || return false
        lo, hi = minimum(d), maximum(d)
        tol = bound_tol * (isfinite(hi - lo) ? hi - lo : 1.0)
        (isfinite(lo) && params[k] - lo < tol) || (isfinite(hi) && hi - params[k] < tol)
    end
    dropped = [k for k in keys(params) if is_free(k) && at_bound(k)]
    free = [k for k in keys(params) if is_free(k) && !(k in dropped)]
    p_free = NamedTuple{Tuple(free)}(Tuple(params[k] for k in free))

    f(q) = data_logl(q) + sum(logpdf(prior[k], q[k]) for k in keys(q))
    h = Newtrinos.hessian(f, p_free)   # NamedTuple of row blocks: stack them to the full matrix
    H = reduce(vcat, [b isa AbstractMatrix ? b : permutedims(b) for b in values(h)])
    A = Symmetric(-(H + H') / 2)
    isposdef(A) || throw(LinearAlgebra.PosDefException(1))   # not a maximum (or a flat direction): cov would have negative variances
    cov = Matrix(inv(A))
    corr = cov ./ sqrt.(diag(cov) * diag(cov)')
    g = ForwardDiff.gradient(f, p_free)
    offset = maximum(abs.(cov * reduce(vcat, [v isa AbstractVector ? v : [v] for v in values(g)])) ./ sqrt.(diag(cov)))

    names = Symbol[]
    for k in free
        v = params[k]
        v isa AbstractVector ? append!(names, [Symbol(k, "[", i, "]") for i in eachindex(v)]) : push!(names, k)
    end
    (; cov, corr, names, dropped, offset)
end

# the Poisson distributions of a forward model leaf (a `Poisson` or a product of them); error for anything else
function _poisson_leaves(d)
    pois = d isa Distributions.Poisson ? [d] :
           d isa Distributions.Product ? vec(d.v) :
           d isa Distributions.ProductDistribution ? vec(d.dists) : []
    (!isempty(pois) && all(p -> p isa Distributions.Poisson, pois)) || error("only Poisson forward models are supported, got $(typeof(d))")
    pois
end

function correlation(result::NewtrinosResult, likelihood; kwargs...)
    bf = bestfit(result)
    ks = Tuple(k for k in keys(bf) if k ∉ (:llh, :log_posterior))
    correlation(likelihood, result.meta["priors"], NamedTuple{ks}(bf); kwargs...)
end

"""
    gof_pvalue(result::NewtrinosResult, experiment::Newtrinos.Experiment) -> (; deviance, ndof, pvalue)
    gof_pvalue(result::NewtrinosResult, experiments::NamedTuple) -> NamedTuple of those, one per experiment

Goodness-of-fit p-value of an experiment at the best fit of `result` (no refit): the Poisson deviance against the saturated model,
`2 * sum(logpdf(Poisson(n), n) - logpdf(model, n))` over all bins with observed counts `n`, compared with a `Chisq(ndof)` distribution
(Wilks): `pvalue = ccdf(Chisq(ndof), deviance)`.

- `ndof` = number of bins - number of free parameters, the free parameters being the non-fixed priors in `result.meta["priors"]` that
  belong to the experiment (`get_params(experiment)`; vector parameters count with their length). This includes the scanned parameters
  (free in the stored priors) and counts prior-constrained nuisance parameters as fully free (conservative).
- Pass a result of a fit of this experiment alone for its own best fit, or the result of a combined fit to test it at the joint best
  fit. For a combined fit the shared parameters make `ndof` approximate.
- The deviance is only defined for Poisson forward models (a `Poisson` product, possibly inside a `NamedTupleDist`); other
  distributions raise an error.
"""
function gof_pvalue(result::NewtrinosResult, experiment::Newtrinos.Experiment)
    bf = bestfit(result)
    ks = keys(get_params(experiment))
    dist = experiment.forward_model(NamedTuple{ks}(bf))

    function deviance(d, x)   # (deviance, number of bins) of one distribution
        if d isa ValueShapes.NamedTupleDist
            parts = [deviance(d[k], x[k]) for k in keys(d)]
            return sum(first, parts), sum(last, parts)
        end
        pois = _poisson_leaves(d)
        n = x isa Number ? [x] : vec(x)
        2 * sum(logpdf(Distributions.Poisson(ni), ni) - logpdf(p, ni) for (p, ni) in zip(pois, n)), length(n)
    end
    dev, nbins = deviance(dist, experiment.assets.observed)

    priors = result.meta["priors"]
    n_free = sum(length(priors[k]) for k in ks if k in keys(priors) && !(priors[k] isa Number || priors[k] isa ValueShapes.ConstValueDist); init = 0)
    ndof = nbins - n_free
    (; deviance = dev, ndof, pvalue = ccdf(Distributions.Chisq(ndof), dev))
end

gof_pvalue(result::NewtrinosResult, experiments::NamedTuple) = map(e -> gof_pvalue(result, e), experiments)
