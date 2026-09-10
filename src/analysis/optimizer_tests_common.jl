# Shared optimizer helpers for optimizer_tests.jl. `include`d both in the main process
# and (via @everywhere) on distributed workers, so it must not reference any variable
# that isn't passed in as an argument.

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
    fit_quality(truth_param, fit_param, prior; eps=1e-8)

Per-roundtrip fit-quality measure over the parameters that were free in `prior`
(fixed/conditioned parameters are excluded, matching `extract_bounds`):

    A = 1/n_params * Σ [ |fit_i - truth_i| / (|truth_i| + |fit_i| + eps) ]^2

Returns `(A, per_param)` where `per_param` is a `Dict{Symbol,Float64}` of each free
parameter's individual squared relative error (before averaging), so per-parameter
quality can be aggregated separately across roundtrips.
"""
function fit_quality(truth_param, fit_param, prior; eps=1e-8)
    per_param = Dict{Symbol, Float64}()
    for k in keys(prior)
        d = prior[k]
        (d isa ValueShapes.ConstValueDist || d isa Number) && continue
        t = Float64(truth_param[k])
        f = Float64(fit_param[k])
        per_param[k] = ((abs(f - t)) / (abs(t) + abs(f) + eps))^2
    end
    A = isempty(per_param) ? NaN : sum(values(per_param)) / length(per_param)
    (A, per_param)
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
    controlled_truth(rng, nominal_params, prior, nsi_param_names; vary_physics=false)

Build an injected-truth NamedTuple starting from `nominal_params` (e.g.
`Newtrinos.get_params(experiments)`), overriding: (1) every free (non-fixed) parameter
in `prior` whose key is in `nsi_param_names`, drawn randomly from its prior
distribution, (2) `:θ₂₃` and `:Δm²₃₁`, drawn randomly from their priors, only if
`vary_physics`, and (3) every parameter `prior` holds **fixed** (`ConstValueDist` /
`Number`, e.g. δCP/θ₁₂/θ₁₃/Δm²₂₁ conditioned in `optimizer_tests.jl`), forced to that
exact fixed value -- `nominal_params` (`Newtrinos.get_params`'s own defaults) can
disagree with the specific constants a hypothesis's `prior` was conditioned to (verified
directly: `p[:δCP]=1.0` vs. the conditioned `0.0`), and `bat_findmode`'s `ExplicitInit`
strictly validates the init point against the prior's constants, raising
`ArgumentError("Cannot set constant value to a different value")` if they differ.

Every other parameter -- nuisances, θ₂₃/Δm²₃₁ (unless varied) -- stays at its nominal
value. This isolates truth variation to the NSI parameter(s) actually being tested by
`prior`'s hypothesis, avoiding prior-pull artifacts from randomly-fluctuated nuisance
parameters landing in low-density prior regions by chance (see `randomize_free_params`
for the same fixed-parameter check used here).
"""
function controlled_truth(rng, nominal_params, prior, nsi_param_names; vary_physics=false)
    overrides = Dict{Symbol, Float64}()
    for k in nsi_param_names
        d = prior[k]
        (d isa ValueShapes.ConstValueDist || d isa Number) && continue
        overrides[k] = rand(rng, d)
    end
    if vary_physics
        for k in (:θ₂₃, :Δm²₃₁)
            d = prior[k]
            (d isa ValueShapes.ConstValueDist || d isa Number) && continue
            overrides[k] = rand(rng, d)
        end
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
    truth_grid(prior, nominal_params, nsi_param_names; n_points=25)

Generate a grid of `n_points` Asimov-truth NamedTuples scanning the free NSI
parameter(s) of `prior`'s hypothesis, matching the dissertation's "recovering injected
NSI hypotheses" test (section 7.6.2): nuisances/θ₂₃/Δm²₃₁ pinned at `nominal_params`
(and every prior-fixed parameter forced to its exact conditioned value, via
`_fix_conditioned` -- same as `controlled_truth`), only the hypothesis's free NSI
parameter(s) varied. `fit_param`/optimization itself is NOT restricted by this
function -- it only controls how the injected truth is built; the subsequent fit is
still free to vary every non-fixed parameter in `prior` (physics + nuisances), exactly
like any other roundtrip.

Auto-detects sampling strategy from which of `nsi_param_names` are free in `prior`:
- Exactly 1 free param that's a lone real-valued NSI param (`:Δ_eμ` or `:Δ_τμ`): 1D
  grid, `n_points` evenly-spaced values across `[minimum(d), maximum(d)]`.
- Exactly 2 free params forming a magnitude/phase pair (`:ε_eμ_abs`+`:δ_eμ`,
  `:ε_eτ_abs`+`:δ_eτ`, or `:ε_μτ_abs`+`:δ_μτ`): 2D grid, `round(sqrt(n_points))` points
  per dimension (e.g. 5x5=25 for `n_points=25`), magnitude x phase.
- More than 2 free params (e.g. `cp_all`'s 8 free NSI params): a full factorial grid is
  combinatorially infeasible (`n_points^8` for an 8D hypothesis at `n_points`/dim), so
  instead draws `n_points` quasi-random samples via a `Sobol.jl` low-discrepancy
  sequence (`SobolSeq`) over all free dimensions simultaneously -- far better space-
  filling coverage of the joint parameter space than either a factorial grid or
  independent-uniform random draws, at a fixed, chosen sample budget.
- 0 free NSI params is not a valid hypothesis for this function and errors.

Returns `(truths, scan_keys, grid_coords)`: `truths::Vector{<:NamedTuple}` (the sampled
truths), `scan_keys::Tuple` (which key(s) were scanned), and `grid_coords` aligned 1:1
with `truths` (`Vector{Float64}` for 1D, `Vector{Tuple{Float64,Float64}}` for 2D,
`Vector{Vector{Float64}}` -- one vector per sample, ordered to match `scan_keys` -- for
the >2D Sobol case).
"""
function truth_grid(prior, nominal_params, nsi_param_names; n_points=25)
    free_nsi = [k for k in nsi_param_names if !(prior[k] isa ValueShapes.ConstValueDist) && !(prior[k] isa Number)]

    magnitude_phase_pairs = Dict(:ε_eμ_abs => :δ_eμ, :ε_eτ_abs => :δ_eτ, :ε_μτ_abs => :δ_μτ)
    standalone_real = (:Δ_eμ, :Δ_τμ)

    if length(free_nsi) == 1 && free_nsi[1] in standalone_real
        k = free_nsi[1]
        d = prior[k]
        lo, hi = Float64(minimum(d)), Float64(maximum(d))
        grid_coords = collect(range(lo, hi, length=n_points))
        truths = [_fix_conditioned(merge(nominal_params, NamedTuple{(k,)}((v,))), prior) for v in grid_coords]
        return truths, (k,), grid_coords
    elseif length(free_nsi) == 2 && haskey(magnitude_phase_pairs, free_nsi[1]) &&
           magnitude_phase_pairs[free_nsi[1]] == free_nsi[2]
        mag_key, phase_key = free_nsi[1], free_nsi[2]
        mag_d, phase_d = prior[mag_key], prior[phase_key]
        n_per_dim = Int(round(sqrt(n_points)))
        mag_vals   = collect(range(Float64(minimum(mag_d)), Float64(maximum(mag_d)), length=n_per_dim))
        phase_vals = collect(range(Float64(minimum(phase_d)), Float64(maximum(phase_d)), length=n_per_dim))
        grid_coords = [(m, ph) for m in mag_vals for ph in phase_vals]
        truths = [_fix_conditioned(merge(nominal_params, NamedTuple{(mag_key,phase_key)}((m,ph))), prior)
                  for (m,ph) in grid_coords]
        return truths, (mag_key, phase_key), grid_coords
    elseif length(free_nsi) > 2
        scan_keys = Tuple(free_nsi)
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
        error("truth_grid: unsupported free-NSI-parameter combination $(free_nsi) -- " *
              "expected exactly one of $(standalone_real), one of the magnitude/phase " *
              "pairs $(collect(magnitude_phase_pairs)), or more than 2 free params for " *
              "Sobol sampling (e.g. cp_all).")
    end
end

"""
    global_seed_search(likelihood, prior, params, strategy, global_iters, rng; n_particles=10)

2-stage seeding: runs a short bounded, gradient-free global search (PSO or simulated
annealing) over the free parameters and returns an updated `params` NamedTuple with the
free entries replaced by the global search's best point, to be used as the starting
point for the local optimizer's polish step. `strategy == "random"` bypasses this and
just draws a fresh random start from the priors via `randomize_free_params`.

`n_particles` sets `Optim.ParticleSwarm`'s swarm size (`strategy == "pso"` only —
`SimulatedAnnealing` is a single-trajectory method with no equivalent knob).

`return_trace=true` additionally returns the best-objective-value-so-far at each
generation (`Optim.f_trace(result)`), for diagnosing how many generations the global
search actually needs before plateauing; default `false` keeps the normal single-value
return for all existing callers.
"""
function global_seed_search(likelihood, prior, params, strategy, global_iters, rng; n_particles=10, return_trace=false)
    if strategy == "random"
        seeded = randomize_free_params(rng, params, prior)
        return return_trace ? (seeded, Float64[]) : seeded
    end

    b = extract_bounds(prior, params)

    function obj(x::AbstractVector)
        all_x = [b.is_free_mask[i] ? x[b.free_x_index[i]] : b.fixed_float[i] for i in 1:length(b.all_keys)]
        full_params = NamedTuple{b.all_keys}(Tuple(all_x))
        llh_val   = logdensityof(likelihood, full_params)
        prior_val = sum(logpdf(b.free_dists[j], x[j]) for j in 1:length(b.free_dists))
        isfinite(llh_val) && isfinite(prior_val) ? -(llh_val + prior_val) : Inf
    end

    galg = strategy == "pso" ? Optim.ParticleSwarm(lower=b.lo_vec, upper=b.hi_vec, n_particles=n_particles) :
                                Optim.SimulatedAnnealing()

    x0 = clamp.(b.x0_vec, b.lo_vec, b.hi_vec)
    result = Optim.optimize(obj, x0, galg, Optim.Options(iterations=global_iters, store_trace=return_trace))
    x_best = Optim.minimizer(result)

    all_x = [b.is_free_mask[i] ? x_best[b.free_x_index[i]] : b.fixed_float[i] for i in 1:length(b.all_keys)]
    seeded = NamedTuple{b.all_keys}(Tuple(all_x))
    return_trace ? (seeded, Optim.f_trace(result)) : seeded
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
function local_find_mle(likelihood, prior, params; fit_method="optim", algorithm=nothing, iterations=nothing, g_tol=nothing, f_tol=nothing, x_tol=nothing, return_trace=false, objective="posterior", ad_backend=nothing)
    if fit_method == "newtrinos"
        result = (Newtrinos.find_mle(likelihood, prior, params)..., missing)
        return return_trace ? (result..., Float64[]) : result
    elseif fit_method == "ext"
        result = (Newtrinos.find_mle_ext(likelihood, prior, params; g_tol=g_tol, maxiters=iterations)..., missing)
        return return_trace ? (result..., Float64[]) : result
    end

    # fit_method == "optim"
    target_prior = objective == "likelihood" ? flat_prior(prior) : prior
    posterior = PosteriorMeasure(likelihood, target_prior)
    n_free = count(!(d isa ValueShapes.ConstValueDist) && !(d isa Number) for d in values(prior))

    @info "Running Optimization (optim, $(typeof(algorithm).name.name)) for point" n_free=n_free iterations=iterations objective=objective

    try
        t0 = time()
        # Newtrinos.find_mle sets this internally (adsel = Newtrinos.select_ad(length(params)));
        # bat_findmode needs an AD backend on the context too, or gradient-based algorithms
        # (LBFGS etc.) fail immediately with "requires automatic differentiation".
        # set_batcontext(ad = Newtrinos.select_ad(n_free)) this caused thread racing / threading conflicts
	ad = ad_backend === nothing ? Newtrinos.select_ad(n_free) : ad_backend # for choosing ad backend when calling local_find_mle
	context = BATContext(ad = ad) # for setting the bat context locally inside bat_findmode and not globally 
        @info "starting optimization"
        res = bat_findmode(posterior, OptimAlg(optalg=algorithm, init=ExplicitInit([params]), maxiters=iterations,
                                                kwargs=(g_tol=g_tol, f_reltol=f_tol, x_reltol=x_tol, store_trace=return_trace)), context) # set bat context locally
        @info "extracting fit results"
        converged = Optim.converged(res.info)
        llh      = logdensityof(likelihood, res.result)
        log_post = logdensityof(posterior, res.result)
        n_iters  = Optim.iterations(res.info)
        @info "Optimization (optim) finished" converged=converged elapsed_s=round(time()-t0, digits=1) iterations=n_iters
        result = (llh, log_post, res.result, converged, n_iters)
        return return_trace ? (result..., Optim.f_trace(res.info.res)) : result
    catch e
	(e isa ArgumentError || e isa AssertionError) || rethrow(e)
        @warn "local_find_mle caught $(typeof(e)), returning NaN" exception=(e, catch_backtrace())
        nan_params = NamedTuple{keys(prior)}(ntuple(_ -> NaN, length(keys(prior))))
        result = (NaN, NaN, nan_params, false, missing)
        return return_trace ? (result..., Float64[]) : result
    end
end

"""
    do_roundtrip(k, work, truth_seeds, data_mode, experiments, cps, seed_strategy,
                 global_iters, start_rngs, fit_method, algorithm, iterations, g_tol, f_tol, x_tol;
                 nseeds=1, pso_particles=10, trace_plots=false, objective="posterior")

Run one truth-injection roundtrip: build the truth-injected likelihood, then run `nseeds`
independent local-optimizer fits (each from its own starting point via
`global_seed_search`) and keep whichever converged to the highest `log_posterior`
(`NaN` treated as `-Inf`) — mirrors `Newtrinos.profile(...; nseeds=N)`'s multi-start
pattern. `nseeds=1` (the default) reduces to a single fit, unchanged from before.

Also computes `fit_quality` (see that function's docstring) comparing the injected truth
against the winning fit's parameters, over the free parameters of this hypothesis.

`trace_plots=true` additionally collects, for every seed `s in 1:nseeds` (not just the
winner), a `(pso_trace, lbfgs_trace)` pair via `global_seed_search`/`local_find_mle`'s
`return_trace`, for convergence-trace plotting. Default `false` skips trace collection
entirely (no `store_trace` overhead) and returns an empty vector for that element.

`objective` is passed straight through to `local_find_mle` (see its docstring): with
`"likelihood"`, `log_post`/`best_log_post`'s seed-selection comparison is the flat-prior
pseudo-posterior, which differs from true MAP posterior only by a constant, so picking
the highest one across seeds still selects the highest-likelihood fit correctly.

Returns `(j, i, truth_param, truth_likelihood, fit_result, posterior, A, per_param, all_seed_traces)`.
"""
function do_roundtrip(k, work, truth_seeds, data_mode, experiments, cps, seed_strategy,
                       global_iters, start_rngs, fit_method, algorithm, iterations, g_tol, f_tol, x_tol;
                       nseeds=1, pso_particles=10, trace_plots=false, objective="posterior")
    j, i = work[k]
    prior_dist = distprod(;cps[j]...)
    truth_param = truth_seeds[k]

    injected_data = data_mode == "toy" ?
        Newtrinos.generate_toy_data(experiments, truth_param) :
        Newtrinos.generate_asimov_data(experiments, truth_param)
    println("generated data")
    truth_likelihood = Newtrinos.generate_likelihood(experiments, injected_data)
    println("generated likelihood")

    best_result = nothing
    best_log_post = -Inf
    all_seed_traces = Tuple{Vector{Float64},Vector{Float64}}[]
    for s in 1:nseeds
        if trace_plots
            start_param, pso_trace = global_seed_search(truth_likelihood, cps[j], truth_param, seed_strategy, global_iters, start_rngs[k]; n_particles=pso_particles, return_trace=true)
        else
            start_param = global_seed_search(truth_likelihood, cps[j], truth_param, seed_strategy, global_iters, start_rngs[k]; n_particles=pso_particles)
            pso_trace = Float64[]
        end
        println("generated start seed $s/$nseeds")
        println("now fitting:")
        result = local_find_mle(truth_likelihood, prior_dist, start_param;
            fit_method=fit_method, algorithm=algorithm, iterations=iterations, g_tol=g_tol, f_tol=f_tol, x_tol=x_tol,
            return_trace=trace_plots, objective=objective)
        if trace_plots
            push!(all_seed_traces, (pso_trace, result[6]))
            result = result[1:5]
        end
        log_post = isnan(result[2]) ? -Inf : result[2]
        if log_post > best_log_post
            best_result = result
            best_log_post = log_post
        end
    end
    fit_result = best_result
    println("finished fitting, now calculating truth posterior")
    posterior = PosteriorMeasure(truth_likelihood, prior_dist)
    A, per_param = fit_quality(truth_param, fit_result[3], cps[j])
    (j, i, truth_param, truth_likelihood, fit_result, posterior, A, per_param, all_seed_traces)
end
