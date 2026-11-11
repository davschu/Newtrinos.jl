# Profile likelihood scan with additional start-value (seeding) options.
#
# `profile_advanced` extends [`profile`](@ref) (unchanged upstream function in analysis_tools.jl) by
#   - `nseeds`:      several random start values per grid point, best fit kept
#   - `seed_result`, `seed_results`: per-point start values interpolated from earlier profile results
#   - `sequential`, `start_from`:    warm-starting each grid point from its already fitted neighbour (BFS order)
# The helpers it uses (`randomize_params`, `_sequential_order`, `make_seed_interpolants`,
# `make_merged_seed_interpolants`) live in analysis_tools.jl.

"""
    _profile_advanced(likelihood, scanpoints, params, cache_dir; map_func=nothing, adsel=...)
        -> NamedTuple

Worker of [`profile_advanced`](@ref): like `Newtrinos._profile`, with the additional per-point seeding
options `nseeds`, `seed_params` and `sequential` (BFS warm-starting from `center_linear`).

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
function _profile_advanced(likelihood, scanpoints, params, cache_dir; map_func=nothing, adsel=_profile_adsel(map_func), starts=[(;)],
                  nseeds=1, seed_params=nothing, sequential=false, center_linear=1)
    # base start values of grid point i: `params`, or the per-point seeds `seed_params` if given
    get_p(i) = isnothing(seed_params) ? params : seed_params[i]

    # one fit per start (start values override the base values `p0`), keep the one with the largest log posterior
    function fit_point(i, p0)
        # start values of parameters that are fixed at this grid point (scanned) do not matter: drop them
        free(st) = NamedTuple(k => v for (k, v) in pairs(st) if !(scanpoints[i][k] isa ValueShapes.ConstValueDist))
        sts = unique(free.(starts))
        fits = [find_mle_cached(likelihood, scanpoints[i], merge(deepcopy(p0), s), cache_dir; adsel = adsel) for s in sts]
        lp = [isnan(f[2]) ? -Inf : f[2] for f in fits]
        fits[argmax(lp)]
    end
    do_work(i) = fit_point(i, get_p(i))

    if sequential
        # Sequential warm-starting: BFS order from center_linear, each point seeded by its parent's optimised result
        grid_shape = size(scanpoints)
        order, parent = _sequential_order(grid_shape, center_linear)
        opt_results_flat = Vector{Any}(undef, prod(grid_shape))

        @info "Sequential profile: $(prod(grid_shape)) scan points, BFS from linear index $center_linear"
        p_idx = ProgressMeter.Progress(prod(grid_shape))
        for i in order
            seed = parent[i] < 0 ? get_p(i) : opt_results_flat[parent[i]][3]
            opt_results_flat[i] = fit_point(i, seed)
            ProgressMeter.next!(p_idx)
        end
        return assemble_profile_results(opt_results_flat, grid_shape)
    end

    if nseeds > 1
        # Multi-seed: expand each scan point into nseeds independent jobs (random start values drawn from the
        # priors), keep the best fit (highest log posterior) per scan point
        n_points = length(scanpoints)
        n_jobs = n_points * nseeds
        job_point = [mod1(j, n_points) for j in 1:n_jobs]
        rng = Random.Xoshiro(42)
        job_params = [randomize_params(rng, get_p(job_point[j]), scanpoints[job_point[j]]) for j in 1:n_jobs]
        run_job(j) = fit_point(job_point[j], job_params[j])

        @info "Running $n_jobs jobs ($n_points scan points × $nseeds seeds)"
        if isnothing(map_func)
            job_results = Vector{Any}(undef, n_jobs)
            @showprogress Threads.@threads for j in 1:n_jobs
                job_results[j] = run_job(j)
            end
        else
            job_results = collect(map_func(run_job, collect(1:n_jobs)))
        end

        best_results = Vector{Any}(undef, n_points)
        for i in 1:n_points
            seed_indices = [(s - 1) * n_points + i for s in 1:nseeds]
            posteriors = [let p = job_results[j][2]; isnan(p) ? -Inf : p end for j in seed_indices]
            best_results[i] = job_results[seed_indices[argmax(posteriors)]]
        end
        return assemble_profile_results(reshape(best_results, size(scanpoints)), size(scanpoints))
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
    profile_advanced(likelihood, priors, vars_to_scan, params;
            cache_dir=nothing, map_func=nothing, adsel=...) -> NewtrinosResult

Run a profile likelihood scan over a parameter grid, like [`profile`](@ref) but with additional
start-value (seeding) options: `nseeds`, `seed_result(s)`, `sequential`/`start_from`.

Creates grid points defined by `vars_to_scan` with [`generate_scanpoints`](@ref).
At each grid point, all nuisance parameters are optimized with `_profile_advanced` , i.e. via [`find_mle_cached`](@ref). 
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
- `nseeds`: with `nseeds > 1` every grid point is fitted from `nseeds` start values drawn at random from the
  priors (see [`randomize_params`](@ref)) and the best fit is kept. Default: 1.
- `seed_result`, `seed_results`: a [`NewtrinosResult`](@ref) (or several) of an earlier profile whose
  smoothed best-fit parameters (see [`make_seed_interpolants`](@ref)) seed each grid point.
- `sequential`: if `true`, scan points are fitted in BFS order from `start_from` (a NamedTuple of axis
  values, e.g. the best fit of a previous run; default: the central grid point), each point seeded by its
  parent's optimised result. Runs locally (`map_func` is ignored).

# Returns
A [`NewtrinosResult`](@ref) NewtrinosResult(axes, values, meta) 
with the scan grid axes, per-point profiling results, and meta data.
"""
function profile_advanced(likelihood, priors, vars_to_scan, params; cache_dir=nothing, map_func=nothing, adsel=_profile_adsel(map_func), starts=[(;)],
                 nseeds=1, seed_result=nothing, seed_results=nothing, sequential=false, start_from=nothing)
    t1 = time()
    # check if there is actually any variable to be profiled over, or if they are all just Numbers
    if all([isa(priors[var], Number) for var in setdiff(keys(priors), keys(vars_to_scan))])
        return scan(likelihood, priors, vars_to_scan, params)
    end

    values, scanpoints = generate_scanpoints(vars_to_scan, priors)
    mesh_flat = vec(_generate_grid(vars_to_scan, priors)[3])   # grid points in the order of `scanpoints`
    if !isnothing(cache_dir)
        if isdir(cache_dir)
            @info "Reusing cache dir `$(cache_dir)`"
        else
            mkdir(cache_dir)
        end
    end
    seed_params = nothing
    if !isnothing(seed_result)
        itps = make_seed_interpolants(seed_result)
        param_keys = collect(keys(itps))
        axis_keys  = keys(seed_result.axes)
        n_axes = length(axis_keys)
        seed_params = [let coords = Tuple(Float64(mesh_flat[i][d]) for d in 1:n_axes)
            merge(params, NamedTuple{Tuple(param_keys)}(itps[k](coords...) for k in param_keys))
        end for i in eachindex(mesh_flat)]
        @info "Using per-point seeds from smoothed result ($(length(param_keys)) params interpolated)"
    elseif !isnothing(seed_results)
        merged_itps = make_merged_seed_interpolants(seed_results, keys(vars_to_scan))
        merged_keys = collect(keys(merged_itps))
        n_axes = length(keys(vars_to_scan))
        seed_params = [let coords = Tuple(Float64(mesh_flat[i][d]) for d in 1:n_axes)
            nt_vals = [merged_itps[k](coords...) for k in merged_keys]
            merge(params, NamedTuple{Tuple(merged_keys)}(nt_vals))
        end for i in eachindex(mesh_flat)]
        @info "Using merged per-point seeds from $(length(seed_results)) files ($(length(merged_keys)) params interpolated)"
    end

    # Determine BFS center for sequential mode
    center_linear = div(prod(size(scanpoints)), 2) + 1  # default: near center
    if sequential
        if !isnothing(start_from)
            axis_keys = collect(keys(vars_to_scan))
            indices = [begin
                ax = values[findfirst(==(k), collect(keys(vars_to_scan)))]
                argmin(abs.(ax .- Float64(start_from[k])))
            end for k in axis_keys]
            center_linear = LinearIndices(size(scanpoints))[CartesianIndex(Tuple(indices))]
            @info "Sequential profile: starting from grid indices $indices (nearest to start_from)"
        else
            @info "Sequential profile: starting from central grid point (linear index $center_linear)"
        end
        isnothing(map_func) || @warn "sequential=true ignores map_func (distributed workers); running locally"
    end

    res = _profile_advanced(likelihood, scanpoints, params, cache_dir; map_func=map_func, adsel=adsel, starts=starts,
                   nseeds=nseeds, seed_params=seed_params, sequential=sequential, center_linear=center_linear)
    t2 = time()
    meta = Dict("task"=> "profile", "priors"=>priors, "vars_to_scan"=>vars_to_scan, "params"=>params, "starts"=>starts, "exec_time"=>t2-t1, "cache_dir"=>cache_dir, "nseeds"=>nseeds, "sequential"=>sequential)
    add_meta!(meta)
    axes = NamedTuple{tuple(keys(vars_to_scan)...)}(values)
    result = NewtrinosResult(axes=axes, values=res, meta=meta)
end
