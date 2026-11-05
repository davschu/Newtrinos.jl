using Newtrinos
using FileIO

"""Save result to JLD2."""
function save_result(result, name)
    FileIO.save(name * ".jld2", Dict("result" => result))
end

"""Plot result and save to PNG. Caller must have `using CairoMakie` in scope."""
function plot_result(result, name, vars_to_scan; title=nothing)
    fig = Figure()
    ax = Axis(fig[1,1])
    plot!(ax, result)
    ax.xlabel = String(collect(keys(vars_to_scan))[1])
    if length(vars_to_scan) == 1
        ax.ylabel = "-2ΔLLH"
    else
        ax.ylabel = String(collect(keys(vars_to_scan))[2])
    end
    if !isnothing(title)
        ax.title = title
    end
    save(name * ".png", fig)
end


"""
    extract_ci(result::Newtrinos.NewtrinosResult, cl::Float64;
               one_sided::Bool = false) -> (lower, upper)

Extract confidence interval boundaries from a 1-D profile likelihood result.

Computes neg2Δllh = -2 * (log_posterior - maximum(log_posterior)) and finds
the x-values at which it crosses `threshold = quantile(Chisq(1), cl)` using
linear interpolation between adjacent grid points.

Returns (lower, upper). For one_sided parameters (abs values bounded at 0),
the best-fit sits at the boundary so only an upper limit is physically
meaningful: returns (0.0, upper).
"""
function extract_ci(result, cl; one_sided=false)
    # Step 1: compute profile statistic
    lp      = result.values.log_posterior         # raw 1-D array
    max_lp  = maximum(lp)
    neg2dlp = -2.0 .* (lp .- max_lp)             # 0 at best-fit, positive everywhere else

    grid    = result.axes[1]                       # 1-D axis vector
    thresh  = quantile(Chisq(1), cl)              # ~1.0 for 68.3%, ~2.706 for 90%

    # Step 2: collect all crossing x-values via linear interpolation
    crossings = Float64[]
    for i in 1:(length(grid) - 1)
        y0, y1 = neg2dlp[i], neg2dlp[i+1]
        x0, x1 = grid[i], grid[i+1]
        # crossing exists when the threshold is bracketed
        if (y0 - thresh) * (y1 - thresh) <= 0
            # linear interpolation: solve y0 + t*(y1-y0) = thresh for t in [0,1]
            t = (thresh - y0) / (y1 - y0)
            push!(crossings, x0 + t * (x1 - x0))
        end
    end

    # Step 3: interpret crossings
    if isempty(crossings)
        # Profile never crosses threshold — sensitivity extends outside the grid
        @warn "No CI crossing found for cl=$cl; grid may be too narrow."
        if one_sided
            return (0.0, grid[end])   # conservative: take grid edge as upper limit
        else
            return (grid[1], grid[end])
        end
    end

    if one_sided
        # For abs parameters: best-fit is at boundary=0, so there is no lower crossing.
        # The interval is [0, first_crossing_above_0].
        upper = minimum(crossings)   # take the first crossing as the upper limit
        return (0.0, upper)
    else
        # Two-sided: take the outermost crossings as lower and upper
        lower = minimum(crossings)
        upper = maximum(crossings)
        return (lower, upper)
    end
end
