# Plots of the `sensitivity.jl --method fisher` results: the static overview figure and the panels it shares with the on-demand explorer
# (nsi/deepcore_8year_analysis/results/sensitivity_results/make_explorer_html.jl). `include`d by sensitivity.jl; the caller loads CairoMakie.
#
#   D = sensitivity_data(results)          results: the Dict that sensitivity.jl saves (or the loaded .jld2); D: plain arrays, parameters x truths
#   static_overview(D; top_n = 25)         one figure
#
# Quantities (all for the tested parameter a and a nuisance j, at one truth):
#   σ_marg, σ_cond         expected uncertainty with all other parameters free / fixed, drawn over the prior width
#   impact / σ_marg(a)     shift of a for a +1 σ_marg(j) move of j (first order, the other parameters re-fitted) in units of σ_marg(a); equal to the correlation of a and j
#   correlation            the correlation matrix of the Fisher covariance at the truth

using Statistics

const C_MARG = RGBf(0x2a / 255, 0x6f / 255, 0xbb / 255)
const C_COND = RGBf(0xd9 / 255, 0x82 / 255, 0x2b / 255)
const DIVERGING = cgrad([C_MARG, RGBf(0.95, 0.95, 0.95), C_COND])    # two hues with a neutral middle
const SEQUENTIAL = cgrad([RGBf(0.95, 0.95, 0.95), C_MARG])

"""
    sensitivity_data(results) -> NamedTuple

Arrange the `fisher` results of `sensitivity.jl` for plotting: all free parameters (`names`, vector parameters expanded) x the truths that have a result
(`truth_index`); a parameter dropped at a truth is NaN there. Fields: `names`, `tested` (names of the tested parameters), `ti` (their positions in `names`),
`tv` (tested parameters x truths: the injected values), `corr` (parameters x parameters x truths), `smarg`, `scond` (parameters x truths), `width`, `truth_index`,
`failed` (indices of the truths without a result), `setup`, `hypothesis`, `information`, `truth_description`.
"""
function sensitivity_data(results)
    fi, truths, st = results["fisher"], results["truths"], results["settings"]
    ok = [i for i in eachindex(fi) if !isnothing(fi[i])]
    isempty(ok) && error("no truth has a result")
    names = Symbol[]
    for i in ok, n in fi[i].names
        n in names || push!(names, n)
    end
    n, K = length(names), length(ok)
    pos = Dict(nm => j for (j, nm) in enumerate(names))
    corr = fill(NaN, n, n, K); smarg = fill(NaN, n, K); scond = fill(NaN, n, K); width = fill(NaN, n)
    for (k, i) in enumerate(ok)
        r = fi[i]
        ix = [pos[nm] for nm in r.names]
        corr[ix, ix, k] = r.corr; smarg[ix, k] = r.σ_marg; scond[ix, k] = r.σ_cond; width[ix] = r.width
    end
    tested = [Symbol(t) for t in st["tested"] if Symbol(t) in names]
    (; names = string.(names), tested, ti = [pos[t] for t in tested], tv = Float64[truths[i][t] for t in tested, i in ok],
       corr, smarg, scond, width, truth_index = ok, failed = [i for i in eachindex(fi) if isnothing(fi[i])],
       setup = st["setup"], hypothesis = st["hypothesis"], information = st["information"], truth_description = get(st, "truth_description", ""))
end

_finite(v) = filter(isfinite, v)
_stat(f, M) = [isempty(_finite(M[i, :])) ? NaN : f(_finite(M[i, :])) for i in axes(M, 1)]   # statistic of every row of M over its finite entries

"""Text for a truth of `D` (index `k` among the shown truths): the tested values for up to 3 tested parameters, else the truth number (`#12`)."""
truth_label(D, k) = length(D.tested) <= 3 ? join(("$(t) = $(round(D.tv[a, k]; sigdigits = 3))" for (a, t) in enumerate(D.tested)), ", ") : "#$(D.truth_index[k])"

# the tested parameters' own σ-panels and impact panels are arranged in rows of up to 4
_panel_cols(nt) = min(nt, 4)

"""The `top_n` nuisances with the largest |impact / σ_marg| on any tested parameter at any truth, strongest first (all of them are in the summary text)."""
function shown_nuisances(D; top_n)
    score(j) = maximum(abs, _finite(vec(D.corr[D.ti, j, :])); init = 0.0)
    first(sort([j for j in eachindex(D.names) if !(j in D.ti)]; by = j -> -score(j)), top_n)
end
selection_note(D, sel) = length(sel) < length(D.names) - length(D.ti) ? "top $(length(sel)) of $(length(D.names) - length(D.ti)) nuisances by largest |impact|" : ""

panel_rows(n) = 24 * n + 170   # pixel height of a panel with n labelled rows

"""σ_marg (solid) and σ_cond (dashed) over the prior width of every tested parameter against its own truth value, in the layout cell `gl`: lines for one tested parameter,
markers for several (the spread at a value is the dependence on the other coordinates). `marker`: an `Observable` truth index to mark (explorer frames)."""
function sigma_truth_panels!(gl, D; marker = nothing)
    nt, K = length(D.ti), length(D.truth_index)
    cols = _panel_cols(nt); ms = K > 40 ? 4 : 7
    for a in 1:nt
        r, c = fldmod1(a, cols)
        ax = Axis(gl[r, c]; title = nt == 1 ? "tested parameter: σ / prior width vs. its truth" : "σ($(D.tested[a])) / prior width",
                  xlabel = "injected truth $(D.tested[a])", ylabel = c == 1 ? "σ($(D.tested[a])) / prior width" : "", yscale = log10)
        x = D.tv[a, :]; w = D.width[D.ti[a]]
        ym, yc = D.smarg[D.ti[a], :] ./ w, D.scond[D.ti[a], :] ./ w
        hlines!(ax, 1.0; color = :gray, linestyle = :dot, linewidth = 1.5)   # σ = prior width: no sensitivity within the prior range
        if nt == 1
            p = sortperm(x)
            scatterlines!(ax, x[p], ym[p]; color = C_MARG, linewidth = 2, markersize = ms, label = "σ_marg (others free)")
            scatterlines!(ax, x[p], yc[p]; color = C_COND, linewidth = 2, markersize = ms, linestyle = :dash, label = "σ_cond (others fixed)")
        else
            scatter!(ax, x, ym; color = C_MARG, markersize = ms, label = "σ_marg")
            scatter!(ax, x, yc; color = C_COND, markersize = ms, marker = :utriangle, label = "σ_cond")
        end
        fin = _finite(vcat(ym, yc))
        isempty(fin) || ylims!(ax, minimum(fin) / 1.5, 12 * maximum(fin))   # head room for the legend
        a == 1 && axislegend(ax; position = :ct, framevisible = false)
        if !isnothing(marker)
            if nt == 1
                vlines!(ax, lift(k -> x[k], marker); color = :black, linestyle = :dot)
            else
                scatter!(ax, lift(k -> Point2f[(x[k], ym[k])], marker); color = :transparent, strokecolor = :black, strokewidth = 2, markersize = ms + 9)
            end
        end
    end
end

"""Dodged horizontal bars of σ_marg and σ_cond over the prior width for the nuisances `order`; `whiskers` (min to max over the truths) are optional."""
function nuisance_bars!(ax, D, order, marg, cond; whiskers = nothing)
    xs = 1:length(order)
    barplot!(ax, xs .- 0.2, marg; width = 0.4, direction = :x, color = C_MARG, label = "σ_marg")
    barplot!(ax, xs .+ 0.2, cond; width = 0.4, direction = :x, color = C_COND, label = "σ_cond")
    if !isnothing(whiskers)
        rangebars!(ax, xs .- 0.2, whiskers.marg_lo, whiskers.marg_hi; direction = :x, color = :black, linewidth = 1, whiskerwidth = 5)
        rangebars!(ax, xs .+ 0.2, whiskers.cond_lo, whiskers.cond_hi; direction = :x, color = :black, linewidth = 1, whiskerwidth = 5)
    end
    ax.yticks = (xs, [j in D.ti ? D.names[j] * " (tested)" : D.names[j] for j in order]); ax.yreversed = true; ax.xlabel = "σ / prior width"
    axislegend(ax; position = :rb, framevisible = false)
end

_nan0(v) = replace(v, NaN => 0.0)

"""impact / σ_marg of the nuisances `sel` on every tested parameter at the truth `k[]` (an `Observable` index), bars coloured by sign, in the layout cell `gl`."""
function impact_bar_panels!(gl, D, sel, k)
    nt = length(D.ti); cols = _panel_cols(nt)
    for a in 1:nt
        r, c = fldmod1(a, cols)
        ax = Axis(gl[r, c]; title = nt == 1 ? "impact of the nuisances on $(D.tested[a]) at the truth" : string(D.tested[a]), xlabel = _impact_xlabel(D, a, nt))
        v = lift(k -> _nan0(D.corr[D.ti[a], sel, k]), k)
        barplot!(ax, 1:length(sel), v; direction = :x, color = lift(v -> [x > 0 ? C_COND : C_MARG for x in v], v))
        ax.yticks = (1:length(sel), D.names[sel]); ax.yreversed = true; xlims!(ax, -1, 1)
        nt > 1 && (ax.xticks = [-1.0, 0.0, 1.0])   # narrow panels: no colliding tick labels
        c == 1 || (ax.yticklabelsvisible = false)
    end
end

"""Heatmap of a correlation matrix `M` of the parameters `keep` with its own colorbar column in the layout cell `gl` (so no label overlaps a neighbour)."""
function matrix_panel!(gl, M, D, keep; title, colormap, colorrange, label)
    n = length(keep)
    ax = Axis(gl[1, 1]; title)
    hm = heatmap!(ax, 1:n, 1:n, M; colormap, colorrange)
    ax.xticks = (1:n, D.names[keep]); ax.yticks = (1:n, D.names[keep]); ax.xticklabelrotation = π / 2; ax.yreversed = true
    Colorbar(gl[1, 2], hm; label)
    colgap!(gl, 12)
    ax
end

_impact_xlabel(D, a, nt) = nt == 1 ? "impact / σ_marg($(D.tested[a]))" : "impact / σ_marg"

"""Per tested parameter the range (min to max over the truths) of the signed impact / σ_marg of every nuisance in `sel`; a bar crossing zero is a sign flip. Cell `gl`."""
function impact_range_panels!(gl, D, sel)
    nt = length(D.ti); cols = _panel_cols(nt)
    for a in 1:nt
        r, c = fldmod1(a, cols)
        imp = D.corr[D.ti[a], sel, :]   # nuisances x truths
        ax = Axis(gl[r, c]; title = nt == 1 ? "impact of the nuisances on $(D.tested[a]) over the truths" : string(D.tested[a]),
                  xlabel = _impact_xlabel(D, a, nt) * ", range over the truths")
        vlines!(ax, 0.0; color = :gray, linestyle = :dot)
        rangebars!(ax, 1:length(sel), _stat(minimum, imp), _stat(maximum, imp); direction = :x, color = C_MARG, linewidth = nt == 1 ? 7 : 5, whiskerwidth = 0)
        ax.yticks = (1:length(sel), D.names[sel]); ax.yreversed = true; xlims!(ax, -1, 1)
        nt > 1 && (ax.xticks = [-1.0, 0.0, 1.0])   # narrow panels: no colliding tick labels
        c == 1 || (ax.yticklabelsvisible = false)
    end
end

"""One truth (`point:`, `bf:`, `nominal`): nothing to resolve over the truths, so the actual numbers: row 1 σ / prior width of the tested parameter(s) and the nuisances,
row 2 the correlation matrix at the truth (signed) | impact / σ_marg of the nuisances on the tested parameter(s)."""
function single_truth_overview(D, sel, keep, note)
    rel(j) = (v = D.smarg[j, 1] / D.width[j]; isfinite(v) ? v : 0.0)
    order = vcat(D.ti, sort(sel; by = j -> -rel(j)))   # the tested parameters first, then the nuisances by σ / width
    h1, h2 = max(420, panel_rows(length(order))), max(700, panel_rows(length(keep)) + 160)
    fig = Figure(size = (1900, h1 + h2 + 80))
    Label(fig[0, 1:2], "$(D.setup), $(D.hypothesis), tested $(join(D.tested, ", ")), truth $(truth_label(D, 1)), $(D.information) information" * (isempty(note) ? "" : " ($note)");
          font = :bold, fontsize = 16)
    ax = Axis(fig[1, 1:2]; title = "σ / prior width at the truth (tested parameter(s) first)")
    nuisance_bars!(ax, D, order, _nan0(D.smarg[order, 1] ./ D.width[order]), _nan0(D.scond[order, 1] ./ D.width[order]))
    gl = GridLayout(); fig[2, 1] = gl
    matrix_panel!(gl, D.corr[keep, keep, 1], D, keep; title = "correlation matrix at the truth", colormap = DIVERGING, colorrange = (-1, 1), label = "correlation")
    gl = GridLayout(); fig[2, 2] = gl
    impact_bar_panels!(gl, D, sel, Observable(1))
    rowsize!(fig.layout, 1, Auto(h1)); rowsize!(fig.layout, 2, Auto(h2))
    colgap!(fig.layout, 50); rowgap!(fig.layout, 30)
    fig
end

"""Overview figure: row 1 σ / prior width (tested parameters vs. their truths | nuisances, mean over the truths with min to max whiskers), row 2 RMS of the
correlation over the truths (sign not shown) | range of the impact / σ_marg of the nuisances. At most `top_n` nuisances are drawn (the ones with the largest |impact|)."""
function static_overview(D; top_n = 25)
    nt, K = length(D.ti), length(D.truth_index)
    sel = shown_nuisances(D; top_n); keep = sort(vcat(D.ti, sel)); note = selection_note(D, sel)
    K == 1 && return single_truth_overview(D, sel, keep, note)
    h1 = max(520, panel_rows(length(sel)), 380 * cld(nt, 4)); h2 = max(700, panel_rows(length(keep)) + 160)
    fig = Figure(size = (1900, h1 + h2 + 80))
    failed = isempty(D.failed) ? "" : "; no result at truth $(join(D.failed, ", "))"
    Label(fig[0, 1:2], "$(D.setup), $(D.hypothesis), tested $(join(D.tested, ", ")): $K truths, $(D.information) information" *
          (isempty(note) ? "" : " ($note)") * failed; font = :bold, fontsize = 16)

    gl = GridLayout(); fig[1, 1] = gl
    sigma_truth_panels!(gl, D)
    wmarg = D.smarg[sel, :] ./ D.width[sel]; wcond = D.scond[sel, :] ./ D.width[sel]
    o = sortperm(_stat(mean, wmarg); rev = true)
    ax = Axis(fig[1, 2]; title = "nuisances: σ / prior width (mean over truths, whiskers min to max)")
    nuisance_bars!(ax, D, sel[o], _stat(mean, wmarg[o, :]), _stat(mean, wcond[o, :]);
                   whiskers = (; marg_lo = _stat(minimum, wmarg[o, :]), marg_hi = _stat(maximum, wmarg[o, :]),
                                 cond_lo = _stat(minimum, wcond[o, :]), cond_hi = _stat(maximum, wcond[o, :])))

    gl = GridLayout(); fig[2, 1] = gl
    rms = [(v = _finite(abs2.(D.corr[a, b, :])); isempty(v) ? NaN : sqrt(mean(v))) for a in keep, b in keep]
    matrix_panel!(gl, rms, D, keep; title = "RMS of the correlation over the truths (sign not shown)", colormap = SEQUENTIAL, colorrange = (0, 1), label = "RMS correlation")
    gl = GridLayout(); fig[2, 2] = gl
    impact_range_panels!(gl, D, sel)

    rowsize!(fig.layout, 1, Auto(h1)); rowsize!(fig.layout, 2, Auto(h2))
    colgap!(fig.layout, 50); rowgap!(fig.layout, 30)
    fig
end
