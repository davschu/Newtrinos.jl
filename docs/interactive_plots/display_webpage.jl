#=
display_webpage.jl

Loads precomputed scan results (see compute_results.jl / results.jld2) and
serves the interactive Bonito dashboard. Run with:

    julia DAVID_project/display_webpage.jl

Everything here is plain data (Dict/NamedTuple/Vector/String) — this script
does not depend on Newtrinos.jl itself.
=#

using Bonito, Observables, WGLMakie, CairoMakie, JLD2, OrderedCollections

# Load precomputed scan results (see compute_results.jl). Everything here is
# plain data (Dict/NamedTuple/Vector/String) — no Newtrinos.jl dependency.
loaded = JLD2.load(joinpath(@__DIR__, "results.jld2"))
params_results = loaded["params_results"]
scan_params     = loaded["scan_params"]     # Vector{String}
cl_labels       = loaded["cl_labels"]
cl_levels_1d    = loaded["cl_levels_1d"]
cl_levels_2d    = loaded["cl_levels_2d"]

# Only the experiment *names* are needed on the display side (as Dict keys for
# iteration) — no configure closures, since all physics work already happened
# in compute_results.jl.
experiment_sets = Dict(name => nothing for name in keys(params_results))

WGLMakie.activate!()

param_names       = scan_params   # already Vector{String}, loaded from results.jld2
params_dropdown   = Dropdown(param_names)
params_dropdown_y = Dropdown(vcat(["None"], param_names))

# Dataset Checkboxes
dayabay_checkbox  = Bonito.Checkbox(true)
minos_checkbox    = Bonito.Checkbox(true)
kamland_checkbox  = Bonito.Checkbox(true)
deepcore_checkbox = Bonito.Checkbox(true)
#orca_checkbox = Bonito.Checkbox(true)
#superk_checkbox = Bonito.Checkbox(true)
#juno_checkbox = Bonito.Checkbox(true)
combined_checkbox = Bonito.Checkbox(true)

# Mass Ordering Checkboxes
NO_checkbox_cb = Bonito.Checkbox(true)
IO_checkbox_cb = Bonito.Checkbox(true)

# Show all confidence levels at once (overrides the slider when checked)
show_all_cl_checkbox = Bonito.Checkbox(false)

# Confidence-level slider (single level shown at a time; order matches cl_labels/cl_levels_1d/cl_levels_2d)
# Disabled (greyed out) while "show all" is checked — `disabled` takes the
# Checkbox's Observable directly, so it stays reactive without rebuilding the widget.
cl_slider = Bonito.Slider(cl_labels; value=cl_labels[1], disabled=show_all_cl_checkbox.value)

# Best-fit marker toggle
bestfit_checkbox = Bonito.Checkbox(true)

experiments = [
    (id = "Daya Bay", color = :blue,   cb = dayabay_checkbox,  results = params_results["Daya Bay"]),
    (id = "MINOS",    color = :orange, cb = minos_checkbox,    results = params_results["MINOS"]),
    (id = "KamLAND",  color = :green,  cb = kamland_checkbox,  results = params_results["KamLAND"]),
    (id = "IceCube Deepcore", color = :yellow, cb = deepcore_checkbox, results = params_results["IceCube Deepcore"]),
    #(id = "KM3NeT/ORCA", color = :purple, cb = orca_checkbox, results = params_results["KM3NeT/ORCA"]),
    #(id = "JUNO (simulation)", color = :purple, cb = juno_checkbox, results = params_results["JUNO (simulation)"]),
    #(id = "Super-K", color = :purple, cb = superk_checkbox, results = params_results["Super-K"]),
    (id = "Combined", color = :red,    cb = combined_checkbox, results = params_results["Combined"]),
]

#Confidence-region summary table helpers: for a given set of parameters, list best fit +
#all 5 CL ranges for every experiment/ordering. The table is filtered to whichever
#parameter(s) are currently on the plot's axes (wired up reactively in the app below).

function cl_interval(xs, Δchi2, level)
    mask = Δchi2 .<= level
    any(mask) || return (NaN, NaN)
    return (minimum(xs[mask]), maximum(xs[mask]))
end

fmt(x) = string(round(x, sigdigits=4))
fmt_range(lo, hi) = isnan(lo) ? "—" : "[$(fmt(lo)), $(fmt(hi))]"

function build_summary_rows(params_to_show)
    rows = []
    for key in params_to_show
        for exp_name in sort(collect(keys(experiment_sets)))
            for ordering in ("NO", "IO")
                r = params_results[exp_name][ordering].results_1d[key]
                xs = r.axes[1]
                lp = r.values.log_posterior
                Δchi2 = -2 .* (lp .- maximum(lp))
                bestfit = xs[argmin(Δchi2)]
                ranges = [cl_interval(xs, Δchi2, lvl) for lvl in cl_levels_1d]
                push!(rows, (param=key, exp=exp_name, ordering=ordering, bestfit=bestfit, ranges=ranges))
            end
        end
    end
    return rows
end

function build_summary_table(params_to_show)
    rows = build_summary_rows(params_to_show)
    header = DOM.tr(
        DOM.th("Parameter"), DOM.th("Experiment"), DOM.th("Ordering"), DOM.th("Best fit"),
        (DOM.th(lbl) for lbl in cl_labels)...
    )
    body = (
        DOM.tr(
            DOM.td(row.param), DOM.td(row.exp), DOM.td(row.ordering), DOM.td(fmt(row.bestfit)),
            (DOM.td(fmt_range(lo, hi)) for (lo, hi) in row.ranges)...
        )
        for row in rows
    )

    Card(
        DOM.div(
            DOM.h3("Confidence-region summary"; class="section-heading"),
            DOM.div(
                DOM.table(DOM.thead(header), DOM.tbody(body...); class="results-table");
                style = Styles("max-height" => "420px", "overflow-y" => "auto")
            )
        );
        width = "100%",
        backgroundcolor = "white",
        shadow_size = "0 1px 3px",
        shadow_color = "rgba(15, 23, 42, 0.12)",
        border_radius = "12px",
        padding = "20px",
    )
end

#Global page styling: fonts, colors, table/slider/section polish.
#Slate accent palette so the Makie plot colors stay the only "loud" colors on the page.
const ACCENT       = "#334155"   # slate-700
const ACCENT_DARK  = "#1e293b"   # slate-800
const ACCENT_LIGHT = "#64748b"   # slate-500
const PAGE_BG       = "#f5f6f8"
const CARD_BG        = "white"
const CARD_SHADOW    = "rgba(15, 23, 42, 0.12)"
const ERROR_RED       = "#b91c1c"

inter_font = DOM.link(
    rel  = "stylesheet",
    href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap",
)

global_stylesheet = Styles(
    CSS("body",
        "font-family" => "'Inter', -apple-system, BlinkMacSystemFont, sans-serif",
        "background-color" => PAGE_BG,
        "color" => ACCENT_DARK,
        "margin" => "0",
    ),
    CSS(".section-heading",
        "font-size" => "0.75rem",
        "font-weight" => "700",
        "letter-spacing" => "0.06em",
        "text-transform" => "uppercase",
        "color" => ACCENT_LIGHT,
        "margin" => "0 0 6px 0",
    ),
    CSS(".results-table",
        "border-collapse" => "collapse",
        "width" => "100%",
        "font-size" => "13px",
    ),
    CSS(".results-table th",
        "position" => "sticky",
        "top" => "0",
        "background-color" => ACCENT,
        "color" => "white",
        "text-align" => "left",
        "padding" => "8px 10px",
        "font-weight" => "600",
    ),
    CSS(".results-table td",
        "padding" => "6px 10px",
        "border-bottom" => "1px solid #e5e7eb",
    ),
    CSS(".results-table tbody tr:nth-child(even)",
        "background-color" => "#f8fafc",
    ),
    CSS(".results-table tbody tr:hover",
        "background-color" => "#eef2f7",
    ),
    CSS("input[type=range]",
        "accent-color" => ACCENT,
        "width" => "100%",
    ),
    CSS("input[type=range]:disabled",
        "opacity" => "0.4",
    ),
    CSS(".download-btn",
        "width" => "100%",
    ),
    CSS(".error-banner",
        "color" => ERROR_RED,
        "background-color" => "#fef2f2",
        "border" => "1px solid #fecaca",
        "border-radius" => "6px",
        "padding" => "8px 10px",
        "font-size" => "0.8rem",
        "white-space" => "pre-wrap",
    ),
    CSS("hr",
        "border" => "none",
        "border-top" => "1px solid #e2e8f0",
        "margin" => "10px 0",
    ),
)

dropdown_style = Styles(
    CSS(
        "border" => "1px solid #cbd5e1",
        "border-radius" => "6px",
        "padding" => "6px 8px",
        "font-size" => "0.9rem",
        "background-color" => "white",
        "color" => ACCENT_DARK,
    ),
    CSS(":focus",
        "outline" => "none",
        "border-color" => ACCENT,
        "box-shadow" => "0 0 0 2px rgba(51, 65, 85, 0.2)",
    ),
)

button_style = Styles(
    CSS(
        "background-color" => ACCENT,
        "color" => "white",
        "border" => "none",
        "border-radius" => "6px",
        "font-weight" => "600",
        "font-size" => "0.85rem",
    ),
    CSS(":hover", "background-color" => ACCENT_DARK),
)

section_label(text) = DOM.div(text; class = "section-heading")

# Triggers a client-side download of `bytes` by clicking a synthetic <a download> link.
function trigger_download(session, bytes::Vector{UInt8}, filename::String, mime::String)
    data_url = Bonito.to_data_url(bytes, mime)
    evaljs(session, js"""
        const a = document.createElement('a');
        a.href = $(data_url);
        a.download = $(filename);
        document.body.appendChild(a);
        a.click();
        a.remove();
    """)
end

app = App() do session
    # Bonito's own message dispatcher only @warns-and-swallows exceptions raised
    # inside Observable `on(...)` callbacks (e.g. a Button click handler), so
    # errors from the download buttons were previously invisible in the notebook.
    # This banner + `run_safely` make them show up in the page itself.
    export_error = Observable("")
    function run_safely(f, label)
        try
            f()
            export_error[] = ""
        catch e
            @error "$label failed" exception = (e, catch_backtrace())
            export_error[] = "$label failed: $(sprint(showerror, e))"
        end
    end

    fig = Figure(size = (1100, 620))
    ax = Axis(fig[1, 1])

    # Resolve which y-parameter to look up: falls back to *some* other parameter
    # when in 1D mode (its data is never shown — the 2D trace is just invisible)
    param_y_lookup = map(params_dropdown.value, params_dropdown_y.value) do px, py
        (py == "None" || py == px) ? first(p for p in param_names if p != px) : py
    end
    is_1d = map((px, py) -> py == "None" || py == px, params_dropdown.value, params_dropdown_y.value)

    # Summary table filtered to whichever parameter(s) are currently on the plot's axes
    # (just x in 1D mode, x and y in 2D mode) — still lists every experiment/ordering.
    summary_table_content = map(params_dropdown.value, params_dropdown_y.value) do px, py
        params_to_show = (py == "None" || py == px) ? [px] : [px, py]
        build_summary_table(params_to_show)
    end

    # Selected CL level(s) (matching what hlines!/contour! expect): a single-element
    # vector normally, or the full level list when "show all" is checked.
    selected_level_1d = map(cl_slider.index, show_all_cl_checkbox.value) do idx, show_all
        show_all ? cl_levels_1d : [cl_levels_1d[idx]]
    end
    selected_level_2d = map(cl_slider.index, show_all_cl_checkbox.value) do idx, show_all
        show_all ? cl_levels_2d : [cl_levels_2d[idx]]
    end

    # Reference Confidence Level line(s) — tracks the slider (or all levels), only meaningful in 1D mode
    hlines!(ax, selected_level_1d, linestyle = :dot, color = :gray, visible = is_1d)

    function get_pair_grid(results_2d, px, py)
        if haskey(results_2d, (px, py))
            r = results_2d[(px, py)]
            return r.axes[1], r.axes[2], r.values.log_posterior
        elseif haskey(results_2d, (py, px))
            r = results_2d[(py, px)]
            return r.axes[2], r.axes[1], permutedims(r.values.log_posterior)
        else
            return nothing
        end
    end

    # Dynamically build line + contour plots for every experiment & mass ordering
    for exp in experiments
        for (ordering, style, ord_cb) in [("NO", :solid, NO_checkbox_cb), ("IO", :dash, IO_checkbox_cb)]

            # --- 1D Δχ² curve ---
            data_1d = map(params_dropdown.value) do key
                param = exp.results[ordering].results_1d[key]
                x = param.axes[1]
                lp = param.values.log_posterior
                Δchi2 = -2 .* (lp .- maximum(lp))
                return (x = x, y = Δchi2)
            end
            vis_1d = map((a, b, c) -> a && b && c, exp.cb.value, ord_cb.value, is_1d)
            lines!(
                ax,
                map(d -> d.x, data_1d),
                map(d -> d.y, data_1d),
                color = exp.color,
                linestyle = style,
                linewidth = 2,
                visible = vis_1d,
                label = "$(exp.id) ($ordering)"
            )

            # --- 1D best-fit marker ---
            bestfit_1d = map(data_1d) do d
                i = argmin(d.y)
                (x = [d.x[i]], y = [d.y[i]])
            end
            vis_1d_star = map((a, b) -> a && b, vis_1d, bestfit_checkbox.value)
            scatter!(
                ax,
                map(b -> b.x, bestfit_1d),
                map(b -> b.y, bestfit_1d),
                marker = ordering == "NO" ? :star5 : :star8,
                markersize = 16,
                color = exp.color,
                strokecolor = :black,
                strokewidth = 1,
                visible = vis_1d_star,
            )

            # --- 2D Δχ² contour ---
            data_2d = map(params_dropdown.value, param_y_lookup) do px, py
                grid = get_pair_grid(exp.results[ordering].results_2d, px, py)
                isnothing(grid) && return (x = [0.0, 1.0], y = [0.0, 1.0], z = zeros(2, 2))
                xs, ys, lp = grid
                return (x = xs, y = ys, z = -2 .* (lp .- maximum(lp)))
            end
            has_2d = map(params_dropdown.value, param_y_lookup) do px, py
                !isnothing(get_pair_grid(exp.results[ordering].results_2d, px, py))
            end
            vis_2d = map((a, b, c, d) -> a && b && !c && d, exp.cb.value, ord_cb.value, is_1d, has_2d)
            contour!(
                ax,
                map(d -> d.x, data_2d),
                map(d -> d.y, data_2d),
                map(d -> d.z, data_2d),
                levels = selected_level_2d,
                color = exp.color,
                linestyle = style,
                visible = vis_2d,
                linewidth = 2,
            )

            # --- 2D best-fit marker ---
            bestfit_2d = map(data_2d) do d
                idx = argmin(d.z)
                (x = [d.x[idx[1]]], y = [d.y[idx[2]]])
            end
            vis_2d_star = map((a, b) -> a && b, vis_2d, bestfit_checkbox.value)
            scatter!(
                ax,
                map(b -> b.x, bestfit_2d),
                map(b -> b.y, bestfit_2d),
                marker = ordering == "NO" ? :star5 : :star8,
                markersize = 16,
                color = exp.color,
                strokecolor = :black,
                strokewidth = 1,
                visible = vis_2d_star,
            )
        end
    end

    # Legend to the right of the plot, outside the axis
    Legend(fig[1, 2], ax, nbanks = 1)

    # Update axis labels & auto-scale view limits on parameter selection change
    onany(params_dropdown.value, params_dropdown_y.value) do px, py
        ax.xlabel = px
        ax.ylabel = (py == "None" || py == px) ? L"\Delta\chi^2" : py
        reset_limits!(ax)
    end
    notify(params_dropdown.value) # Initialize axis bounds on first load

    # Re-skin the dropdowns for display: Dropdown's CSS is fixed at construction
    # time (baked into its `.style` field), so we wrap the *same* underlying
    # options/option_index Observables from `params_dropdown`/`params_dropdown_y`
    # in freshly-styled Dropdown instances — both views stay in sync since they
    # share the same Observables, but only the styled one is placed in the DOM.
    styled_params_dropdown = Dropdown(
        params_dropdown.options; index = params_dropdown.option_index, style = dropdown_style,
    )
    styled_params_dropdown_y = Dropdown(
        params_dropdown_y.options; index = params_dropdown_y.option_index, style = dropdown_style,
    )

    # --- Download: current plot as vector PDF (Legend included) ---
    # FileIO's `save` infers output format from the file extension, so it needs
    # a real path (not a bare IOBuffer) — write to a tempfile, read the bytes back.
    download_pdf_btn = Bonito.Button("⬇ Download plot (PDF)"; style = button_style, class = "download-btn")
    on(download_pdf_btn.value) do _
        run_safely("PDF export") do
            CairoMakie.activate!()
            path = tempname() * ".pdf"
            try
                save(path, fig)
                trigger_download(session, read(path), "newtrinos_plot.pdf", "application/pdf")
            finally
                WGLMakie.activate!()
                isfile(path) && rm(path)
            end
        end
    end

    # --- Download: currently-selected scan data as JLD2 (all experiments & orderings) ---
    # JLD2.jldsave also requires a real filesystem path (memory-mapped I/O), same tempfile pattern.
    download_jld2_btn = Bonito.Button("⬇ Download data (JLD2)"; style = button_style, class = "download-btn")
    on(download_jld2_btn.value) do _
        run_safely("JLD2 export") do
            px = params_dropdown.value[]
            py = params_dropdown_y.value[]
            params_to_show = (py == "None" || py == px) ? [px] : [px, py]

            export_data = Dict(
                exp_name => Dict(
                    ordering => (
                        results_1d = Dict(k => params_results[exp_name][ordering].results_1d[k] for k in params_to_show),
                        results_2d = Dict(
                            kv for kv in params_results[exp_name][ordering].results_2d
                            if kv[1][1] in params_to_show && kv[1][2] in params_to_show
                        ),
                    )
                    for ordering in ("NO", "IO")
                )
                for exp_name in keys(experiment_sets)
            )

            path = tempname() * ".jld2"
            try
                JLD2.jldsave(path; export_data)
                filename = "newtrinos_scan_$(join(params_to_show, "_")).jld2"
                trigger_download(session, read(path), filename, "application/octet-stream")
            finally
                isfile(path) && rm(path)
            end
        end
    end

    # UI Controls Sidebar Layout — Export moved below the other controls, per feedback
    controls = Card(
        Col(
            section_label("Axes"),
            DOM.div(DOM.b("Parameter (x): "), styled_params_dropdown),
            DOM.div(DOM.b("Parameter (y): "), styled_params_dropdown_y),
            DOM.hr(),
            section_label("Confidence Level"),
            DOM.div(DOM.b("Level: "), cl_slider.value),
            cl_slider,
            DOM.div(show_all_cl_checkbox, " Show all confidence levels"),
            DOM.div(bestfit_checkbox, " Show best-fit ★"),
            DOM.hr(),
            section_label("Datasets"),
            DOM.div(dayabay_checkbox,  " Daya Bay"),
            DOM.div(minos_checkbox,    " MINOS"),
            DOM.div(kamland_checkbox,  " KamLAND"),
            DOM.div(deepcore_checkbox,  " IceCube Deepcore"),
            #DOM.div(orca_checkbox,  " KM3NeT/ORCA"),
            #DOM.div(juno_checkbox,  " JUNO"),
            #DOM.div(superk_checkbox,  " Super Kamiokande"),
            DOM.div(combined_checkbox, " Combined Fit"),
            DOM.hr(),
            section_label("Mass Ordering"),
            DOM.div(NO_checkbox_cb, " Normal (NO) — solid"),
            DOM.div(IO_checkbox_cb, " Inverted (IO) — dashed"),
            DOM.hr(),
            section_label("Export"),
            download_pdf_btn,
            download_jld2_btn,
            DOM.div(
                export_error;
                class = "error-banner",
                style = map(m -> isempty(m) ? "display:none" : "display:block", export_error),
            ),
            ;
            gap = "8px"
        );
        width = "300px",
        backgroundcolor = CARD_BG,
        shadow_size = "0 1px 3px",
        shadow_color = CARD_SHADOW,
        border_radius = "12px",
        padding = "20px",
    )

    plot_card = Card(
        fig;
        width = "100%",
        backgroundcolor = CARD_BG,
        shadow_size = "0 1px 3px",
        shadow_color = CARD_SHADOW,
        border_radius = "12px",
        padding = "16px",
    )

    header = DOM.div(
        DOM.h1(
            "Newtrinos.jl — Global Fit Results";
            style = Styles("font-size" => "1.6rem", "font-weight" => "700", "margin" => "0 0 4px 0", "color" => ACCENT_DARK),
        ),
        DOM.p(
            "Interactive Δχ² scans and confidence regions across experiments.";
            style = Styles("font-size" => "0.95rem", "color" => ACCENT_LIGHT, "margin" => "0"),
        );
        style = Styles("margin-bottom" => "18px"),
    )

    return DOM.div(
        inter_font,
        global_stylesheet,
        header,
        Col(
            Row(plot_card, controls; columns = "1fr 300px"),
            summary_table_content;
            gap = "18px"
        );
        style = Styles(
            "max-width" => "1600px",
            "margin" => "0 auto",
            "padding" => "24px",
        ),
    )
end

server = Server(app, "127.0.0.1", 8081)
print(server)

# `julia display_webpage.jl` exits as soon as the last statement runs, which
# would immediately kill the server started above. Block here so the process
# (and the server) stays alive until you stop it with Ctrl+C.
println("\nPress Ctrl+C to stop the server.")
wait(Condition())
