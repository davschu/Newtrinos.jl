using Newtrinos
using FileIO

"""Configure experiments using their built-in defaults."""
function configure_experiments(experiment_list)
    pairs = (Symbol(lowercase(exp)) => getproperty(getproperty(Newtrinos, Symbol(lowercase(exp))), :configure)() for exp in experiment_list)
    return (; pairs...)
end

"""Configure experiments with a shared physics override."""
function configure_experiments(experiment_list, physics)
    pairs = (Symbol(lowercase(exp)) => getproperty(getproperty(Newtrinos, Symbol(lowercase(exp))), :configure)(physics) for exp in experiment_list)
    return (; pairs...)
end

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
