module sns_flux

using LinearAlgebra
using Distributions
using CSV
using DataFrames
using ..Newtrinos

@kwdef struct SNSFlux <: Newtrinos.Physics
    params::NamedTuple
    priors::NamedTuple
    assets::NamedTuple
    flux::Function
end
const gf = 1.1663787e-11
const me = 0.510998
const mmu = 105.6
const mtau = 1776.86
const mpi = 139.57
const sw2 = 0.231
const gu = 1/2 - 2*2/3*sw2
const gd = -(1/2) + 2*1/3*sw2
const alph = 1/137
const ep = (mpi^2 - mmu^2) / (2*mpi)

function configure(; 
        exposure, distance, 
        ecut= mmu/2., E_bin_width=0.5, tcut=6000.0, time_bin_width=100.0,
        use_data=true, flux_folder=joinpath(@__DIR__, "..", "experiments", "coherent", "coherent_2020", "csi", "snsFlux2D_CSV"), 
        beam_power=1.4, proton_energy=1.0)
    # Call get_assets and assign the result to `assets`
    assets = get_assets(; use_data, exposure, distance, beam_power, proton_energy, flux_folder, ecut, E_bin_width, tcut, time_bin_width)

    # Return the configured SNSFlux object
    return SNSFlux(
        params = get_params(use_data),
        priors = get_priors(use_data),
        assets = assets,
        flux = get_flux(use_data, assets),
    )
end

function get_params(use_data)
    if use_data
        (
            norm_mu = 1.0,
            norm_e = 1.0,
            norm_mu_bar = 1.0,
            norm_e_bar = 1.0,
        )
    else
        (
            sns_nu_per_POT = 0.09,
        )
    end
end

function get_priors(use_data)
    if use_data
        (
            norm_mu = truncated(Normal(1.0, 0.1), 0.7, 1.3),
            norm_e = truncated(Normal(1.0, 0.1), 0.7, 1.3),
            norm_mu_bar = truncated(Normal(1.0, 0.1), 0.7, 1.3),
            norm_e_bar = truncated(Normal(1.0, 0.1), 0.7, 1.3),
        )
    else
        (
            sns_nu_per_POT = truncated(Normal(0.09, 0.009), 0.05, 0.15),
        )
    end
end

function get_assets(; use_data, exposure, distance,
                        beam_power, proton_energy,
                        flux_folder, ecut, E_bin_width, tcut, time_bin_width)
    if use_data
        #@info "Loading SNS Flux data"
        # Function to read and rebin flux data from CSV files
        function read_flux_data(file_path, ecut, tcut, time_bin_width)
            # Read CSV; first column = energy, remaining columns = flux values
            data = CSV.read(file_path, DataFrame; header=true, normalizenames=false)

            # Energy centers (first column)
            energy_centers = data[:, 1]

            # Time centers: column headers (strings) converted to numeric type matching energy_centers
            time_centers = parse.(eltype(energy_centers), names(data)[2:end])

            # Apply energy and time cuts
            energy_mask = energy_centers .<= ecut
            time_mask   = time_centers   .<= tcut

            # Filter centers
            energy_centers = energy_centers[energy_mask]
            time_centers   = time_centers[time_mask]

            # Extract weight submatrix
            weights = Matrix(data[:, 2:end])
            weights = weights[energy_mask, time_mask]

            # Compute bin widths
            dE = diff(energy_centers)
            dE = vcat(dE, dE[end])
            dt = diff(time_centers)
            dt = vcat(dt, dt[end])

            # Scale weights by bin widths (area-normalized)
            weights .*= dE .* dt'

            # --- Rebin time dimension if requested ---
            if time_bin_width !== nothing
                # Construct new time bin edges
                tmin, tmax = first(time_centers), last(time_centers)
                new_edges = collect(tmin:time_bin_width:tmax)
                new_centers = (new_edges[1:end-1] .+ new_edges[2:end]) ./ 2

                # Initialize rebinned weight matrix
                rebinned_weights = zeros(size(weights, 1), length(new_centers))

                # Assign each original time bin to a new time bin index
                bin_index = searchsortedlast.(Ref(new_edges), time_centers)

                for (i, bi) in enumerate(bin_index)
                    if bi < length(new_edges)
                        rebinned_weights[:, bi] .+= weights[:, i]
                    end
                end

                # Replace original time data
                time_centers = new_centers
                weights = rebinned_weights

                # Update dt for new bins
                dt = fill(time_bin_width, length(time_centers))
            end

            return (energy_centers, time_centers, weights)
        end
        s_per_gwh = 3.6e6 / beam_power  # seconds per GWh
        distance = distance * 5.07e10 # in MeV^-1
        eta = exposure * s_per_gwh / (4 * pi * distance^2)
        # Define file paths for the flux components
        file_path_mu = joinpath(flux_folder, "convolved_energy_time_of_nu_mu.csv")
        file_path_e = joinpath(flux_folder, "convolved_energy_time_of_nu_e.csv")
        file_path_mu_bar = joinpath(flux_folder, "convolved_energy_time_of_anti_nu_mu.csv")
        file_path_e_bar = joinpath(flux_folder, "convolved_energy_time_of_anti_nu_e.csv")

        # Load flux data for each component
        E_mu, T_mu, flux_mu = read_flux_data(file_path_mu, ecut, tcut, time_bin_width)
        E_e, T_e, flux_e = read_flux_data(file_path_e, ecut, tcut, time_bin_width)
        E_mu_bar, T_mu_bar, flux_mu_bar = read_flux_data(file_path_mu_bar, ecut, tcut, time_bin_width)
        E_e_bar, T_e_bar, flux_e_bar = read_flux_data(file_path_e_bar, ecut, tcut, time_bin_width)

        # Store the loaded data as assets
        assets = (;
            E = E_mu,
            T = T_mu,
            flux_mu = eta * flux_mu,
            flux_e = eta * flux_e,
            flux_mu_bar = eta * flux_mu_bar,
            flux_e_bar = eta * flux_e_bar,
        )
    else
        #@info "Priming SNS Flux functions"
        # Uniform grid in E
        npts = Int(round((ecut - 0.5)/E_bin_width))
        E = collect(range(0.5, stop=ecut, length=max(npts, 2)))
        distance = distance * 5.07e10 # in MeV^-1
        POT_per_GWhr = (3.6e12 / (proton_energy * 1.602e-10))
        eta = POT_per_GWhr * exposure / (4 * π * (distance)^2)
        E0=ep
        E_mu, flux_mu = flux_nu_mu(E, E0, eta, E_bin_width)
        E_e, flux_e = flux_nu_e(E, eta, ecut, E_bin_width)
        E_mu_bar, flux_mu_bar = flux_nu_mu_bar(E, eta, ecut, E_bin_width)
        assets = (;
            E = E_mu,
            flux_mu,
            flux_e,
            flux_mu_bar,
        )
    end
end


# Returns (E, flux) arrays for nu_mu
function flux_nu_mu(E, E0, eta, bin_width)
    dE = diff(E)
    @assert all(isapprox.(dE, dE[1]; atol=1e-8))
    avg_dE = mean(dE)
    dE = fill(avg_dE, length(E))

    # Element type that can hold eta (possibly Dual)
    T = promote_type(eltype(E), typeof(eta))
    flux = zeros(T, length(E))

    mask = isapprox.(E, E0; atol=avg_dE/2)
    if !any(mask)
        return (E, flux)
    end

    sigma = avg_dE
    weights = exp.(-0.5 .* ((E[mask] .- E0) ./ sigma).^2)
    norm = sum(weights .* dE[mask])
    if norm != 0
        weights ./= norm
    end
    flux[mask] .= eta .* weights .* dE[mask]
    return (E, flux)
end

function flux_nu_e(E, eta, Emax, bin_width)
    # Element type that can hold eta (possibly Dual)
    T = promote_type(eltype(E), typeof(eta))
    flux = zeros(T, length(E))
    mask = (E .>= 0) .& (E .<= Emax)
    flux[mask] .= eta .* 192 .* (E[mask].^2 ./ mmu^3) .* (0.5 .- E[mask] ./ mmu)
    return (E, flux .* T(bin_width))
end

function flux_nu_mu_bar(E, eta, Emax, bin_width)
    # Element type that can hold eta (possibly Dual)
    T = promote_type(eltype(E), typeof(eta))
    flux = zeros(T, length(E))
    mask = (E .>= 0) .& (E .<= Emax)
    flux[mask] .= eta .* 64 .* (E[mask].^2 ./ mmu^3) .* (0.75 .- E[mask] ./ mmu)
    return (E, flux .* T(bin_width))
end

function get_flux(use_data, assets)
    if use_data
        # Extract fluxes and energy/time centers from assets
        E = assets.E
        T = assets.T
        flux_mu = assets.flux_mu
        flux_e = assets.flux_e
        flux_mu_bar = assets.flux_mu_bar
        flux_e_bar = assets.flux_e_bar

        # Return the closure directly
        return function (params)
            # Extract normalization parameters
            norm_mu = params.norm_mu
            norm_e = params.norm_e
            norm_mu_bar = params.norm_mu_bar
            norm_e_bar = params.norm_e_bar

            # Apply normalizations
            flux_mu_scaled = flux_mu .* norm_mu
            flux_e_scaled = flux_e .* norm_e
            flux_mu_bar_scaled = flux_mu_bar .* norm_mu_bar
            flux_e_bar_scaled = flux_e_bar .* norm_e_bar

            # Compute total flux
            total_flux = flux_mu_scaled .+ flux_e_scaled .+ flux_mu_bar_scaled .+ flux_e_bar_scaled

            return (;
                E,  # Energy grid
                T,  # Time grid
                total_flux,
                flux_e = flux_e_scaled,
                flux_mu = flux_mu_scaled,
                flux_mu_bar = flux_mu_bar_scaled,
                flux_e_bar = flux_e_bar_scaled,
            )
        end
    else
        # Extract fluxes and energy centers from assets
        E = assets.E
        flux_mu = assets.flux_mu
        flux_e = assets.flux_e
        flux_mu_bar = assets.flux_mu_bar

        # Return the closure directly
        return function (params)
            # Use params for nu_per_POT
            nu_per_POT_val = params.sns_nu_per_POT

            # Scale fluxes by nu_per_POT
            flux_mu_scaled = flux_mu .* nu_per_POT_val
            flux_e_scaled = flux_e .* nu_per_POT_val
            flux_mu_bar_scaled = flux_mu_bar .* nu_per_POT_val

            # Compute total flux
            total_flux = flux_mu_scaled .+ flux_e_scaled .+ flux_mu_bar_scaled

            return (;
                E,  # Energy grid
                total_flux,
                flux_e = flux_e_scaled,
                flux_mu = flux_mu_scaled,
                flux_mu_bar = flux_mu_bar_scaled,
            )
        end
    end
end

end