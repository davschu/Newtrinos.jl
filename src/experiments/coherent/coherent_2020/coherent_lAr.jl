module coherent_lAr

using DataFrames
using CSV
using Distributions
using Interpolations
using LinearAlgebra
using Statistics
using DataStructures
using BAT
using CairoMakie
using Logging
using StatsBase
using ..Helpers
import ..Newtrinos

@kwdef struct COHERENT_LAR <: Newtrinos.Experiment
    physics::NamedTuple
    params::NamedTuple
    priors::NamedTuple
    assets::NamedTuple
    forward_model::Function
    plot::Function
end

function configure(physics)
    assets = get_assets(physics)
    # Dynamically build params/priors for cevns_xsec using isotope list (with Rn_nom)
    cevns_params, cevns_priors = Newtrinos.cevns_xsec.build_params_and_priors(assets.isotopes)
    # Reconfigure cevns_xsec with correct params/priors
    cevns_xsec = Newtrinos.cevns_xsec.configure(cevns_params, cevns_priors)
    # Update physics NamedTuple with new cevns_xsec
    physics = (;physics.sns_flux, cevns_xsec)
    return COHERENT_LAR(
        physics = physics,
        params = get_params(assets.ss_bkg_nom, assets.pbrn_nom, assets.delbrn_nom),
        priors = get_priors(assets.ss_bkg_nom, assets.pbrn_nom, assets.delbrn_nom),
        assets = assets,
        forward_model = get_forward_model(physics, assets),
        plot = get_plot(physics, assets)
    )
end

function get_params(ss_bkg_nom, pbrn_nom, delbrn_nom)
    params = (
        coherent_lar_qfa_a = 0.246,  # QF polynomial coefficients
        coherent_lar_qfa_b = 0.00078,
        coherent_lar_mass = 24.4,  # kg
        pbrn_norm= pbrn_nom,  # Normalization factor for BRN
        delbrn_norm= delbrn_nom,  # Normalization factor for delBRN
        ss_bkg_norm= ss_bkg_nom,  # Normalization factor for SS background
        )
end

# TODO!
function get_priors(ss_bkg_nom, pbrn_nom, delbrn_nom)
    priors = (
        coherent_lar_qfa_a = Normal(0.246, 0.006),
        coherent_lar_qfa_b = Normal(0.00078, 0.00009),
        coherent_lar_mass = truncated(Normal(24.4, 0.61), 0.0, 24.4 + 3 * 0.61),
        pbrn_norm= truncated(Normal(pbrn_nom, 0.32 * pbrn_nom), 0.0, pbrn_nom + 3 * 0.32 * pbrn_nom),  # Normalization factor for BRN
        delbrn_norm= truncated(Normal(delbrn_nom, 1.0 * delbrn_nom), 0.0, delbrn_nom + 3 * 1.0 * delbrn_nom),  # Normalization factor for delBRN
        ss_bkg_norm= truncated(Normal(ss_bkg_nom, 0.008 * ss_bkg_nom), 0.0, ss_bkg_nom + 3 * 0.008 * ss_bkg_nom),  # Normalization factor for SS background
        )
end

function get_assets(physics, datadir = @__DIR__)
    @info "Loading coherent lAr data"


    er_edges = collect(3:0.5:300) # keVnr
    isotopes = [
        (fraction=1.0, mass=37.3e3, Z=18, N=22, Rn_key=:Rn_Ar, Rn_nom=4.1039) # Ar-37
    ] # List of isotopes with [fraction, Nuclear mass (GeV), Z, N=A-Z, Rn_key]
    Nt = (1.0/0.039948) * 6.023e+23
    
    resolution = 0.58  # a/Eee and b*Eee
    # Import Data
    eff_data = CSV.read(joinpath(datadir, "lAr/CENNS10AnlAEfficiency.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: keVee, keVnr, efficiency
    pbrn = CSV.read(joinpath(datadir, "lAr/brnpdf.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: keVee, F90, t, counts/bin/6.12GWh
    delbrn = CSV.read(joinpath(datadir, "lAr/delbrnpdf.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: keVee, F90, t, counts/bin/6.12GWh
    ss_bkg = CSV.read(joinpath(datadir, "lAr/bkgpdf.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: keVee, F90, t, counts/bin/6.12GWh
    observed_df = CSV.read(joinpath(datadir, "lAr/datanobkgsub.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: keVee, F90, t, counts/bin/6.12GWh
    
    # Define bin centers as the average of consecutive edges
    er_centers = collect((er_edges[1:end-1] + er_edges[2:end]) ./ 2) # Reconstruct bin centers from edges

    # Match data: output bin centers 5, 15, 25, ..., 55
    out_width = 10.0
    out_edges = collect(0.0:out_width:120.0)         # Edges: [0, 10, 20, ..., 120]
    out_centers = midpoints(out_edges)              # Centers: [5, 15, 25, ..., 115]
    
    # Bin observed data into out_centers using Helpers.rebin
    observed = Helpers.rebin(observed_df, out_edges; var_col=1, count_col=4)

    # Bin background PDFs ONCE
    pbrn_binned = Helpers.rebin(pbrn, out_edges; var_col=1, count_col=4)
    delbrn_binned = Helpers.rebin(delbrn, out_edges; var_col=1, count_col=4)
    ss_bkg_binned = Helpers.rebin(ss_bkg, out_edges; var_col=1, count_col=4)
    # Get initial nominal value for Bkg normalizations
    ss_bkg_nom = sum(ss_bkg_binned)
    @info "Initial SS background normalization: $ss_bkg_nom"
    pbrn_nom = sum(pbrn_binned)
    @info "Initial BRN background normalization: $pbrn_nom"
    delbrn_nom = sum(delbrn_binned)
    @info "Initial delBRN background normalization: $delbrn_nom"
    distance = 27.5 # m
    exposure = 6.12 # GWh


    assets = (;
        observed,
        er_edges,
        er_centers,
        out_centers,
        isotopes,
        Nt,
        resolution,
        eff_data,
        pbrn_binned,
        pbrn_nom,
        delbrn_binned,
        delbrn_nom,
        ss_bkg_binned,
        ss_bkg_nom,
        distance,
        exposure,
    )
end

function qf(er_centers, params)
    a = params.coherent_lar_qfa_a
    b = params.coherent_lar_qfa_b
    vals = @. (a + b * er_centers) * er_centers
    z = vals isa AbstractArray ? zero(eltype(vals)) : zero(vals)
    return max.(vals, z)  # clamp to >= 0 with typed zero
end

# AD-safe piecewise-linear interpolation for efficiency
function eff(keVee, assets)
    x = assets.eff_data[:, 1]  # Float64 knots
    y = assets.eff_data[:, 3]  # Float64 values
    T = eltype(keVee)
    vals = similar(keVee, promote_type(T, Float64))
    xmin, xmax = first(x), last(x)
    for i in eachindex(keVee)
        u = keVee[i]
        if u <= xmin
            vals[i] = zero(u)                        # typed zero
        elseif u >= xmax
            vals[i] = one(u) * y[end]               # promote to T
        else
            # Find rightmost x[j] ≤ u
            j = findlast(x .<= u)
            if j === nothing
                vals[i] = zero(u)
            elseif j == length(x)
                vals[i] = one(u) * y[end]
            else
                x0 = x[j]; x1 = x[j + 1]
                y0 = y[j]; y1 = y[j + 1]
                t = (u - x0) / (x1 - x0)
                vals[i] = one(u) * y0 + t * (y1 - y0)  # linear interp, returns T
            end
        end
    end
    return vals
end

function sigma_keVee(keVee, assets)
    a = assets.resolution
    return @. a * sqrt(keVee)  # keVee ≥ 0 due to qf clamp; works with Duals
end

function construct_response_matrix(params, assets)
    er_centers = assets.er_centers
    out_centers = assets.out_centers
    n_out = length(out_centers)
    n_er = length(er_centers)

    # Precompute
    keVee = qf(er_centers, params)
    sigma = sigma_keVee(keVee, assets)
    eff_vals = eff(keVee, assets)

    # First column to set element type
    col1 = begin
        d = Distributions.Normal(keVee[1], sigma[1])
        v = Distributions.pdf.(d, out_centers)
        s = sum(v)
        s == 0 ? zero.(v) : (v ./ s)
    end
    A = Array{eltype(col1)}(undef, n_out, n_er)
    A[:, 1] = eff_vals[1] .* col1

    for j in 2:n_er
        d = Distributions.Normal(keVee[j], sigma[j])
        v = Distributions.pdf.(d, out_centers)
        s = sum(v)
        v = s == 0 ? zero.(v) : (v ./ s)
        A[:, j] = eff_vals[j] .* v
    end

    # Ensure non-negative with typed zero
    T = eltype(A)
    return max.(A, zero(T))
end

function build_rate_matrix(er_centers, enu_centers, nupar, physics, params, Rn_key)
    physics.cevns_xsec.diff_xsec_lar(er_centers, enu_centers, params, nupar, Rn_key)
end

function get_expected(params, physics, assets)
    response_matrix = construct_response_matrix(params, assets)
    flux = physics.sns_flux.flux(exposure = assets.exposure, distance = assets.distance, params = params)

    dNdEr_all = nothing
    for iso in assets.isotopes
        nupar = [iso.fraction, iso.mass, iso.Z, iso.N]
        rate_matrix = build_rate_matrix(
            assets.er_centers .* 1e-3,  # MeV
            flux.E,                     # MeV
            nupar,
            physics,
            params,
            iso.Rn_key,
        )
        dNdEr = vec(sum(rate_matrix .* permutedims(flux.total_flux), dims = 2))
        if dNdEr_all === nothing
            dNdEr_all = zero.(dNdEr)  # picks Dual when AD is active
        end
        dNdEr_all .+= iso.fraction .* dNdEr
    end

    dEr = diff(assets.er_edges .* 1e-3)  # MeV
    int_rate = params.coherent_lar_mass .* assets.Nt .* dNdEr_all .* dEr
    predicted_counts = response_matrix * int_rate
    return predicted_counts
end

function get_backgrounds(params, assets)
    scale_template(template, norm) =
        sum(template) > 0 ?
            norm .* (template ./ sum(template)) :
            fill(zero(norm), length(template))
    pbrn = scale_template(assets.pbrn_binned, params.pbrn_norm)
    delbrn = scale_template(assets.delbrn_binned, params.delbrn_norm)
    ss_bkg = scale_template(assets.ss_bkg_binned, params.ss_bkg_norm)
    return (pbrn, delbrn, ss_bkg)
end

function get_forward_model(physics, assets)
    function forward_model(params)
        signal = get_expected(params, physics, assets)
        bkg_pbrn, bkg_delbrn, bkg_ss_bkg = get_backgrounds(params, assets)
        total_bkg = bkg_pbrn .+ bkg_delbrn .+ bkg_ss_bkg
        exp_events = signal .+ total_bkg
        distprod(Poisson.(exp_events))
    end
end

function get_plot(physics, assets)
    function plot(params, data=assets.observed)
        nothing
    end
end


end
