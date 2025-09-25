module coherent_csi

using DataFrames
using CSV
using Distributions
using SpecialFunctions
using LinearAlgebra
using Statistics
using DataStructures
using BAT
using CairoMakie
using Logging
using StatsBase
using ..Helpers
import ..Newtrinos

@kwdef struct COHERENT_CSI <: Newtrinos.Experiment
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
    return COHERENT_CSI(
        physics = physics,
        params = get_params(assets.ss_bkg_nom, assets.brn_nom, assets.nin_nom),
        priors = get_priors(assets.ss_bkg_nom, assets.brn_nom, assets.nin_nom),
        assets = assets,
        forward_model = get_forward_model(physics, assets),
        plot = get_plot(physics, assets)
    )
end

function get_params(ss_bkg_nom, brn_nom, nin_nom)
    params = (
        coherent_csi_eff_a = 1.32045,
        coherent_csi_eff_b = 0.285979,
        coherent_csi_eff_c = 10.8646,
        coherent_csi_eff_d = -0.333322,
        coherent_csi_qfa_a = 0.0554628,  # QF polynomial coefficients
        coherent_csi_qfa_b = 4.30681,
        coherent_csi_qfa_c = -111.707,
        coherent_csi_qfa_d = 840.384,
        brn_norm= brn_nom,  # Normalization factor for BRN
        nin_norm= nin_nom,  # Normalization factor for NIN
        ss_bkg_norm= ss_bkg_nom,  # Normalization factor for SS background
        )
end

# TODO!
function get_priors(ss_bkg_nom, brn_nom, nin_nom)
    priors = (
        coherent_csi_eff_a = Normal(1.32045, 0.02),
        coherent_csi_eff_b = Normal(0.285979, 0.0006),
        coherent_csi_eff_c = Normal(10.8646, 1.),
        coherent_csi_eff_d = Normal(-0.333322, 0.03),
        coherent_csi_qfa_a = Normal(0.0554628, 0.0059),
        coherent_csi_qfa_b = Normal(4.30681, 0.79),
        coherent_csi_qfa_c = Normal(-111.707, 26.15),
        coherent_csi_qfa_d = Normal(840.384, 244.82),
        brn_norm= truncated(Normal(brn_nom, 0.25 * brn_nom), 0.0, brn_nom + 3 * 0.25 * brn_nom),  # Normalization factor for BRN
        nin_norm= truncated(Normal(nin_nom, 0.36 * nin_nom), 0.0, nin_nom + 3 * 0.36 * nin_nom),  # Normalization factor for NIN
        ss_bkg_norm= truncated(Normal(ss_bkg_nom, 0.021 * ss_bkg_nom), 0.0, ss_bkg_nom + 3 * 0.021 * ss_bkg_nom),  # Normalization factor for SS background
        )
end

function get_assets(physics, datadir = @__DIR__)
    @info "Loading coherent csi data"


    er_edges = LinRange(3, 200, Int((200-3)/0.5)) # keV
    isotopes = [
        (fraction=0.49, mass=123.8e3, Z=55, N=78, Rn_key=:Rn_Cs, Rn_nom=5.7242),  # Cs-133
        (fraction=0.51, mass=118.21e3, Z=53, N=74, Rn_key=:Rn_I, Rn_nom=5.7242) # I-127
    ] # List of isotopes with [fraction, Nuclear mass (GeV), Z, N=A-Z, Rn_key]
    Nt = 2 * (14.6/0.25981) * 6.023e+23
    light_yield = 13.35  # PE/keVee
    resolution = [0.0749, 9.56]  # a/Eee and b*Eee
    # Import Data
    brnPE_df = CSV.read(joinpath(datadir, "csi/brnPE.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: PE, BRN PDF
    ninPE_df = CSV.read(joinpath(datadir, "csi/ninPE.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: PE, NIN PDF
    ssBkg_df = CSV.read(joinpath(datadir, "csi/dataBeamOnAC.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: PE, timestamp
    observed_df = CSV.read(joinpath(datadir, "csi/dataBeamOnC.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: PE, timestamp

    # Reconstruct bin edges from centers
    er_centers = midpoints(er_edges)

    pe_width = 5.0
    out_edges = collect(2.5:pe_width:202.5)  # PE bin edges: [2.5, 7.5, 12.5, ..., 202.5]
    out_centers = midpoints(out_edges) # Bin centers: [5, 10, 15, ..., 200]

    # For event lists (PE only, e.g. ssBkg, observed), use Helpers.bin
    ssBkg = Helpers.bin(ssBkg_df, out_edges; col=1)
    observed = Helpers.bin(observed_df, out_edges; col=1)

    # For PDFs (PE, count), use Helpers.rebin
    brnPE = Helpers.rebin(brnPE_df, out_edges; var_col=1, count_col=2)
    ninPE = Helpers.rebin(ninPE_df, out_edges; var_col=1, count_col=2)

    # Get initial nominal value for Bkg normalizations
    ss_bkg_nom = sum(ssBkg)
    @info "Initial SS background normalization: $ss_bkg_nom"
    brn_nom = sum(brnPE)
    @info "Initial BRN background normalization: $brn_nom"
    nin_nom = sum(ninPE)
    @info "Initial NIN background normalization: $nin_nom"

    distance = 19.3 # m
    exposure = 13.99 # GWh


    assets = (;
        observed,
        er_edges,
        er_centers,
        out_edges,
        out_centers,
        isotopes,
        Nt,
        light_yield,
        resolution,
        brnPE,
        brn_nom,
        ninPE,
        nin_nom,
        ssBkg,
        ss_bkg_nom,
        distance,
        exposure,
    )
end

# QF: accepts scalar or array, returns same shape, type-generic
function qf(er, params)
    a = params.coherent_csi_qfa_a
    b = params.coherent_csi_qfa_b
    c = params.coherent_csi_qfa_c
    d = params.coherent_csi_qfa_d
    x = er .* 1e-3                         # MeV
    vals = (a .* x .+ b .* x.^2 .+ c .* x.^3 .+ d .* x.^4) .* 1e3  # keVee
    z = vals isa AbstractArray ? zero(eltype(vals)) : zero(vals)
    return max.(vals, z)
end

# Efficiency: accepts scalar or array of PE, type-generic
function eff(pe, params)
    a = params.coherent_csi_eff_a
    b = params.coherent_csi_eff_b
    c = params.coherent_csi_eff_c
    d = params.coherent_csi_eff_d
    vals = @. a / (1 + exp(-b * (pe - c))) + d
    z = vals isa AbstractArray ? zero(eltype(vals)) : zero(vals)
    return max.(vals, z)
end

# Gamma PDF (k-θ parameterization), AD-friendly
@generated function _typed_one(::Type{T}) where {T}
    :(one(T))
end
@inline function _gamma_pdf(x, k, θ)
    # x, k, θ should be of a common Real-like type (possibly Dual)
    # Use log-pdf for numerical stability and AD-compatibility
    # Assumes x > 0 (our PE bin edges are > 0)
    return exp((k - 1) * log(x) - x / θ - k * log(θ) - loggamma(k))
end

"""
Compute gamma-smearing probabilities per PE bin (Simpson integration), AD-safe.

Arguments:
- Eee: recoil energy in keVee (can be Dual)
- pe_centers: PE bin centers (vector)
- pe_edges: PE bin edges (vector, length = length(pe_centers)+1)
- resolution: [a, b] with width parameters
- light_yield: keVee -> PE

Returns:
- probs: vector of probabilities per PE bin (same length as pe_centers)
"""
function gamma_pdf_integrated_over_bins(Eee, pe_centers, pe_edges, resolution, light_yield)
    n_bins = length(pe_centers)
    T = promote_type(eltype(pe_centers), typeof(Eee))
    probs = zeros(T, n_bins)

    # If Eee == 0, there is no signal; avoid divisions by Eee
    if iszero(Eee)
        return probs
    end

    a = resolution[1] / Eee
    b = resolution[2] * Eee
    k = one(T) + b
    θ = inv(a * (one(T) + b))

    # Integrate Gamma(k, θ) over each [lo, hi] via Simpson's rule
    for i in eachindex(pe_centers)
        lo = T(pe_edges[i])
        hi = T(pe_edges[i + 1])
        mid = (lo + hi) / 2
        probs[i] = (hi - lo) / 6 * (_gamma_pdf(lo, k, θ) + 4 * _gamma_pdf(mid, k, θ) + _gamma_pdf(hi, k, θ))
    end

    s = sum(probs)
    if !iszero(s)
        probs ./= s
    end
    return probs
end

# Single ER-bin response column, AD-safe
function response_matrix_per_er_bin(keVnr, params, assets)
    keVee = qf(keVnr, params)                    # scalar (possibly Dual)
    weights = gamma_pdf_integrated_over_bins(keVee, assets.out_centers, assets.out_edges,
                                             assets.resolution, assets.light_yield)
    s = sum(weights)
    if iszero(s)
        return weights  # already zeros of correct type
    end
    eff_vals = eff(assets.out_centers, params)
    return weights .* eff_vals
end

# Full response matrix, AD-safe
function construct_response_matrix(params, assets)
    n_out = length(assets.out_centers)
    n_er = length(assets.er_centers)

    first_col = response_matrix_per_er_bin(first(assets.er_centers), params, assets)
    Tcol = eltype(first_col)
    A = Array{Tcol}(undef, n_out, n_er)
    A[:, 1] = max.(first_col, zero(eltype(first_col)))

    for j in 2:n_er
        col = response_matrix_per_er_bin(assets.er_centers[j], params, assets)
        A[:, j] = length(col) == n_out ? max.(col, zero(eltype(col))) : fill(zero(Tcol), n_out)
    end
    return A
end

function build_rate_matrix(er_centers, enu_centers, nupar, physics, params, Rn_key)
    physics.cevns_xsec.diff_xsec_csi(er_centers, enu_centers, params, nupar, Rn_key)
end


function get_expected(params, physics, assets)
    response_matrix = construct_response_matrix(params, assets)

    flux = physics.sns_flux.flux(exposure = assets.exposure, distance = assets.distance, params = params)

    function dNdEr(iso)
        nupar = [iso.fraction, iso.mass, iso.Z, iso.N]
        rate_matrix = build_rate_matrix(
            assets.er_centers .* 1e-3,  # MeV
            flux.E,                     # MeV
            nupar,
            physics,
            params,
            iso.Rn_key,
        )
        return iso.fraction .* vec(sum(rate_matrix .* permutedims(flux.total_flux), dims = 2))
    end
    
    dNdEr_all = sum([dNdEr(iso) for iso in assets.isotopes])
    dEr = diff(assets.er_edges .* 1e-3)           # MeV
    int_rate = assets.Nt .* dNdEr_all .* dEr
    predicted_counts = response_matrix * int_rate
    return predicted_counts
end

function get_backgrounds(params, assets)
    scale_template(template, norm) = sum(template) > 0 ?
        norm .* (template ./ sum(template)) :
        fill(zero(norm), length(template))
    brnPE = scale_template(assets.brnPE, params.brn_norm)
    ninPE = scale_template(assets.ninPE, params.nin_norm)
    ssBkg = scale_template(assets.ssBkg, params.ss_bkg_norm)
    return (brnPE, ninPE, ssBkg)
end

function get_forward_model(physics, assets)
    function forward_model(params)
        signal = get_expected(params, physics, assets)
        brnPE, ninPE, ssBkg = get_backgrounds(params, assets)
        total_bkg = brnPE .+ ninPE .+ ssBkg
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
