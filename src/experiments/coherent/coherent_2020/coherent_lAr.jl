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
        params = get_params(),
        priors = get_priors(),
        assets = assets,
        forward_model = get_forward_model(physics, assets),
        plot = get_plot(physics, assets)
    )
end

function get_params()
    params = (
        coherent_lar_qfa_a = 0.246,  # QF polynomial coefficients
        coherent_lar_qfa_b = 0.00078,
        coherent_lar_mass = 24.4,  # kg
        brn_norm= 497.0,  # Normalization factor for BRN
        delbrn_norm= 33.0,  # Normalization factor for delBRN
        ss_bkg_norm= 3154.0,  # Normalization factor for SS background
        )
end

# TODO!
function get_priors()
    priors = (
        coherent_lar_qfa_a = Normal(0.246, 0.006),
        coherent_lar_qfa_b = Normal(0.00078, 0.00009),
        coherent_lar_mass = truncated(Normal(24.4, 0.61), 0.0, Inf),
        brn_norm= truncated(Normal(497.0, 160.0), 0.0, Inf),  # Normalization factor for BRN
        delbrn_norm= truncated(Normal(33.0, 33.0), 0.0, Inf),  # Normalization factor for delBRN
        ss_bkg_norm= truncated(Normal(3154.0, 25.0), 0.0, Inf),  # Normalization factor for SS background
        )
end

function get_assets(physics, datadir = @__DIR__)
    @info "Loading coherent lAr data"


    er_edges = collect(3:0.5:100) # keVnr
    isotopes = [
        (fraction=1.0, mass=37.3e3, Z=18, N=22, Rn_key=:Rn_Ar, Rn_nom=4.1039) # Ar-37
    ] # List of isotopes with [fraction, Nuclear mass (GeV), Z, N=A-Z, Rn_key]
    Nt = (1.0/0.039948) * 6.023e+23
    
    resolution = 0.58  # a/Eee and b*Eee
    # Import Data
    eff_data = CSV.read(joinpath(datadir, "lAr/CENNS10AnlAEfficiency.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: keVee, keVnr, efficiency
    brn = CSV.read(joinpath(datadir, "lAr/brnpdf.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: keVee, F90, t, counts/bin/6.12GWh
    delbrn = CSV.read(joinpath(datadir, "lAr/delbrnpdf.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: keVee, F90, t, counts/bin/6.12GWh
    ss_bkg = CSV.read(joinpath(datadir, "lAr/bkgpdf.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: keVee, F90, t, counts/bin/6.12GWh
    observed_df = CSV.read(joinpath(datadir, "lAr/datanobkgsub.txt"), DataFrame, comment="#", header=false, delim=' ') # columns: keVee, F90, t, counts/bin/6.12GWh
    
    # Define bin centers as the average of consecutive edges
    er_centers = collect((er_edges[1:end-1] + er_edges[2:end]) ./ 2) # Reconstruct bin centers from edges

    # Match data: output bin centers 5, 15, 25, ..., 55
    out_width = 10.0
    out_edges = collect(0.0:out_width:60.0)         # Edges: [0, 10, 20, ..., 60]
    out_centers = midpoints(out_edges)              # Centers: [5, 15, 25, ..., 55]
    
    # Bin observed data into out_centers using Helpers.rebin
    observed = Helpers.rebin(observed_df, out_edges; var_col=1, count_col=4)

    # Bin background PDFs ONCE
    brn_binned = Helpers.rebin(brn, out_edges; var_col=1, count_col=4)
    delbrn_binned = Helpers.rebin(delbrn, out_edges; var_col=1, count_col=4)
    ss_bkg_binned = Helpers.rebin(ss_bkg, out_edges; var_col=1, count_col=4)

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
        brn_binned,
        delbrn_binned,
        ss_bkg_binned,
        distance,
        exposure,
    )
end

function qf(er_centers, params)
    a = params.coherent_lar_qfa_a
    b = params.coherent_lar_qfa_b
    return @. (a + b * er_centers) * er_centers  # Convert to keVee
end

function eff(keVee, assets)
    x = assets.eff_data[:, 1]
    y = assets.eff_data[:, 3]
    itp = LinearInterpolation(x, y, extrapolation_bc=Flat())
    vals = similar(keVee)
    for i in eachindex(keVee)
        if keVee[i] < minimum(x)
            vals[i] = 0.0
        elseif keVee[i] > maximum(x)
            vals[i] = y[end]
        else
            vals[i] = itp(keVee[i])
        end
    end
    return vals
end

function sigma_keVee(keVee, assets)
    a = assets.resolution
    return @. a * sqrt(keVee)  # keVee resolution  
end

function construct_response_matrix(params, assets)
    er_centers = assets.er_centers
    out_centers = assets.out_centers
    n_out = length(out_centers)
    n_er = length(er_centers)

    # Compute keVee and sigma_keVee for all er_centers
    keVee = qf(er_centers, params)
    sigma = sigma_keVee(keVee, assets)
    eff_vals = eff(keVee, assets)

    # Build the response matrix: each column is the PDF for one er_center
    A = zeros(n_out, n_er)
    for j in 1:n_er
        # Debug: print type and value if error is likely
        try
            if keVee[j] <= 0 || sigma[j] <= 0 || eff_vals[j] <= 0
                continue
            end
        catch err
            println("DEBUG: j=", j, " keVee[j]=", keVee[j], " type=", typeof(keVee[j]))
            println("DEBUG: sigma[j]=", sigma[j], " type=", typeof(sigma[j]))
            println("DEBUG: eff_vals[j]=", eff_vals[j], " type=", typeof(eff_vals[j]))
            rethrow(err)
        end
        pdf_vals = pdf.(Normal(keVee[j], sigma[j]), out_centers)
        if sum(pdf_vals) > 0
            pdf_vals ./= sum(pdf_vals)
        end
        A[:, j] .= eff_vals[j] .* pdf_vals
    end
    A[A .< 0] .= 0
    return A
end

function build_rate_matrix(er_centers, enu_centers, nupar, physics, params, Rn_key)
    physics.cevns_xsec.diff_xsec_lar(er_centers, enu_centers, params, nupar, Rn_key)
end


function get_expected(params, physics, assets)
    response_matrix = construct_response_matrix(params, assets)
    flux = physics.sns_flux.flux(exposure=assets.exposure, distance=assets.distance, params=params)
    dNdEr_all = zeros(eltype(assets.er_centers), size(assets.er_centers))
    for iso in assets.isotopes
        nupar = [iso.fraction, iso.mass, iso.Z, iso.N]
        rate_matrix = build_rate_matrix(
            assets.er_centers * 1e-3,  # convert to MeV
            flux.E,         # Neutrino energy centers (MeV)
            nupar,
            physics,
            params,
            iso.Rn_key,
        )
        dNdEr = sum(rate_matrix .* flux.total_flux', dims=2)
        dNdEr = dropdims(dNdEr, dims=2)
        dNdEr_all .+= iso.fraction * dNdEr
    end
    int_rate = params.coherent_lar_mass * assets.Nt * dNdEr_all .* diff(assets.er_edges * 1e-3)
    predicted_counts = response_matrix * int_rate
    return predicted_counts
end

function get_backgrounds(params, assets)
    scale_template(template, norm) = sum(template) > 0 ? norm * (template / sum(template)) : zeros(length(template))
    brn = scale_template(assets.brn_binned, params.brn_norm)
    delbrn = scale_template(assets.delbrn_binned, params.delbrn_norm)
    ss_bkg = scale_template(assets.ss_bkg_binned, params.ss_bkg_norm)
    return (brn, delbrn, ss_bkg)
end

function get_forward_model(physics, assets)
    function forward_model(params)
        # Signal prediction
        signal = get_expected(params, physics, assets)
        # Backgrounds
        bkg_brn, bkg_delbrn, bkg_ss_bkg = get_backgrounds(params, assets)
        total_bkg = bkg_brn .+ bkg_delbrn .+ bkg_ss_bkg
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
