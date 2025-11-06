module cevns_xsec

using LinearAlgebra
using Distributions
using SpecialFunctions
using ..Newtrinos

@kwdef struct CevnsXsec <: Newtrinos.Physics
    params::NamedTuple
    priors::NamedTuple
    diff_xsec::Function
end

# Only keep the dynamic configure (isotope_keys) and the (params, priors) version
function configure(isotopes, er_centers, enu_centers)
    # Build assets from isotopes
    assets = get_assets(isotopes, er_centers, enu_centers)
    params, priors = build_params_and_priors(isotopes)
    CevnsXsec(
        params = params,
        priors = priors,
        diff_xsec = get_diff_xsec(assets),
    )
end

function configure(params::NamedTuple, priors::NamedTuple)
    CevnsXsec(
        params = params,
        priors = priors,
        diff_xsec_lar = get_diff_xsec_lar(),
        diff_xsec_csi = get_diff_xsec_csi(),
    )
end

# Dynamic parameter/prior builder for isotope-specific Rn keys, using isotope list
function build_params_and_priors(isotopes)
    param_dict = Dict(
        :cevns_xsec_a => 0.0,
        :cevns_xsec_b => 0.0,
        :cevns_xsec_c => 0.0,
        :cevns_xsec_d => 0.0,
    )
    prior_dict = Dict{Symbol, Distributions.Distribution}(
        :cevns_xsec_a => Uniform(-2, 2),
        :cevns_xsec_b => Uniform(-2000.0, 2000.0),
        :cevns_xsec_c => Uniform(-3, 3),
        :cevns_xsec_d => Uniform(-1e6, 1e6),
    )
    for iso in isotopes
        param_dict[iso.Rn_key] = iso.Rn_nom
        prior_dict[iso.Rn_key] = truncated(Normal(iso.Rn_nom, 1), 0.0, iso.Rn_nom + 3 * 1)
    end
    return (NamedTuple(param_dict), NamedTuple(prior_dict))
end

const gf=1.1663787e-11
const me = 0.510998
const mmu = 105.6
const mtau=1776.86
const mpi = 139.57
const sw2 = 0.231
const gu = 1/2 - 2*2/3*sw2
const gd = -(1/2) + 2*1/3*sw2
const alph = 1/137
const ep= (mpi^2-mmu^2)/(2*mpi)

function get_assets(isotopes, er_centers, enu_centers)
    # Extract isotope data into a structured format
    @info "Configuring CEvNS cross-section assets"
    isotope_data = Dict(iso.Rn_key => (
        mass = iso.mass,
        Z = iso.Z,
        N = iso.N,
        fraction = iso.fraction,
        Rn_nom = iso.Rn_nom
    ) for iso in isotopes)

    # Return assets as a NamedTuple
    return (
        isotopes = isotope_data,  # Dictionary of isotope data keyed by Rn_key
        er_centers = er_centers,
        enu_centers = enu_centers
    )
end

# Helm-like nuclear form factor squared, generic in type
function ffsq(er, mn, rn)
    r0 = rn / 197.326963
    arg = 2 * mn * er
    q = sqrt(max(arg, zero(arg)))                 # typed zero
    j1 = sphericalbesselj(1, q * r0)
    denom = q * r0
    # Use short-circuiting branch to avoid evaluating the division when denom == 0
    ratio = iszero(denom) ? one(j1) : (3 * j1) / denom
    exp_factor = exp(-((q * (0.9 / 197.326963))^2) / 2)
    return (ratio * exp_factor)^2
end

# Vectorized differential cross section dσ/dEr (n_er × n_enu), AD-safe
function ds(er, enu, params, nupar, Rn_key)
    mN = nupar[1]
    Z  = nupar[2]
    N  = nupar[3]
    rn = params[Rn_key]

    # Per-Er prefactor (n_er,)
    C1d = (gf^2 / (4 * pi)) .* ffsq.(er, mN, rn)

    qwsq = (N - (1 - 4 * sw2) * Z)^2

    # SM deformation parameters as scalars; no tiny length-4 arrays needed
    c1 = qwsq * (params.cevns_xsec_a + one(params.cevns_xsec_a))
    c2 = qwsq * (params.cevns_xsec_b - one(params.cevns_xsec_b))
    c3 = qwsq * (params.cevns_xsec_c - one(params.cevns_xsec_c))
    c4 = qwsq * (params.cevns_xsec_d + one(params.cevns_xsec_d))

    # Grids
    er_grid  = reshape(er, :, 1)          # (n_er, 1)
    enu_grid = reshape(enu, 1, :)         # (1, n_enu)

    base2 = er_grid ./ enu_grid
    base3 = mN .* er_grid ./ (2 .* enu_grid.^2)
    base4 = (er_grid.^2) ./ (enu_grid.^2)

    xf = c1 .+ c2 .* base2 .+ c3 .* base3 .+ c4 .* base4
    zxf = zero(eltype(xf))
    heav = max.(xf, zxf)                  # type-stable clamp

    # Broadcast C1d over columns, keep element type generic
    res = (C1d .* mN) .* heav
    return res
end

function get_diff_xsec_lar()
    function diff_xsec_lar(er_centers, enu_centers, params, nupar, Rn_key)
        ds(er_centers, enu_centers, params, nupar, Rn_key)
    end
end

function get_diff_xsec_csi()
    function diff_xsec_csi(er_centers, enu_centers, params, nupar, Rn_key)
        ds(er_centers, enu_centers, params, nupar, Rn_key)
    end
end

function get_diff_xsec(assets)
    # Extract assets
    er_centers = assets.er_centers
    enu_centers = assets.enu_centers
    isotopes = assets.isotopes

    # Return a callable function that computes the differential cross-section
    return function (params)
        # Determine the element type of params (e.g., from the first element)
        param_type = eltype(params[:cevns_xsec_a])

        # Compute cross-section for each isotope and store in a dictionary
        xsec_dict = Dict{Symbol, Matrix{param_type}}()
        for (Rn_key, iso) in isotopes
            mass = iso.mass
            Z = iso.Z
            N = iso.N

            # Call ds() for the current isotope
            xsec_dict[Rn_key] = ds(er_centers, enu_centers, params, (mass, Z, N), Rn_key)
        end

        return xsec_dict  # Dictionary of cross-section matrices keyed by Rn_key
    end
end
end # module cevns_xsec