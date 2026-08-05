module earth_layers

using CSV, DataFrames
using StatsBase
using StaticArrays, ArraysOfArrays, StructArrays
using DataStructures

using ..Newtrinos
export configure
export PREM
export PREM12

const datadir = @__DIR__

"""
    DensityModel

Abstract type for Earth density profile models.

Each subtype provides a recipe for dividing the Earth into concentric shells of constant
density. Currently the only implementation is [`PREM`](@ref).
"""
abstract type DensityModel end

"""
    PREM <: DensityModel

Preliminary Reference Earth Model (Dziewonski & Anderson, 1981).

Reads the tabulated PREM density profile from `PREM_1s.csv` and groups radial shells into
zones defined by density boundaries. The proton-to-nucleon fraction is assumed constant
across all layers.

# Fields
- `zones::Array{Float64} = [0., 4., 7.5, 12.5, 13.1]`: density boundaries [g/cm³]
  defining the constant-density zones. Adjacent PREM rows whose density falls within the
  same bin are averaged into a single [`Newtrinos.osc.Layer`](@ref).
- `p_fractions::Float64 = 0.5`: proton number fraction ``Y_p = N_p / (N_p + N_n)``.
- `atm_heihgt::Float64 = 20.`: atmospheric shell thickness [km] added above the Earth's
  surface (density = 0).
"""
@kwdef struct PREM <: DensityModel
    zones::Array{Float64} = [0., 4., 7.5, 12.5, 13.1]
    p_fractions::Float64 = 0.5
    atm_heihgt::Float64 = 20.
end

"""
    PREM12 <: DensityModel

12-layer PREM approximation with radius-defined shell boundaries and separate core/mantle
electron fractions, as used in IceCube's NSI analysis (Dziewonski & Anderson, 1981; PRD 104,
072006 (2021), Sec. II.B).

Unlike [`PREM`](@ref) (which bins PREM rows into density-threshold zones and applies a single
proton/electron fraction to the whole Earth), this model bins by fixed radius boundaries and
distinguishes the core's electron fraction from the mantle's, matching the paper's stated
``Y_e^c = 0.466`` (inner + outer core) and ``Y_e^m = 0.496`` (mantle).

# Fields
- `r_boundaries::Array{Float64}`: 13 radius boundaries [km] defining 12 concentric shells,
  from the Earth's center (0) to its surface (6371).
- `r_cmb::Float64 = 3480.`: core-mantle boundary radius [km]; shells with outer radius
  `<= r_cmb` use `Y_e_core`, all others use `Y_e_mantle`.
- `Y_e_core::Float64 = 0.466`: relative electron-to-nucleon number density in the core.
- `Y_e_mantle::Float64 = 0.496`: relative electron-to-nucleon number density in the mantle.
- `atm_heihgt::Float64 = 20.`: atmospheric shell thickness [km] added above the Earth's
  surface (density = 0).
"""
@kwdef struct PREM12 <: DensityModel
    r_boundaries::Array{Float64} = [0.0, 1221.5, 2000.0, 3000.0, 3480.0,
                                     3700.0, 4200.0, 4700.0, 5200.0, 5700.0,
                                     5970.0, 6250.0, 6371.0]
    r_cmb::Float64 = 3480.0
    Y_e_core::Float64 = 0.466
    Y_e_mantle::Float64 = 0.496
    atm_heihgt::Float64 = 20.
end

"""
    EarthLayers <: Newtrinos.Physics

Configured Earth density model, returned by [`configure`](@ref).

Provides the layer structure and path-computation functions needed by the oscillation
module's matter-effect calculations (see [`Newtrinos.osc.SI`](@ref)).

# Fields
- `cfg::DensityModel`: the density model used to build this module.
- `params::NamedTuple`: oscillation parameters.
- `priors::NamedTuple`: prior distributions.
- `compute_layers::Function`: closure `compute_layers() -> StructVector{Layer}` returning
  the concentric density shells.
- `compute_paths::Function`: function
  `compute_paths(cz, layers; r_detector) -> VectorOfVectors{Path}` computing the layer
  traversal for each cosine-zenith value.
"""
@kwdef struct EarthLayers <: Newtrinos.Physics
    cfg::DensityModel
    params::NamedTuple
    priors::NamedTuple
    compute_layers::Function
    compute_paths::Function
end

"""
    configure(cfg::DensityModel=PREM()) -> EarthLayers

Create a fully configured Earth density physics module.

# Arguments
- `cfg::DensityModel`: density model (defaults to [`PREM`](@ref)).

# Returns
An [`EarthLayers`](@ref) instance with `compute_layers` and `compute_paths` closures and empty parameters and priors.

# Examples
```julia
using Newtrinos

earth = Newtrinos.earth_layers.configure()
layers = earth.compute_layers()
paths = earth.compute_paths([-1.0, -0.5, 0.0], layers)
```
"""
function configure(cfg::DensityModel=PREM())
    EarthLayers(
        cfg=cfg,
        params = (;),
        priors = (;),
        compute_layers = get_compute_layers(cfg),
        compute_paths = compute_paths
        )
end


"""
    get_compute_layers(cfg::PREM) -> Function

Construct a closure that builds the Earth layer structure from the PREM density profile.

The returned function `compute_layers()` reads `PREM_1s.csv`, bins rows by the density
boundaries in `cfg.zones`, and returns a `StructVector{Layer}` with one entry per zone
plus an atmospheric shell.

# Arguments
- `cfg::PREM`: PREM configuration with zone boundaries, proton fraction, and atmosphere height.

# Returns
A zero-argument closure `compute_layers() -> StructVector{Layer}`.
"""
function get_compute_layers(cfg::PREM)
    function compute_layers()
        
        PREM = CSV.read(joinpath(datadir, "PREM_1s.csv"), DataFrame, header=["radius","depth","density","Vpv","Vph","Vsv","Vsh","eta","Q-mu","Q-kappa"])
        # density boundaries to define the constant density zones
        
        radii = Float64[]
        ave_densities = Float64[]
        
        push!(radii, 6371+cfg.atm_heihgt)
        push!(ave_densities, 0.)
        
        for i in 1:length(cfg.zones)-1
            mask = (PREM.density .< cfg.zones[i+1]) .& (PREM.density .>= cfg.zones[i])
            push!(radii, maximum(PREM.radius[mask]))
            push!(ave_densities, mean(PREM.density[mask]))
        end
        
        layers = StructArray{Newtrinos.Layer}((radii, ave_densities .* cfg.p_fractions, ave_densities .* (1 .- cfg.p_fractions)))
    end
end

"""
    get_compute_layers(cfg::PREM12) -> Function

Construct a closure that builds the Earth layer structure from the PREM density profile,
using [`PREM12`](@ref)'s fixed radius boundaries and core/mantle-dependent electron fraction.

The returned function `compute_layers()` reads `PREM_1s.csv`, bins rows by the radius
boundaries in `cfg.r_boundaries` (iterating outer shell to inner, matching the radius-descending
convention [`compute_paths`](@ref) expects), and returns a `StructVector{Layer}` with one entry
per shell plus an atmospheric shell.

# Arguments
- `cfg::PREM12`: PREM12 configuration with radius boundaries, core-mantle boundary, electron
  fractions, and atmosphere height.

# Returns
A zero-argument closure `compute_layers() -> StructVector{Layer}`.
"""
function get_compute_layers(cfg::PREM12)
    function compute_layers()
        PREM = CSV.read(joinpath(datadir, "PREM_1s.csv"), DataFrame, header=["radius","depth","density","Vpv","Vph","Vsv","Vsh","eta","Q-mu","Q-kappa"])

        radii = Float64[]
        p_densities = Float64[]
        n_densities = Float64[]

        push!(radii, 6371+cfg.atm_heihgt)
        push!(p_densities, 0.)
        push!(n_densities, 0.)

        for i in length(cfg.r_boundaries)-1:-1:1
            r_inner, r_outer = cfg.r_boundaries[i], cfg.r_boundaries[i+1]
            mask = (PREM.radius .>= r_inner) .& (PREM.radius .< r_outer)
            !any(mask) && continue
            ρ = mean(PREM.density[mask])
            Y_e = r_outer <= cfg.r_cmb ? cfg.Y_e_core : cfg.Y_e_mantle
            push!(radii, maximum(PREM.radius[mask]))
            push!(p_densities, ρ * Y_e)
            push!(n_densities, ρ * (1 - Y_e))
        end

        layers = StructArray{Newtrinos.Layer}((radii, p_densities, n_densities))
    end
end

"""
    ray_circle_path_length(r, y, cz) -> Real

Compute the chord length of a ray through a circle of radius `r`.

The ray originates at radial distance `y` from the Earth's centre with direction
cosine `cz` (cosine of the zenith angle). Returns zero if the ray does not intersect
the circle or if the chord is shorter than 1 km (numerical noise filter).

# Arguments
- `r`: radius of the spherical shell [km].
- `y`: radial position of the detector [km].
- `cz`: cosine of the zenith angle (−1 = vertically upgoing through the core).

# Returns
Path length through the shell [km], or zero if no intersection.
"""
function ray_circle_path_length(r, y, cz)
    # Compute the discriminant
    disc = r^2 - y^2 + (y * cz)^2
    T = typeof(disc)

    if disc < 0
        return zero(T)  # No intersection
    end

    sqrt_disc = sqrt(disc)

    # Compute intersection points
    t1 = - y * cz - sqrt_disc
    t2 = - y * cz + sqrt_disc

    # Compute path length, ensuring we only count positive t-values
    L = max(zero(T), t2 - max(zero(T), t1))

    if L < 1
        return zero(T)
    end
    L
end

# ToDo: could probably skip layers smaller than few km and "absorb" those into the next outer layer

"""
    compute_paths(cz::Number, layers, r_detector) -> StructArray{Path}
    compute_paths(cz::AbstractArray, layers; r_detector=6369) -> VectorOfVectors{Path}

Compute the sequence of [`Newtrinos.osc.Path`](@ref) segments a neutrino traverses through the Earth.

For a single cosine-zenith value, determines which [`Newtrinos.osc.Layer`](@ref) shells are crossed
using [`ray_circle_path_length`](@ref), then builds an ordered list of segments. Layers
below the detector are traversed twice (entry and exit), while layers above the detector
are traversed once.

The array method broadcasts over multiple cosine-zenith values and returns a
`VectorOfVectors{Path}`.

# Arguments
- `cz`: cosine of the zenith angle (scalar or array). −1 is vertically upgoing.
- `layers::StructVector{Layer}`: Earth density layers from `compute_layers()`.
- `r_detector::Real`: radial position of the detector [km] (default 6369, approximate
  IceCube depth).

# Returns
- Scalar method: `StructArray{Path}` with `length` and `layer_idx` columns.
- Array method: `VectorOfVectors{Path}`, one `Vector{Path}` per zenith angle.
"""
function compute_paths(cz::Number, layers, r_detector)
    radii = layers.radius
    intersections = ray_circle_path_length.(radii, r_detector, cz)
    for i in 1:length(intersections) - 1
        intersections[i] -= intersections[i+1]
    end
    mask = intersections .> 0.
    rs = radii[mask]
    intersections = intersections[mask]

    n_layers_outside = sum(radii .>= r_detector)

    n_layers = 2 * (length(intersections) - n_layers_outside) + n_layers_outside

    lengths_traversed = zeros(n_layers)
    layer_idx_traversed = zeros(Int, n_layers)

    for i in 1:length(intersections)
        if (i < n_layers_outside) | (i == length(intersections))
            lengths_traversed[i] = intersections[i]
            layer_idx_traversed[i] = i
        elseif i == n_layers_outside
            #len_det = -cz * (rs[i] - r_detector) # straight line approximation  
            len_det = r_detector * cz + sqrt(rs[i]^2 - r_detector^2 + (r_detector*cz)^2) # improved more exact result 
            inter = intersections[i] - len_det
            lengths_traversed[i] = inter/2 + len_det
            layer_idx_traversed[i] = i
            lengths_traversed[end-i+n_layers_outside] = inter/2
            layer_idx_traversed[end-i+n_layers_outside] = i              
        else
            lengths_traversed[i] = intersections[i]/2
            layer_idx_traversed[i] = i
            lengths_traversed[end-i+n_layers_outside] = intersections[i]/2
            layer_idx_traversed[end-i+n_layers_outside] = i
        end
    end

    la = StructArray{Newtrinos.Path}((lengths_traversed, layer_idx_traversed))
    
end

function compute_paths(cz::AbstractArray, layers; r_detector = 6369)
    VectorOfVectors{Newtrinos.Path}(compute_paths.(cz, Ref(layers), r_detector));
end

end