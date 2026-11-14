### ANALYSIS SETUP FILE FOR NSI ANALYSIS OF THE DEEPCORE 8-YEAR DATA SET ###

using Distributions
using Distributions
using DataStructures
using Accessors
using Newtrinos

function setup(; ftype = Float64)

    ### PHYSICS AND EXPERIMENTS ###

    osc = Newtrinos.osc.configure(
        Newtrinos.osc.OscillationConfig(
            flavour     = Newtrinos.osc.ThreeFlavour(),
            interaction = Newtrinos.osc.NSI_Standard()
        )
    )
    atm_flux     = Newtrinos.atm_flux.configure()
    earth_layers = Newtrinos.earth_layers.configure(Newtrinos.earth_layers.PREM12())
    xsec         = Newtrinos.xsec.configure()
    physics      = (; osc, atm_flux, earth_layers, xsec)
    experiments  = (; deepcore = Newtrinos.deepcore_8y.configure(physics))

    p = Newtrinos.get_params(experiments)
    priors = Newtrinos.get_priors(experiments)

    # CP-conserving approximation, fix standard oscillation params that are not relevant to the NSI analysis 
    conditional_vars = Dict(:δCP=>0.0, :θ₁₂=>ftype(33.62 * π/180), :θ₁₃=>ftype(8.54  * π/180), :Δm²₂₁=>ftype(7.40e-5))
    priors = Newtrinos.condition(priors, conditional_vars, p)

    # fix priors of NSI and theta_23
    @reset priors.θ₂₃   = Uniform(ftype(30*π/180), ftype(60*π/180))
    @reset priors.Δ_eμ     = Uniform(-ftype(5),   ftype(5))
    @reset priors.Δ_τμ     = Uniform(-ftype(0.1), ftype(0.1))
    @reset priors.ε_eμ_abs = Uniform(ftype(0),    ftype(0.3))
    @reset priors.ε_eτ_abs = Uniform(ftype(0),    ftype(0.35))
    @reset priors.ε_μτ_abs = Uniform(ftype(0),    ftype(0.07))
    @reset priors.δ_eμ     = Uniform(ftype(0),    ftype(2π))  
    @reset priors.δ_eτ     = Uniform(ftype(0),    ftype(2π))
    @reset priors.δ_μτ     = Uniform(ftype(0),    ftype(2π))

    # hypotheses: every hypothesis frees the listed NSI parameters and pins all other NSI parameters to 0.

    nsi_names = (:Δ_eμ, :Δ_τμ, :ε_eμ_abs, :ε_eτ_abs, :ε_μτ_abs, :δ_eμ, :δ_eτ, :δ_μτ)

    function nsi_hypothesis(tested)
        pinned = Dict{Symbol, Float64}(k => 0.0 for k in nsi_names if !(k in tested))
        prior = isempty(pinned) ? deepcopy(priors) : Newtrinos.condition(priors, pinned, p)
        (; prior, null = NamedTuple{Tuple(tested)}(ntuple(_ -> 0.0, length(tested))))
    end

    hypotheses = OrderedDict(   # ordered: the scripts list/run hypotheses in this order
        "eμ"    => nsi_hypothesis((:ε_eμ_abs, :δ_eμ)),
        "eτ"    => nsi_hypothesis((:ε_eτ_abs, :δ_eτ)),
        "μτ"    => nsi_hypothesis((:ε_μτ_abs, :δ_μτ)),
        "ee_μμ" => nsi_hypothesis((:Δ_eμ,)),
        "ττ_μμ" => nsi_hypothesis((:Δ_τμ,)),
        "all"   => nsi_hypothesis(nsi_names), # all NSI parameters free
    )

    #return the setup as a named tuple of name, experiments, params, priors, hypotheses
    (; name = "deepcore_8y", experiments, params = p, priors, hypotheses)
end
