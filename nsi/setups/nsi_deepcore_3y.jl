# Analysis setup: NSI (Standard parametrization) with the DeepCore 3-year high-statistics sample B.
#
# A setup file describes the PROBLEM (physics, data, priors, hypotheses) and nothing about HOW it is solved:
# optimizer settings (PSO, LBFGS, seeds, ...) stay command-line arguments of the scripts. Scripts take
# `--setup <this file>`, `include` it and call `setup(; ftype)`. To analyse something else (another dataset,
# another physics model, e.g. the CP phase) write a new file with the same `setup` function; the scripts
# themselves do not change.
#
# `setup` returns a NamedTuple with
#   name          short label of the setup; used for the results folder and stored in every saved result
#   experiments   NamedTuple of configured experiments (all built from the same `physics`)
#   params        start values of all parameters (Newtrinos.get_params(experiments))
#   priors        priors of all parameters (before any hypothesis is applied)
#   hypotheses    OrderedDict  name => (; prior, null)   -- see below
#
# ---------------------------------------------------------------------------------------------------------
# HOW TO SPECIFY A HYPOTHESIS
#
# A hypothesis is one question of the form "is the effect X present, and can the fit recover it?".
# It consists of two things:
#
#   prior   the prior in which the parameters of X are FREE and every parameter that does not belong to the
#           analysis is fixed. Build it from `priors` with `Newtrinos.condition(priors, Dict(:name => value, ...), params)`:
#           every parameter listed there becomes a constant, every other parameter keeps its prior.
#   null    a NamedTuple with the TESTED parameters (the parameters of X) and their values when X is absent,
#           e.g. (ε_eμ_abs = 0.0, δ_eμ = 0.0) for "no eμ non-standard interaction".
#
# From these two the scripts derive everything else:
#   * tested parameters       = keys(null); they are the injected/scanned truth values
#   * null model              = condition(prior, null, params), i.e. the same prior with X switched off
#                               (used e.g. by the `smseed` stage of `--strategy`)
#   * scan type               = number of tested parameters: 1 -> 1D grid, 2 -> 2D grid, more -> Sobol sample
#   * nuisance/other parameters = all remaining FREE parameters of `prior`; they are fitted in every roundtrip
#                               and are the ones redrawn from their priors by `--fluctuate-nuisances`
# The order of the names in `null` is the order of the scan axes (for a 2D hypothesis: first name = x axis).
#
# Example for a different analysis, the CP phase (one tested parameter, scanned in 1D):
#     hypotheses = Dict("δCP" => (; prior = priors_with_δCP_free, null = (δCP = 0.0,)))
#
# ---------------------------------------------------------------------------------------------------------
# HOW TO FIX NUISANCE PARAMETERS OR CHANGE THEIR PRIORS
#
#   * fix a parameter (it then stays at that value in every fit and truth):
#         priors = Newtrinos.condition(priors, Dict(:atm_flux_updown_sigma => 0.0), params)
#   * change a prior (the parameter stays free):
#         @reset priors.deepcore_lifetime = Uniform(ftype(0), ftype(3.8))
#   Do this on `priors` BEFORE building the hypotheses: every hypothesis is derived from `priors`.
# ---------------------------------------------------------------------------------------------------------

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
    experiments  = (; deepcore = Newtrinos.deepcore_3y.configure(physics))

    ### PARAMETERS AND PRIORS ###

    p = Newtrinos.get_params(experiments)
    priors = Newtrinos.get_priors(experiments)

    # Standard oscillation parameters fixed to the paper values, plus the two atmospheric-flux systematics
    # that are not fitted in these analyses.
    conditional_vars = Dict(:δCP=>0.0, :θ₁₂=>ftype(33.62 * π/180), :θ₁₃=>ftype(8.54  * π/180), :Δm²₂₁=>ftype(7.40e-5),
                            :atm_flux_updown_sigma=>ftype(0.0), :atm_flux_nuenuebar_sigma=>ftype(0.0))
    priors = Newtrinos.condition(priors, conditional_vars, p)
    @reset priors.θ₂₃   = Uniform(ftype(30*π/180), ftype(60*π/180))
    @reset priors.Δm²₃₁ = Uniform(ftype(0.93e-3 + 7.40e-5), ftype(3.93e-3 + 7.40e-5))  # paper gives a range for Δm²₂₃; flat prior on the sum (NO, Δm²₂₁ fixed)

    # Detector-systematic prior ranges as in the paper (Table III)
    @reset priors.deepcore_lifetime        = Uniform(ftype(0), ftype(3.8))   # lifetime: [0, 3.8] years
    @reset priors.deepcore_ice_scattering  = Truncated(Normal(ftype(1), ftype(0.1)), ftype(0.9), ftype(1.1))   # truncated at ±1σ
    @reset priors.deepcore_ice_absorption  = Truncated(Normal(ftype(1), ftype(0.1)), ftype(0.9), ftype(1.1))
    @reset priors.deepcore_opt_eff_lateral = Truncated(Normal(ftype(0), ftype(1.)), ftype(-2), ftype(2.5))      # asymmetric: +2.5σ/-2σ

    # NSI parameter bounds, Standard parametrization
    @reset priors.Δ_eμ     = Uniform(-ftype(5),   ftype(5))
    @reset priors.Δ_τμ     = Uniform(-ftype(0.1), ftype(0.1))
    @reset priors.ε_eμ_abs = Uniform(ftype(0),    ftype(0.3))
    @reset priors.ε_eτ_abs = Uniform(ftype(0),    ftype(0.35))
    @reset priors.ε_μτ_abs = Uniform(ftype(0),    ftype(0.07))
    @reset priors.δ_eμ     = Uniform(ftype(0),    ftype(2π))   # [0, 2π] so that a general grid covers the whole phase
    @reset priors.δ_eτ     = Uniform(ftype(0),    ftype(2π))
    @reset priors.δ_μτ     = Uniform(ftype(0),    ftype(2π))

    ### HYPOTHESES ###
    # Every hypothesis frees the listed NSI parameters and pins all other NSI parameters to 0 (their null value).
    # The tested names are listed magnitude before phase; for "all" in the order of `nsi_names`.

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
        "all"   => nsi_hypothesis(nsi_names),
    )

    (; name = "deepcore_3y", experiments, params = p, priors, hypotheses)
end
