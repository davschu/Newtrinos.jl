# Throwaway timing probe: measures ONE realistic (k=1, s=1) PSO+LBFGS task in isolation,
# at the real job's settings, to find out whether chaining PSO straight into LBFGS (as
# roundtrips.jl does per seed-attempt) costs meaningfully more than the ~4.4h PSO-alone /
# ~1h truth-seeded-LBFGS figures from test_pso.jl -- which never chains the two stages.
#
# Deliberately single-threaded / no grid-point or seed-attempt parallelism: this measures
# one task's own wall-clock cost, not how it behaves under 75-way contention (that's a
# separate question, already addressed earlier). Run with `julia --threads=1` or similar;
# thread count here doesn't matter since nothing in this script is parallelized.
#
# Usage: julia --project=.. timing_probe.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Distributions
using DensityInterface
using BAT
using DataStructures
using MeasureBase
using ADTypes
using Newtrinos
using Accessors
using Random
import Optim
import ForwardDiff
using ValueShapes

include(joinpath(@__DIR__, "..", "..", "src", "analysis", "optimizer_tests_common.jl"))

### Settings -- matching the cancelled job (228104) exactly ###
const HYPOTHESIS = "μτ"
const PSO_PARTICLES = 500
const GLOBAL_ITERATIONS = 200
const LBFGS_ITERATIONS = 2000
const G_TOL = 1e-6
const STEP_SIZE = 1.0
const FLUCTUATE_NUISANCES = true
const OBJECTIVE = "posterior"

### PHYSICS CONFIG (verbatim from roundtrips.jl:244-256) ###

function configure_physics()
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
    (; deepcore = Newtrinos.deepcore.configure(physics))
end

println("Configuring physics...")
experiments = configure_physics()

ftype = Float64
p = Newtrinos.get_params(experiments)
priors = Newtrinos.get_priors(experiments)
conditional_vars = Dict(:δCP=>0.0, :θ₁₂=>ftype(33.62 * π/180), :θ₁₃=>ftype(8.54  * π/180), :Δm²₂₁=>ftype(7.40e-5),
                         :atm_flux_updown_sigma=>ftype(0.0), :atm_flux_nuenuebar_sigma=>ftype(0.0))
priors = Newtrinos.condition(priors, conditional_vars, p)
@reset priors.θ₂₃   = Uniform(ftype(30*π/180), ftype(60*π/180))
@reset priors.Δm²₃₁ = Uniform(ftype(0.93e-3 + 7.40e-5), ftype(3.93e-3 + 7.40e-5))

@reset priors.deepcore_lifetime        = Uniform(ftype(0), ftype(3.8))
@reset priors.deepcore_ice_scattering  = Truncated(Normal(ftype(1), ftype(0.1)), ftype(0.9), ftype(1.1))
@reset priors.deepcore_ice_absorption  = Truncated(Normal(ftype(1), ftype(0.1)), ftype(0.9), ftype(1.1))
@reset priors.deepcore_opt_eff_lateral = Truncated(Normal(ftype(0), ftype(1.)), ftype(-2), ftype(2.5))

@reset priors.Δ_eμ     = Uniform(-ftype(5),   ftype(5))
@reset priors.Δ_τμ     = Uniform(-ftype(0.1), ftype(0.1))
@reset priors.ε_eμ_abs = Uniform(ftype(0),    ftype(0.3))
@reset priors.ε_eτ_abs = Uniform(ftype(0),    ftype(0.35))
@reset priors.ε_μτ_abs = Uniform(ftype(0),    ftype(0.07))
@reset priors.δ_eμ     = Uniform(ftype(0),    ftype(2π))
@reset priors.δ_eτ     = Uniform(ftype(0),    ftype(2π))
@reset priors.δ_μτ     = Uniform(ftype(0),    ftype(2π))

cp_μτ = Newtrinos.condition(priors, Dict(:Δ_eμ => 0.0, :Δ_τμ => 0.0, :ε_eμ_abs => 0.0, :ε_eτ_abs => 0.0, :δ_eμ => 0.0, :δ_eτ => 0.0), p)

cp_hyp = cp_μτ
nsi_param_names = (:Δ_eμ, :Δ_τμ, :ε_eμ_abs, :ε_eτ_abs, :ε_μτ_abs, :δ_eμ, :δ_eτ, :δ_μτ)

nuisance_names = (
    :atm_flux_delta_spectral_index, :atm_flux_nuenuebar_sigma, :atm_flux_nuenumu_sigma,
    :atm_flux_numunumubar_sigma, :atm_flux_updown_sigma, :atm_flux_uphorizontal_sigma,
    :deepcore_atm_muon_scale, :deepcore_ice_absorption, :deepcore_ice_scattering,
    :deepcore_lifetime, :deepcore_opt_eff_headon, :deepcore_opt_eff_lateral,
    :deepcore_opt_eff_overall, :nc_norm, :nutau_cc_norm,
)
physics_names = (:θ₂₃, :Δm²₃₁)
fluctuated_names = (nuisance_names..., physics_names...)

### ONE GRID POINT: use the hypothesis's own nominal/truth point (NSI params at their
### midpoint via truth_grid with n_points=1 is awkward for a 2D grid -- just build the
### truth directly at nominal NSI values, matching "truth = 0 shift" for simplicity) ###

truths, scan_keys, grid_coords = truth_grid(cp_hyp, p, nsi_param_names; n_points=4)
k = 1  # first of the 2x2 grid for μτ
truth_param_nominal = truths[k]

println("="^80)
println("Settings: hypothesis=$HYPOTHESIS pso_particles=$PSO_PARTICLES global_iterations=$GLOBAL_ITERATIONS " *
        "lbfgs_iterations=$LBFGS_ITERATIONS g_tol=$G_TOL fluctuate_nuisances=$FLUCTUATE_NUISANCES")
println("Grid point coord: $(grid_coords[k]), scan_keys=$scan_keys")

seed_rng = Random.Xoshiro(1000 + k)
truth_param = FLUCTUATE_NUISANCES ?
    fluctuate_nuisances(seed_rng, truth_param_nominal, cp_hyp, fluctuated_names) : truth_param_nominal

println("Generating Asimov data + likelihood...")
t_setup0 = time()
injected_data = Newtrinos.generate_asimov_data(experiments, truth_param)
truth_likelihood = Newtrinos.generate_likelihood(experiments, injected_data)
prior_dist = distprod(;cp_hyp...)
truth_posterior = PosteriorMeasure(truth_likelihood, prior_dist)
truth_log_post = logdensityof(truth_posterior, truth_param)
setup_elapsed = time() - t_setup0
println("Setup done in $(round(setup_elapsed, digits=1))s. truth_log_post=$(round(truth_log_post, digits=4))")

algorithm = make_algorithm("lbfgs", STEP_SIZE)
attempt_rng = Random.Xoshiro(1000 * k + 1)

### STAGE 1: PSO global seed search ###
println("="^80)
println("Starting PSO global seed search (pso_particles=$PSO_PARTICLES, global_iterations=$GLOBAL_ITERATIONS)...")
t_pso0 = time()
start_param = global_seed_search(truth_likelihood, cp_hyp, truth_param, "pso", GLOBAL_ITERATIONS, attempt_rng;
                                  n_particles=PSO_PARTICLES, objective=OBJECTIVE)
pso_elapsed = time() - t_pso0
pso_log_post = logdensityof(truth_posterior, start_param)
println("PSO done in $(round(pso_elapsed, digits=1))s ($(round(pso_elapsed/3600, digits=2))h). " *
        "PSO-result log_post=$(round(pso_log_post, digits=4)) (truth_log_post=$(round(truth_log_post, digits=4)))")

### STAGE 2: LBFGS polish from PSO's result (NOT from truth) ###
println("="^80)
println("Starting LBFGS polish from PSO's seed (iterations=$LBFGS_ITERATIONS, g_tol=$G_TOL)...")
t_lbfgs0 = time()
llh, log_post, fit_param, converged, n_iters = local_find_mle(truth_likelihood, prior_dist, start_param;
    fit_method="optim", algorithm=algorithm, iterations=LBFGS_ITERATIONS, g_tol=G_TOL, f_tol=0.0, x_tol=0.0,
    ad_backend=ADTypes.AutoForwardDiff(), objective=OBJECTIVE)
lbfgs_elapsed = time() - t_lbfgs0
println("LBFGS done in $(round(lbfgs_elapsed, digits=1))s ($(round(lbfgs_elapsed/3600, digits=2))h). " *
        "converged=$converged n_iters=$n_iters log_post=$(round(log_post, digits=4))")

println("="^80)
println("SUMMARY (one (k=1,s=1) task, chained PSO->LBFGS, no concurrency):")
println("  setup:  $(round(setup_elapsed, digits=1))s")
println("  PSO:    $(round(pso_elapsed, digits=1))s ($(round(pso_elapsed/3600, digits=2))h)")
println("  LBFGS:  $(round(lbfgs_elapsed, digits=1))s ($(round(lbfgs_elapsed/3600, digits=2))h)")
println("  total:  $(round(setup_elapsed+pso_elapsed+lbfgs_elapsed, digits=1))s ($(round((setup_elapsed+pso_elapsed+lbfgs_elapsed)/3600, digits=2))h)")
println("  final Δχ²_recovery = $(round(2*(log_post - truth_log_post), digits=4))")
