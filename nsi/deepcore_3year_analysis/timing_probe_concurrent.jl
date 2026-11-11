# Concurrency-contention probe: runs N copies of the EXACT same PSO+LBFGS task as
# timing_probe.jl, concurrently via Threads.@threads, to see whether per-task wall-clock
# time degrades under real concurrent load compared to the isolated single-task baseline
# (timing_probe.jl measured 3.58h PSO + 139.6s LBFGS = 3.62h total, single-threaded).
#
# Each of the N tasks is independent (own RNG, own truth-injected likelihood built from
# the SAME grid point but a different seed index so Asimov/fluctuation draws differ
# slightly -- doesn't matter for this probe, only the TIMING matters, not the physics
# result). If per-task time stays close to ~3.6h even at N-way concurrency, contention is
# not the explanation and something else is going on. If it balloons to several times
# that, contention is confirmed as (at least part of) the explanation for the 24h timeout.
#
# Usage: julia --project=.. --threads=N timing_probe_concurrent.jl N [pso_particles] [global_iterations] [lbfgs_iterations]
# Lighter settings (e.g. pso_particles=50 global_iterations=30 lbfgs_iterations=300) give a
# much faster contention signal -- run once with N_CONCURRENT=1 to get a light isolated
# baseline, then again with the real N to compare per-task slowdown at the same settings.

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
using Base.Threads
using Statistics

include(joinpath(@__DIR__, "..", "..", "src", "analysis", "optimizer_tests_common.jl"))

const N_CONCURRENT = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 8

### Settings -- default to the cancelled job (228104) exactly; override via extra args for
### a much faster/cheaper contention signal, e.g.:
###   julia --threads=16 timing_probe_concurrent.jl 16 50 30 300
const PSO_PARTICLES = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 500
const GLOBAL_ITERATIONS = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 200
const LBFGS_ITERATIONS = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 2000
const G_TOL = 1e-6
const STEP_SIZE = 1.0
const FLUCTUATE_NUISANCES = true
const OBJECTIVE = "posterior"

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

truths, scan_keys, grid_coords = truth_grid(cp_hyp, p, nsi_param_names; n_points=4)
truth_param_nominal = truths[1]
algorithm = make_algorithm("lbfgs", STEP_SIZE)

println("="^80)
println("Running $N_CONCURRENT concurrent PSO+LBFGS tasks (pso_particles=$PSO_PARTICLES, " *
        "global_iterations=$GLOBAL_ITERATIONS, lbfgs_iterations=$LBFGS_ITERATIONS) on $(Threads.nthreads()) threads...")

task_elapsed = Vector{NTuple{3,Float64}}(undef, N_CONCURRENT)  # (setup, pso, lbfgs)

Threads.@threads for i in 1:N_CONCURRENT
    seed_rng = Random.Xoshiro(1000 + i)
    truth_param = FLUCTUATE_NUISANCES ?
        fluctuate_nuisances(seed_rng, truth_param_nominal, cp_hyp, fluctuated_names) : truth_param_nominal

    t_setup0 = time()
    injected_data = Newtrinos.generate_asimov_data(experiments, truth_param)
    truth_likelihood = Newtrinos.generate_likelihood(experiments, injected_data)
    prior_dist = distprod(;cp_hyp...)
    truth_posterior = PosteriorMeasure(truth_likelihood, prior_dist)
    truth_log_post = logdensityof(truth_posterior, truth_param)
    setup_elapsed = time() - t_setup0

    attempt_rng = Random.Xoshiro(1000 * i + 1)

    t_pso0 = time()
    start_param = global_seed_search(truth_likelihood, cp_hyp, truth_param, "pso", GLOBAL_ITERATIONS, attempt_rng;
                                      n_particles=PSO_PARTICLES, objective=OBJECTIVE)
    pso_elapsed = time() - t_pso0

    t_lbfgs0 = time()
    llh, log_post, fit_param, converged, n_iters = local_find_mle(truth_likelihood, prior_dist, start_param;
        fit_method="optim", algorithm=algorithm, iterations=LBFGS_ITERATIONS, g_tol=G_TOL, f_tol=0.0, x_tol=0.0,
        ad_backend=ADTypes.AutoForwardDiff(), objective=OBJECTIVE)
    lbfgs_elapsed = time() - t_lbfgs0

    task_elapsed[i] = (setup_elapsed, pso_elapsed, lbfgs_elapsed)
    println("  [task $i/$N_CONCURRENT] setup=$(round(setup_elapsed,digits=1))s " *
            "pso=$(round(pso_elapsed,digits=1))s ($(round(pso_elapsed/3600,digits=2))h) " *
            "lbfgs=$(round(lbfgs_elapsed,digits=1))s ($(round(lbfgs_elapsed/3600,digits=2))h) " *
            "converged=$converged n_iters=$n_iters")
end

println("="^80)
println("SUMMARY ($N_CONCURRENT concurrent tasks on $(Threads.nthreads()) threads, " *
        "pso_particles=$PSO_PARTICLES, global_iterations=$GLOBAL_ITERATIONS, lbfgs_iterations=$LBFGS_ITERATIONS):")
pso_times = [t[2] for t in task_elapsed]
lbfgs_times = [t[3] for t in task_elapsed]
println("  PSO:   mean=$(round(mean(pso_times),digits=2))s min=$(round(minimum(pso_times),digits=2))s max=$(round(maximum(pso_times),digits=2))s")
println("  LBFGS: mean=$(round(mean(lbfgs_times),digits=2))s min=$(round(minimum(lbfgs_times),digits=2))s max=$(round(maximum(lbfgs_times),digits=2))s")
println("  (compare against a separate N_CONCURRENT=1 run at the SAME pso_particles/global_iterations/lbfgs_iterations for the isolated baseline)")
