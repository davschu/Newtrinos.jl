# Hybrid contention probe: N_WORKERS separate Julia processes (own GC each, like
# timing_probe_distributed.jl), each internally running N_THREADS_PER_WORKER concurrent
# tasks via Threads.@threads (like timing_probe_concurrent.jl) -- total concurrency =
# N_WORKERS * N_THREADS_PER_WORKER. Tests whether this hybrid shape lands at a better
# contention sweet spot than pure Distributed (many single-threaded workers) or pure
# Threads (one process, many threads) at the same total concurrency.
#
# Usage: julia --project=.. -p N_WORKERS timing_probe_hybrid.jl N_WORKERS N_THREADS_PER_WORKER [pso_particles] [global_iterations] [lbfgs_iterations]
# Each worker is started with --threads=N_THREADS_PER_WORKER via addprocs' exeflags.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Distributed

const N_WORKERS = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 4
const N_THREADS_PER_WORKER = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
const PSO_PARTICLES = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 50
const GLOBAL_ITERATIONS = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 30
const LBFGS_ITERATIONS = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 300
const N_CONCURRENT = N_WORKERS * N_THREADS_PER_WORKER

println("nworkers() = $(nworkers()) (expected $N_WORKERS -- launch with `-p $N_WORKERS`), " *
        "each with $N_THREADS_PER_WORKER threads -> total concurrency $N_CONCURRENT")

@everywhere begin
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
end

# Separate @everywhere block: `@reset` (Accessors macro) must not be expanded in the same
# atomically-`eval`'d block as `using Accessors` above (world-age ordering) -- same fix as
# timing_probe_distributed.jl.
@everywhere begin
    include(joinpath(@__DIR__, "..", "..", "src", "analysis", "optimizer_tests_common.jl"))

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

    function run_one_task(i::Int, pso_particles::Int, global_iterations::Int, lbfgs_iterations::Int,
                           experiments, cp_hyp, p, truth_param_nominal, fluctuated_names)
        algorithm = make_algorithm("lbfgs", 1.0)

        seed_rng = Random.Xoshiro(1000 + i)
        truth_param = fluctuate_nuisances(seed_rng, truth_param_nominal, cp_hyp, fluctuated_names)

        t_setup0 = time()
        injected_data = Newtrinos.generate_asimov_data(experiments, truth_param)
        truth_likelihood = Newtrinos.generate_likelihood(experiments, injected_data)
        prior_dist = distprod(;cp_hyp...)
        setup_elapsed = time() - t_setup0

        attempt_rng = Random.Xoshiro(1000 * i + 1)

        t_pso0 = time()
        start_param = global_seed_search(truth_likelihood, cp_hyp, truth_param, "pso", global_iterations, attempt_rng;
                                          n_particles=pso_particles, objective="posterior")
        pso_elapsed = time() - t_pso0

        t_lbfgs0 = time()
        llh, log_post, fit_param, converged, n_iters = local_find_mle(truth_likelihood, prior_dist, start_param;
            fit_method="optim", algorithm=algorithm, iterations=lbfgs_iterations, g_tol=1e-6, f_tol=0.0, x_tol=0.0,
            ad_backend=ADTypes.AutoForwardDiff(), objective="posterior")
        lbfgs_elapsed = time() - t_lbfgs0

        return (i=i, setup=setup_elapsed, pso=pso_elapsed, lbfgs=lbfgs_elapsed, converged=converged, n_iters=n_iters)
    end

    # Runs on ONE worker: builds physics once (shared read-only across this worker's own
    # threads, same pattern as roundtrips.jl sharing `experiments` across grid points),
    # then runs `n_tasks_here` of run_one_task CONCURRENTLY via Threads.@threads using this
    # worker's own --threads=N_THREADS_PER_WORKER. task_ids are GLOBAL indices (unique
    # across all workers) so each task still gets a distinct, reproducible RNG seed.
    function run_worker_batch(task_ids::Vector{Int}, pso_particles::Int, global_iterations::Int, lbfgs_iterations::Int)
        ftype = Float64
        experiments = configure_physics()
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

        n_here = length(task_ids)
        batch_results = Vector{Any}(undef, n_here)
        Threads.@threads for j in 1:n_here
            batch_results[j] = run_one_task(task_ids[j], pso_particles, global_iterations, lbfgs_iterations,
                                             experiments, cp_hyp, p, truth_param_nominal, fluctuated_names)
        end
        return batch_results
    end
end

println("="^80)
println("Running $N_CONCURRENT PSO+LBFGS tasks on $N_WORKERS workers x $N_THREADS_PER_WORKER threads each " *
        "(pso_particles=$PSO_PARTICLES, global_iterations=$GLOBAL_ITERATIONS, lbfgs_iterations=$LBFGS_ITERATIONS)...")

# Partition 1:N_CONCURRENT into N_WORKERS contiguous chunks of N_THREADS_PER_WORKER task ids each.
chunks = [collect((w-1)*N_THREADS_PER_WORKER+1 : w*N_THREADS_PER_WORKER) for w in 1:N_WORKERS]

t0 = time()
batch_results = pmap(chunk -> run_worker_batch(chunk, PSO_PARTICLES, GLOBAL_ITERATIONS, LBFGS_ITERATIONS), chunks)
total_wall = time() - t0
results = reduce(vcat, batch_results)

for r in results
    println("  [task $(r.i)/$N_CONCURRENT] setup=$(round(r.setup,digits=1))s " *
            "pso=$(round(r.pso,digits=1))s ($(round(r.pso/3600,digits=2))h) " *
            "lbfgs=$(round(r.lbfgs,digits=1))s ($(round(r.lbfgs/3600,digits=2))h) " *
            "converged=$(r.converged) n_iters=$(r.n_iters)")
end

println("="^80)
println("SUMMARY ($N_WORKERS workers x $N_THREADS_PER_WORKER threads = $N_CONCURRENT total, " *
        "pso_particles=$PSO_PARTICLES, global_iterations=$GLOBAL_ITERATIONS, lbfgs_iterations=$LBFGS_ITERATIONS):")
pso_times = [r.pso for r in results]
lbfgs_times = [r.lbfgs for r in results]
using Statistics
println("  PSO:   mean=$(round(mean(pso_times),digits=2))s min=$(round(minimum(pso_times),digits=2))s max=$(round(maximum(pso_times),digits=2))s")
println("  LBFGS: mean=$(round(mean(lbfgs_times),digits=2))s min=$(round(minimum(lbfgs_times),digits=2))s max=$(round(maximum(lbfgs_times),digits=2))s")
println("  total wall clock: $(round(total_wall,digits=1))s")
