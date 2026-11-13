# ARCHIVED (reference only). Setup-driven version of this script, replaced by src/analysis/optimizer_tests.jl
# (compare optimizers on the same problems: --optimizers SPEC...). It still calls `generate_seed`, which was removed
# when strategies (find_mle_advanced) were introduced, so it no longer runs against the current helpers; the exact scripts
# that ran are available via the git tag `nsi-3y-scripts` and the frozen copy in this folder.
#
# Isolated head-to-head comparison of Optim.jl's PSO (the global-seed-search algorithm
# currently used in nsi_scan.jl/test_pso.jl) against Metaheuristics.jl's PSO, DE, and ABC,
# on the SAME objective, bounds, and truth points, under an EQUAL total function-evaluation
# budget (f_calls_limit for Metaheuristics == n_particles*global_iterations for Optim).
#
# This is a diagnostic-only script to decide whether Metaheuristics.jl is worth adopting --
# it does NOT touch nsi_scan.jl, test_pso.jl, roundtrips.jl, or optimizer_tests_common.jl.
#
# Hypothesis tested: one hypothesis of the setup file (default "eμ": 2 tested parameters) -- the
# simplest case, fast enough to sanity-check trace shapes before considering a harder (e.g. all-parameter) case.
#
# For each of N_TRUTH independently drawn truth points, every algorithm is run ONCE with
# the SAME total evaluation budget, and its best-so-far objective value is recorded against
# cumulative function evaluations (not generation index, since algorithms take different
# numbers of evaluations per generation -- evaluations is the fair common x-axis).
#
# A single truth-seeded LBFGS reference fit per truth point (exactly the Phase-1 pattern in
# test_pso.jl) gives the "actual optimum" each algorithm is compared against.
#
# The N_TRUTH loop below is `Threads.@threads`-parallelized (each truth point is fully
# independent: its own RNG, own Asimov data/likelihood, own result slot) -- run with
# `julia --threads=N` to use it; per-thread progress lines may interleave in stdout.
#
# Usage: julia --threads=4 nsi/compare_pso_metaheuristics.jl [setup file [hypothesis name [results dir]]]
#   defaults: nsi/setups/nsi_deepcore_3y.jl, "eμ", results/<setup name>/pso_metaheuristics
#
# Output (in the results dir): compare_<hypothesis>_<n>points.{jld2,pdf}

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Distributions
using DensityInterface
using BAT
using MeasureBase
using ADTypes
using Newtrinos
using FileIO
using Accessors
using Random
using ValueShapes
using CairoMakie
using Statistics
import Optim
import Metaheuristics

include(joinpath(@__DIR__, "..", "..", "src", "analysis", "optimizer_tests_common.jl"))

### SETTINGS ###

const N_TRUTH      = 4      # independent random truth points
const N_PARTICLES  = 100     # population size, shared across all 4 algorithms
const GLOBAL_ITERS = 100     # Optim PSO generation budget -- defines the shared eval budget
const FEVALS       = N_PARTICLES * GLOBAL_ITERS
const SEED         = 1234
const SETUP_FILE   = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "setups", "nsi_deepcore_3y.jl")
const HYPOTHESIS   = length(ARGS) >= 2 ? ARGS[2] : "eμ"

# LBFGS settings for the truth-seeded reference fit (mirrors test_pso.jl's REF_LBFGS_*).
const REF_LBFGS_ITERATIONS = 1000
const REF_LBFGS_GTOL       = 1e-6
const REF_LBFGS_STEP_SIZE  = 1.0
const REF_LBFGS_FTOL       = 0.0
const REF_LBFGS_XTOL       = 0.0

### SETUP (physics, experiments, priors, hypotheses: defined by the setup file) ###

include(abspath(SETUP_FILE))
S = finalize_setup(setup(; ftype = Float64))
hyp = hypothesis_info(S, HYPOTHESIS)

experiments = S.experiments
p = S.params
cp_hyp = hyp.prior                  # prior of the hypothesis under test
tested_names = hyp.tested           # tested parameters = keys of the hypothesis' `null`
scan_keys = collect(tested_names)
const RESULTS_DIR = length(ARGS) >= 3 ? ARGS[3] : joinpath("results", S.name, "pso_metaheuristics")

prior_dist = distprod(;cp_hyp...)
algorithm  = make_algorithm("lbfgs", REF_LBFGS_STEP_SIZE)

println("Hypothesis '$HYPOTHESIS': scan_keys=$(Tuple(scan_keys)), n_truth=$N_TRUTH, " *
        "n_particles=$N_PARTICLES, global_iters=$GLOBAL_ITERS, total_fevals=$FEVALS")

### OBJECTIVE (same posterior objective generate_seed uses, built directly here so it
### can be handed to both Optim.optimize and Metaheuristics.optimize without going through
### generate_seed's own Optim-specific call) ###

function make_objective(likelihood, b)
    function obj(x::AbstractVector)
        all_x = [b.is_free_mask[i] ? x[b.free_x_index[i]] : b.fixed_float[i] for i in 1:length(b.all_keys)]
        full_params = NamedTuple{b.all_keys}(Tuple(all_x))
        llh_val = logdensityof(likelihood, full_params)
        prior_val = sum(logpdf(b.free_dists[j], x[j]) for j in 1:length(b.free_dists))
        isfinite(llh_val) && isfinite(prior_val) ? -(llh_val + prior_val) : Inf
    end
end

### PER-TRUTH RUNS ###

struct AlgoRun
    name::String
    cumulative_evals::Vector{Int}
    best_so_far::Vector{Float64}
    final_obj::Float64
    elapsed::Float64
end

truth_results = Vector{NamedTuple}(undef, N_TRUTH)

println("="^80)
println("Running $N_TRUTH truth points x 4 algorithms (Optim-PSO, MH-PSO, MH-DE, MH-ABC)...")

Threads.@threads for t in 1:N_TRUTH
    truth_rng = Random.Xoshiro(SEED + t)
    truth_param = controlled_truth(truth_rng, p, cp_hyp, tested_names)

    injected_data    = Newtrinos.generate_asimov_data(experiments, truth_param)
    truth_likelihood = Newtrinos.generate_likelihood(experiments, injected_data)

    t0 = time()
    ref_llh, ref_log_post, ref_fit_param, ref_converged, ref_n_iters = local_find_mle(
        truth_likelihood, prior_dist, truth_param;
        fit_method="optim", algorithm=algorithm, iterations=REF_LBFGS_ITERATIONS,
        g_tol=REF_LBFGS_GTOL, f_tol=REF_LBFGS_FTOL, x_tol=REF_LBFGS_XTOL, objective="posterior",
        ad_backend=ADTypes.AutoForwardDiff())
    ref_elapsed = time() - t0
    ref_neglogpost = -ref_log_post

    println("  [truth $t/$N_TRUTH] ref_log_post=$(round(ref_log_post,digits=4)) " *
            "ref_converged=$ref_converged ref_elapsed=$(round(ref_elapsed,digits=1))s")

    b = extract_bounds(cp_hyp, p)
    obj = make_objective(truth_likelihood, b)

    algo_runs = Vector{AlgoRun}(undef, 4)

    # 1) Optim.jl PSO -- reuse generate_seed as-is, no reimplementation.
    rng1 = Random.Xoshiro(SEED + 1_000_000 + t)
    t0 = time()
    _, f_trace = generate_seed(truth_likelihood, cp_hyp, p, "pso", GLOBAL_ITERS, rng1;
                                     n_particles=N_PARTICLES, return_trace=true, objective="posterior")
    elapsed1 = time() - t0
    best_so_far1 = accumulate(min, f_trace)
    cum_evals1 = N_PARTICLES .* collect(1:length(f_trace))
    algo_runs[1] = AlgoRun("Optim-PSO", cum_evals1, best_so_far1, best_so_far1[end], elapsed1)

    # 2-4) Metaheuristics.jl PSO / DE / ABC, all capped at the SAME total f_calls_limit.
    # Note: ABC's evals-per-generation also depends on its Ne/No (employed/onlooker bee
    # counts, both default div(N+1,2)) in addition to N, so its generation count at a given
    # f_calls_limit need not match PSO/DE's -- f_calls_limit still equalizes the TOTAL
    # evaluation budget across all three, which is what this comparison controls for.
    bounds = permutedims(hcat(b.lo_vec, b.hi_vec))
    # `iterations` is set to a large value (not left at Options' default of 1) so that
    # f_calls_limit -- not a premature iteration cap -- is the binding stopping criterion;
    # Metaheuristics.jl stops as soon as ANY configured limit is hit.
    mh_opts(sd) = Metaheuristics.Options(f_calls_limit=FEVALS, iterations=10_000, store_convergence=true, seed=sd)
    mh_algos = [
        ("MH-PSO", Metaheuristics.PSO(N=N_PARTICLES, options=mh_opts(SEED+2_000_000+t))),
        ("MH-DE",  Metaheuristics.DE(N=N_PARTICLES,  options=mh_opts(SEED+3_000_000+t))),
        ("MH-ABC", Metaheuristics.ABC(N=N_PARTICLES, options=mh_opts(SEED+4_000_000+t))),
    ]
    for (i, (name, alg)) in enumerate(mh_algos)
        t0 = time()
        state = Metaheuristics.optimize(obj, bounds, alg)
        elapsed = time() - t0
        n_fes, fxs = Metaheuristics.convergence(state)
        best_so_far = accumulate(min, fxs)
        algo_runs[1 + i] = AlgoRun(name, collect(Int, n_fes), best_so_far, best_so_far[end], elapsed)
    end

    for r in algo_runs
        dchi2_final = 2 * (r.final_obj - ref_neglogpost)
        println("    $(rpad(r.name,10)) final_Δχ²_vs_reference=$(round(dchi2_final,digits=4)) " *
                 "elapsed=$(round(r.elapsed,digits=1))s n_evals_recorded=$(length(r.cumulative_evals))")
    end

    truth_results[t] = (truth_param=truth_param, ref_llh=ref_llh, ref_log_post=ref_log_post,
                         ref_fit_param=ref_fit_param, ref_converged=ref_converged,
                         ref_n_iters=ref_n_iters, ref_elapsed=ref_elapsed, algo_runs=algo_runs)
end

### SAVE RESULTS (.jld2, before plotting) ###

mkpath(RESULTS_DIR)
output_tag = "$(HYPOTHESIS)_$(N_TRUTH)points"

algo_names = [r.name for r in truth_results[1].algo_runs]

FileIO.save(joinpath(RESULTS_DIR, "compare_$(output_tag).jld2"), Dict(
    "hypothesis"  => HYPOTHESIS,
    "setup"       => S.name,
    "scan_keys"   => collect(scan_keys),
    "n_truth"     => N_TRUTH,
    "algo_names"  => algo_names,

    "n_particles"  => N_PARTICLES,
    "global_iters" => GLOBAL_ITERS,
    "total_fevals" => FEVALS,

    "truth_params"   => [tr.truth_param for tr in truth_results],
    "ref_llh"        => [tr.ref_llh for tr in truth_results],
    "ref_log_post"   => [tr.ref_log_post for tr in truth_results],
    "ref_fit_params" => [tr.ref_fit_param for tr in truth_results],
    "ref_converged"  => [tr.ref_converged for tr in truth_results],
    "ref_n_iters"    => [tr.ref_n_iters for tr in truth_results],
    "ref_elapsed"    => [tr.ref_elapsed for tr in truth_results],

    # Per (truth, algorithm) traces, keyed flat as "<algo>_cumulative_evals" / "..._best_so_far"
    # with one entry per truth point (same ordering as truth_params above), plus final
    # objective/elapsed -- enough to fully reconstruct the plot and summary table below
    # without re-running any optimizer.
    (vcat(
        [["$(name)_cumulative_evals" => [tr.algo_runs[i].cumulative_evals for tr in truth_results],
          "$(name)_best_so_far"      => [tr.algo_runs[i].best_so_far for tr in truth_results],
          "$(name)_final_obj"        => [tr.algo_runs[i].final_obj for tr in truth_results],
          "$(name)_elapsed"          => [tr.algo_runs[i].elapsed for tr in truth_results]]
         for (i, name) in enumerate(algo_names)]...)...),

    "settings" => Dict(
        "hypothesis"   => HYPOTHESIS,
        "n_truth"      => N_TRUTH,
        "n_particles"  => N_PARTICLES,
        "global_iters" => GLOBAL_ITERS,
        "total_fevals" => FEVALS,
        "seed"         => SEED,
        "ref_lbfgs_iterations" => REF_LBFGS_ITERATIONS,
        "ref_lbfgs_gtol"       => REF_LBFGS_GTOL,
        "ref_lbfgs_step_size"  => REF_LBFGS_STEP_SIZE,
    ),
))
println("Saved results to $(joinpath(RESULTS_DIR, "compare_$(output_tag).jld2"))")

### SUMMARY TABLE (stdout) ###

println("="^80)
println("Summary: mean final Δχ²_vs_reference and mean wall-clock time per algorithm")
for (i, name) in enumerate(algo_names)
    dchi2s = [2 * (tr.algo_runs[i].final_obj - (-tr.ref_log_post)) for tr in truth_results]
    elapsed = [tr.algo_runs[i].elapsed for tr in truth_results]
    println("  $(rpad(name,10)) mean_Δχ²=$(round(mean(dchi2s),digits=4)) " *
            "std_Δχ²=$(round(std(dchi2s),digits=4)) mean_elapsed=$(round(mean(elapsed),digits=2))s")
end

### PLOT: best-so-far objective vs. cumulative function evaluations, one panel per truth point ###

n_cols = min(N_TRUTH, 2)
n_rows = cld(N_TRUTH, n_cols)
fig = Figure(size=(600 * n_cols, 450 * n_rows))

algo_colors = Dict(zip(algo_names, Makie.wong_colors()[1:length(algo_names)]))

for t in 1:N_TRUTH
    row = (t - 1) ÷ n_cols + 1
    col = (t - 1) % n_cols + 1
    tr = truth_results[t]
    ref_neglogpost = -tr.ref_log_post

    # pseudolog10: Δχ²_vs_reference is SIGNED (negative = an algorithm's trajectory beat the
    # LBFGS reference optimum, possible since the reference is only a single local fit) --
    # same convention as test_pso.jl's trace plot, see that file for the full rationale.
    ax = Axis(fig[row, col], xlabel="cumulative function evaluations",
              ylabel="Δχ²_vs_reference = 2×(best_so_far - ref_neglogpost)",
              yscale=Makie.pseudolog10, title="truth point $t")

    for r in tr.algo_runs
        Δχ²_vs_reference = 2 .* (r.best_so_far .- ref_neglogpost)
        lines!(ax, r.cumulative_evals, Δχ²_vs_reference; color=algo_colors[r.name], linewidth=2, label=r.name)
    end
    hlines!(ax, [0.0]; color=:black, linestyle=:dash, linewidth=2)
    t == 1 && axislegend(ax; position=:rt, framevisible=false)
end

Label(fig[0, 1:n_cols], "PSO/DE/ABC convergence comparison -- hypothesis=$HYPOTHESIS, n_truth=$N_TRUTH, " *
      "n_particles=$N_PARTICLES, total_fevals=$FEVALS", fontsize=16, font=:bold)

resize_to_layout!(fig)
save(joinpath(RESULTS_DIR, "compare_$(output_tag).pdf"), fig; pt_per_unit=1.0)
println("Saved plot to $(joinpath(RESULTS_DIR, "compare_$(output_tag).pdf"))")
