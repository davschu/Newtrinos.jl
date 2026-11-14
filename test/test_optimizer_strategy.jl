# Tests of the optimizer registry and of strategies (chains of stages) in src/analysis/optimizer_tests_common.jl,
# on a toy likelihood with a known optimum -- no experiment data needed.

module OptimizerStrategyTests

using Test
using Random
using Distributions
using DensityInterface
using BAT
using MeasureBase
using ValueShapes
using ADTypes
using Newtrinos
using Accessors
using DataStructures
import Optim
using FileIO

include(joinpath(@__DIR__, "..", "src", "analysis", "optimizer_tests_common.jl"))

# Gaussian likelihood around (a, b, c) = (0.3, -0.2, 1.0) (llh = -0.1 at the optimum because of the -0.1 d term, d fixed at 1)
struct ToyLikelihood end
DensityInterface.DensityKind(::ToyLikelihood) = DensityInterface.IsDensity()
DensityInterface.logdensityof(::ToyLikelihood, p) = -0.5 * ((p.a - 0.3)^2 / 0.1^2 + (p.b + 0.2)^2 / 0.2^2 + (p.c - 1.0)^2 / 0.5^2) - 0.1 * p.d

const PRIOR = (a = Uniform(-1.0, 1.0), b = Uniform(-1.0, 1.0), c = Uniform(0.0, 3.0), d = 1.0)
const NULL_PRIOR = (a = Uniform(-1.0, 1.0), b = 0.0, c = Uniform(0.0, 3.0), d = 1.0)   # tested parameter b pinned at its null value
const START = (a = 0.0, b = 0.0, c = 1.5, d = 1.0)
const OPTIMUM_LLH = -0.1

@testset "optimizer registry and strategies" begin

    @testset "spec parsing" begin
        sp = parse_optimizer_spec("mh-de:N=30,F=0.5")
        @test sp.name == "mh-de" && sp.settings == (N = 30, F = 0.5) && sp.label == "mh-de:N=30,F=0.5"
        @test parse_optimizer_spec("lbfgsb").settings == (;)
        @test parse_optimizer_spec("lbfgsb:iterations=50,reltol=1e-9").settings == (iterations = 50, reltol = 1e-9)
        @test parse_optimizer_spec("lbfgsb:ad=forward").settings.ad === :forward
        @test parse_optimizer_spec("mh-de:evals=500").settings.evals == 500
        @test parse_optimizer_spec("manual:a=0.2,b=-1").settings == (a = 0.2, b = -1)    # a manual stage takes any parameter name
        @test_throws ErrorException parse_optimizer_spec("nope")
        @test_throws ErrorException parse_optimizer_spec("lbfgsb:foo=1")
        @test_throws ErrorException parse_optimizer_spec("lbfgsb:iterations")
    end

    @testset "strategy order and aliases" begin
        @test [s.name for s in normalize_strategy(["lbfgsb"])] == ["lbfgsb"]
        @test [s.name for s in normalize_strategy(["mh-de", "lbfgsb", "lbfgsb"])] == ["mh-de", "lbfgsb", "lbfgsb"]
        @test [s.name for s in normalize_strategy(["random", "lbfgsb"])] == ["random", "lbfgsb"]
        @test [(s.name, s.null) for s in normalize_strategy(["smseed", "lbfgsb"])] == [("lbfgsb", true), ("lbfgsb", false)]
        @test [s.name for s in normalize_strategy(["mle", "lbfgsb"])] == ["random", "lbfgsb", "lbfgsb"]
        # seed, alias and global stages only first; a global stage after a local one (or another global) is rejected
        @test_throws ErrorException normalize_strategy(["lbfgsb", "mh-de"])
        @test_throws ErrorException normalize_strategy(["optim-pso", "lbfgsb", "optim-pso"])
        @test_throws ErrorException normalize_strategy(["mh-de", "mh-pso"])
        @test_throws ErrorException normalize_strategy(["lbfgsb", "random"])
        @test_throws ErrorException normalize_strategy(["random", "smseed"])
        @test_throws ErrorException normalize_strategy(["lbfgsb", "manual:a=1.0"])
        @test_throws ErrorException normalize_strategy(String[])
    end

    @testset "single local stage" begin
        r = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["lbfgsb"], rng = Random.Xoshiro(1))
        @test r.llh ≈ OPTIMUM_LLH atol = 1e-4
        @test r.converged === true
        @test r.params.d == 1.0                     # the fixed parameter stays fixed
        @test length(r.stages) == 1 && r.seed_params == START && r.strategy == "lbfgsb"
        @test r.n_calls > 0 && r.time > 0
    end

    @testset "a single global stage is not polished" begin
        r = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["mh-de:N=20,evals=300"], rng = Random.Xoshiro(1))
        @test ismissing(r.converged)                 # a global search has no convergence flag
        @test r.n_calls == 300                       # `evals` is the stopping rule
        @test length(r.stages) == 1
        @test isfinite(r.llh) && r.llh <= OPTIMUM_LLH + 1e-9
        # the harness default budget applies when a global stage gives no `evals`, an explicit `evals` beats it
        @test find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["mh-de:N=20"], budget = 120, rng = Random.Xoshiro(1)).n_calls == 120
        @test find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["mh-de:N=20,evals=80"], budget = 120, rng = Random.Xoshiro(1)).n_calls == 80
        r = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["optim-pso:particles=10,iters=15"], rng = Random.Xoshiro(1))
        @test ismissing(r.converged) && length(r.stages) == 1
    end

    @testset "chains" begin
        r = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["mh-de:N=20,evals=400", "lbfgsb"], rng = Random.Xoshiro(1))
        @test length(r.stages) == 2 && r.stages[1].label == "mh-de:N=20,evals=400" && r.stages[2].label == "lbfgsb"
        @test r.llh ≈ OPTIMUM_LLH atol = 1e-4
        @test r.converged === true
        @test r.seed_params == r.stages[1].params    # the last stage started from the global stage's result
        @test r.time ≈ sum(s.time for s in r.stages)
        @test r.strategy == "mh-de:N=20,evals=400 lbfgsb"
        # settings of a stage are used (a single iteration cannot converge)
        r1 = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["lbfgsb:iterations=1"], rng = Random.Xoshiro(1))
        @test r1.converged === false && r1.n_iters <= 2
        # `local_defaults` fill in settings the stage does not give, an explicit stage setting wins, unknown keys are ignored
        r2 = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["lbfgsb"], rng = Random.Xoshiro(1), local_defaults = (; iterations = 1, nonexistent = 5))
        @test r2.converged === false
        r3 = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["lbfgsb:iterations=1000"], rng = Random.Xoshiro(1), local_defaults = (; iterations = 1))
        @test r3.converged === true
    end

    @testset "seed stages" begin
        # random: a reproducible draw from the prior, then the fit
        a = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["random", "lbfgsb"], rng = Random.Xoshiro(7))
        b = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["random", "lbfgsb"], rng = Random.Xoshiro(7))
        @test a.params == b.params && a.seed_params == b.seed_params && a.seed_params != START
        @test a.llh ≈ OPTIMUM_LLH atol = 1e-4
        # manual: given parameters replace the start, unlisted ones keep their value, a pinned parameter keeps its pinned value
        m = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["manual:a=0.25,b=-0.15", "lbfgsb"], rng = Random.Xoshiro(1))
        @test m.seed_params.a == 0.25 && m.seed_params.b == -0.15 && m.seed_params.c == START.c
        @test m.llh ≈ OPTIMUM_LLH atol = 1e-4
        @test find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["manual:d=5.0", "lbfgsb"], rng = Random.Xoshiro(1)).params.d == 1.0
        @test_throws ErrorException find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["manual:nonexistent=1.0", "lbfgsb"], rng = Random.Xoshiro(1))
        # manual through the Julia API with a full NamedTuple
        api = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = [manual_stage((; a = 0.2, b = 0.1, c = 2.0)), "lbfgsb"], rng = Random.Xoshiro(1))
        @test api.seed_params == (a = 0.2, b = 0.1, c = 2.0, d = 1.0)
        # smseed: a fit with the tested parameter pinned at its null value, then the full fit
        s = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["smseed", "lbfgsb"], null_prior = NULL_PRIOR, rng = Random.Xoshiro(1))
        @test s.seed_params.b == 0.0                       # the tested parameter is still at its null value ...
        @test s.seed_params.a ≈ 0.3 atol = 1e-3            # ... and the nuisance is tuned
        @test s.params.b ≈ -0.2 atol = 1e-3
        @test s.llh ≈ OPTIMUM_LLH atol = 1e-4
        @test_throws ErrorException find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["smseed", "lbfgsb"], rng = Random.Xoshiro(1))   # needs null_prior
        # mle: random start + local fit as the seed of the next stage
        e = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["mle", "lbfgsb:iterations=50"], rng = Random.Xoshiro(2))
        @test length(e.stages) == 3 && e.strategy == "mle lbfgsb:iterations=50"
        @test e.seed_params.a ≈ 0.3 atol = 1e-3
    end

    @testset "several attempts" begin
        multi = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["random", "lbfgsb"], nseeds = 4, rng = Random.Xoshiro(3))
        @test multi.llh ≈ OPTIMUM_LLH atol = 1e-4
        # explicit starts: one attempt per start, the best is kept
        st = [merge(START, (; a = -0.9)), merge(START, (; a = 0.5))]
        r = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["lbfgsb"], starts = st, rng = Random.Xoshiro(1))
        @test r.llh ≈ OPTIMUM_LLH atol = 1e-4
    end

    @testset "strategy arguments, cost and trace" begin
        a = parse_strategy_arg("optim-pso:particles=10 lbfgsb @3")
        @test a.stages == ["optim-pso:particles=10", "lbfgsb"] && a.nseeds == 3 && a.label == "optim-pso:particles=10 lbfgsb @3"
        b = parse_strategy_arg("lbfgsb"; nseeds = 2)
        @test b.stages == ["lbfgsb"] && b.nseeds == 2 && b.label == "lbfgsb @2"
        @test parse_strategy_arg("lbfgsb @1").label == "lbfgsb"
        @test_throws ErrorException parse_strategy_arg("lbfgsb mh-de")        # order is checked
        @test_throws ErrorException parse_strategy_arg("lbfgsb @0")
        # cost: summed over all attempts and stages, the trace runs through all of them
        one = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["mh-de:N=20,evals=100", "lbfgsb"], nseeds = 1, rng = Random.Xoshiro(1))
        three = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["mh-de:N=20,evals=100", "lbfgsb"], nseeds = 3, rng = Random.Xoshiro(1))
        @test one.n_attempts == 1 && three.n_attempts == 3
        @test one.total_calls == sum(s.n_calls for s in one.stages) && one.total_calls > 100
        @test three.total_calls > 2 * one.total_calls && three.total_time >= three.time   # (time() is coarse on Windows: tiny fits can read 0)
        @test length(three.trace_x) == length(three.trace_y) && issorted(three.trace_x)
        @test issorted(three.trace_y; rev = true)                             # best objective so far never gets worse
        @test three.trace_y[end] ≈ -three.log_post atol = 1e-6                # ends at the kept fit's objective value
        @test three.trace_x[end] <= three.trace_span ≈ three.total_calls
        # a global-only strategy: the trace is the global search's
        g = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["mh-de:N=20,evals=100"], rng = Random.Xoshiro(1))
        @test g.total_calls == 100 && length(g.trace_x) == 100 && g.trace_x[end] == 100
        # the stage of every trace point: a chain has both stages, in order; a single stage only stage 1
        @test length(three.trace_stage) == length(three.trace_x) && issorted(three.trace_stage[1:length(one.trace_stage)])
        @test Set(three.trace_stage) == Set([1, 2]) && all(==(1), g.trace_stage)
    end

    @testset "profile_advanced" begin
        ref_llh(b) = -0.5 * ((b + 0.2) / 0.2)^2 - 0.1       # llh of the toy likelihood with b fixed and a, c at their optimum
        vars = OrderedDict(:b => 5)
        r = profile_advanced(ToyLikelihood(), PRIOR, vars, START; strategy = ["lbfgsb"], seed = 5)
        @test r isa NewtrinosResult && size(r.values.llh) == (5,)
        @test r.axes.b == r.values.b
        @test r.values.llh ≈ ref_llh.(r.axes.b) atol = 1e-3
        @test r.meta["task"] == "profile" && r.meta["strategy"] == ["lbfgsb"] && r.meta["nseeds"] == 1
        runs = r.meta["optimizer_runs"]
        @test size(runs["converged"]) == (5,) && all(==(true), runs["converged"]) && all(>(0), runs["n_calls"])
        @test length(runs["trace_x"][1]) == length(runs["trace_y"][1]) == length(runs["trace_stage"][1]) > 0
        @test runs["stages"][1][1].label == "lbfgsb"
        # a hoisted seed stage is computed once and shared; nseeds attempts per point; a global first stage runs per point
        rr = profile_advanced(ToyLikelihood(), PRIOR, vars, START; strategy = ["random", "lbfgsb"], nseeds = 2, seed = 5)
        @test rr.values.llh ≈ ref_llh.(rr.axes.b) atol = 1e-3
        @test all(==(2), rr.meta["optimizer_runs"]["n_attempts"])
        rs = profile_advanced(ToyLikelihood(), PRIOR, vars, START; strategy = ["smseed", "lbfgsb"], null_prior = NULL_PRIOR, seed = 5)
        @test rs.values.llh ≈ ref_llh.(rs.axes.b) atol = 1e-3
        @test_throws ErrorException profile_advanced(ToyLikelihood(), PRIOR, vars, START; strategy = ["random"])   # seed only
        rg = profile_advanced(ToyLikelihood(), PRIOR, vars, START; strategy = ["mh-de:N=20,evals=300", "lbfgsb"], seed = 5)
        @test rg.values.llh ≈ ref_llh.(rg.axes.b) atol = 1e-3
        # same seed: same result (per-point RNG streams, independent of scheduling)
        again = profile_advanced(ToyLikelihood(), PRIOR, vars, START; strategy = ["random", "lbfgsb"], nseeds = 2, seed = 5)
        @test again.values.log_posterior == rr.values.log_posterior
        # warm starting: sequential (BFS from the nearest grid point) and from an earlier result; warmed points run only local stages
        rq = profile_advanced(ToyLikelihood(), PRIOR, vars, START; strategy = ["mh-de:N=20,evals=300", "lbfgsb"], sequential = true, start_from = (b = -0.2,), seed = 5)
        @test rq.meta["sequential"] && rq.meta["warm_strategy"] == ["lbfgsb"]
        @test rq.values.llh ≈ ref_llh.(rq.axes.b) atol = 1e-3
        @test count(s -> length(s) == 2, rq.meta["optimizer_runs"]["stages"]) == 1      # only the root point ran the global stage
        rw = profile_advanced(ToyLikelihood(), PRIOR, vars, START; strategy = ["lbfgsb"], seed_result = r, seed = 5)
        @test rw.values.llh ≈ ref_llh.(rw.axes.b) atol = 1e-3
        # explicit start overrides
        ro = profile_advanced(ToyLikelihood(), PRIOR, vars, START; strategy = ["lbfgsb"], starts = [(a = -0.9,), (a = 0.9,)], seed = 5)
        @test all(==(2), ro.meta["optimizer_runs"]["n_attempts"])
        @test ro.values.llh ≈ ref_llh.(ro.axes.b) atol = 1e-3
        # a 2-D grid has grid-shaped arrays
        r2 = profile_advanced(ToyLikelihood(), PRIOR, OrderedDict(:a => 3, :b => 4), START; strategy = ["lbfgsb"], seed = 5)
        @test size(r2.values.llh) == (3, 4) && size(r2.meta["optimizer_runs"]["time"]) == (3, 4)
        # the result can be saved and loaded again
        mktempdir() do dir
            f = joinpath(dir, "res.jld2")
            FileIO.save(f, "result", r)
            back = FileIO.load(f)["result"]
            @test back isa NewtrinosResult && back.values.llh == r.values.llh && back.meta["optimizer_runs"]["converged"] == runs["converged"]
        end
    end

    @testset "run_optimizer and registry" begin
        for spec in ("lbfgsb", "optim-lbfgs", "optim-neldermead", "optim-pso:particles=10,iters=10", "optim-sa:iters=200", "mh-pso:N=20,evals=200",
                     "mh-abc:N=20,evals=200", "bbo:evals=200", "bbo-xnes:evals=200", "bbo-dxnes:evals=200")
            r = run_optimizer(spec, ToyLikelihood(), PRIOR, START; rng = Random.Xoshiro(1))
            @test isfinite(r.llh) && r.llh <= OPTIMUM_LLH + 1e-9
            @test r.label == spec && r.n_calls > 0
        end
        # a seed or alias stage is not an optimizer
        @test_throws ErrorException run_optimizer("random", ToyLikelihood(), PRIOR, START)
        # registering a new optimizer makes it usable in strategies
        register_optimizer!("toy-grid", :global, (counter, b, cfg, budget, rng) -> (foreach(x -> counter([x, x, x]), range(0.0, 0.9; length = 5)); (; converged = missing, n_iters = 5, native = Float64[])),
                            (; evals = nothing))
        r = find_mle_advanced(ToyLikelihood(), PRIOR, START; strategy = ["toy-grid", "lbfgsb"], rng = Random.Xoshiro(1))
        @test r.llh ≈ OPTIMUM_LLH atol = 1e-4
        @test r.stages[1].n_calls == 5
    end
end

end # module
