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
