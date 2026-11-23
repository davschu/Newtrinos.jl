# Tests of the truth selection of `sensitivity.jl --truth` (`truth_from_spec`, `truth_grid(ranges = ...)` in src/analysis/optimizer_tests_common.jl)
# and of the data arrangement and figures in src/analysis/sensitivity_plots.jl, on mock priors and mock fisher results -- no experiment data needed.

module SensitivityTests

using Test
using Random
using LinearAlgebra
using Distributions
using DensityInterface
using BAT
using MeasureBase
using ValueShapes
using ADTypes
using Newtrinos
using Accessors
using DataStructures
using FileIO
using CairoMakie
import Optim
import Sobol

include(joinpath(@__DIR__, "..", "src", "analysis", "optimizer_tests_common.jl"))
include(joinpath(@__DIR__, "..", "src", "analysis", "sensitivity_plots.jl"))

const PRIOR = (a = Uniform(-1.0, 1.0), b = Uniform(-1.0, 1.0), c = Uniform(0.0, 3.0), d = 1.0)   # d is held fixed
const NOMINAL = (a = 0.0, b = 0.5, c = 1.5, d = 7.0)   # d differs from its pinned value on purpose

@testset "truth selection" begin

    @testset "nominal and grid" begin
        t = truth_from_spec("nominal", PRIOR, NOMINAL, (:a,))
        @test t.kind == :nominal && length(t.truths) == 1 && t.truths[1].d == 1.0 && t.truths[1].a == 0.0   # the pinned parameter is forced
        g = truth_from_spec("grid", PRIOR, NOMINAL, (:a,); n_points = 5)
        @test g.kind == :grid && [x.a for x in g.truths] ≈ [-1, -0.5, 0, 0.5, 1] && g.scan_keys == (:a,)
        gi = truth_from_spec("grid", PRIOR, NOMINAL, (:a,); n_points = 4, interior = true)
        @test [x.a for x in gi.truths] ≈ [-0.75, -0.25, 0.25, 0.75]
    end

    @testset "range" begin
        r = truth_from_spec("range:0.2,0.6", PRIOR, NOMINAL, (:a,); n_points = 5)
        @test r.kind == :range && [x.a for x in r.truths] ≈ range(0.2, 0.6; length = 5) && all(x -> x.b == 0.5 && x.d == 1.0, r.truths)
        @test occursin("a in [0.2, 0.6]", r.description)
        # two tested parameters: only the listed one is restricted, a 2D grid with round(sqrt(n)) points per dimension
        r2 = truth_from_spec("range:a=0.1,0.3", PRIOR, NOMINAL, (:a, :b); n_points = 9)
        @test length(r2.truths) == 9 && all(x -> 0.1 <= x.a <= 0.3, r2.truths) && extrema(x.b for x in r2.truths) == (-1.0, 1.0)
        # three tested parameters: Sobol points inside the given intervals
        r3 = truth_from_spec("range:a=0,0.5;c=1,2", PRIOR, NOMINAL, (:a, :b, :c); n_points = 16)
        @test length(r3.truths) == 16 && all(x -> 0 <= x.a <= 0.5 && 1 <= x.c <= 2 && -1 <= x.b <= 1, r3.truths)
        # cell-centred values stay off the interval bounds
        ri = truth_from_spec("range:0,1", PRIOR, NOMINAL, (:a,); n_points = 4, interior = true)
        @test [x.a for x in ri.truths] ≈ [0.125, 0.375, 0.625, 0.875]
        # truth_grid itself
        tg = first(truth_grid(PRIOR, NOMINAL, (:a,); n_points = 3, ranges = Dict(:a => (0.0, 0.5))))
        @test [x.a for x in tg] ≈ [0, 0.25, 0.5]
    end

    @testset "point" begin
        p = truth_from_spec("point:a=0.3, c=2.0", PRIOR, NOMINAL, (:a,))
        @test p.kind == :point && length(p.truths) == 1
        t = p.truths[1]
        @test t.a == 0.3 && t.c == 2.0 && t.b == 0.5 && t.d == 1.0
        @test occursin("a = 0.3", p.description)
        @test truth_from_spec("point:b=-0.5", PRIOR, NOMINAL, (:a,)).truths[1].b == -0.5   # a nuisance
    end

    @testset "best fit of a result file" begin
        mktempdir() do dir
            res = NewtrinosResult(axes = (a = [-0.5, 0.0, 0.5],),
                                  values = (log_posterior = [-3.0, -1.0, -2.0], llh = [-3.0, -1.0, -2.0], b = [0.1, 0.2, 0.3], c = [1.0, 1.1, 1.2]), meta = Dict())
            f = joinpath(dir, "res.jld2")
            FileIO.save(f, "result", res)
            bf = @test_logs (:warn, r"does not contain d") truth_from_spec("bf:$f", PRIOR, NOMINAL, (:a,))   # d is not in the result
            @test bf.kind == :bf && length(bf.truths) == 1
            t = bf.truths[1]
            @test t.a == 0.0 && t.b == 0.2 && t.c == 1.1 && t.d == 1.0
            g = joinpath(dir, "other.jld2")
            FileIO.save(g, "x", 1)
            @test_throws "has no key \"result\"" truth_from_spec("bf:$g", PRIOR, NOMINAL, (:a,))
            @test_throws "not found" truth_from_spec("bf:$(joinpath(dir, "nope.jld2"))", PRIOR, NOMINAL, (:a,))
        end
    end

    @testset "errors" begin
        err(spec, tested = (:a,)) = truth_from_spec(spec, PRIOR, NOMINAL, tested)
        @test_throws "must be nominal" err("everywhere")
        @test_throws "takes no argument" err("grid:3")
        @test_throws "needs an argument" err("point:")
        @test_throws "outside its prior range" err("range:0.5,1.5")
        @test_throws "lo < hi" err("range:0.5,0.2")
        @test_throws "expected lo,hi" err("range:0.5")
        @test_throws "exactly one free tested parameter" err("range:0,1", (:a, :b))
        @test_throws "not a free tested parameter" err("range:c=0,1")
        @test_throws "given twice" err("range:a=0,1;a=0,0.5")
        @test_throws "outside its prior range" err("point:a=2")
        @test_throws "is not a parameter" err("point:zz=1")
        @test_throws "fixed in the prior" err("point:d=2")
        @test_throws "finite number" err("point:a=abc")
        @test_throws "expected name=value" err("point:a")
        @test_throws "given twice" err("point:a=0.1,a=0.2")
    end
end

# mock of the `fisher` results of sensitivity.jl: nt tested parameters, n_nuis nuisances, K truths (each with a random valid correlation matrix)
function mock_results(; nt = 1, n_nuis = 6, K = 5, failed = Int[], rng = Random.Xoshiro(3))
    tested = [Symbol("t", i) for i in 1:nt]
    names = vcat(tested, [Symbol("n", i) for i in 1:n_nuis])
    n = length(names)
    truths = [NamedTuple{Tuple(tested)}(Tuple(randn(rng, nt))) for _ in 1:K]
    fisher = Any[]
    for k in 1:K
        k in failed && (push!(fisher, nothing); continue)
        A = randn(rng, n, 3); C = A * A' + I
        σ = exp.(randn(rng, n) .- 1)
        cov = C .* (σ * σ') ./ (sqrt.(diag(C)) * sqrt.(diag(C))')
        σ_marg = sqrt.(diag(cov)); corr = cov ./ (σ_marg * σ_marg')
        push!(fisher, (; names, cov, corr, σ_marg, σ_cond = 0.5 .* σ_marg, width = fill(2.0, n), impact = zeros(nt, n), σ_fixed = zeros(nt, n),
                       tested_idx = collect(1:nt), offset = 0.0, dropped = Symbol[]))
    end
    Dict("settings" => Dict("setup" => "mock", "hypothesis" => "h", "tested" => string.(tested), "information" => "expected", "truth_description" => "mock"),
         "truths" => truths, "fisher" => fisher)
end

@testset "sensitivity plots" begin
    @testset "data arrangement" begin
        res = mock_results(nt = 1, n_nuis = 6, K = 5, failed = [3])
        D = sensitivity_data(res)
        @test D.failed == [3] && D.truth_index == [1, 2, 4, 5] && size(D.tv) == (1, 4)
        @test D.names == ["t1", "n1", "n2", "n3", "n4", "n5", "n6"] && D.ti == [1] && D.tested == [:t1]
        @test size(D.corr) == (7, 7, 4) && size(D.smarg) == (7, 4) && all(D.width .== 2.0)
        @test all(diag(D.corr[:, :, 1]) .≈ 1)
        sel = shown_nuisances(D; top_n = 3)
        @test length(sel) == 3 && !(1 in sel)
        score(j) = maximum(abs, D.corr[1, j, :])
        @test issorted(score.(sel); rev = true)
        @test selection_note(D, sel) != "" && selection_note(D, shown_nuisances(D; top_n = 10)) == ""
        @test truth_label(D, 1) == "t1 = $(round(D.tv[1, 1]; sigdigits = 3))"
        @test_throws ErrorException sensitivity_data(mock_results(K = 2, failed = [1, 2]))
    end

    @testset "figures" begin
        mktempdir() do dir
            # one tested parameter, many nuisances and truths: only the top_n nuisances are drawn
            D = sensitivity_data(mock_results(nt = 1, n_nuis = 40, K = 50, failed = [7]))
            fig = static_overview(D; top_n = 10)
            @test fig isa Figure
            save(joinpath(dir, "one.pdf"), fig)
            @test filesize(joinpath(dir, "one.pdf")) > 0
            # several tested parameters (small multiples)
            D3 = sensitivity_data(mock_results(nt = 6, n_nuis = 8, K = 12))
            fig3 = static_overview(D3; top_n = 5)
            save(joinpath(dir, "six.pdf"), fig3)
            @test filesize(joinpath(dir, "six.pdf")) > 0
            # a single truth: the actual numbers (signed correlation matrix, impact bars, tested parameter among the sigma bars), not ranges over truths
            D1 = sensitivity_data(mock_results(nt = 1, n_nuis = 4, K = 1))
            fig1 = static_overview(D1)
            save(joinpath(dir, "single.pdf"), fig1)
            @test filesize(joinpath(dir, "single.pdf")) > 0
            axes1 = [a for a in fig1.content if a isa Axis]
            titles1 = [string(a.title[]) for a in axes1]
            @test "correlation matrix at the truth" in titles1
            impact_ax = axes1[findfirst(t -> startswith(t, "impact"), titles1)]
            @test any(p -> p isa Makie.BarPlot, impact_ax.scene.plots)   # not empty (a range bar of zero length is invisible)
            sigma_ax = axes1[findfirst(t -> startswith(t, "σ / prior width"), titles1)]
            @test any(l -> occursin("(tested)", l), sigma_ax.yticks[][2])
            # several tested parameters at one truth
            D3 = sensitivity_data(mock_results(nt = 3, n_nuis = 6, K = 1))
            save(joinpath(dir, "single3.pdf"), static_overview(D3))
            @test filesize(joinpath(dir, "single3.pdf")) > 0
        end
    end
end

end # module
