using DataStructures

@testset "Analysis Tools" begin

    @testset "sort_nt" begin
        nt = (c=3, a=1, b=2)
        sorted = Newtrinos.sort_nt(nt)
        @test keys(sorted) == (:a, :b, :c)
        @test sorted.a == 1
        @test sorted.b == 2
        @test sorted.c == 3

        # Already sorted
        nt2 = (a=1, b=2)
        @test Newtrinos.sort_nt(nt2) == nt2

        # Single element
        nt1 = (x=42,)
        @test Newtrinos.sort_nt(nt1) == nt1
    end

    @testset "safe_merge" begin
        a = (x=1, y=2)
        b = (z=3,)
        merged = Newtrinos.safe_merge(a, b)
        @test haskey(merged, :x)
        @test haskey(merged, :y)
        @test haskey(merged, :z)
        @test merged.x == 1
        @test merged.z == 3

        # Duplicate keys with same values should work
        c = (x=1,)
        merged2 = Newtrinos.safe_merge(a, c)
        @test merged2.x == 1

        # Duplicate keys with different values should error
        d = (x=99,)
        @test_throws ErrorException Newtrinos.safe_merge(a, d)

        # Result should be sorted by key
        e = (z=1, a=2)
        f = (m=3,)
        merged3 = Newtrinos.safe_merge(e, f)
        @test keys(merged3) == (:a, :m, :z)
    end

    @testset "NewtrinosResult" begin
        axes = (x=[1.0, 2.0, 3.0],)
        values = (llh=[-10.0, -5.0, -8.0], log_posterior=[-10.0, -5.0, -8.0])
        result = Newtrinos.NewtrinosResult(axes=axes, values=values)

        @test result.axes == axes
        @test result.values == values
        @test result.meta isa Dict

        # bestfit should find the maximum log_posterior
        bf = Newtrinos.bestfit(result)
        @test bf.log_posterior == -5.0
        @test bf.x == 2.0  # axis value at best fit
    end

    @testset "get_params and get_priors" begin
        osc = Newtrinos.osc.configure()

        # get_params should work with physics objects
        params = Newtrinos.get_params(osc)
        @test params isa NamedTuple
        @test haskey(params, :θ₁₂)

        # get_priors should work with physics objects
        priors = Newtrinos.get_priors(osc)
        @test priors isa NamedTuple
        @test haskey(priors, :θ₁₂)

        # Params and priors should have same keys
        @test keys(params) == keys(priors)
    end
end
