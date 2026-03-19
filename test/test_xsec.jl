@testset "Cross Sections" begin

    @testset "SimpleScaling configuration" begin
        xs = Newtrinos.xsec.configure()
        @test xs isa Newtrinos.xsec.Xsec
        @test haskey(xs.params, :nc_norm)
        @test haskey(xs.params, :nutau_cc_norm)
        @test xs.params.nc_norm ≈ 1.0
        @test xs.params.nutau_cc_norm ≈ 1.0
    end

    @testset "SimpleScaling scale function" begin
        xs = Newtrinos.xsec.configure()

        # NC interaction returns nc_norm
        @test xs.scale(:numu, :NC, xs.params) == xs.params.nc_norm
        @test xs.scale(:nue, :NC, xs.params) == xs.params.nc_norm

        # nutau CC returns nutau_cc_norm
        @test xs.scale(:nutau, :CC, xs.params) == xs.params.nutau_cc_norm

        # Other CC returns 1.0 (type-stable with params)
        result = xs.scale(:numu, :CC, xs.params)
        @test result ≈ 1.0
        @test typeof(result) == typeof(xs.params.nc_norm)

        # With modified params
        mod_params = (nc_norm=1.5, nutau_cc_norm=0.8)
        @test xs.scale(:nue, :NC, mod_params) == 1.5
        @test xs.scale(:nutau, :CC, mod_params) == 0.8
    end
end
