using Distributions
using Test
using Newtrinos

@testset "Cross Sections" begin

    @testset "SimpleScaling scale function" begin
        xs = Newtrinos.xsec.configure()

        # NC interaction returns nc_norm
        @test xs.scale(:numu, :NC, xs.params) == xs.params.nc_norm
        @test xs.scale(:nue, :NC, xs.params) == xs.params.nc_norm

        # nutau CC returns nutau_cc_norm
        @test xs.scale(:nutau, :CC, xs.params) == xs.params.nutau_cc_norm

        # Other CC returns 1.0 (type-stable with params)
        @test xs.scale(:numu, :CC, xs.params) ≈ 1.0

        # With modified params
        mod_params = (nc_norm=1.5, nutau_cc_norm=0.8)
        @test xs.scale(:nue, :NC, mod_params) == 1.5
        @test xs.scale(:nutau, :CC, mod_params) == 0.8
    end

    @testset "Default params within prior support" begin
        for (name, cfg) in [("SimpleScaling", Newtrinos.xsec.SimpleScaling()),
                            ("Differential_H2O", Newtrinos.xsec.Differential_H2O())]
            xs = Newtrinos.xsec.configure(cfg)
            for key in keys(xs.params)
                @test Distributions.insupport(xs.priors[key], xs.params[key])
            end
        end
    end

    @testset "Differential_H2O scale NC/CC" begin
        xs = Newtrinos.xsec.configure(Newtrinos.xsec.Differential_H2O())
        E = [1.0, 5.0, 10.0]

        # NC returns scalar nc_norm regardless of flavor or anti flag
        @test xs.scale(E, :numu, :NC, false, xs.params) == 1.0
        @test xs.scale(E, :nue, :NC, false, xs.params) == 1.0
        @test xs.scale(E, :nutau, :NC, true, xs.params) == 1.0

        # Modified nc_norm
        mod_params = merge(xs.params, (nc_norm=1.3,))
        @test xs.scale(E, :nue, :NC, false, mod_params) == 1.3

        #CC ratios
        # With all norms=1, ratios sum to 1 so result ≈ 1.0
        result = xs.scale(E, :numu, :CC, false, xs.params)
        @test all(result .≈ 1.0)
        result_anti = xs.scale(E, :numu, :CC, true, xs.params)
        @test all(result_anti .≈ 1.0)

        # Result length matches input
        @test length(result) == length(E)

        # Non-negativity
        @test all(result .>= 0)
    end

    @testset "Differential_H2O nutau CC multiplier" begin
        xs = Newtrinos.xsec.configure(Newtrinos.xsec.Differential_H2O())
        E = [0.5, 1.0, 5.0, 10.0]

        # nutau CC = numu CC * nutau_cc_norm
        mod_params = merge(xs.params, (nutau_cc_norm=0.5,))
        result_nutau = xs.scale(E, :nutau, :CC, false, mod_params)
        result_numu = xs.scale(E, :numu, :CC, false, mod_params)
        @test all(result_nutau .≈ result_numu .* 0.5)

        # Non-nutau flavors unaffected by nutau_cc_norm
        result_nue = xs.scale(E, :nue, :CC, false, mod_params)
        @test all(result_nue .≈ result_numu)
    end

    @testset "Differential_H2O channel weighting" begin
        xs = Newtrinos.xsec.configure(Newtrinos.xsec.Differential_H2O())
        E = [0.5, 1.0, 2.0, 5.0, 10.0]

        # All CC norms = k → result ≈ k
        k = 1.7
        uniform_params = merge(xs.params, (cc1p1h_norm=k, cc2p2h_norm=k, cc1pi_norm=k, ccother_norm=k, ccdis_norm=k))
        result = xs.scale(E, :numu, :CC, false, uniform_params)
        @test all(result .≈ k)

        # Doubling one channel norm → result > 1.0 at higher energies
        mod_params = merge(xs.params, (ccdis_norm=2.0,))
        result = xs.scale([2.0, 5.0, 10.0, 20.0], :numu, :CC, false, mod_params)
        @test all(result .> 1.0)

        # Anti vs non-anti differ with non-uniform norms
        r_nu = xs.scale(E, :numu, :CC, false, mod_params)
        r_anti = xs.scale(E, :numu, :CC, true, mod_params)
        @test !all(r_nu .≈ r_anti)
    end
end
