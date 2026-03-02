using Newtrinos
using Test
using Distributions
using LinearAlgebra
using StaticArrays

@testset "osc.jl" begin
    @testset "params and priors" begin 
    #test get_params and get_priors for different model configurations
        configs= [
            Newtrinos.osc.ThreeFlavour(), 
            Newtrinos.osc.ThreeFlavourXYCP(), 
            Newtrinos.osc.Sterile(), 
            Newtrinos.osc.ADD()]#=, 
            Newtrinos.osc.Darkdim_Lambda(), 
            Newtrinos.osc.Darkdim_Masses(), 
            Newtrinos.osc.Darkdim_cas()]=#

        for cfg in configs
            @test Newtrinos.osc.get_params(cfg) isa NamedTuple
            @test Newtrinos.osc.get_priors(cfg) isa NamedTuple
            @test keys(Newtrinos.osc.get_params(cfg)) == keys(Newtrinos.osc.get_priors(cfg))
        end
        
        #test if the parameters are correctly calculated
        #eventually: test against PDG live values?
        #basic parameters
        angles=(θ₁₂ = 0.58725, θ₁₃ = 0.14543, θ₂₃ = 0.85563)
        NO_masses=(Δm²₂₁ = 7.53e-5, Δm²₃₁ = 2.4e-3 + 7.53e-5)
        IO_masses=(Δm²₂₁ = 7.53e-5, Δm²₃₁ = -(2.4e-3 - 7.53e-5))
        δCP = (δCP = 1.,)
        #define expected parameter sets
        ThreeFlavour_Params_NO=merge(angles, NO_masses, δCP)
        ThreeFlavour_Params_IO=merge(angles, IO_masses, δCP)
        #these are all for NO because they they dont affect m31
        ThreeFlavourXYCP_Params=merge(angles, NO_masses, (δCPshell = [1., 0.],))
        Sterile_Params = merge(ThreeFlavour_Params_NO, (Δm²₄₁ = 1., θ₁₄ = 0.1, θ₂₄ = 0.1, θ₃₄ = 0.1,))
        ADD_Params = merge(ThreeFlavour_Params_NO, (m₀ = 0.01, ADD_radius =1e-2,))
        Darkdim_Lambda_Params = merge(angles, δCP, (Darkdim_radius=0.1, ca1=1e-5, ca2=1e-5, ca3=1e-5, λ₁= 1., λ₂ = 1., λ₃ = 1.,))
        Darkdim_Masses_Params = merge(ThreeFlavour_Params_NO, (m₀ = 0.01, Darkdim_radius=0.1,  λ₁= 1., λ₂ = 1., λ₃ = 1.,))
        Darkdim_cas_Params = merge(ThreeFlavour_Params_NO, (m₀ = 0.01, Darkdim_radius=0.1, ca1=1e-5, ca2=1e-5, ca3=1e-5,))

        #tests
        #isapprox cannot handle NamedTuples -> loop over keys and test each parameter separately, with an appropriate tolerance
        for key in keys(ThreeFlavour_Params_NO)
            @test getfield(Newtrinos.osc.get_params(Newtrinos.osc.ThreeFlavour(ordering=:NO)), key) ≈ getfield(ThreeFlavour_Params_NO, key) atol=1e-4
        end
        for key in keys(ThreeFlavour_Params_IO)
            @test getfield(Newtrinos.osc.get_params(Newtrinos.osc.ThreeFlavour(ordering=:IO)), key) ≈ getfield(ThreeFlavour_Params_IO, key) atol=1e-4
        end
        for key in keys(ThreeFlavourXYCP_Params)
            @test getfield(Newtrinos.osc.get_params(Newtrinos.osc.ThreeFlavourXYCP()), key) ≈ getfield(ThreeFlavourXYCP_Params, key) atol=1e-4
        end
        for key in keys(Sterile_Params)
            @test getfield(Newtrinos.osc.get_params(Newtrinos.osc.Sterile()), key) ≈ getfield(Sterile_Params, key) atol=1e-4
        end
        for key in keys(ADD_Params)
            @test getfield(Newtrinos.osc.get_params(Newtrinos.osc.ADD()), key) ≈ getfield(ADD_Params, key) atol=1e-4
        end
        #darkdim not yet exported from osc file
        #=for key in keys(Darkdim_Lambda_Params)
            @test getfield(Newtrinos.osc.get_params(Newtrinos.osc.Darkdim_Lambda()), key) ≈ getfield(Darkdim_Lambda_Params, key) atol=1e-4
        end
        for key in keys(Darkdim_Masses_Params)
            @test getfield(Newtrinos.osc.get_params(Newtrinos.osc.Darkdim_Masses()), key) ≈ getfield(Darkdim_Masses_Params, key) atol=1e-4
        end
        for key in keys(Darkdim_cas_Params)
            @test getfield(Newtrinos.osc.get_params(Newtrinos.osc.Darkdim_cas()), key) ≈ getfield(Darkdim_cas_Params, key) atol=1e-4
        end=#
        
        
        #test if the priors are correctly calculated
        
        priors_ThreeFlavour_NO = Newtrinos.osc.get_priors(Newtrinos.osc.ThreeFlavour(ordering=:NO))
        priors_ThreeFlavour_IO = Newtrinos.osc.get_priors(Newtrinos.osc.ThreeFlavour(ordering=:IO))
        priors_ThreeFlavourXYCP = Newtrinos.osc.get_priors(Newtrinos.osc.ThreeFlavourXYCP())
        priors_Sterile = Newtrinos.osc.get_priors(Newtrinos.osc.Sterile())
        priors_ADD = Newtrinos.osc.get_priors(Newtrinos.osc.ADD())
        #priors_Darkdim_Lambda = Newtrinos.osc.get_priors(Newtrinos.osc.Darkdim_Lambda())
        #priors_Darkdim_Masses = Newtrinos.osc.get_priors(Newtrinos.osc.Darkdim_Masses())
        #priors_Darkdim_cas = Newtrinos.osc.get_priors(Newtrinos.osc.Darkdim_cas())

        #ThreeFlavour
        #TODO: test if param is contained in prior, with @test insupport(prior_distr,param)
        for key in keys(priors_ThreeFlavour_NO)
            @test getfield(priors_ThreeFlavour_NO, key) isa Uniform
        end
        @test getfield(priors_ThreeFlavour_NO, :Δm²₃₁) != getfield(priors_ThreeFlavour_IO, :Δm²₃₁)
        #ThreeFlavourXYCP
        @test getfield(priors_ThreeFlavourXYCP, :δCPshell) isa MvNormal
        #Sterile
        for key in [:θ₁₄, :θ₂₄, :θ₃₄, :Δm²₄₁]
            @test getfield(priors_Sterile, key) isa Uniform
        end
        #ADD
        @test getfield(priors_ADD, :m₀ ) isa LogUniform
        @test getfield(priors_ADD, :ADD_radius) isa LogUniform
        #Darkdim_Lambda
        #=@test getfield(priors_Darkdim_Lambda, :Darkdim_radius) isa LogUniform
        @test !haskey(priors_Darkdim_Lambda, :Δm²₂₁)
        @test !haskey(priors_Darkdim_Lambda, :Δm²₃₁)
        for key in [:ca1, :ca2, :ca3, :λ₁, :λ₂, :λ₃]
            @test getfield(priors_Darkdim_Lambda, key) isa Uniform
        end
        #Darkdim_Masses
        @test getfield(priors_Darkdim_Masses, :m₀) isa LogUniform
        @test getfield(priors_Darkdim_Masses, :Darkdim_radius) isa LogUniform
        for key in [:ca1, :ca2, :ca3]
            @test getfield(priors_Darkdim_Masses, key) isa Uniform
        end
        #Darkdim_cas
        @test getfield(priors_Darkdim_cas, :m₀) isa LogUniform
        @test getfield(priors_Darkdim_cas, :Darkdim_radius) isa LogUniform
        for key in [:λ₁, :λ₂, :λ₃]
            @test getfield(priors_Darkdim_cas, key) isa Uniform
        end=#
    end

    @testset "oscillation functions" begin
        #define test parameter sets
        test_params_1 = (θ₁₂ = 0.0, θ₁₃ = 0.0, θ₂₃ = 0.0, δCP = 0.0)
        test_params_2 = (θ₁₂ = 0.6, θ₁₃ = 0.2, θ₂₃ = 0.7, δCP = 0.4)
        test_params_3 = (θ₁₂ = 0.4, θ₁₃ = 0.15, θ₂₃ = 0.9, δCP = -0.5)

        #test get_PMNS
        for params in [test_params_1,test_params_2, test_params_3]
            s12, c12 = sin(params.θ₁₂), cos(params.θ₁₂)
            s13, c13 = sin(params.θ₁₃), cos(params.θ₁₃)
            s23, c23 = sin(params.θ₂₃), cos(params.θ₂₃)
            cp = cis(params.δCP) # exp(i * δ)
            cp_conj = conj(cp)   # exp(-i * δ)
            #PMNS standard parametrization
            expected_PMNS = [
            c13*c12                     c13*s12                     s13*cp_conj
            -c23*s12-s23*c12*s13*cp     c23*c12-s23*s12*s13*cp      s23*c13
            s23*s12-c23*s13*c12*cp      -s23*c12-s12*s13*c23*cp     c23*c13   
            ] 
            @test Newtrinos.osc.get_PMNS(params) ≈ expected_PMNS atol=1e-6
            @test Newtrinos.osc.get_PMNS(params) isa SMatrix{3,3}
        end

        #test get_abs_masses()
        @test Newtrinos.osc.get_abs_masses((m₀=1.5, Δm²₂₁=2.3, Δm²₃₁=3.1)) isa Tuple
        @test collect(Newtrinos.osc.get_abs_masses((m₀=1.5, Δm²₂₁=2.3, Δm²₃₁=3.1))) ≈ [1.5, 2.1330729, 2.3130067] atol=1e-6
        @test collect(Newtrinos.osc.get_abs_masses((m₀=0.4, Δm²₂₁=1.2, Δm²₃₁=-2.4))) ≈ [1.6, 1.9390719, 0.4] atol=1e-6
        
        #test osc_kernel()
        U = @SMatrix [cos(π/4) sin(π/4) 0; -sin(π/4) cos(π/4) 0; 0 0 1]
        H = @SVector [0.0, 2.5e-5, 1] 
        e, l, σₑ = 1.0, 100.0, 2       
        #test simple kernel
        result_simple = Newtrinos.osc.osc_kernel(U, H, e, l)
        @test size(result_simple) == (3, 3)
        @test result_simple' * result_simple ≈ I atol=1e-12
        #@test check result_simple values
        #test osc_kernel with low pass filter
        result_lowpass=Newtrinos.osc.osc_kernel(U, H, e, l, σₑ)
        @test length(result_lowpass) == 2
        @test size(result_lowpass[1]) == (3, 3)
        @test size(result_lowpass[2]) == (3,)
        @test result_lowpass' * result_lowpass ≈ I atol=1e-12 
        #test matrix elements against results from rigorous calculation
        u11,u22,u33, u12,u13,u21,u23,u31,u32 = U[1,1], U[2,2], U[3,3], U[1,2], U[1,3], U[2,1], U[2,3], U[3,1], U[3,2]
        phi_simple = -F_units * 1im * (l / e) .* H
        phi_decay = - 2 * abs.(-1im*phi_simple) * σₑ^2
        phi_lowpass = phi_simple + phi_decay
        
        for phi in [phi_simple, phi_lowpass] 
            expected_kernel_matrix= [
                u11^2*exp(phi[1])+u12^2*exp(phi[2])+u13^2*exp(phi[3])           u11*u21*exp(phi[1])+u12*u22*exp(phi[2])+u13*u23*exp(phi[3])     u11*u31*exp(phi[1])+u12*u32*exp(phi[2])+u13*u33*exp(phi[3])
                u21*u11*exp(phi[1])+u22*u12*exp(phi[2])+u23*u13*exp(phi[3])     u21^2*exp(phi[1])+u22^2*exp(phi[2])+u23^2*exp(phi[3])           u21*u31*exp(phi[1])+u22*u32*exp(phi[2])+u23*u33*exp(phi[3])
                u31*u11*exp(phi[1])+u32*u12*exp(phi[2])+u33*u13*exp(phi[3])     u31*u21*exp(phi[1])+u32*u22*exp(phi[2])+u33*u23*exp(phi[3])     u31^2*exp(phi[1])+u32^2*exp(phi[2])+u33^2*exp(phi[3])]
            if phi == phi_simple
                @test collect(Newtrinos.osc.osc_kernel(U, H, e, l)) ≈ expected_kernel_matrix atol=5e-6
            elseif phi == phi_lowpass
                @test collect(Newtrinos.osc.osc_kernel(U, H, e, l, σₑ)[1]) ≈ expected_kernel_matrix atol=5e-6
                @test Newtrinos.osc.osc_kernel(U, H, e, l, σₑ)[2] ≈ exp.(phi_decay) atol=5e-6
            end
        end

        #test compute_matter_matrices()
        #TODO: understand first and then test!
        
        #=@test compute_matter_matrices()#two variants
        @test osc_reduce() #two variants
        @test matter_osc_per_e() #two variants
        @test select() #two variants
        @test propagate() #five variants
        @test get_osc_prob() 
        @test osc_prob() #two variants=#
    end 

    @testset "get matrices" begin
        #test get_matrices for different model configurations
        #@test get_matrices() 
    end
end
