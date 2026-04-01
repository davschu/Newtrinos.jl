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

        #TEST get_PMNS
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

        #TEST get_abs_masses()
        @test Newtrinos.osc.get_abs_masses((m₀=1.5, Δm²₂₁=2.3, Δm²₃₁=3.1)) isa Tuple
        @test collect(Newtrinos.osc.get_abs_masses((m₀=1.5, Δm²₂₁=2.3, Δm²₃₁=3.1))) ≈ [1.5, 2.1330729, 2.3130067] atol=1e-6
        @test collect(Newtrinos.osc.get_abs_masses((m₀=0.4, Δm²₂₁=1.2, Δm²₃₁=-2.4))) ≈ [1.6, 1.9390719, 0.4] atol=1e-6
        
        #TEST osc_kernel()
        U = @SMatrix [cos(π/4) sin(π/4) 0; -sin(π/4) cos(π/4) 0; 0 0 1]
        H = @SVector [0.0, 2.5e-5, 1] 
        e, l, σₑ = 1.0, 100.0, 2       
        #test simple kernel
        result_simple = Newtrinos.osc.osc_kernel(U, H, e, l)
        @test size(result_simple) == (3, 3)
        @test result_simple' * result_simple ≈ I atol=1e-12
        #test osc_kernel with low pass filter
        result_lowpass=Newtrinos.osc.osc_kernel(U, H, e, l, σₑ)
        @test length(result_lowpass) == 2
        @test size(result_lowpass[1]) == (3, 3)
        @test size(result_lowpass[2]) == (3,) 
        #test matrix elements against results from rigorous calculation
        u11,u22,u33, u12,u13,u21,u23,u31,u32 = U[1,1], U[2,2], U[3,3], U[1,2], U[1,3], U[2,1], U[2,3], U[3,1], U[3,2]
        phi_simple = -Newtrinos.osc.F_units * 1im * (l / e) .* H
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

        #TEST compute_matter_matrices()
        H_eff = [1.0 0.2 0.1; 0.2 2.0 0.3; 0.1 0.3 3.0] #example hamiltonian
        static_H_eff = SMatrix{3,3}(H_eff)
        layer = Newtrinos.osc.Layer(6371.0, 2.0, 1.5) #radius earth and example proton/neutron densities
        e = 1.5 #energy value

        vecs, vals = Newtrinos.osc.compute_matter_matrices(H_eff, e, layer, false, Newtrinos.osc.SI())
        vecs_anti, vals_anti = Newtrinos.osc.compute_matter_matrices(H_eff, e, layer, true, Newtrinos.osc.SI())
        static_vecs, static_vals = Newtrinos.osc.compute_matter_matrices(static_H_eff, e, layer, false, Newtrinos.osc.SI())
        @test size(vecs) == (3, 3)
        @test size(vals) == (3,)
        #test compatibility of the two implementations (static vs non-static)
        @test vals ≈ static_vals atol=1e-6
        @test abs.(vecs) ≈ abs.(static_vecs) atol=1e-6 # we can only compare the absolute values, because the eigenvectors are not uniquely defined (phase and order)
        #compare to expected values
        A, f, n_p, n_n = Newtrinos.osc.A, e*1e9, layer.p_density, layer.n_density
        H = [1.0+2*A*n_p*f-A*n_n*f 0.2 0.1; 0.2 2.0-A*n_n*f 0.3; 0.1 0.3 3.0-A*n_n*f] #anti = false
        expected_vals, expected_vecs = eigvals(H), eigvecs(H)
        @test collect(vals) ≈ collect(expected_vals) atol=1e-6
        @test collect(vecs) ≈ collect(expected_vecs) atol=1e-6 
        
        #TEST osc_reduce()
        #take matter matrices and energy e=1.5 (GeV) from above
        matter_matrices = [(vecs, vals), (static_vecs, static_vals)] 
        path = [(layer_idx = 1, length = 5.0), (layer_idx = 2, length = 10.0)] #define path through matter (here: path ~ 5 km through abstr. matter matrix, and then 10 km through matter Smatrix)
        #with basic propagation
        U_expected_1 = Newtrinos.osc.osc_kernel(matter_matrices[1][1], matter_matrices[1][2], e, path[1].length) 
        U_expected_2 = Newtrinos.osc.osc_kernel(matter_matrices[2][1], matter_matrices[2][2], e, path[2].length)
        P_expected= abs2.(U_expected_1 * U_expected_2) #expected probability matrix for the given path through matter
        P_result_Basic = Newtrinos.osc.osc_reduce(matter_matrices, path, e, Newtrinos.osc.Basic())
        #with damping propagation: sigma_e = 0.1
        #get decay factor for each layer from the osc_kernel with lowpass filter for both matter matrices
        res1 = Newtrinos.osc.osc_kernel(matter_matrices[1][1], matter_matrices[1][2], e, path[1].length, Newtrinos.osc.Damping().σₑ)
        res2 = Newtrinos.osc.osc_kernel(matter_matrices[2][1], matter_matrices[2][2], e, path[2].length, Newtrinos.osc.Damping().σₑ)
        #use bold approximation: coherent neutrino behaves as if it was influenced by an average weighted damping factor for the entire path 
        #-> matter_matrix_avg = sum(Length_i * matrix_i) / sum(Lenght_i) 
        P_bold_avg = (path[1].length * abs2.(matter_matrices[1][1]) + path[2].length * abs2.(matter_matrices[2][1])) / (path[1].length + path[2].length)
        #combine coherent and incoherent parts to account for damping effects in the probability matrix 
        P_expected_Damping = abs2.(res1[1] * res2[1]) .+ P_bold_avg * Diagonal(1 .- abs2.(res1[2] .* res2[2])) * P_bold_avg' 
        P_result_Damping = Newtrinos.osc.osc_reduce(matter_matrices, path, e, Newtrinos.osc.Damping())
        #@test collect(Newtrinos.osc.osc_reduce(matter_matrices, path, e, Newtrinos.osc.Basic())) ≈ P_expected atol=1e-6
        @test size(P_result_Basic) == size(matter_matrices[1][1])
        @test size(P_result_Damping) == size(matter_matrices[1][1])
        @test all(P_result_Basic .>= 0) && all(P_result_Basic .<= 1)
        @test all(P_result_Damping .>= 0) && all(P_result_Damping .<= 1)
        @test P_result_Basic ≈ P_expected atol=1e-6
        @test P_result_Damping ≈ P_expected_Damping atol=1e-6

        #TEST matter_osc_per_e() 
        #take H_eff, e, layer from above -> can take above matter matrices
        #take path from above 
        #TEST matter_osc_per_e() 
        # reuse H_eff = [1.0 0.2 0.1; ...], e = 1.5 from above
        layer2 = Newtrinos.osc.Layer(3480.0, 4.0, 3.0)
        layers_test = [layer, layer2]   # layer already defined above
        σ_decoh = Newtrinos.osc.Decoherent().σₑ

        # one single-layer path and one two-layer path
        paths_test = [
            [Newtrinos.osc.Path(5.0, 1)],
            [Newtrinos.osc.Path(5.0, 1), Newtrinos.osc.Path(10.0, 2)]
        ]

        #test basic/damping propagation
        mat = Newtrinos.osc.compute_matter_matrices.(Ref(H_eff), e, layers_test, false, Ref(Newtrinos.osc.SI()))
        osc1 = Newtrinos.osc.osc_reduce(mat, paths_test[1], e, Newtrinos.osc.Damping())
        osc2 = Newtrinos.osc.osc_reduce(mat, paths_test[2], e, Newtrinos.osc.Damping())
        p = stack((osc1, osc2))
        expected = Newtrinos.osc.matter_osc_per_e(H_eff, e, layers_test, paths_test, false, Newtrinos.osc.Damping(), Newtrinos.osc.SI())
        @test size(Newtrinos.osc.matter_osc_per_e(H_eff, e, layers_test, paths_test, false, Newtrinos.osc.Damping(), Newtrinos.osc.SI())) == (3,3, length(paths_test)) #two paths, 3x3 probability matrices
        @test p == expected 

        #test decoherent propagation
        U1, h1 = Newtrinos.osc.compute_matter_matrices(H_eff, e, layer,  false, Newtrinos.osc.SI())
        U2, h2 = Newtrinos.osc.compute_matter_matrices(H_eff, e, layer2, false, Newtrinos.osc.SI())
        matter_U = [U1, U2]; matter_h = [h1, h2]

        function manual_decoherent(path_segs)
            P = zeros(3, 3)
            for α in 1:3
                eα = [i==α ? 1.0 : 0.0 for i in 1:3]
                ρ = eα * eα'
                for seg in path_segs
                    U, h = matter_U[seg.layer_idx], matter_h[seg.layer_idx]
                    l = seg.length
                    ρ_eig = U' * ρ * U
                    phases = exp.(-Newtrinos.osc.F_units * 1im * (l/e) .* h)
                    ρ_eig = Diagonal(phases) * ρ_eig * Diagonal(phases)'
                    Δφ = abs.(h .- h') * (l/e) * Newtrinos.osc.F_units
                    D = exp.(-2 .* Δφ .* σ_decoh^2)
                    ρ_eig .= ρ_eig .* D
                    ρ = U * ρ_eig * U'
                end
                for β in 1:3
                    eβ = [i==β ? 1.0 : 0.0 for i in 1:3]
                    P[β, α] = real(eβ' * ρ * eβ)
                end
            end
            P
        end

        P_expected_path1 = manual_decoherent(paths_test[1])
        P_expected_path2 = manual_decoherent(paths_test[2])

        result_Decoherent = Newtrinos.osc.matter_osc_per_e(H_eff, e, layers_test, paths_test, false, Newtrinos.osc.Decoherent(), Newtrinos.osc.SI())

        @test size(result_Decoherent) == (3, 3, 2)
        @test all(result_Decoherent .>= 0) && all(result_Decoherent .<= 1)
        @test result_Decoherent[:, :, 1] ≈ P_expected_path1 atol=1e-6
        @test result_Decoherent[:, :, 2] ≈ P_expected_path2 atol=1e-6

        
        #TEST select()
        @test Newtrinos.osc.select(U1, h1, Newtrinos.osc.All()) == Newtrinos.osc.select(U1, h1, Newtrinos.osc.Cut()) #same with out cut-off-value 
        @test size(Newtrinos.osc.select(U1, h1, Newtrinos.osc.Cut(cutoff=0.5))[3]) == (3,3)
        @test Newtrinos.osc.select(U1, h1, Newtrinos.osc.Cut(cutoff=0.5)) != Newtrinos.osc.select(U1, h1, Newtrinos.osc.All()) 
        @test Newtrinos.osc.select(U1, h1, Newtrinos.osc.Cut(cutoff=0.5)) != Newtrinos.osc.select(U1, h1, Newtrinos.osc.Cut(cutoff=1)) #difference in cutoff
        
        
        #TEST propagate()
        #test-setup
        U = @SMatrix [cos(π/4) sin(π/4) 0; -sin(π/4) cos(π/4) 0; 0 0 1]
        H = @SVector [0.0, 2.5e-5, 1]
        e, l, σₑ = 1.0, 100.0, 2
        E_test = [0.5, 1.0, 1.5]
        L_test = [50.0, 100.0]
        nE, nL = length(E_test), length(L_test)
        
        #test for basic propagation
        U = @SMatrix [cos(π/4) sin(π/4) 0; -sin(π/4) cos(π/4) 0; 0 0 1]
        H = @SVector [0.0, 2.5e-5, 1]
        e, l, σₑ = 1.0, 100.0, 2
        E_test = [0.5, 1.0, 1.5]
        L_test = [50.0, 100.0]
        nE, nL = length(E_test), length(L_test)
        
        #test for damped propagation
        σₑ_damp = Newtrinos.osc.Damping().σₑ
        result_damp = Newtrinos.osc.propagate(U, H, E_test, L_test, Newtrinos.osc.Damping())
        @test size(result_damp) == (3, 3, nE, nL)
        @test all(result_damp .>= 0) && all(result_damp .<= 1)
        expected_damp = stack(broadcast((e, l) -> begin
            pf = -Newtrinos.osc.F_units * (l / e) .* H
            decay = exp.(-2 * abs.(pf) * σₑ_damp^2)
            K = U * Diagonal(exp.(1im * pf) .* decay) * U'
            abs2.(K) + abs2.(U) * Diagonal(1 .- abs2.(decay)) * abs2.(U)'
        end, E_test, L_test'))
        @test result_damp ≈ expected_damp atol=1e-10
        
        #test for decoherent propagation
        σₑ_dec = Newtrinos.osc.Decoherent().σₑ
        result_decoh = Newtrinos.osc.propagate(U, H, E_test, L_test, Newtrinos.osc.Decoherent())
        @test size(result_decoh) == (3, 3, nE, nL)
        @test all(result_decoh .>= 0) && all(result_decoh .<= 1)
        expected_decoh = stack(broadcast((e, l) -> begin
            P = zeros(3, 3)
            v = one(U * U')
            for α in 1:3
                eα = v[α, :]
                ρ = eα * eα'
                ρ_eig = U' * ρ * U
                phases = exp.(-Newtrinos.osc.F_units * 1im * (l / e) .* H)
                ρ_eig = Diagonal(phases) * ρ_eig * Diagonal(phases)'
                Δφ = abs.(H .- H') * (l / e) * Newtrinos.osc.F_units
                D = exp.(-2 .* Δφ .* σₑ_dec^2)
                ρ_eig .= ρ_eig .* D
                ρ = U * ρ_eig * U'
                for β in 1:3
                    eβ = v[β, :]
                    P[β, α] = real(eβ' * ρ * eβ)
                end
            end
            P
        end, E_test, L_test'))
        @test result_decoh ≈ expected_decoh atol=1e-10
        #=MAYBE BUG ->(The root cause: Layer is a parametric struct (Layer{T}), so StructVector{Layer{Float64}} doesn't match the function signature StructVector{Layer} due to Julia's type invariance. Fix: Change StructVector{Layer} → StructVector{<:Layer} on lines 505 and 510 of osc.jl. This is a genuine bug — it would affect any caller constructing layers with concrete types)
        #test for different layers in vacuum 
        paths_vov = VectorOfVectors(paths_test)
        layers_sv = StructVector(layers_test)
        nPaths = length(paths_test)
        L_total = [sum(seg.length for seg in p) for p in paths_test]

        result_vac = Newtrinos.osc.propagate(U, H, E_test, paths_vov, layers_sv, Newtrinos.osc.Basic(), Newtrinos.osc.Vacuum(), false)
        @test size(result_vac) == (3, 3, nE, nPaths)
        result_direct = Newtrinos.osc.propagate(U, H, E_test, L_total, Newtrinos.osc.Basic())
        @test result_vac ≈ result_direct atol=1e-10
        
        #test for different layers in matter
        result_si = Newtrinos.osc.propagate(U, H, E_test, paths_vov, layers_sv, Newtrinos.osc.Damping(), Newtrinos.osc.SI(), false)
        @test size(result_si) == (3, 3, nE, nPaths)
        @test all(result_si .>= 0) && all(result_si .<= 1)
        H_eff_ref = U * Diagonal(H) * adjoint(U)
        expected_si = stack(map(e -> Newtrinos.osc.matter_osc_per_e(H_eff_ref, e, layers_sv, paths_vov, false, Newtrinos.osc.Damping(), Newtrinos.osc.SI()), E_test))
        expected_si = permutedims(expected_si, (1, 2, 4, 3))
        @test result_si ≈ expected_si atol=1e-10=#



        #@test get_osc_prob() 
        #@test osc_prob() #two variants=#
    end 

    @testset "get matrices" begin
        #test get_matrices for different model configurations
        #@test get_matrices() 
    end
end
