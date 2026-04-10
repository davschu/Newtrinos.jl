using Newtrinos
using Test

@testset "Newtrinos.jl" begin
    include("test_helpers.jl")
    include("test_osc.jl")
    include("test_earth_layers.jl")
    include("test_xsec.jl")
    include("test_analysis.jl")
    include("test_autodiff.jl")
    include("test_regression.jl")
    include("osc_unit_tests.jl")
    include("test_snsflux.jl")
end
