using SafeTestsets
using Test

#! format: off
@safetestset "Aqua" begin @time include("aqua.jl") end
@safetestset "Utilities" begin @time include("utilities.jl") end
@safetestset "Parameter tests" begin @time include("parameter_tests.jl") end
@safetestset "Coupler Compatibility" begin @time include("coupler_compatibility.jl") end
#! format: on

nothing
