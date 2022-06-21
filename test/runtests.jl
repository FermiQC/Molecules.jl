using Molecules
using Test
using StaticArrays

@testset "Molecules" begin
    include("test_parsing.jl")
    include("test_transformations.jl")
    include("test_properties.jl")
    include("test_sym.jl")
    include("test_point_groups.jl")
    include("test_character_tables.jl")
end