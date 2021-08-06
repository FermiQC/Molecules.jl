module Molecules

using Unitful
using PeriodicTable

struct Atom{T}
    Z::Int16
    A::Int16
    xyz::Vector{T}
end

# Alias for a vector of atoms
Molecule = Vector{A} where A <: Atom


# Include rest of the code
# include("ParseAtoms.jl")
# include("Manipulations.jl")
# include("Properties.jl")
# include("RotationalEnergy.jl")
# include("PointGroup.jl")

end # module
