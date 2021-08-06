module Molecules

using Unitful
using PeriodicTable

import PhysicalConstants.CODATA2018: a_0

bohr_to_angstrom = convert(Float64, a_0 / 1u"Å")

struct Atom{I,F}
    Z::I
    A::F
    xyz::Vector{F}
end

# Alias for a vector of atoms
Molecule = Vector{A} where A <: Atom

include("Parse.jl")
# include("Manipulations.jl")
# include("Properties.jl")
# include("RotationalEnergy.jl")
# include("PointGroup.jl")

end # module
