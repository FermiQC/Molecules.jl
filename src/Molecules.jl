module Molecules

using Unitful
using PeriodicTable
using LinearAlgebra
using Formatting

import PhysicalConstants.CODATA2018: a_0

export Atom, Molecule

bohr_to_angstrom = convert(Float64, a_0 / 1u"Å")

struct Atom{F}
    Z
    mass
    xyz::Vector{F}
end

# For checking if two atoms are the same
function Base.:(==)(A1::Atom, A2::Atom)
    A1.Z    == A2.Z    && 
    A1.mass == A2.mass &&
    A1.xyz  == A2.xyz
end

# Alias for a vector of atoms
Molecule = Vector{A} where A <: Atom

# Tolerance for equivalence evaluations
tol = 1E-5

include("Parse.jl")
include("Transformations.jl")
#include("ZMAT.jl")
include("Properties.jl")
# include("RotationalEnergy.jl")
include("Symmetry/Symmetry.jl")
# include("Symmetry/Symmetry.jl")
end # module
