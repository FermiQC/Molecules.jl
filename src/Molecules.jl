module Molecules

using Unitful
using PeriodicTable
using LinearAlgebra
using Formatting

import PhysicalConstants.CODATA2018: a_0

export Atom, Molecule

# Tolerance for equivalence evaluations
const tol = 1E-5

# Conversion factor 
const bohr_to_angstrom = convert(Float64, a_0 / 1u"Å")
const angstrom_to_bohr = 1.0/bohr_to_angstrom

struct Atom
    Z
    mass
    xyz
end

struct Molecule
    atoms::Vector{Atom}
    charge::Int
    multiplicity::Int
    Nα::Int
    Nβ::Int
end

# Default charge = 0
Molecule(atoms::Vector{T}) where T <: Atom = Molecule(atoms, 0) 
# Default multiplicity = singlet or doublet 
function Molecule(atoms::Vector{T}, charge::Int) where T <: Atom 
    # Compute number of electrons
    nelec = sum(a.Z for a in atoms) - charge
    Molecule(atoms, charge, (nelec % 2) + 1)
end
function Molecule(atoms::Vector{T}, charge::Int, multiplicity::Int) where T <: Atom

    if multiplicity < 1
        throw(DomainError("multiplicity input ($multiplicity)", "Multiplicity value cannot be less than one"))
    end

    # Compute number of electrons
    nelec = -charge
    for i in eachindex(atoms)
        nelec += atoms[i].Z
    end

    # If the number of electrons turns out negative returns an error
    if nelec ≤ 0
        throw(DomainError("negative number of electrons", "Invalid charge ($charge) for given molecule."))
    end

    # Mult = 2Ms + 1 thus the number of unpaired electrons (taken as α) is Mult-1 = 2Ms
    αexcess = multiplicity - 1

    # The number of β electrons must be equal the number of doubly occupied orbitals (Ndo)
    # Ndo = (nelec - n_of_unpaired_elec)/2 this must be an integer
    if isodd(nelec - αexcess)
        throw(DomainError("charge and multiplicity", "Invalid charge ($charge) and multiplicity ($multiplicity) combination."))
    end

    Nβ = (nelec - αexcess) ÷ 2
    Nα = nelec - Nβ

    return Molecule(atoms, charge, multiplicity, Nα, Nβ)
end

# Check if two atoms are the same
Base.:(==)(A1::Atom, A2::Atom) = A1.Z == A2.Z && A1.mass == A2.mass && A1.xyz == A2.xyz

include("Parse.jl")
include("Transformations.jl")
include("Properties.jl")
include("Symmetry/Symmetry.jl")
include("Misc.jl")
end # module
