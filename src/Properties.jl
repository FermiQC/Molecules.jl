"""
    Molecules.symbol(A::Atom)

Returns a String representing the atomic symbol of the atom.
"""
function symbol(A::Atom)
    return elements[A.Z].symbol
end

"""
    Molecules.nuclear_repulsion(A::Atom, B::Atom)

Returns nuclear respulsion energy between atoms A and B in atomic units.
"""
function nuclear_repulsion(A::Atom, B::Atom)
    return (A.Z*B.Z)/(√((A.xyz.-B.xyz)⋅(A.xyz.-B.xyz))/bohr_to_angstrom)
end

"""
    Molecules.nuclear_repulsion(atoms::Vector{<:Atom}) where A <: Atom

Returns the total nuclear respulsion energy of a group of atoms. 
""" 
function nuclear_repulsion(atoms::Vector{<:Atom})
    N = length(atoms)
    E = 0.0
    for i in 1:N
        for j in (i+1):N
            E += nuclear_repulsion(atoms[i], atoms[j])
        end
    end
    return E
end
nuclear_repulsion(M::Molecule) = nuclear_repulsion(M.atoms)

"""
    Molecules.∇nuclear_repulsion(atoms::Vector{<:Atom}, i) where A <: Atom

Returns the derivative of the nuclear repulsion w.r.t the center i. Units Eₕ/Å
""" 
function ∇nuclear_repulsion(atoms::Vector{<:Atom}, i)
    E = zeros(3)
    At = atoms[i]

    for B in atoms
        B == At ? continue : nothing

        D = √((At.xyz.-B.xyz)⋅(At.xyz.-B.xyz)) * angstrom_to_bohr

        E -= (B.Z / D^3) .* (At.xyz - B.xyz) * angstrom_to_bohr
    end

    return At.Z .* E 
end
∇nuclear_repulsion(M::Molecule, i) = ∇nuclear_repulsion(M.atoms, i)

"""
    Molecules.center_of_mass(atoms::Vector{<:Atom}, i) where A <: Atom

Returns the center of mass of a group if atoms.
""" 
function center_of_mass(atoms::Vector{<:Atom})

    # Compute the mass-weighted XYZ
    cm_xyz = sum(a.xyz * a.mass for a in atoms)

    # Compute total mass
    M = sum(a.mass for a in atoms)

    # Return center of mass
    return cm_xyz ./ M
end
center_of_mass(M::Molecule) = center_of_mass(M.atoms)

"""
    Molecules.nuclear_dipole(atoms::Molecule, o = (0.0, 0.0, 0.0))

Returns the nuclear dipole moment w.r.t the origin o.
"""
function nuclear_dipole(atoms::Vector{<:Atom}, o = [0.0, 0.0, 0.0]) 
    charges = [a.Z for a in atoms]
    r = [a.xyz .- o for a in atoms]
    return sum(r .* charges)
end
nuclear_dipole(M::Molecule, o = [0.0, 0.0, 0.0]) = nuclear_dipole(M.atoms, o)