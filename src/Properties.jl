"""
    Molecules.symbol(A::Atom)

Returns a String representing the atomic symbol of the atom.
"""
function symbol(A::Atom)
    return elements[A.Z].symbol
end

"""
    Molecules.nuclear_repulsion(A::Atom, B::Atom)

Returns the repulsion energy between atoms A and B in Atomic units.
"""
function nuclear_repulsion(A::Atom, B::Atom)
    return (A.Z*B.Z)/(√((A.xyz.-B.xyz)⋅(A.xyz.-B.xyz))/bohr_to_angstrom)
end

"""
    Molecules.nuclear_repulsion(atoms::Molecule)

Returns the total nuclear respulsion energy of a group of atoms. 
"""
function nuclear_repulsion(atoms::Molecule)
    N = length(atoms)
    E = 0.0
    for i in 1:N
        for j in (i+1):N
            E += nuclear_repulsion(atoms[i], atoms[j])
        end
    end
    return E
end

function ∇nuclear_repulsion(atoms::Molecule, i)
    E = zeros(3)

    A = atoms[i]

    for B in atoms
        B == A ? continue : nothing

        D = √((A.xyz.-B.xyz)⋅(A.xyz.-B.xyz))/bohr_to_angstrom

        E -= (B.Z / D^3) .* (A.xyz - B.xyz)/bohr_to_angstrom
    end

    return A.Z .* E 
end

"""
    Molecules.center_of_mass(atoms::Molecule)

Returns the center of mass of a group if atoms.
"""
function center_of_mass(atoms::Molecule)

    # Compute the mass-weighted XYZ
    cm_xyz = sum(a.xyz * a.mass for a in atoms)

    # Compute total mass
    M = sum(a.mass for a in atoms)

    # Return center of mass
    return cm_xyz ./ M
end


"""
    Molecules.nuclear_dipole(atoms::Molecule, o = (0.0, 0.0, 0.0))

Returns the nuclear dipole moment
"""
function nuclear_dipole(atoms::Molecule, o = [0.0, 0.0, 0.0])

    charges = [a.Z for a in atoms]
    r = [a.xyz .- o for a in atoms]
    return sum(r .* charges)
end