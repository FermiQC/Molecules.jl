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

Returns the derivative of the nuclear repulsion w.r.t the center i. Units Eₕ/a₀
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
    Molecules.∇2nuclear_repulsion(atoms::Vector{<:Atom}, i, j) where A <: Atom

Returns the second derivative of the nuclear repulsion energy w.r.t. the
centers (atoms) i and j, as a 3x3 matrix in units Eₕ/a₀². 
"""
function ∇2nuclear_repulsion(atoms::Vector{<:Atom}, i, j)
    H = zeros(3, 3)
    Ai = atoms[i]
    Aj = atoms[j]

    if i == j
        for B in atoms
            B == Ai && continue

            r = (Ai.xyz .- B.xyz) .* angstrom_to_bohr
            R = √(r ⋅ r)

            for l in 1:3, k in 1:3
                H[k, l] -= Ai.Z * B.Z * ((k == l ? 1.0/R^3 : 0.0) - 3*r[k]*r[l]/R^5)
            end
        end
    else
        r = (Ai.xyz .- Aj.xyz) .* angstrom_to_bohr
        R = √(r ⋅ r)

        for l in 1:3, k in 1:3
            H[k, l] = Ai.Z * Aj.Z * ((k == l ? 1.0/R^3 : 0.0) - 3*r[k]*r[l]/R^5)
        end
    end

    return H
end
∇2nuclear_repulsion(M::Molecule, i, j) = ∇2nuclear_repulsion(M.atoms, i, j)

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
