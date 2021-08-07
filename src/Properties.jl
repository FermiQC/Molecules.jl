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