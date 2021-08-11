"""
    Molecules.translate!(atoms::Molecule, r::Vector{T})

Translate atoms into position r.
"""
function translate!(atoms::Molecule, r::Vector{T}) where T
    for a in atoms
        for i = 1:3
            a.xyz[i] -= r[i]
        end
    end
end

"""
    Molecules.translate(atoms::Molecule)

Return a copy of atoms translated into position r.
"""
function translate(atoms::Molecule, r::Vector{T}) where T
    newatoms = deepcopy(atoms)
    translate!(newatoms, r)
    return newatoms
end