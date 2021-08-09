export centerXYZ, calcCOM

function calcCOM(Molecule)
    len = size(Molecule)[1]
    suum = zeros(Float64, (3,))
    mass = 0.0
    for i = 1:len
        mass += Molecule[i].mass
        suum += Molecule[i].mass * Molecule[i].xyz
    end
    suum /= mass
    return suum
end

function centerXYZ(Molecule)
    com = calcCOM(Molecule)
    len = size(Molecule)[1]
    out = []
    for i = 1:len
        xyz = Molecule[i].xyz - com
        push!(out, Atom(Molecule[i].Z, Molecule[i].mass, xyz))
    end
    return(out)
end

