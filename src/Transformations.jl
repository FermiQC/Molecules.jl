"""
    Molecules.translate!(atoms::Vector{<:Atom}, r::Vector{T})

Translate atoms into position r.
"""
function translate!(atoms::Vector{<:Atom}, r)
    for a in 1:length(atoms)
        atoms[a] = translate(atoms[a], r)
    end
end

"""
    Molecules.translate(atoms::Vector{<:Atom})

Return a copy of atoms translated into position r.
"""
function translate(atoms::Vector{<:Atom}, r)
    newatoms = deepcopy(atoms)
    translate!(newatoms, r)
    return newatoms
end

"""
    Molecules.translate(atom::A, r) where A <: Atom

Return a new `Atom` translated into position r.
"""
function translate(atom::A, r) where A <: Atom
    return Atom(atom.Z, atom.mass, atom.xyz - r)
end

"""
    Molecules.Cn(v::Vector, n::Int)

Returns an n rotation matrix about vector v
"""
function Cn(v, n::Integer)
    θ = 2*π / n
    return rotation_matrix(v, θ)
end

"""
    Molecules.rotation_matrix(v, θ)

Returns a rotation matrix by θ about a rotation axis v
"""
function rotation_matrix(V, θ)
    cosθ = cos(θ)
    sinθ = sin(θ)
    v = normalize(V)
    a = [1,2,3]
    O = zeros(eltype(v), (3,3))
    O .+= 1 - cosθ
    for i = 1:3, j = 1:3
        if i == j
            O[i,i] *= v[i]^2
            O[i,i] += cosθ
        else
            O[i,j] *= v[i]*v[j]
            b = [i,j]
            C = [i for i in a if i ∉ b][1]
            if i < j
                O[i,j] += (-1)^(i+j) * v[C]*sinθ
            else
                O[i,j] += (-1)^(i+j-1) * v[C]*sinθ
            end
        end
    end
    return SMatrix{3,3}(O)
end

function rotation_matrixd(v, θ)
    r = deg2rad(θ)
    return rotation_matrix(v, r)
end

"""
    Molecules.reflection_matrix(v::Vector)

Returns a reflection matrix through the plane with normal vector v
"""
function reflection_matrix(V)
    O = zeros(eltype(V), (3,3))
    v = normalize(V)
    for i = 1:3, j = i:3
        if i == j
            O[i,i] = 1 - 2*v[i]^2
        else
            O[i,j] = -2 * v[i] * v[j]
            O[j,i] = O[i,j]
        end
    end
    return SMatrix{3,3}(O)
end

function σ(v)
    return reflection_matrix(v)
end

"""
    Molecules.Sn(v, n)

Returns Sn improper rotation about vector v
"""    
function Sn(v, n)
    return Cn(v, n) * σ(v)
end

"""
    Molecules.inversion_matrix()

Returns a diagonal matrix with -1 along the diagonal
"""
function inversion_matrix()
    return SMatrix{3,3}(
        [-1.0  0.0 0.0;
          0.0 -1.0 0.0;
          0.0  0.0 -1.0]
    )
end

function i()
    return inversion_matrix()
end
    
"""
    Molecules.transform!(A::Vector{<:Atom}, O)

Transforms xyz coordinates of each atom in A by some Matrix O
"""
function transform!(atoms::Vector{<:Atom}, O)
    for a in 1:length(atoms)
        atoms[a] = transform(atoms[a], O)
    end
end

"""
    Molecules.transform(A::Vector{<:Atom}, O)

Retruns a list of atoms with transformed xyz coordinates by some Matrix O
"""
function transform(atoms::Vector{<:Atom}, O)
    newatoms = deepcopy(atoms)
    transform!(newatoms, O)
    return newatoms
end

function transform(atom::A, O) where A <: Atom
    return Atom(atom.Z, atom.mass, O * atom.xyz)
end

"""
    Molecules.issame(A::Vector{<:Atom}, B::Vector{<:Atom})

Checks to see if two arrays of atoms are equivalent under permutation
Returns true or false
"""
function isequivalent(A::Vector{<:Atom}, B::Vector{<:Atom})
    len = length(A)
    if size(atom_map(A,B), 1) == len
        return true
    end
    return false
end

function atom_map(A::Vector{<:Atom}, B::Vector{<:Atom})
    h = []
    len = length(A)
    #println(A, "\n", B)
    for i = 1:len
        for j = 1:len
            if A[i].mass == B[j].mass
                zs = broadcast(abs, A[i].xyz - B[j].xyz)
                c =  isapprox(zs, [0.0,0.0,0.0], atol=tol)
                #println(zs, " : ", c)
                #if isapprox(broadcast(abs, A[i].xyz - B[j].xyz), [0.0,0.0,0.0], rtol=1E-5)
                if c
                    push!(h, j)
                    break
                end
            end
        end
    end
    return h
end

