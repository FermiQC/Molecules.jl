using Molecules
using Molecules.PointGroup
using LinearAlgebra

export SeaPlane

path = joinpath(@__DIR__, "../../test/xyz/benzene.xyz")
mol = Molecules.parse_file(path)
D = buildD(mol)
sea = findSEA(D,5)

function CheckSea(sea)
    z = length(sea)
    SeaSet = 0 
    for i in 1:z
        if length(sea[i]) > 1
            SeaSet = sea[i]
        end
    end
    return SeaSet
end
function IsEqual(A::Vector{Float64},B::Vector{Float64},epsilon::Float64=1.0e-6)
    result = abs(dot(A,B)) 
    return abs( result - 1  ) <= epsilon

end
function SeaPlane(SeaSet)
    SeaSet = CheckSea(sea)
    A = 0
    B = 0
    R = []
    Normal = 0
    list_norm = []
    for a in 1:length(SeaSet) 
        for b in 1 + a:length(SeaSet) 
            A = SeaSet[a]
            B = SeaSet[b]
            rx = mol[A].xyz[1] - mol[B].xyz[1]
            ry = mol[A].xyz[2] - mol[B].xyz[2]
            rz = mol[A].xyz[3] - mol[B].xyz[3]
            R = [rx, ry, rz]
            normR = (norm(R))
            Normal = R / normR
            list_norm = append!(list_norm, [Normal])
            
        end
    end
    return list_norm; 
end
SeaSet = CheckSea(sea)
axes = SeaPlane(SeaSet)
#axes = [...]
if size(axes, 1) < 1
    return nothing
end

unique_vectors = [axes[1]]
for i in axes
    check = true
    for j in unique_vectors
        if IsEqual(i,j)
            check = false
            break
        end
    end
    if check
        push!(unique_vectors, i)
    end
end
println(unique_vectors)