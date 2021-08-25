
export SeaPlane, is_there_sigma

#path = joinpath(@__DIR__, "../../test/xyz/benzene.xyz")
#mol = Molecules.parse_file(path)
#D = buildD(mol)
#sea = findSEA(D,5)

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

function IsEqual(A::Vector{AbstractFloat}, B::Vector{AbstractFloat})
    result = abs(dot(A,B)) 
    return isapprox(result, 1.0, atol = 1.0E-6)
end

function SeaPlane(mol, SeaSet)
#    SeaSet = CheckSea(sea)
    A = 0
    B = 0
    R = []
    Normal = 0
    list_norm = []
    for a = 1:length(SeaSet) 
        for b = a+1:length(SeaSet) 
#            A = SeaSet[a]
#            B = SeaSet[b]
#            rx = mol[A].xyz[1] - mol[B].xyz[1]
#            ry = mol[A].xyz[2] - mol[B].xyz[2]
#            rz = mol[A].xyz[3] - mol[B].xyz[3]
#            R = [rx, ry, rz]
            R = mol[SeaSet[a]].xyz - mol[SeaSet[b]].xyz
#            normR = norm(R)
            Normal = R ./ norm(R)
            list_norm = append!(list_norm, [Normal])
            
        end
    end
    return list_norm; 
end

#SeaSet = CheckSea(sea)
#axes = SeaPlane(SeaSet)

function extract_unique_vectors(axes)
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
    return unique_vectors
end

function is_there_sigma(mol, SEAs, paxis)
    SeaSet = CheckSea(sea)
    axes = SeaPlane(SeaSet)
    u_axes = extract_unique_vectors(axes)
    for i in u_axes
        d = abs(i⋅paxis)
        if !isapprox(d, 1.0, atol = 1E-5)
            return true
        end
    end
    return false
end
