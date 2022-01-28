export buildD, findSEA, checkSEA

struct DistanceMatrix
    atoms::Vector{Any}
    distances::AbstractArray
end

mutable struct SEA{I,F}
    label::String
    set::Vector{I}
    paxis::Vector{F}
end

function distance(A::Atom, B::Atom)
    a = A.xyz
    b = B.xyz
    return (((a[1]-b[1])^2) + ((a[2]-b[2])^2) + ((a[3]-b[3])^2))^0.5
end

function buildD(mol::Molecule)
    len = size(mol)[1]
    atoms = []
    arr = zeros(Float64, (len,len))
    for i = 1:len
        push!(atoms, mol[i])
    end
    for i = 1:len, j = i+1:len
        arr[i,j] = distance(mol[i], mol[j])
        arr[j,i] = arr[i,j]
    end
    return DistanceMatrix(atoms, arr)
end

function checkSEA(A::Vector{Float64}, B::Vector{Float64}, δ)::Bool
    a = sort(A)
    b = sort(B)
    z = a - b
    for i in z
        if abs(i) < 10.0^(-δ)
            continue
        else
            return false
        end
    end
    return true
end

function findSEA(D::DistanceMatrix, δ::Int)
    len = size(D.distances)[1]
    out = Tuple[]
    for i = 1:len, j = i+1:len
        if cmp(D.atoms[i].mass, D.atoms[j].mass) == 0
            c = checkSEA(D.distances[:,i], D.distances[:,j], δ)
            if c
                push!(out, (i,j))
            end
        end
    end
    nonos = Int[]
    biggerun = SEA[]
    for i = 1:len
        if i in nonos
            continue
        else
            biggun = Int[]
            push!(biggun,i)
        end

        for k in out
            if i in k
                if i == k[1]
                    push!(biggun, k[2])
                    push!(nonos, k[2])
                else
                    push!(biggun, k[1])
                    push!(nonos, k[1])
                end
            end
        end
    sea = SEA("None", biggun, [0.0, 0.0, 0.0])
    push!(biggerun, sea)
    end

    return biggerun
end


