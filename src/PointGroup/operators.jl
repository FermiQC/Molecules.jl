export Cn, σ, Sn, i, transform, issame

function Cn(v, n)
    θ = 2*π / n
    cosθ = cos(θ)
    sinθ = sin(θ)
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
    return O
end

function σ(v)
    O = zeros(eltype(v), (3,3))
    for i = 1:3, j = i:3
        if i == j
            O[i,i] = 1 - 2*v[i]^2
        else
            O[i,j] = -2 * v[i] * v[j]
            O[j,i] = O[i,j]
        end
    end
    return O
end

function Sn(v, n)
    return Cn(v, n) * σ(v)
end

function i()
    a = zeros(Int8, (3,3))
    for i = 1:3
        a[i,i] = -1
    end
    return a
end

function transform(A, O)
    len = size(A)[1]
    out = []
    for i = 1:len
        newxyz = O * A[i].xyz
        push!(out, Atom(A[i].Z, A[i].mass, newxyz))
    end
    return out
end

function issame(A, B)
    h = []
    len = size(A)[1]
    for i = 1:len
        for j = 1:len
            if A[i].mass == B[j].mass
                if isapprox(A[i].xyz, B[j].xyz, rtol=1E-5)
                    push!(h, i)
                    break
                end
            end
        end
    end
    return size(h)[1] == len
end