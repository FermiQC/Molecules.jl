export Cn, σ, Sn, i

function Cn(v, n)
    θ = 2*π / n
    cosθ = cos(θ)
    sinθ = sin(θ)
    a = [1,2,3]
    O = zeros(Float64, (3,3))
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
    O = zeros(Float64, (3,3))
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