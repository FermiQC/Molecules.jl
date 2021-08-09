using LinearAlgebra

export calcmoit, eigenmoit

function calcmoit(Atoms)
    I = zeros(Float64, (3,3))
    len = size(Atoms)[1]
    for i = 1:3, j = 1:3
        if i == j
            a = i % 3 + 1
            b = (i+1) % 3 + 1
            suum = 0.0
            for k = 1:len
                suum += Atoms[k].mass * (Atoms[k].xyz[a]^2 + Atoms[k].xyz[b]^2)
            end
            I[i,i] = suum
        else
            suum = 0.0
            for k = 1:len
                suum -= Atoms[k].mass * Atoms[k].xyz[i] * Atoms[k].xyz[j]
            end
            I[i,j] = suum
            I[j,i] = I[i,j]
        end
    end
    return I
end

function eigenmoit(moit)
    A = LinearAlgebra.Symmetric(moit)
    eval, evec = LinearAlgebra.eigen(A)
    return (eval, evec)
end