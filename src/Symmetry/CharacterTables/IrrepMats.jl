c3 = cos(π/3)
s3 = sin(π/3)
c5 = cos(π/5)
s5 = sin(π/5)
c25 = cos(2π/5)
s25 = sin(2π/5)
irrm_C3v = Dict(
    "A1" => [1,1,1,1,1,1],
    "A2" => [1,1,1,-1,-1,-1],
    "E"  => [[1 0;0 1], [-c3 -s3; s3 -c3], [-c3 s3; -s3 -c3], [1 0;0 -1], [-c3 -s3; -s3 c3], [-c3 s3; s3 c3]]
)
irrm_C4v = Dict(
    "A1" => [1,1,1,1,1,1,1,1],
    "A2" => [1,1,1,1,-1,-1,-1,-1],
    "B1" => [1,-1,-1,1,1,1,-1,-1],
    "B2" => [1,-1,-1,1,-1,-1,1,1],
    "E"  => [[1 0;0 1],[0 1; -1 0],[0 -1; 1 0],[-1 0; 0 -1],[1 0; 0 -1],[-1 0; 0 1],[0 1; 1 0],[0 -1; -1 0]]
)
irrm_C5v = Dict(
    "A1" => [1,1,1,1,1,1,1,1,1,1],
    "A2" => [1,1,1,1,1,-1,-1,-1,-1,-1],
    "E1" => [[1 0;0 1],[c5 s5; -s5 c5],[c5 -s5; s5 c5],[c25 s25; -s25 c25],[c25 -s25; s25 c25],[c25 s25; c25 s25],[c25 s25; c25 s25],[c25 s25; c25 s25],[c25 s25; c25 s25],[1 0;0 -1]],
)
function irrep_things()
    return irrm_C4v
end

function mtable_check(irrm, mtable)
    l = size(mtable)[1]
    for i = 1:l, j = 1:l
        if mtable[i,j] ∈ multifly(irrm, i, j)
            continue
        else
            return false
        end
    end
    return true
end

function multifly(irrm, a, b)
    l = length(irrm)
    out = []
    errl = []
    for i = 1:l
        r = irrm[a]*irrm[b]
        push!(errl, r)
        if isapprox(irrm[i], r, atol = 1e-10)
            push!(out, i)
        end
    end
    return out
    return 0
end

function goat_chk(irrm)
    g = length(irrm)
    d = size(irrm[1])[1]
    gc = zeros(d,d,d,d)
    for i = 1:d
        for j = 1:d
            for l = 1:d
                for m = 1:d
                    gc[i,j,l,m] = goat_chk(irrm, i, j, l, m)
                end
            end
        end
    end
    return gc
end

function goat_chk(irrm, i, j, l, m)
    a = 0.0
    for r = 1:length(irrm)
        a += irrm[r][i,l]*irrm[r][j,m]
    end
    return a
end