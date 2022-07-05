c3 = cos(π/3)
s3 = sin(π/3)
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

function irrep_things()
    return irrm_C4v
end

function mtable_from_irrm(irrm)
    l = length(irrm)
    mtab = zeros(l,l)
    for i = 1:l
        for j = 1:l
            mtab[i,j] = multifly(irrm, i, j)
        end
    end
    return mtab
end

function multifly(irrm, a, b)
    l = length(irrm)
    errl = []
    for i = 1:l
        r = irrm[a]*irrm[b]
        push!(errl, r)
        println(irrm[i])
        println(isapprox(irrm[i], irrm[a]*irrm[b]))
        if isapprox(irrm[i], r, atol = 1e-10)
            return i
        end
    end
    return errl[1]
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