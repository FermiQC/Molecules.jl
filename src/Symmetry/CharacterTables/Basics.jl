struct Symel
    symbol
    rrep
end

struct Chartable
    name
    symels
    irreps
    characters
end

struct PG
    str
    family
    n
    subfamily
end

function Base.:(==)(A::Symel, B::Symel)
    if sum(A.rrep .- B.rrep) < tol
        return true
    else
        return false
    end
end