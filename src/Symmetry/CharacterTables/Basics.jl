struct Symel
    symbol
    rrep
end

struct Chartable
    name
    irreps
    classes
    class_orders
    characters
end

struct SymText
    pg
    symels
    ctab
    class_map
    atom_map
end

struct PG
    str
    family
    n
    subfamily
end

function reduce(n, i)
    g = gcd(n, i)
    return div(n, g), div(i, g)
end

function gcd(A, B)
    # A quick implementation of the Euclid algorithm for finding the greatest common divisor
    a = max(A,B)
    b = min(A,B)
    if a == 0
        return b
    elseif b == 0
        return a
    else
        r = a % b
        gcd(b, r)
    end
end

function Base.:(==)(A::Symel, B::Symel)
    if sum(A.rrep .- B.rrep) < tol
        return true
    else
        return false
    end
end
using Formatting
function string_repr(ctab::Chartable)
    l = length(ctab.classes) + 1
    longstrang = "-"^(l*10)*"\n"
    out = longstrang
    out *= "$(ctab.name) Character Table:\n"*(" "^15)
    for i = 1:length(ctab.classes)
        out *= "$(format("{:10s}",ctab.classes[i]))"
    end
    out *= "\n"
    for i = 1:length(ctab.irreps)
        out *= "$(format("{:10s}",ctab.irreps[i]))"
        for j = 1:length(ctab.characters[i,:])
            out *= "$(format("{:>10.3f}",ctab.characters[i,j]))"
        end
        out *= "\n"
    end
    out *= longstrang
    return out
end

#function show(io::IO, ::MIME"text/plain", symels::Vector{Symel})
#    print(io, string_repr(symels))
#end

function show(io::IO, ::MIME"text/plain", ctab::Chartable)
    print(io, string_repr(ctab))
end
