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
    mult_table
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

function string_repr(symels::Vector{Symel})
    out = ""
    for s in symels
        out *= "$(format("{:>8s}",s.symbol))"*": $(s.rrep)\n"
    end
    return out
end

function string_repr(symtext::SymText)
    out = "SymText for the $(symtext.pg) Point Group\n"
    out *= string_repr(symtext.ctab)
    out *= "\nSymmetry Elements (Symels)\n"
    out *= string_repr(symtext.symels)
    out *= "\nMultiplication table (row×column)\n"
    h = length(symtext.symels)
    out *= " "^(6)*"|"
    for i = 1:h
        out *= "$(format("{:>5d}",symtext.mult_table[1,i]))"
    end
    longstrang = "-"^6*"|"*"-"^(5*h)
    out *= "\n"*longstrang*"\n"
    for i = 1:h
        out *= "$(format("{:>5d}",symtext.mult_table[i,1])) |"
        for j = 1:h
            out *= "$(format("{:>5d}",symtext.mult_table[i,j]))"
        end
        out *= "\n"
    end
    out *= "\nClass Map\n$(symtext.class_map)\n\nAtom Map\n"
    c,d = size(symtext.atom_map)
    out *= " "^(9)*"|"
    for s = 1:d
        out *= "$(format("{:>8s}",symtext.symels[s].symbol))"
    end
    longstrang2 = "-"^(9)*"|"*"-"^(8*d)
    out *= "\n"*longstrang2*"\n"
    for a = 1:c
        out *= "$(format("{:>8d} |",a))"
        for b = 1:d
            out *= "$(format("{:>8d}",symtext.atom_map[a,b]))"
        end
        out *= "\n"
    end
    return out
end

function show(io::IO, ::MIME"text/plain", symels::Vector{Symel})
    print(io, string_repr(symels))
end

function show(io::IO, ::MIME"text/plain", ctab::Chartable)
    print(io, string_repr(ctab))
end

function show(io::IO, ::MIME"text/plain", symtext::SymText)
    print(io, string_repr(symtext))
end

function getindex(ctab::Chartable, irrep::String)
    return ctab.characters[findall(x->x==irrep, ctab.irreps)[1],:]
end