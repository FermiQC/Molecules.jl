# Not using this anymore
struct MTable
    symels
    table
end

struct Group
    idxs
    elements
    mtable
    classes
    order
    ctable
end

function multifly(symels::Vector{Symel}, A::Symel, B::Symel)
    Crrep = A.rrep * B.rrep
    for (i,g) in enumerate(symels)
        if isapprox(Crrep, g.rrep)
            return i,g
        end
    end
    throw(ArgumentError("No match found for Symels $(A.symbol) and $(B.symbol)!"))
end

function build_mult_table(symels)
    h = length(symels)
    t = zeros(Int64,h,h)
    for (i,a) in enumerate(symels)
        for (j,b) in enumerate(symels)
            t[i,j] = multifly(symels, a, b)[1]
        end
    end
    return t
end

function find_inv(mtable, g)
    for i = 1:size(mtable)[1]
        if mtable[i,g] == 1
            return i
        end
    end
end

function find_classes(mtable)
    h = length(mtable.symels)
    assigned_elements = []
    classes = []
    for a = 1:h
        if a ∉ assigned_elements
            results = []
            for b = 1:h
                if a != b
                    binv = find_inv(mtable, b)
                    r1 = mtable.table[binv,a]
                    r = mtable.table[r1,b]
                    push!(results, r)
                end
            end
            push!(classes,unique(results))
            for i in unique(results)
                push!(assigned_elements, i)
            end
        end
    end
    return classes
end

function build_regular_repr(mtab)
    l = size(mtab.table)[1]
    reg_reprs = []
    mtab_rarr = zeros(l,l)
    for i = 1:l
        i_inv = find_inv(mtab, i)
        mtab_rarr[i,:] .= mtab.table[i_inv,:]
    end
    for r = 1:l
        reg_repr = zeros(Int64,l,l)
        for i = 1:l
            for j = 1:l
                if mtab_rarr[i,j] == r
                    reg_repr[i,j] = 1
                end
            end
        end
        push!(reg_reprs, reg_repr)
    end
    return(reg_reprs)
end

function jordanize(m)
    evals, evecs = eigen(m)
    println(evals)
end