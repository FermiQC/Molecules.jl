
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

function multifly(symels, A::Symel, B::Symel)
    Crrep = A.rrep * B.rrep
    for (i,g) in enumerate(symels)
        s = sum(abs.(g.rrep - Crrep))
        if abs(s) < tol
            return i,g
        end
    end
    println(i,g)
    throw(ArgumentError("No match found!"))
end

function build_mult_table(symels)
    h = length(symels)
    t = zeros(Int64,h,h)
    for (i,a) in enumerate(symels)
        for (j,b) in enumerate(symels)
            t[i,j] = multifly(symels, a, b)[1]
        end
    end
    return MTable(symels, t)
end

function find_inv(mtable, g)
    for i = 1:length(mtable.symels)
        if mtable.table[i,g] == 1
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
