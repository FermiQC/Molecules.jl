export find_rotation_sets

struct rotation_element{F}
    axis::Vector{F}
    order::Int
end

function Base.:(==)(A::rotation_element, B::rotation_element)
    isapprox(A.axis, B.axis, rtol=1E-5) &&
    A.order == B.order
end

function rotation_set_intersection(rotation_set)
    out = deepcopy(rotation_set[1])
    len = size(rotation_set)[1]
    if len > 2
        for i = 1:size(rotation_set)[1]
            intersect!(out, rotation_set[i])
        end
    end
    return out
end

function find_rotation_sets(mol::Molecule, SEAs)
    out = []
    tol = 1E-5
    for sea in SEAs
        outinloop = []
        len = size(sea.set)[1]
        #sea = [mol[i] for i in seaidx]
        if len < 2
            sea.label = "Single Atom"
        elseif len == 2
            sea.label = "Linear"
            sea.paxis = mol[sea.set[1]].xyz - mol[sea.set[2]].xyz
        else
            moit = eigenmoit(calcmoit([mol[i] for i in sea.set]))
            Ia, Ib, Ic = moit[1]
            Iav, Ibv, Icv = [moit[2][:,idx] for idx = 1:3]
            if isapprox(Ia, Ib, rtol=tol) && isapprox(Ia, Ic, rtol=tol)
                sea.label = "Spherical"
            elseif isapprox(Ia+Ib, Ic, rtol=tol)
                paxis = Icv
                if isapprox(Ia, Ib, rtol=tol)
                    sea.label = "Regular Polygon"
                    sea.paxis = paxis
                    for i = 2:len
                        if isfactor(len, i)
                            re = rotation_element(paxis, i)
                            push!(outinloop, re)
                        end
                    end
                else
                    sea.label = "Irregular Polygon"
                    sea.paxis = paxis
                    for i = 2:len-1
                        if isfactor(len, i)
                            re = rotation_element(paxis, i)
                            push!(outinloop, re)
                        end
                    end
                end
            else
                if Ia != Ib != Ic
                    sea.label = "Asymmetric Rotor"
                    for i in [Iav, Ibv, Icv]
                        re = rotation_element(i, 2)
                        push!(outinloop, re)
                    end
                else
                    if Ia == Ib
                        paxis = Icv
                        sea.label = "Oblate Symmetric Top"
                        sea.paxis = paxis
                    else
                        paxis = Iav
                        sea.label = "Prolate Symmetric Top"
                        sea.paxis = paxis
                    end
                    k = div(len, 2)
                    for i = 2:k
                        if isfactor(k, i)
                            re = rotation_element(paxis, i)
                            push!(outinloop, re)
                        end
                    end
                end
            end
            push!(out, outinloop)
        end
    end
    return out
end

function isfactor(n, a)
    if n % a == 0
        return true
    else
        return false
    end
end

function find_rotations(mol::Molecule, rotation_set)
    molmoit = eigenmoit(calcmoit(mol))
    if molmoit[1][1] == 0.0 && molmoit[1][2] == molmoit[1][3]
        paxis = normalize(mol[1].xyz)
        re = rotation_element(paxis, 0)
        return [re]
    end
    rsi = rotation_set_intersection(rotation_set)
    out = []
    for i in rsi
        rmat = Molecules.Cn(i.axis, i.order)
        molB = Molecules.transform(mol, rmat)
        c = Molecules.isequivalent(molB, mol)
        if c
            push!(out, i)
        end
    end
    return out
end

function find_c2(mol::Molecule, SEAs)
    len = size(mol)[1]
    for sea in SEAs
        a = c2a(mol, sea)
        if a !== nothing
            return a 
        else
            b = c2b(mol, sea)
            if b !== nothing
                return b
            else

                if sea.label == "linear"
                    for sea2 in SEAs
                        if sea == sea2
                            continue
                        elseif sea2.label == "linear"
                            c = c2c(mol, sea, sea2)
                            if c !== nothing
                                return c
                            end
                        end
                    end
                end     
            end
        end    
    end
    return nothing
end

function is_there_ortho_c2(mol::Molecule, cn_axis::Vector{T}, SEAs) where T
    len = size(mol)[1]
    for sea in SEAs
        if c2a(mol, cn_axis, sea)
            return true
        elseif c2b(mol, cn_axis, sea)
            return true
        elseif sea.label == "linear"
            for sea2 in SEAs
                if sea == sea2
                    continue
                elseif sea2.label == "linear"
                    if c2c(mol, cn_axis, sea, sea2)
                        return true
                    end
                end
            end
        end         
    end
    return false
end

function c2a(mol, cn_axis, sea)
    len = size(sea.set)[1]
    for i = 1:len-1, j = i+1:len
        midpoint = normalize(mol[sea.set[i]].xyz + mol[sea.set[j]].xyz)
        d = abs(midpoint ⋅ cn_axis)
        if midpoint == [0.0, 0.0, 0.0]
            continue
        elseif isapprox(d, 1.0, rtol=1E-5)
            continue
        else
            c2 = Molecules.Cn(midpoint, 2)
            molB = Molecules.transform(mol, c2)
            check = Molecules.isequivalent(mol, molB)
            if check
                return true
            else
                continue
            end
        end
    end
    return false
end

function c2b(mol, cn_axis, sea)
    len = size(sea.set)[1]
    for i = 1:len
        c2_axis = normalize(mol[sea.set[i]].xyz)
        d = c2_axis ⋅ cn_axis
        if isapprox(d, 1.0, rtol=1E-5)
            continue
        else
            c2 = Molecules.Cn(c2_axis, 2)
            molB = Molecules.transform(mol, c2)
            check = Molecules.isequivalent(mol, molB)
            if check
                return true
            else
                continue
            end
        end
    end
    return false
end

function c2c(mol, cn_axis, sea1, sea2)
    rij = mol[sea1.set[1]].xyz - mol[sea1.set[2]].xyz
    rkl = mol[sea2.set[1]].xyz - mol[sea2.set[2]].xyz
    c2_axis = cross(rij, rkl)
    d = c2_axis ⋅ cn_axis
    if isapprox(d, 1.0, rtol=1E-5)
        return false
    else
        c2 = Molecules.Cn(c2_axis, 2)
        molB = Molecules.transform(mol, c2)
        return Molecules.isequivalent(mol, molB)
    end
end

function c2a(mol, sea)
    len = size(sea.set)[1]
    for i = 1:len-1, j = i+1:len
        midpoint = normalize(mol[sea.set[i]].xyz + mol[sea.set[j]].xyz)
        if midpoint == [0.0, 0.0, 0.0]
            continue
        else
            c2 = Molecules.Cn(midpoint, 2)
            molB = Molecules.transform(mol, c2)
            check = Molecules.isequivalent(mol, molB)
            if check
                return midpoint
            else
                continue
            end
        end
    end
    return nothing
end

function c2b(mol, sea)
    len = size(sea.set)[1]
    for i = 1:len
        c2_axis = normalize(mol[sea.set[i]].xyz)
        c2 = Molecules.Cn(c2_axis, 2)
        molB = Molecules.transform(mol, c2)
        check = Molecules.isequivalent(mol, molB)
        if check
            return c2_axis
        else
            continue
        end
    end
    return nothing
end

function c2c(mol, sea1, sea2)
    rij = mol[sea1.set[1]].xyz - mol[sea1.set[2]].xyz
    rkl = mol[sea2.set[1]].xyz - mol[sea2.set[2]].xyz
    c2_axis = cross(rij, rkl)
    c2 = Molecules.Cn(c2_axis, 2)
    molB = Molecules.transform(mol, c2)
    if Molecules.isequivalent(mol, molB)
        return c2_axis
    else
        return nothing
    end
end

function highest_ordered_axis(rotations::Vector{Any})
    len = size(rotations)[1]
    ns = []
    for i = 1:len
        push!(ns, rotations[i].order)
    end
    return sort(ns)[len]
end

function is_there_sigmah(mol::Molecule, paxis)
    σh = Molecules.reflection_matrix(paxis)
    molB = Molecules.transform(mol, σh)
    return Molecules.isequivalent(mol, molB)
end

function num_C2(mol::Molecule, SEAs)
    axes = []
    for sea in SEAs
        a = all_c2a(mol, sea)
        if a !== nothing
            for i in a
                push!(axes, i)
            end
        end

        b = all_c2b(mol, sea)
        if b !== nothing
            for i in b
                push!(axes, i)
            end
        end
    end
    if size(axes)[1] < 1
        return nothing
    end
    unique_axes = [axes[1]]
    for i in axes
        check = true
        for j in unique_axes
            if issame_axis(i,j)
                check = false
                break
            end
        end
        if check
            push!(unique_axes, i)
        end
    end
    return size(unique_axes)[1]
end

function all_c2a(mol, sea)
    out = []
    len = size(sea.set)[1]
    for i = 1:len-1, j = i+1:len
        midpoint = normalize(mol[sea.set[i]].xyz + mol[sea.set[j]].xyz)
        if midpoint == [0.0, 0.0, 0.0]
            continue
        else
            c2 = Molecules.Cn(midpoint, 2)
            molB = Molecules.transform(mol, c2)
            check = Molecules.isequivalent(mol, molB)
            if check
                push!(out, midpoint)
            else
                continue
            end
        end
    end
    if size(out)[1] < 1
        return nothing
    end
    return out
end

function all_c2b(mol, sea)
    out = []
    len = size(sea.set)[1]
    for i = 1:len
        c2_axis = normalize(mol[sea.set[i]].xyz)
        c2 = Molecules.Cn(c2_axis, 2)
        molB = Molecules.transform(mol, c2)
        check = Molecules.isequivalent(mol, molB)
        if check
            push!(out, c2_axis)
        else
            continue
        end
    end
    if size(out)[1] < 1
        return nothing
    end
    return out
end

function issame_axis(a, b)
    d = abs(a ⋅ b)
    return isapprox(d, 1.0, rtol=1E-5)
end

function is_there_sigmav(mol, SEAs, paxis)
    axes = []
    for sea in SEAs
        len = size(sea.set,1)
        if len < 2
            continue
        end
        A = sea.set[1]
        for i = 2:len
            B = sea.set[i]
            n = normalize(mol[A].xyz - mol[B].xyz)
            σ = Molecules.reflection_matrix(n)
            molB = Molecules.transform(mol, σ)
            check = Molecules.isequivalent(mol, molB)
            if check
                push!(axes, n)
            else
                continue
            end
        end
    end
    if size(axes,1) < 1
        if mol_is_planar(mol)
            return true
        else
            return false
        end
    end
    unique_axes = [axes[1]]
    for i in axes
        check = true
        for j in unique_axes
            if issame_axis(i,j)
                check = false
                break
            end
            if check
                push!(unique_axes, i)
            end
        end
    end

    for i in unique_axes
        d = abs(i ⋅ paxis)
        if isapprox(d, 1.0, atol=1E-5)
            continue
        else
            return true
        end
    end
    return false
end

function mol_is_planar(mol)
    t = eltype(mol[1].xyz)
    len = size(mol,1)
    mat = zeros(Float64, (3,len))
    for i = 1:len
        mat[:,i] = mol[i].xyz
    end
    rank = LinearAlgebra.rank(mat)
    if rank < 3
        return true
    end
    return false
end