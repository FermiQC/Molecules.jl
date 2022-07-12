export find_rotation_sets

struct rotation_element{F}
    axis::Vector{F}
    order::Int
end

function Base.:(==)(A::rotation_element, B::rotation_element)
    issame_axis(A.axis, B.axis) &&
    A.order == B.order
end

function issame_axis(a, b)
    """
    Checks if two axes are equivalent
    |a⋅b| = 1 assuming ||a|| = ||b|| = 1
    """
    A = normalize(a)
    B = normalize(b)
    d = abs(A ⋅ B)
    return isapprox(d, 1.0, atol=tol)
end

function rotation_set_intersection(rotation_set)
    """
    Finds the intersection of the elements in the sets
    within 'rotation_set'
    """
    out = deepcopy(rotation_set[1])
    len = size(rotation_set)[1]
    if len > 2
        for i = 1:size(rotation_set)[1]
            intersect!(out, rotation_set[i])
        end
    end
    return out
end

function find_rotation_sets(mol::Vector{<:Atom}, SEAs)
    """
    Loop through SEAs to label each set, and find principal axis and potential rotation orders
    as outlined in the flowchart in DOI: 10.1002/jcc.23493

    Returns a list of lists with each SEA's potential rotations
    """
    out_all_SEAs = []
    for sea in SEAs
        out_per_SEA = []
        len = size(sea.set)[1]
        if len < 2
            sea.label = "Single Atom"
        elseif len == 2
            sea.label = "Linear"
            sea.paxis = mol[sea.set[1]].xyz - mol[sea.set[2]].xyz
        else
            moit = eigenmoit(calcmoit([mol[i] for i in sea.set]))
            Ia, Ib, Ic = moit[1]
            Iav, Ibv, Icv = [moit[2][:,idx] for idx = 1:3]
            if isapprox(Ia, Ib, atol=tol) && isapprox(Ia, Ic, atol=tol)
                sea.label = "Spherical"
            elseif isapprox(Ia+Ib, Ic, atol=tol)
                paxis = Icv
                if isapprox(Ia, Ib, atol=tol)
                    sea.label = "Regular Polygon"
                    sea.paxis = paxis
                    for i = 2:len
                        if isfactor(len, i)
                            re = rotation_element(paxis, i)
                            push!(out_per_SEA, re)
                        end
                    end
                else
                    sea.label = "Irregular Polygon"
                    sea.paxis = paxis
                    for i = 2:len-1
                        if isfactor(len, i)
                            re = rotation_element(paxis, i)
                            push!(out_per_SEA, re)
                        end
                    end
                end
            else
                if !(isapprox(Ia, Ib, atol=tol) || isapprox(Ib, Ic, atol = tol))
                    sea.label = "Asymmetric Rotor"
                    for i in [Iav, Ibv, Icv]
                        re = rotation_element(i, 2)
                        push!(out_per_SEA, re)
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
                            push!(out_per_SEA, re)
                        end
                    end
                end
            end
            push!(out_all_SEAs, out_per_SEA)
        end
    end
    return out_all_SEAs
end

function isfactor(n, a)
    """
    Simple function, is a a factor of n? Return true/false
    """
    if n % a == 0
        return true
    else
        return false
    end
end

function find_rotations(mol::Vector{<:Atom}, rotation_set)
    """
    Find all actual rotation elements for the Molecule
    based off of potential rotation elements from 'find_rotation_sets'
    """
    molmoit = eigenmoit(calcmoit(mol))
    if molmoit[1][1] == 0.0 && molmoit[1][2] == molmoit[1][3]
        # This block executes if the Molecule is linear
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

function find_c2(mol::Vector{<:Atom}, SEAs)
    """
    Find C2 rotations as outlined in the paper DOI: 10.1002/jcc.23493
    Uses 'c2a', 'c2b', and 'c2c' to determine the presence of at least
    one C2 axis.
    """
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
                if sea.label == "Linear"
                    for sea2 in SEAs
                        if sea == sea2
                            continue
                        elseif sea2.label == "Linear"
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

function is_there_ortho_c2(mol::Vector{<:Atom}, cn_axis, SEAs) where T
    """
    Finds C2 axes orthogonal to principal axis 'paxis' of Molecule
    When passing a principal axis to 'c2a', 'c2b', and 'c2c', returns
    true/false if there exists at least one C2 orthogonal to principal axis
    """
    len = size(mol)[1]
    for sea in SEAs
        c2b_chk = c2b(mol, cn_axis, sea)
        if c2b_chk !== nothing
            return true, c2b_chk
            #return true
        else c2a(mol, cn_axis, sea)
            c2a_chk = c2a(mol, cn_axis, sea)
            if c2a_chk !== nothing
                return true, c2a_chk
                #return true
            elseif sea.label == "Linear"
                for sea2 in SEAs
                    if sea == sea2
                        continue
                    elseif sea2.label == "Linear"
                        c2c_chk = c2c(mol, cn_axis, sea, sea2)
                        if c2c_chk !== nothing
                            return true, c2c_chk
                            #return true
                        end
                    end
                end
            end         
        end
    end
    return false, nothing
end

function num_C2(mol::Vector{<:Atom}, SEAs)
    """
    Counts number of unique C2 axes. Particularly useful for
    finding high symmetry groups. Uses C2 algorithms a and b
    but unlike previous functions, finds all unique C2 axes and 
    returns the count.
    """
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
    return size(unique_axes)[1], unique_axes
end

function c2a(mol, sea)
    """
    Constructs a line going from the COM to the midpoint of every pair
    of atoms in a set of SEAs (assuming that midpoint is not the COM).
    That line is a potential C2 axis
    """
    len = size(sea.set)[1]
    for i = 1:len-1, j = i+1:len
        midpoint = mol[sea.set[i]].xyz + mol[sea.set[j]].xyz
        if midpoint == [0.0, 0.0, 0.0]
            continue
        else
            c2 = Molecules.Cn(midpoint, 2)
            molB = Molecules.transform(mol, c2)
            check = Molecules.isequivalent(mol, molB)
            if check
                return normalize(midpoint)
            else
                continue
            end
        end
    end
    return nothing
end

function c2b(mol, sea)
    """
    Constructs a line from the COM to each atom in a set of SEAs.
    That line is a potential C2 axis.
    """
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
    """
    Given two linear sets of SEAs (ij and kl), a potential C2 axis
    is the vector r_C2 = (i-j) x (k-l)
    """
    rij = mol[sea1.set[1]].xyz - mol[sea1.set[2]].xyz
    rkl = mol[sea2.set[1]].xyz - mol[sea2.set[2]].xyz
    c2_axis = normalize(cross(rij, rkl))
    c2 = Molecules.Cn(c2_axis, 2)
    molB = Molecules.transform(mol, c2)
    if Molecules.isequivalent(mol, molB)
        return c2_axis
    else
        return nothing
    end
end

function c2a(mol, cn_axis, sea)
    len = size(sea.set)[1]
    for i = 1:len-1, j = i+1:len
        midpoint = mol[sea.set[i]].xyz + mol[sea.set[j]].xyz
        if midpoint == [0.0, 0.0, 0.0]
            continue
        elseif issame_axis(midpoint, cn_axis)
            continue
        else
            c2 = Molecules.Cn(midpoint, 2)
            molB = Molecules.transform(mol, c2)
            check = Molecules.isequivalent(mol, molB)
            if check
                #return true
                return normalize(midpoint)
            else
                continue
            end
        end
    end
    return nothing
    #return false
end

function c2b(mol, cn_axis, sea)
    len = size(sea.set)[1]
    for i = 1:len
        c2_axis = normalize(mol[sea.set[i]].xyz)
        if issame_axis(c2_axis, cn_axis)
            continue
        else
            c2 = Molecules.Cn(c2_axis, 2)
            molB = Molecules.transform(mol, c2)
            check = Molecules.isequivalent(mol, molB)
            if check
                #return true
                return c2_axis
            else
                continue
            end
        end
    end
    return nothing
    #return false
end

function c2c(mol, cn_axis, sea1, sea2)
    rij = mol[sea1.set[1]].xyz - mol[sea1.set[2]].xyz
    rkl = mol[sea2.set[1]].xyz - mol[sea2.set[2]].xyz
    c2_axis = normalize(cross(rij, rkl))
    if issame_axis(c2_axis, cn_axis)
        return nothing
        #return false
    else
        c2 = Molecules.Cn(c2_axis, 2)
        molB = Molecules.transform(mol, c2)
        if Molecules.isequivalent(mol, molB)
            return c2_axis
        end
        #return Molecules.isequivalent(mol, molB)
    end
end

function all_c2a(mol, sea)
    out = []
    len = size(sea.set)[1]
    for i = 1:len-1, j = i+1:len
        midpoint = mol[sea.set[i]].xyz + mol[sea.set[j]].xyz
        if midpoint == [0.0, 0.0, 0.0]
            continue
        else
            c2 = Molecules.Cn(midpoint, 2)
            molB = Molecules.transform(mol, c2)
            check = Molecules.isequivalent(mol, molB)
            if check
                push!(out, normalize(midpoint))
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

function highest_ordered_axis(rotations::Vector{Any})
    # Returns the largest element of the set of rotation orders
    len = size(rotations)[1]
    ns = []
    for i = 1:len
        push!(ns, rotations[i].order)
    end
    return sort(ns)[len]
end

function is_there_σh(mol, paxis)
    """
    Checks for reflection across plane with normal vector
    equal to principal rotation axis
    """
    σh = Molecules.reflection_matrix(paxis)
    molB = Molecules.transform(mol, σh)
    return Molecules.isequivalent(mol, molB)
end

function is_there_σv(mol, SEAs, paxis)
    """
    Checks for the presence of reflection planes that
    contain the principal c2_axis

    Returns true/false and normal vector to a symmetry plane/nothing.
    """
    axes = []
    for sea in SEAs
        len = size(sea.set,1)
        if len < 2
            continue
        end
        A = sea.set[1]
        for i = 2:len
            B = sea.set[i]
            #n = normalize(mol[A].xyz - mol[B].xyz)
            n = mol[A].xyz - mol[B].xyz
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
            return true, planar_mol_axis(mol)
        else
            return false, nothing
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
        end
        if check
            push!(unique_axes, i)
        end
    end
    for i in unique_axes
        if issame_axis(i, paxis)
            continue
        else
            return true, i
        end
    end
    return false, nothing
end

function mol_is_planar(mol)
    """
    Checks if Molecule is planar
    """
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

function planar_mol_axis(mol)
    "Finds normal vector of molecule plane, 
    assuming molecule is planar as checked by mol_is_planar"
    len = size(mol,1)
    for i = 1:len
        for j = i:len
            a = normalize(mol[i].xyz)
            b = normalize(mol[j].xyz)
            chk = dot(a, b)
            if !isapprox(chk, 1.0, atol=tol)
                return normalize(cross(a,b))
            end
        end
    end
    return nothing
end

function find_C3s_for_Ih(mol)
    """
    Finds the twenty C3 axes for an Ih point group so the paxis and saxis can be defined
    """
    c3_axes = []
    for i in mol, j in mol, k in mol
        if i != j && i != k
            rij = i.xyz - j.xyz
            rjk = j.xyz - k.xyz
            rik = i.xyz - k.xyz
            nij = norm(rij)
            njk = norm(rjk)
            nik = norm(rik)
            if isapprox(nij, njk, atol=tol) && isapprox(nij, nik, atol=tol)
                c3_axis = normalize(cross(rij, rjk))
                c3 = Molecules.Cn(c3_axis, 3)
                molB = Molecules.transform(mol, c3)
                if Molecules.isequivalent(mol, molB)
                    push!(c3_axes, c3_axis)
                end
            end
        end
    end
    unique_axes = [c3_axes[1]]
    for i in c3_axes
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
    chk = size(unique_axes)[1]
    if chk != 10
        throw(ErrorException("Unexpected number of C3 axes for Ih point group: Found $(chk) unique C3 axes"))
    end
    return unique_axes
end

function find_C4s_for_Oh(mol)
    """
    Finds the three C4 axes for an Oh point group so the paxis and saxis can be defined
    """
    c4_axes = []
    for i in mol, j in mol, k in mol, l in mol
        if i != j && k != l && i != k
            rij = i.xyz - j.xyz
            rjk = j.xyz - k.xyz
            rkl = k.xyz - l.xyz
            ril = i.xyz - l.xyz
            nij = norm(rij)
            njk = norm(rjk)
            nkl = norm(rkl)
            nil = norm(ril)
            if isapprox(nij, njk, atol=tol) && isapprox(nkl, nil, atol=tol) && isapprox(nij, nkl, atol=tol)
                c4_axis = normalize(cross(rij, rjk))
                c4 = Molecules.Cn(c4_axis, 4)
                molB = Molecules.transform(mol, c4)
                if Molecules.isequivalent(mol, molB)
                    push!(c4_axes, c4_axis)
                end
            end
        end
    end
    unique_axes = [c4_axes[1]]
    for i in c4_axes
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
    chk = size(unique_axes)[1]
    if chk != 3
        throw(ErrorException("Unexpected number of C4 axes for Oh point group: Found $(chk) unique C4 axes"))
    end
    return unique_axes
end