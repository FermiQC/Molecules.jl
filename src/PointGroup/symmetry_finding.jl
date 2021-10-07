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

function find_rotation_sets(mol::Molecule, SEAs)
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

function find_rotations(mol::Molecule, rotation_set)
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

function find_c2(mol::Molecule, SEAs)
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

function is_there_ortho_c2(mol::Molecule, cn_axis::Vector{T}, SEAs) where T
    """
    Finds C2 axes orthogonal to principal axis 'paxis' of Molecule
    When passing a principal axis to 'c2a', 'c2b', and 'c2c', returns
    true/false if there exists at least one C2 orthogonal to principal axis
    """
    len = size(mol)[1]
    for sea in SEAs
        if c2a(mol, cn_axis, sea)
            return true
        elseif c2b(mol, cn_axis, sea)
            return true
        elseif sea.label == "Linear"
            for sea2 in SEAs
                if sea == sea2
                    continue
                elseif sea2.label == "Linear"
                    if c2c(mol, cn_axis, sea, sea2)
                        return true
                    end
                end
            end
        end         
    end
    return false
end

function num_C2(mol::Molecule, SEAs)
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
    return size(unique_axes)[1]
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
                return midpoint
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
        c2_axis = mol[sea.set[i]].xyz
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
    c2_axis = cross(rij, rkl)
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
        c2_axis = mol[sea.set[i]].xyz
        if issame_axis(c2_axis, cn_axis)
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
    if issame_axis(c2_axis, cn_axis)
        return false
    else
        c2 = Molecules.Cn(c2_axis, 2)
        molB = Molecules.transform(mol, c2)
        return Molecules.isequivalent(mol, molB)
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
        c2_axis = mol[sea.set[i]].xyz
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

function is_there_sigmah(mol::Molecule, paxis)
    """
    Checks for reflection across plane with normal vector
    equal to principal rotation axis
    """
    σh = Molecules.reflection_matrix(paxis)
    molB = Molecules.transform(mol, σh)
    return Molecules.isequivalent(mol, molB)
end

function is_there_sigmav(mol, SEAs, paxis)
    """
    Checks for the presence of reflection planes that
    contain the principal c2_axis
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
            return true
        else
            return false
        end
    end
    unique_axes = [axes[1]]
    # This part is VERY inefficient!!!
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
        if issame_axis(i, paxis)
            continue
        else
            return true
        end
    end
    return false
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