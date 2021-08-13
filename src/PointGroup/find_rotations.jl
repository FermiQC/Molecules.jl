export find_rotation_sets

struct rotation_element
    axis::Vector{Any}
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

function find_rotation_sets(mol::Molecule, SEAs::AbstractArray)
    out = []
    for seaidx in SEAs
        outinloop = []
        len = size(seaidx)[1]
        sea = [mol[i] for i in seaidx]
        if len < 2
            println("single atom")
        elseif len == 2
            println("linear")
        else
            moit = eigenmoit(calcmoit(sea))
            Ia, Ib, Ic = moit[1]
            Iav, Ibv, Icv = [moit[2][:,idx] for idx = 1:3]
            if Ia + Ib == Ic
                paxis = Icv
                if Ia == Ib
                    for i = 2:len
                        if isfactor(len, i)
                            re = rotation_element(paxis, i)
                            push!(outinloop, re)
                        end
                    end
                else
                    for i = 2:len-1
                        if isfactor(len, i)
                            re = rotation_element(paxis, i)
                            push!(outinloop, re)
                        end
                    end
                end
            else
                if Ia != Ib != Ic
                    for i in [Iav, Ibv, Icv]
                        re = rotation_element(i, 2)
                        push!(outinloop, re)
                    end
                else
                    if Ia == Ib
                        paxis = Icv
                    else
                        paxis = Iav
                    end
                    k = len / 2
                    for i = 2:k
                        if isfactor(k, i)
                            re = rotation_element(paxis, i)
                            push!(outinloop, re)
                        end
                    end
                end
            end
        end
        push!(out, outinloop)
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
    println(rotation_set_intersection(rotation_set))
end