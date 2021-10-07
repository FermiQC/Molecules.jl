export find_point_group
import Molecules: center_of_mass, translate


function find_point_group(mol::Molecule)
    # Center Molecule on COM
    mol = translate(mol, center_of_mass(mol))
    # Calculate moment of inertia tensor (moit) eigenvalues
    Ia_mol, Ib_mol, Ic_mol = eigenmoit(calcmoit(mol))[1]

    if Ia_mol == 0.0
        # Linear Molecule
        if Molecules.isequivalent(Molecules.transform(mol, Molecules.inversion_matrix()), mol)
            return("Dinfh")
        else
            return("Cinfv")
        end
    elseif isapprox(Ia_mol, Ib_mol, atol=tol) && isapprox(Ia_mol, Ic_mol, atol = tol)
        # Molecule with high symmetry
        D = buildD(mol)
        SEAs = findSEA(D, 5)
        n = num_C2(mol, SEAs)
        invertable = Molecules.isequivalent(Molecules.transform(mol, Molecules.inversion_matrix()), mol)
        if n == 15
            return("Ih")
        elseif n == 9
            if invertable
                return("Oh")
            else
                return("O")
            end
        elseif n == 3
            if invertable
                return("Th")
            else
                return("Td")
            end
        end
    else
        # Build distance matrix and find symmetry equivalent atoms
        D = buildD(mol)
        SEAs = findSEA(D, 5)
        # Find sets of rotation elements
        rot_set = find_rotation_sets(mol, SEAs)
        if size(rot_set)[1] < 1
            # Path followed if no rotation elements found by 'find_rotation_sets'
            # Therefore, maximum rotation order of 2
            c2 = find_c2(mol, SEAs)
            if c2 !== nothing
                c2_ortho = is_there_ortho_c2(mol, c2, SEAs)
                σh = is_there_sigmah(mol, c2)
                if c2_ortho
                    if σh
                        return("D2h")
                    else
                        σv = is_there_sigmav(mol, SEAs, c2)
                        if σv
                            return("D2d")
                        else
                            return("D2")
                        end
                    end
                else
                    if σh
                        return("C2h")
                    else
                        σv = is_there_sigmav(mol, SEAs, c2)
                        if σv
                            return("C2v")
                        else
                            S4 = Molecules.Sn(c2, 4)
                            molB = Molecules.transform(mol, S4)
                            check = Molecules.isequivalent(mol, molB)
                            if check
                                return("S4")
                            else
                                return("C2")
                            end
                        end
                    end
                end
            else
                if Molecules.isequivalent(Molecules.transform(mol, Molecules.inversion_matrix()), mol)
                    return("Ci")
                else
                    σv = is_there_sigmav(mol, SEAs, [0.0,0.0,0.0])
                    if σv
                        return("Cs")
                    else
                        return("C1")
                    end
                end
            end
        else
            # Path followed if rotation elements found by 'find_rotation_sets'
            rots = find_rotations(mol, rot_set)
            Cn = highest_ordered_axis(rots)
            paxis = rots[1].axis
            q1 = is_there_ortho_c2(mol, paxis, SEAs)
            q2 = is_there_sigmah(mol, paxis)
            q3 = is_there_sigmav(mol, SEAs, paxis)
            if q1
                if q2
                    return("D"*string(Cn)*"h")
                else
                    if q3
                        return("D"*string(Cn)*"d")
                    else
                        return("D"*string(Cn))
                    end
                end
            else
                if q2
                    return("C"*string(Cn)*"h")
                else
                    if q3
                        return("C"*string(Cn)*"v")
                        return("C",Cn,"v")
                    else
                        S2n = Molecules.Sn(paxis, Cn*2)
                        molB = Molecules.transform(mol, S2n)
                        check = Molecules.isequivalent(mol, molB)
                        if check
                        return("S"*string(Cn*2))
                            return("S",Cn*2)
                        else
                        return("C"*string(Cn))
                            return("C",Cn)
                        end
                    end
                end
            end
        end
    end
end