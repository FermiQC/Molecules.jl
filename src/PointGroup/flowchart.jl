export find_point_group
import Molecules: center_of_mass, translate


function find_point_group(mol::Molecule)
    mol = translate(mol, center_of_mass(mol))
    Ia_mol, Ib_mol, Ic_mol = eigenmoit(calcmoit(mol))[1]
    #println(Ia_mol, "   ", Ib_mol, "   ", Ic_mol)
    if Ia_mol == 0.0
        if Molecules.isequivalent(Molecules.transform(mol, Molecules.inversion_matrix()), mol)
            println("D∞h")
        else
            println("C∞h")
        end
    elseif isapprox(Ia_mol, Ib_mol, rtol=1E-5) && isapprox(Ia_mol, Ic_mol, rtol = 1E-5)
        D = buildD(mol)
        SEAs = findSEA(D, 5)
        println("High sym") # Keep going
        n = num_C2(mol, SEAs)
        println(n)
        invertable = Molecules.isequivalent(Molecules.transform(mol, Molecules.inversion_matrix()), mol)
        if n == 15
            println("Ih")
        elseif n == 9
            if invertable
                println("Oh")
            else
                println("O")
            end
        elseif n == 3
            if invertable
                println("Th")
            else
                println("Td")
            end
        else
            println("ERROR!",n)
        end
    else
        D = buildD(mol)
        SEAs = findSEA(D, 5)
        rot_set = find_rotation_sets(mol, SEAs)

        if size(rot_set)[1] < 1
            c2 = find_c2(mol, SEAs)
            if c2 !== nothing
                c2_ortho = is_there_ortho_c2(mol, c2, SEAs)
                σh = is_there_sigmah(mol, c2)
                if c2_ortho
                    if σh
                        println("D2h")
                    else
                        println("D2") # Keep going
                        σv = is_there_sigmav(mol, SEAs, c2)
                        if σv
                            println("D2d")
                        else
                            println("D2")
                        end
                    end
                else
                    if σh
                        println("C2h")
                    else
                        σv = is_there_sigmav(mol, SEAs, c2)
                        if σv
                            println("C2v")
                        else
                            S4 = Molecules.Sn(c2, 4)
                            molB = Molecules.transform(mol, S4)
                            check = Molecules.isequivalent(mol, molB)
                            if check
                                println("S4")
                            else
                                println("C2")
                            end
                        end
                    end
                end
            else
                if Molecules.isequivalent(Molecules.transform(mol, Molecules.inversion_matrix()), mol)
                    println("Ci")
                else
                    σv = is_there_sigmav(mol, SEAs, [0.0,0.0,0.0])
                    if σv
                        println("Cs")
                    else
                        println("C1")
                    end
                end
            end
        else
            rots = find_rotations(mol, rot_set)
            Cn = highest_ordered_axis(rots)
            paxis = rots[1].axis
            q1 = is_there_ortho_c2(mol, paxis, SEAs)
            q2 = is_there_sigmah(mol, paxis)
            q3 = is_there_sigmav(mol, SEAs, paxis)
            if q1
                if q2
                    println("D",Cn,"h")
                else
                    if q3
                        println("D",Cn,"d")
                    else
                        println("D",Cn)
                    end
                end
            else
                if q2
                    println("C",Cn,"h")
                else
                    if q3
                        println("C",Cn,"v")
                    else
                        S2n = Molecules.Sn(paxis, Cn*2)
                        molB = Molecules.transform(mol, S2n)
                        check = Molecules.isequivalent(mol, molB)
                        if check
                            println("S",Cn*2)
                        else
                            println("C",Cn)
                        end
                    end
                end
            end
        end
    end
end