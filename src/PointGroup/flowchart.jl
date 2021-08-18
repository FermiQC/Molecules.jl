export find_point_group
import Molecules: center_of_mass, translate


function find_point_group(mol::Molecule)
    mol = translate(mol, center_of_mass(mol))
    Ia_mol, Ib_mol, Ic_mol = eigenmoit(calcmoit(mol))[1]
    if Ia_mol == 0.0
        println("Linear")
    elseif Ia_mol == Ib_mol == Ic_mol
        println("High sym")
    else
        D = buildD(mol)
        SEAs = findSEA(D, 5)
        rot_set = find_rotation_sets(mol, SEAs)
        
        if size(rot_set)[1] < 1
            println("No Cn s.t. n>2")
        else
            rots = find_rotations(mol, rot_set)
            Cn = highest_ordered_axis(rots)
            paxis = rots[1].axis
            q1 = is_there_ortho_c2(mol, paxis, SEAs)
            q2 = is_there_sigmah(mol, paxis)
            if q1
                if q2
                    println("D",Cn,"h")
                else
                    println("D",Cn)
                end
            else
                if q2
                    println("C",Cn,"h")
                else
                    println("C",Cn)
                end
            end
        end
    end
end