export find_point_group
import Molecules: center_of_mass, translate

"""
    Molecules.Symmetry.find_point_group(mol::Vector{Atoms})

    Determines point group of mol by searching for identifying symmetry operations.
    Returns a point group string ("D3h", "C5v", "Oh", etc.), a primary axis, and a secondary axis.
    The primary and secondary axes define vectors we can rotate generated symmetry elements onto.
    This rotation is carried out by the function Molecules.Symmetry.CharacterTables.rotate_symels_to_mol().
    
    The procedure for detecting symmetry is similar to that employed in doi:10.1002/jcc.23493.
"""
function find_point_group(mol::Molecule)
    # Translate Molecule to COM origin
    mol = translate(mol, center_of_mass(mol))
    # Calculate moment of inertia tensor (moit) eigenvalues
    Ia_mol, Ib_mol, Ic_mol = eigenmoit(calcmoit(mol))[1]
    paxis = [0;0;0]
    saxis = [0;0;0]
    if Ia_mol == 0.0
        # Linear Molecule
        if Molecules.isequivalent(Molecules.transform(mol, Molecules.inversion_matrix()), mol)
            pg = "Dinfh"
        else
            pg = "Cinfv"
        end
    elseif isapprox(Ia_mol, Ib_mol, atol=tol) && isapprox(Ia_mol, Ic_mol, atol = tol)
        # Molecule with high symmetry, count number of C2 operations and look for i operation to determine PG
        D = buildD(mol)
        SEAs = findSEA(D, 5)
        n = num_C2(mol, SEAs)
        invertable = Molecules.isequivalent(Molecules.transform(mol, Molecules.inversion_matrix()), mol)
        if n == 15
            pg = "Ih"
        elseif n == 9
            if invertable
                pg = "Oh"
            else
                pg = "O"
            end
        elseif n == 3
            if invertable
                pg = "Th"
            else
                pg = "Td"
            end
        end
    else
        # Build distance matrix and find symmetry equivalent atoms
        D = buildD(mol)
        SEAs = findSEA(D, 5)
        # Find sets of rotation elements, may return nothing if there is only a rotation of order 2
        rot_set = find_rotation_sets(mol, SEAs)
        if size(rot_set)[1] < 1
            # Path followed if no rotation elements found by 'find_rotation_sets'
            # Therefore, maximum rotation order of 2
            c2 = find_c2(mol, SEAs)
            if c2 !== nothing
                paxis = c2
                c2_ortho_chk, c2_ortho = is_there_ortho_c2(mol, c2, SEAs)
                σh = is_there_σh(mol, c2)
                if c2_ortho_chk
                    saxis = c2_ortho
                    if σh
                        pg = "D2h"
                    else
                        σv_chk, σv = is_there_σv(mol, SEAs, c2)
                        if σv_chk
                            pg = "D2d"
                        else
                            pg = "D2"
                        end
                    end
                else
                    if σh
                        pg = "C2h"
                    else
                        σv_chk, σv = is_there_σv(mol, SEAs, c2)
                        if σv_chk
                            if σv !== nothing
                                saxis = normalize(σv)
                            end
                            pg = "C2v"
                        else
                            S4 = Molecules.Sn(c2, 4)
                            molB = Molecules.transform(mol, S4)
                            check = Molecules.isequivalent(mol, molB)
                            if check
                                pg = "S4"
                            else
                                pg = "C2"
                            end
                        end
                    end
                end
            else
                if Molecules.isequivalent(Molecules.transform(mol, Molecules.inversion_matrix()), mol)
                    pg = "Ci"
                else
                    σv_chk, σv = is_there_σv(mol, SEAs, [0.0,0.0,0.0])
                    if σv_chk
                        if σv !== nothing
                            paxis = σv
                        end
                        pg = "Cs"
                    else
                        pg = "C1"
                    end
                end
            end
        else
            # Path followed if rotation elements found by 'find_rotation_sets'
            rots = find_rotations(mol, rot_set)
            Cn = highest_ordered_axis(rots)
            paxis = rots[1].axis
            q1, c2_ortho = is_there_ortho_c2(mol, paxis, SEAs) # Check for C2 orthogonal to paxis
            q2 = is_there_σh(mol, paxis) # Check for reflection plane orthogonal to paxis
            q3, σv = is_there_σv(mol, SEAs, paxis) # Check for reflection planes containing paxis
            if q1
                saxis = c2_ortho
                if q2
                    pg = "D"*string(Cn)*"h"
                else
                    if q3
                        pg = "D"*string(Cn)*"d"
                    else
                        pg = "D"*string(Cn)
                    end
                end
            else
                if q2
                    pg = "C"*string(Cn)*"h"
                else
                    if q3
                        if σv !== nothing
                            saxis = normalize(cross(paxis, σv))
                        end
                        pg = "C"*string(Cn)*"v"
                    else
                        S2n = Molecules.Sn(paxis, Cn*2)
                        molB = Molecules.transform(mol, S2n)
                        check = Molecules.isequivalent(mol, molB)
                        if check
                            pg = "S"*string(Cn*2)
                        else
                            pg = "C"*string(Cn)
                        end
                    end
                end
            end
        end
    end
    return pg, paxis, saxis
end