@testset "Properties" begin
    
    # Fake atoms
    atoms = [
        Atom(1, 1.008, [5.0, 5.0, 0.0])
        Atom(2, 1.008, [5.0, -5.0, 0.0])
        Atom(3, 1.008, [-5.0, -5.0, 0.0])
        Atom(4, 1.008, [-5.0, 5.0, 0.0])
    ]
    mol = Molecule(atoms)

    @test Molecules.center_of_mass(mol) == [0.0, 0.0, 0.0]
    @test Molecules.nuclear_repulsion(mol) ≈ 1.6816285798739845
    @test Molecules.nuclear_dipole(mol) == [-20.0, 0.0, 0.0]
    @test Molecules.symbol(atoms[3]) == "Li"

    # Test nuclear repulsion gradient
    findif = zeros(3)
    h = 1e-8
    ap = deepcopy(atoms)
    am = deepcopy(atoms)
    r = [SVector(-h,0.0,0.0), SVector(0.0,-h,0.0), SVector(0.0,0.0,-h)]
    for i = 1:length(atoms)
        for k = 1:3 
            # Translate XYZ coordinates
            ap[i] = Molecules.translate(ap[i], r[k]) 
            am[i] = Molecules.translate(am[i], -r[k]) 
            findif[k] = (Molecules.nuclear_repulsion(ap) - Molecules.nuclear_repulsion(am)) / (2h)
            # Reverse transformation
            ap[i] = Molecules.translate(ap[i], -r[k]) 
            am[i] = Molecules.translate(am[i], r[k]) 
        end
        s =  sum((findif .- Molecules.∇nuclear_repulsion(mol, i) * Molecules.angstrom_to_bohr).^2)
        @test √(s/3) < 1e-7
    end

    # Test nuclear repulsion Hessian (finite difference of the gradient just
    # validated above, same translate/undo pattern, same angstrom_to_bohr
    # correction since ∇nuclear_repulsion's raw output is per-bohr already)
    for i = 1:length(atoms), j = 1:length(atoms)
        findif2 = zeros(3, 3)
        ap = deepcopy(atoms)
        am = deepcopy(atoms)
        for k = 1:3
            ap[j] = Molecules.translate(ap[j], r[k])
            am[j] = Molecules.translate(am[j], -r[k])
            findif2[:, k] .= (Molecules.∇nuclear_repulsion(ap, i) .- Molecules.∇nuclear_repulsion(am, i)) ./ (2h)
            ap[j] = Molecules.translate(ap[j], -r[k])
            am[j] = Molecules.translate(am[j], r[k])
        end
        s = sum((findif2 .- Molecules.∇2nuclear_repulsion(mol, i, j) .* Molecules.angstrom_to_bohr).^2)
        @test √(s/9) < 1e-6
    end
end