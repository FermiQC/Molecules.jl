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
    for i = 1:length(atoms)
        for k = 1:3 
            ap[i].xyz[k] += h 
            am[i].xyz[k] -= h 
            findif[k] = (Molecules.nuclear_repulsion(ap) - Molecules.nuclear_repulsion(am)) / (2h)
            ap[i].xyz[k] -= h 
            am[i].xyz[k] += h 
        end
        s =  sum((findif .- Molecules.∇nuclear_repulsion(mol, i) * Molecules.angstrom_to_bohr).^2)
        @test √(s/3) < 1e-7
    end
end