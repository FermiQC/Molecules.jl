@testset "Properties" begin
    
    # Fake atoms
    atoms = [
        Atom(1, 1.008, [5.0, 5.0, 0.0])
        Atom(2, 1.008, [5.0, -5.0, 0.0])
        Atom(3, 1.008, [-5.0, -5.0, 0.0])
        Atom(4, 1.008, [-5.0, 5.0, 0.0])
    ]

    @test Molecules.center_of_mass(atoms) == [0.0, 0.0, 0.0]
    @test Molecules.nuclear_repulsion(atoms) ≈ 1.6816285798739845
end