import Molecules: translate

@testset "Transformations" begin
    atoms = [
        Atom(1, 1.008, [1.0, 2.0, 3.0])
    ]

    natoms = translate(atoms, [1.0, 2.0, 3.0])

    @test natoms[1].xyz == [0.0, 0.0, 0.0]
end