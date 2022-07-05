@testset "Parsing" begin
    # Nice one
    molstringA = """
    C   -2.131551124300    2.286168823700    0.000000000000
    H   -1.061551124300    2.286168823700    0.000000000000
    H   -2.488213906200    1.408104616400    0.496683911300
    H   -2.488218762100    2.295059432700   -1.008766153900
    H   -2.488220057000    3.155340844300    0.512081313000"""

    # A bit messy
    molstringB = """
    C      -2.1315511243  2.2861688237    .0000000000                 
    H    -1.0615511243      2.2861688237     0.
    H    -2.4882139062   1.4081046164  0.4966839113                 
    H      -2.4882187621    2.2950594327     -1.0087661539                 
    H  -2.4882200570  3.1553408443      .5120813130"""

    # Messy with atomic number and masses
    molstringC = """
    6 12.011     -2.1315511243  2.2861688237    .0000000000                 
    H    -1.0615511243      2.2861688237     0.
    H    -2.4882139062   1.4081046164  0.4966839113                 
    1 1.008     -2.4882187621    2.2950594327     -1.0087661539                 
    H  -2.4882200570  3.1553408443      .5120813130"""

    atoms = [
        Atom(6, 12.011, [-2.1315511243,     2.2861688237,     0.0000000000])
        Atom(1,  1.008, [-1.0615511243,     2.2861688237,     0.0000000000])
        Atom(1,  1.008, [-2.4882139062,     1.4081046164,     0.4966839113])
        Atom(1,  1.008, [-2.4882187621,     2.2950594327,    -1.0087661539])
        Atom(1,  1.008, [-2.4882200570,     3.1553408443,     0.5120813130])
    ]

    @test Molecules.parse_string(molstringA) == atoms
    @test Molecules.parse_string(molstringB) == atoms
    @test Molecules.parse_string(molstringC) == atoms

    atoms = [
        Atom(8, 15.999, [1.2091536548,      1.7664118189,     -0.0171613972])
        Atom(1, 1.008,  [2.1984800075,      1.7977100627,      0.0121161719])
        Atom(1, 1.008,  [0.9197881882,      2.4580185570,      0.6297938832])
    ]

    file_path = joinpath(@__DIR__, "xyz/water.xyz")
    @test Molecules.parse_file(file_path) == atoms

    mol = Molecule(molstringC)
    @test Molecules.get_xyz(mol) |> chomp == molstringA 

    # Test parsing in bohr
    atoms = Molecules.parse_string("""
    H 0.0 0.0 0.0
    H 1.0 0.0 0.0
    """; unit=:bohr)
    @test Molecules.nuclear_repulsion(atoms) == 1.0

    # Error due invalid unit
    @test_throws ArgumentError Molecules.parse_string("""
    H 0.0 0.0 0.0"""; unit=:borabora)

    # Test molecule object
    mol = Molecule(molstringA)
    @test mol.Nα == mol.Nβ == length(mol.atoms) == 5
    @test_throws DomainError Molecule(molstringA, 11, 2)
    @test_throws DomainError Molecule(molstringA, 0, 0)
    @test_throws DomainError Molecule(molstringA, 1, 1)

    # Test erros during parsing
    @test_throws ArgumentError Molecule("Some nonsense")
    @test_throws ArgumentError Molecule("X 1.0 1.0 1.0")
    @test_throws ArgumentError Molecule("X 5.0 1.0 1.0 1.0")
end