import Molecules.PointGroup: findSEA, buildD, find_rotation_sets, checkSEA, find_rotations, is_there_ortho_c2, find_c2, find_point_group
import Molecules: center_of_mass, translate
using Molecules
using Molecules.PointGroup

r2 = 2.0^0.5
r3 = (3.0^0.5)/2.0
Square = [Atom(1, 1.0, [-0.5, 0.0, 0.0]), Atom(1, 1.0, [0.5, 0.0, 0.0]), Atom(1, 1.0, [0.0, 0.5, 0.0]), Atom(1, 1.0, [0.0, -0.5, 0.0])]
Triangle = [Atom(1, 1.0, [0.0, 0.0, 0.0]), Atom(1, 1.0, [1.0, 0.0, 0.0]), Atom(1, 1.0, [0.5, (3.0^0.5)/2.0, 0.0])]
Dodecagon = [Atom(1, 1.0, [r3, 0.5, 0.0]), Atom(1, 1.0, [0.5, r3, 0.0]), Atom(1, 1.0, [0.0, 1.0, 0.0]), Atom(1, 1.0, [-0.5, r3, 0.0]), Atom(1, 1.0, [-r3, 0.5, 0.0]), Atom(1, 1.0, [-1.0, 0.0, 0.0]), Atom(1, 1.0, [-r3, -0.5, 0.0]), Atom(1, 1.0, [-0.5, -r3, 0.0]), Atom(1, 1.0, [0.0, -1.0, 0.0]), Atom(1, 1.0, [0.5, -r3, 0.0]), Atom(1, 1.0, [r3, -0.5, 0.0]), Atom(1, 1.0, [1.0, 0.0, 0.0])]
Line = [Atom(1, 1.0, [-0.5, 0.0, 0.0]), Atom(1, 1.0, [0.5, 0.0, 0.0])]
RecPrism = [Atom(1, 1.0, [1.0, 1.0, 1.0]), Atom(1, 1.0, [-1.0, 1.0, 1.0]), Atom(1, 1.0, [1.0, -1.0, 1.0]), Atom(1, 1.0, [-1.0, -1.0, 1.0]), Atom(1, 1.0, [1.0, 1.0, -2.0]), Atom(1, 1.0, [-1.0, 1.0, -2.0]), Atom(1, 1.0, [1.0, -1.0, -2.0]), Atom(1, 1.0, [-1.0, -1.0, -2.0])]
RecAnti = [Atom(1, 1.0, [1.0, 1.0, 1.0]), Atom(1, 1.0, [1.0, -1.0, 1.0]), Atom(1, 1.0, [-1.0, 1.0, 1.0]), Atom(1, 1.0, [-1.0, -1.0, 1.0]), Atom(1, 1.0, [r2, 0.0, -2.0]), Atom(1, 1.0, [0.0, r2, -2.0]), Atom(1, 1.0, [-r2, 0.0, -2.0]), Atom(1, 1.0, [0.0, -r2, -2.0])]
RecPrismO = [Atom(1, 1.0, [1.0, 1.0, 1.0]), Atom(1, 1.0, [-1.0, 1.0, 1.0]), Atom(1, 1.0, [1.0, -1.0, 1.0]), Atom(1, 1.0, [-1.0, -1.0, 1.0]), Atom(8, 15.99, [1.0, 1.0, -2.0]), Atom(8, 15.99, [-1.0, 1.0, -2.0]), Atom(8, 15.99, [1.0, -1.0, -2.0]), Atom(8, 15.99, [-1.0, -1.0, -2.0])]

names = ["C1", "Ci", "Cs", "Cs_nonplanar", "C2", "C3", "S4", "C2v", "C3v", "C2h", "C3h", "D2", "D3", "D2d", "D4d", "D3h", "D4h", "D5h", "D6h", "D12h", "D100h", 
         "Cinfv", "Dinfh", "Th", "Td", "cube", "octahedron", "dodecahedron", "icosahedron"]
pgs = ["C1", "Ci", "Cs", "Cs", "C2", "C3", "S4", "C2v", "C3v", "C2h", "C3h", "D2", "D3", "D2d", "D4d", "D3h", "D4h", "D5h", "D6h", "D12h", "D100h", 
       "Cinfv", "Dinfh", "Th", "Td", "Oh", "Oh", "Ih", "Ih"]
#names = ["D100h"]
#pgs = ["D100h"]
@testset "Point Groups" begin
    for i = 1:size(names)[1]
        path = joinpath(@__DIR__, "sxyz/"*names[i]*".xyz")
        mol = Molecules.parse_file(path)
        s = find_point_group(mol)
        @test s == pgs[i]
    end
end
#out = []
#for i = eachindex(mols)
#    path = joinpath(@__DIR__, "xyz/"*mols[i]*".xyz")
#    mol = Molecules.parse_file(path)
#    
#    D = buildD(mol)
#
#    sea = findSEA(D, 5)
#    push!(out, sea)
#
#end
#println(out)