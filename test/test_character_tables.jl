using Molecules.Symmetry.CharacterTables

rt2 = sqrt(2.0)
rt3 = sqrt(3.0)
twocos72 = 2*cosd(72.0)
twocos144 = 2*cosd(144.0)

include("pgs/Cn.jl")
include("pgs/Cnv.jl")
include("pgs/Cnh.jl")
include("pgs/Sn.jl")
include("pgs/Dn.jl")
include("pgs/Dnd.jl")
include("pgs/Dnh.jl")
allofem = ["C1", "C2",  "C3",  "C4",  "C5",  "C6",
                 "C2v", "C3v", "C4v", "C5v", "C6v",
                 "C2h", "C3h", "C4h", "C5h", "C6h",
                 "S4",  "S6",  "S8",  "S10",
                 "D2",  "D3",  "D4",  "D5",  "D6",
                 "D2d", "D3d", "D4d", "D5d", "D6d",
                 "D2h", "D3h", "D4h", "D5h", "D6h"]

@testset "Symmetry" begin
    include("test_symel_gen.jl")
    include("test_ctab_gen.jl")
end
