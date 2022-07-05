module CharacterTables
import Base
import Base: show
using LinearAlgebra
using Molecules
using Molecules:i,Cn,Sn,σ
export Symel
tol = 1E-5

include("Basics.jl")
include("PointGroupGenerators.jl")
include("CharacterTableGenerators.jl")
include("MultiplicationTable.jl")
include("IrrepMats.jl")
include("Main.jl")
end