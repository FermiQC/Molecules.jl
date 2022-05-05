module CharacterTables
import Base
using LinearAlgebra
using Molecules
using Molecules:i,Cn,Sn,σ
export Symel
tol = 1E-5

include("Basics.jl")
include("PointGroupGenerators.jl")
include("CharacterTableGenerators.jl")
include("Main.jl")

end