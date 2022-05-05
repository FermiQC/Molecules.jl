module Symmetry

using Molecules

tol = 1E-5

include("sea.jl")
include("moit.jl")
include("symmetry_finding.jl")
include("flowchart.jl")
include("CharacterTables/CharacterTables.jl")

end

