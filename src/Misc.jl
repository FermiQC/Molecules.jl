import Base: show

# Returns a nicely formatted string with all the molecule's information
string_repr(M::Molecule) = """
    Molecule:

    $(Molecules.get_xyz(M.atoms))

    Charge: $(M.charge)
    Multiplicity: $(M.multiplicity)
    Nuclear repulsion: $(@sprintf("%15.10f", Molecules.nuclear_repulsion(M)))"""

string_repr(A::Atom) = """
    $(elements[A.Z].name) Z = $(A.Z) M = $(@sprintf("%6.4f", A.mass))
    Center: $(@sprintf("%15.10f  %15.10f  %15.10f", A.xyz...))"""

# Pretty printing
function show(io::IO, ::MIME"text/plain", X::Union{Atom, Molecule})
    print(io, string_repr(X))
end
