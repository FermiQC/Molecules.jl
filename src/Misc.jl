import Base: show

# Returns a nicely formatted string with all the molecule's information
function string_repr(M::Molecule)
    out = ""
    out = out*format("Molecule:\n\n")
    out = out*format(Molecules.get_xyz(M.atoms))
    out = out*format("\n")
    out = out*format("\nCharge: {}   ", M.charge)
    out = out*format("Multiplicity: {}   \n", M.multiplicity)
    out = out*format("Nuclear repulsion: {:15.10f}", Molecules.nuclear_repulsion(M))
    return out
end

function string_repr(A::Atom)
    out = elements[A.Z].name * " Z = $(A.Z) M = "
    out = out*format("{:6.4f}", A.mass)
    out = out*format("\n Center: {:15.10f}  {:15.10f}  {:15.10f}", A.xyz...)
end

# Pretty printing
function show(io::IO, ::MIME"text/plain", X::T) where T <: Union{Atom, Molecule}
    print(io, string_repr(X))
end