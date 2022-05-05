"""
    Molecules.function parse_file(file::String; unit=:angstrom, F=Float64)

Reads a xyz file and return a vector o `Atom` objects. Units can be indicated through the keyword argument `unit`.
Coordinates are stored as `F` (`Float64` by default) where `F` is keyword arguments.

See also: parse_string
"""
function parse_file(file::String; unit=:angstrom, F=Float64)
    molstring = read(file, String)
    parse_string(molstring::String; unit=unit, F=F)
end

"""
    Molecules.function parse_string(molstring::String; unit=:angstrom, F=Float64)

Reads a String representing a XYZ file and return a vector o `Atom` objects. Units can be indicated through the keyword argument `unit`.
Coordinates are stored as `F` (`Float64` by default) where `F` is keyword arguments.

See also: parse_file
"""
function parse_string(molstring::String; unit=:angstrom, F=Float64)

    # Get a list of Atom objects from String
    if unit == :bohr
        conv = bohr_to_angstrom
    elseif unit == :angstrom
        conv = 1.0
    else
        throw(ArgumentError("unsupported unit requested for XYZ parsing: $unit"))
    end

    # Regex for parsing
    # Match a floating point, e.g. XYZ coordinate or mass, with at least one blank space on the left of it
    re_float = r"\s+([+-]?(?:\d+[.]?\d*|\d*[.]?\d+))"
    # Match an atom id (String for atomic symbol or integer for atomic number) with one or more blank spaces on the left of it
    re_id = r"\s*(\d+|\w{1,2})"

    atoms = Atom{F}[]

    # Possible line formats
    # ID Mass X Y Z
    r1 = re_id * re_float * re_float * re_float * re_float
    # ID X Y Z
    r2 = re_id * re_float * re_float * re_float

    line_num = 1
    for line in split(strip(molstring), "\n")

        if contains(line, r1)
            m = match(r1, line)
            # id is taken as a Symbol if a String is given ("H" → :H), otherwise convert to number ("6" → 6)
            id = occursin(r"\d+", m.captures[1]) ? parse(Int64, m.captures[1]) : Symbol(m.captures[1])

            # Try to parse Mass string into a Float
            try
                mass = parse(F, m.captures[2])
            catch ArgumentError
                throw(ArgumentError("Failed to process data in line $line_num:\n $(line)"))
            end
            str_xyz = m.captures[3:5]

        elseif contains(line, r2)
            m = match(r2, line)
            # id is taken as a Symbol if a String is given ("H" → :H), otherwise convert to number ("6" → 6)
            id = occursin(r"\d+", m.captures[1]) ? parse(Int64, m.captures[1]) : Symbol(m.captures[1])

            # Get mass from PeriodicTable
            mass = convert(F, elements[id].atomic_mass / 1u"u")
            str_xyz = m.captures[2:4]
        else
            throw(ArgumentError("Failed to process data in line $line_num:\n $(line)"))
        end

        Z = elements[id].number

        # Convert String vector to Float vector
        xyz = zeros(F, 3)
        try
            xyz .= parse.(F, str_xyz).*conv
        catch ArgumentError
            throw(ArgumentError("Failed to process XYZ coordinates in line $line_num:\n $(m[2:4])"))
        end

        push!(atoms, Atom(Z, mass, xyz))
        line_num += 1
    end

    return atoms
end

"""
    get_xyz(M::Molecule)

Returns a XYZ string in angstrom for the given Molecule.
"""
function get_xyz(M::Molecule)
    molstring = ""
    for A in M
        molstring *= format("{}   {: 15.12f}   {: 15.12f}   {: 15.12f}\n", Molecules.symbol(A), A.xyz...)
    end
    return molstring
end