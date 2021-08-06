function parse_file(file::String; unit=:angstrom, F=Float64, I=Int16)
    molstring = read(file, String)
    parse_string(molstring::String; unit=unit, F=F, I=I)
end

function parse_string(molstring::String; unit=:angstrom, F=Float64, I=Int16)

    # Get a list of Atom objects from String
    if unit == :bohr
        conv = bohr_to_angstrom
    elseif unit == :angstrom
        conv = 1.0
    else
        throw(ArgumentError("unsupported unit requested for XYZ parsing: $unit"))
    end

    # Regex for parsing
    # Match a floating point, e.g. XYZ coordinate or mass, with at least one blank space in front of it
    re_float = r"\s+([+-]?(?:\d+[.]?\d*|\d*[.]?\d+))"
    # Match an atom id (String for atomic symbol or integer for atomic number) with one or more blank spaces in front of it
    re_id = r"\s*(\d*|\w{1,2})"

    atoms = Atom{I,F}[]

    # Possible line formats
    # ID Mass X Y Z
    r1 = re_id * re_float * re_float * re_float * re_float
    # ID X Y Z
    r2 = re_id * re_float * re_float * re_float

    line_num = 1
    for line in split(strip(molstring), "\n")

        if contains(line, r1)
            m = match(r1, line)
            id = Symbol(m.captures[1])

            # Try to parse Mass string into a Float
            try
                mass = parse(F, m.captures[2])
            catch ArgumentError
                throw(ArgumentError("Failed to process data in line $line_num:\n $(line)"))
            end
            str_xyz = m.captures[3:5]

        elseif contains(line, r2)
            m = match(r2, line)
            id = Symbol(m.captures[1])

            # Get mass from PeriodicTable
            mass = convert(F, elements[id].atomic_mass / 1u"u")
            str_xyz = m.captures[2:4]
        else
            throw(ArgumentError("Failed to process data in line $line_num:\n $(line)"))
        end

        Z = I(elements[Symbol(id)].number)

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