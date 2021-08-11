function xyz_from_zmat(zmat_string::String)
    # Match an atom id (String for atomic symbol or integer for atomic number) with one or more blank spaces on the left of it
    re_id = r"\s*(\d+|\w{1,2})"

    # Match an integer (mapping to an atom) followed by a symbol or value (bond length, angle or dihedral). e.g. "1 R" or "2 106.5"
    re_map = r"\s*(\d+)\s+?([^\s]+)"

    # Match a definition e.g. "R = 1.024"
    re_def = r"(\w+)\s*[=]\s*(\d+)"

    line_num = 0
    entries = []
    definitions = []
    for line in split(strip(zmat_string), "\n")

        line_num += 1

        if isempty(strip(line, [' ']))
            continue
        end

        # Required format for first line:  "ATOM_ID"
        if line_num == 1
            m = occursin(re_id, line) ? match(re_id, line) : throw(ArgumentError("Failed to process first ZMAT entry:\n $(line)"))
            push!(entries, m.captures...)

        # Required format for second line:  "ATOM_ID INT VALUE"
        # e.g. "H 1 R" "1 1 1.025"
        elseif line_num == 2
            re_2line = re_id*re_map
            m = occursin(re_2line, line) ? match(re_2line, line) : throw(ArgumentError("Failed to process second ZMAT entry:\n $(line)"))
            println(m)
            push!(entries, m.captures...)
        
        # Required format for third line:  "ATOM_ID INT VALUE INT VALUE"
        # e.g. "H 1 R 2 A" "1 1 1.025 2 120"
        elseif line_num == 3
            re_3line = re_id*re_map*re_map
            m = occursin(re_3line, line) ? match(re_3line, line) : throw(ArgumentError("Failed to process third ZMAT entry:\n $(line)"))
            push!(entries, m.captures...)

        # Remaining lines are expected to be "ATOM_ID INT VALUE INT VALUE INT VALUE" or "VAR = VALUE"
        else
            re_line = re_id*re_map*re_map*re_map
            if occursin(re_def, line)
                m = match(re_def, line)
                push!(definitions, m.captures...)
            elseif occursin(re_line, line)
                m = match(re_line, line)
                push!(entries, m.captures...)
            else
                throw(ArgumentError("Failed to process line $line_num on ZMAT:\n $(line)"))
            end
        end
    end

    line_format = FormatExpr("{:<3}   {}")
    out = ""
    return entries, definitions
end