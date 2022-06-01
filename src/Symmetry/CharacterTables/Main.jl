
function pg_to_symels(PG)
    pg = parse_pg_str(PG)
    symels = Vector{Symel}([Symel("E", [1 0 0; 0 1 0; 0 0 1])])
    σh = [1 0 0; 0 1 0; 0 0 -1]
    if pg.family == "C"
        if pg.subfamily == "h"
            push!(symels, Symel("sigmah", σh)) # sigmah
            if pg.n % 2 == 0
                push!(symels, Symel("i", i()))
            end
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n)
            symels = vcat(symels, cns, sns)
        elseif pg.subfamily == "v"
            cns = generate_Cn(pg.n)
            if pg.n % 2 == 0
                n = div(pg.n, 2)
                σds = generate_σd(n)
            else
                n = pg.n
                σds = Vector{Symel}([])
            end
            σvs = generate_σv(n)
            symels = vcat(symels, cns, σvs, σds)
        elseif pg.subfamily == "s"
            push!(symels, Symel("sigmah", σh))
        elseif pg.subfamily == "i"
            push!(symels, Symel("i", i()))
        elseif isnothing(pg.subfamily)
            cns = generate_Cn(pg.n)
            symels = vcat(symels, cns)
        else
            throw(ArgumentError("Unidentified Point Group!"))
        end
    elseif pg.family == "D"
        if pg.subfamily == "h"
            push!(symels, Symel("sigmah", σh))
            if pg.n % 2 == 0
                push!(symels, Symel("i", i()))
                n = pg.n >> 1
                σds = generate_σd(n)
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = vcat(c2ps, c2pps)
            else
                n = pg.n
                σds = Vector{Symel}([])
                c2s = generate_C2p(pg.n)
            end
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n)
            σvs = generate_σv(n)
            #c2s = generate_C2(pg.n)
            symels = vcat(symels, cns, sns, σvs, σds, c2s)
        elseif pg.subfamily == "d"
            if pg.n % 2 == 0
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = vcat(c2ps, c2pps)
            else
                c2s = generate_C2p(pg.n)
                push!(symels, Symel("i", i()))
            end
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n * 2, true)
            σds = generate_σd(pg.n)
            symels = vcat(symels, cns, sns, σds, c2s)
        elseif isnothing(pg.subfamily)
            cns = generate_Cn(pg.n)
            if pg.n % 2 == 0
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = vcat(c2ps, c2pps)
            else
                c2s = generate_C2p(pg.n)
            end
            symels = vcat(symels, cns, c2s)
        else
            throw(ArgumentError("Oh shit, the trout population"))
        end
    elseif pg.family == "S"
        if isnothing(pg.subfamily) & (pg.n % 2 == 0)
            n = div(pg.n, 2)
            if n % 2 != 0
                push!(symels, Symel("i", i()))
            end
            cns = generate_Cn(n)
            sns = generate_Sn(pg.n, true)
            symels = vcat(symels, cns, sns)
        else
            throw(ArgumentError("Oh shit, the trout population"))
        end
    else
        if pg.family == "T"
            if pg.subfamily == "h"
                Ths = generate_Th()
                symels = vcat(symels, Ths)
            elseif pg.subfamily == "d"
                Tds = generate_Td()
                symels = vcat(symels,Tds)
            else
                Ts = generate_T()
                symels = vcat(symels,Ts)
            end
        elseif pg.family == "O"
            if pg.subfamily == "h"
                Ohs = generate_Oh()
                symels = vcat(symels, Ohs)
            else
                Os = generate_O()
                symels = vcat(symels, Os)
            end
        elseif pg.family == "I"
            if pg.subfamily == "h"
                Ihs = generate_Ih()
                symels = vcat(symels, Ihs)
            else
                Is = generate_I()
                symels = vcat(symels, Is)
            end
        else
            throw(ArgumentError("Unidentified Point Group!"))
        end
    end
    return symels
end

function parse_pg_str(s)
    re = r"([A-Z]+)(\d+)?([a-z]+)?"
    m = match(re, s)
    family, n, subfamily = m.captures
    if !isnothing(n)
        n = parse(Int, n)
    end
    if !isnothing(subfamily)
        subfamily = string(subfamily)
    end
    family = string(family)
    return PG(s, family, n, subfamily)
end

function pg_to_chartab(PG)
    pg = parse_pg_str(PG)
    irreps = []
    if pg.family == "C"
        if pg.subfamily == "v"
            irreps, classes, chars = Cnv_irr(pg.n)
        elseif pg.subfamily == "h"
            irreps, classes, chars = Cnh_irr(pg.n)
        else
            irreps, classes, chars = Cn_irrmat(pg.n)
        end
    elseif pg.family == "D"
        if pg.subfamily == "d"
            irreps, classes, chars = Dnd_irr(pg.n)
        elseif pg.subfamily == "h"
            irreps, classes, chars = Dnh_irr(pg.n)
        else
            irreps, classes, chars = Dn_irr(pg.n)
        end
    elseif pg.family == "S"
        irreps, classes, chars = Sn_irr(pg.n)
    else
        cp3 = cos(pi/3)
        pr5 = 0.5*(1.0+sqrt(5.0))
        mr5 = 0.5*(1.0-sqrt(5.0))
        if pg.family == "T"
            if pg.subfamily == "h"
                irreps, classes, chars = (["Ag","Au","Eg","Eu","Tg","Tu"],
                 [1,2,3,4,5,6],
                 [1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0;
                  1.0  1.0  1.0  1.0 -1.0 -1.0 -1.0 -1.0;
                  2.0  cp3  cp3  2.0  2.0  cp3  cp3  1.0;
                  2.0  cp3  cp3  2.0 -2.0 -cp3 -cp3 -1.0;
                  3.0  0.0  0.0 -1.0  1.0  0.0  0.0 -1.0;
                  3.0  0.0  0.0 -1.0 -1.0  0.0  0.0  1.0])
            elseif pg.subfamily == "d"
                irreps, classes, chars = (["A1","A2","E","T1","T2"],
                 [1,2,3,4,5],
                 [1.0  1.0  1.0  1.0  1.0;
                  1.0  1.0  1.0 -1.0 -1.0;
                  2.0 -1.0  2.0  0.0  0.0;
                  3.0  1.0 -1.0  1.0 -1.0;
                  3.0 -1.0 -1.0 -1.0  1.0])
            else
                irreps, classes, chars = (["A","E","T"],
                 [1,2,3],
                 [1.0  1.0  1.0  1.0;
                  2.0  cp3  cp3  2.0;
                  3.0  0.0  0.0 -1.0])
            end
        elseif pg.family == "O"
            if pg.subfamily == "h"
                irreps, classes, chars = (["A1g","A2g","Eg","T1g","T2g","A1u","A2u","Eu","T1u","T2u"],
                 [1,2,3,4,5,6,7,8,9,10],
                 [1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0;
                  1.0  1.0 -1.0 -1.0  1.0  1.0 -1.0  1.0  1.0 -1.0;
                  2.0 -1.0  0.0  0.0  2.0  2.0  0.0 -1.0  2.0  0.0;
                  3.0  0.0 -1.0  1.0 -1.0  3.0  1.0  0.0 -1.0 -1.0;
                  3.0  0.0  1.0 -1.0 -1.0  3.0 -1.0  0.0 -1.0  1.0;
                  1.0  1.0  1.0  1.0  1.0 -1.0 -1.0 -1.0 -1.0 -1.0;
                  1.0  1.0 -1.0 -1.0  1.0 -1.0  1.0 -1.0 -1.0  1.0;
                  2.0 -1.0  0.0  0.0  2.0 -2.0  0.0  1.0 -2.0  0.0;
                  3.0  0.0 -1.0  1.0 -1.0 -3.0 -1.0  0.0  1.0  1.0;
                  3.0  0.0  1.0 -1.0 -1.0 -3.0  1.0  0.0  1.0 -1.0])
            else
                irreps, classes, chars = (["A1","A2","E","T1","T2"],
                 [1,2,3,4,5],
                 [1.0  1.0  1.0  1.0  1.0;
                  1.0 -1.0  1.0  1.0 -1.0;
                  2.0  0.0  2.0 -1.0  0.0;
                  3.0  1.0 -1.0  0.0 -1.0;
                  3.0 -1.0 -1.0  0.0  1.0])
            end
        elseif pg.family == "I"
            if pg.subfamily == "h"
                irreps, classes, chars = (["Ag","T1g","T2g","Gg","Hg","Au","T1u","T2u","Gu","Hu"],
                 [1,2,3,4,5,6,7,8,9,10],
                 [1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0;
                  3.0  pr5  mr5  0.0 -1.0  3.0  mr5  pr5  0.0 -1.0;
                  3.0  mr5  pr5  0.0 -1.0  3.0  pr5  mr5  0.0 -1.0;
                  4.0 -1.0 -1.0  1.0  0.0  4.0 -1.0 -1.0  1.0  0.0;
                  5.0  0.0  0.0 -1.0  1.0  5.0  0.0  0.0 -1.0  1.0;
                  1.0  1.0  1.0  1.0  1.0 -1.0 -1.0 -1.0 -1.0 -1.0;
                  3.0  pr5  mr5  0.0 -1.0 -3.0 -mr5 -pr5  0.0  1.0;
                  3.0  mr5  pr5  0.0 -1.0 -3.0 -pr5 -mr5  0.0  1.0;
                  4.0 -1.0 -1.0  1.0  0.0 -4.0  1.0  1.0 -1.0  0.0;
                  5.0  0.0  0.0 -1.0  1.0 -5.0  0.0  0.0  1.0 -1.0])
            else
                irreps, classes, chars = (["A","T1","T2","G","H"],
                 [1,2,3,4,5],
                 [1.0  1.0  1.0  1.0  1.0;
                  3.0  pr5  mr5  0.0 -1.0;
                  3.0  mr5  pr5  0.0 -1.0;
                  4.0 -1.0 -1.0  1.0  0.0;
                  5.0  0.0  0.0 -1.0  1.0])
            end
        else
            throw(ArgumentError("Unrecognized Point Group"))
        end
    end
    return Chartable(PG, irreps, classes, chars)
end

function generate_symel_to_class_map(symels::Vector{Symel}, ctab::Chartable)
    pg = parse_pg_str(ctab.name)
    display(ctab)
    for s in symels
        println(s.symbol, "  ", s.rrep)
    end
    ncls = length(ctab.classes)
    nsymel = length(symels)
    class_map = zeros(Int64, nsymel)
    class_map[1] = 1 # E is always first
    if pg.family == "C"
        if pg.subfamily == "h"
            if pg.n % 2 == 0
                class_map[4:pg.n+2] .= 2:pg.n # C_n
                class_map[3] = pg.n+1 # i
                class_map[2] = pg.n+(pg.n>>1)+1 # σh
                for i = pg.n+3:2*pg.n # S_n
                    if i > 3*(pg.n>>1)+1
                        class_map[i] = i - (pg.n>>1)
                    else
                        class_map[i] = i + (pg.n>>1)-1
                    end
                end
            else
                for i = 2:pg.n+1 # C_n
                    class_map[i] = i-1
                end
                class_map[2] = pg.n+1 # σh
                for i = pg.n+2:2*pg.n # S_n
                    class_map[i] = i
                end
            end
        elseif pg.subfamily == "v"
            # The last class is σv (and then σd if n is even), and the last symels are also these!
            cn_class_map!(class_map, pg.n, 0)
            if pg.n % 2 == 0
                class_map[end-pg.n+1:end-(pg.n>>1)] .= ncls-1
                class_map[end-(pg.n>>1)+1:end] .= ncls
            else
                class_map[end-pg.n+1:end] .= ncls
            end
        else
            class_map[2:end] .= 2:nsymel
        end
    elseif pg.family == "S"
        if pg.n % 4 == 0
            for i = 2:pg.n
                if i <= pg.n>>1
                    class_map[i] = 2*i-1
                else
                    class_map[i] = 2*(i-pg.n>>1)
                end
            end
        else
            class_map[2] = (pg.n>>1)+1 # i
            class_map[3:(pg.n>>1)+1] .= 2:pg.n>>1 # C_n
            for i = (pg.n>>1)+2:pg.n # S_n
                println(i)
                if i > (pg.n>>1)+(pg.n>>2)+1
                    class_map[i] = i-((pg.n>>1)-1)>>1
                else
                    class_map[i] = i + ((pg.n>>1)-1)>>1
                end
            end
        end
    elseif pg.family == "D"
        if pg.subfamily == "h"
            if pg.n % 2 == 0
                class_map[2] = ncls-2
                class_map[3] = (ncls>>1)+1
                cn_class_map!(class_map, pg.n, 2)
            else
                class_map[2] = (ncls>>1)+1
                cn_class_map!(class_map, pg.n, 1)
            end
        elseif pg.subfamily == "d"
        else
        end
    else
    end
    return class_map
end

function cn_class_map!(class_map, n, offset)
    for i = 2:n
        if i > (n>>1)+1
            class_map[i+offset] = n-i+2
        else
            class_map[i+offset] = i
        end
    end
end

function rotate_symels_to_mol(symels, paxis, saxis)
    φ, θ, χ = get_euler_angles(paxis, saxis)
    dc = dc_mat(φ, θ, χ)
    dct = deepcopy(dc)
    dct = transpose!(dct, dc)
    new_symels = Vector{Symel}([])
    for s in symels
        push!(new_symels, Symel(s.symbol, dct*s.rrep*dc))
    end
    return new_symels
end

function get_euler_angles(paxis, saxis)
    x = [1.0;0.0;0.0]
    y = [0.0;1.0;0.0]
    z = [0.0;0.0;1.0]
    ynew = paxis×saxis
    zxp = normalize!(z×paxis)
    if isnan(zxp[1])
        φ = 0.0
    else
        xproj = zxp⋅x
        if xproj <= 0
            φ = acos(y⋅zxp)
        else
            φ = 2*π-acos(y⋅zxp)
        end
    end
    rφ = Molecules.rotation_matrix(z,φ)
    yN = rφ*y
    xN = rφ*x
    θ = acos(z⋅paxis)
    rθ = Molecules.rotation_matrix(yN,θ)
    x2N = rθ*xN
    Z = rθ*z
    χ = acos(x2N⋅saxis)
    rχ = Molecules.rotation_matrix(Z,χ)
    X = rχ*x2N
    Y = rχ*yN
    #println("Euler Check")
    #println(Z⋅paxis)
    #println(X⋅saxis)
    #println(ynew⋅Y)
    return φ, θ, χ
end

function dc_mat(φ, θ, χ)
    sp = sin(φ)
    cp = cos(φ)
    st = sin(θ)
    ct = cos(θ)
    sc = sin(χ)
    cc = cos(χ)
    direction_cosine = [cp*ct*cc-sp*sc sp*ct*cc+cp*sc -cc*st; -cp*ct*sc-sp*cc -sp*ct*sc+cp*cc sc*st; st*cp st*sp ct]
    return direction_cosine
end

function get_atom_mapping(mol, symels)
    "symels after transformation"
    amap = zeros(Int, length(mol), length(symels))
    for (a, atom) in enumerate(mol)
        for (s, symel) in enumerate(symels)
            w = where_you_go(mol, atom, symel)
            if w !== nothing
                amap[a,s] = w
            else
                throw(ErrorException("Atom $(atom) not mapped to another atom under symel $(symel)"))
            end
        end
    end
    return amap
end

function where_you_go(mol, atom, symel)
    ratom = Atom(atom.Z, atom.mass, symel.rrep*atom.xyz)
    len = size(mol,1)
    for i = 1:len
        if isapprox(mol[i].xyz, ratom.xyz, atol=tol)
            return i
        end
    end
    return nothing
end

function get_salcs(mol, symels, ctab)
    display(ctab)
    println(typeof(symels))
end

function do_things(fn)
    mol = Molecules.parse_file(fn)
    mol = Molecules.translate(mol, Molecules.center_of_mass(mol))
    pg, paxis, saxis = Molecules.Symmetry.find_point_group(mol)
    symels = pg_to_symels(pg)
    symels = rotate_symels_to_mol(symels, paxis, saxis)
    ctab = pg_to_chartab(pg)
    class_map = generate_symel_to_class_map(symels, ctab)
    println(class_map)
    #amap = get_atom_mapping(mol, symels)
    #println(amap)
    #salcs = get_salcs(mol, symels, ctab)
end

function do_things2(pg)
    symels = pg_to_symels(pg)
    ctab = pg_to_chartab(pg)
    class_map = generate_symel_to_class_map(symels, ctab)
    println(class_map)
end