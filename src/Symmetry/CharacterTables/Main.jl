
function pg_to_symels(PG)
    pg = parse_pg_str(PG)
    symels = [Symel("E", [1 0 0; 0 1 0; 0 0 1])]
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
                σds = []
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
                n = div(pg.n, 2)
                σds = generate_σd(n)
                c2s = generate_C2(pg.n)
            else
                n = pg.n
                σds = []
                c2s = generate_C2(pg.n)
            end
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n)
            σvs = generate_σv(n)
            #c2s = generate_C2(pg.n)
            symels = vcat(symels, cns, sns, σvs, σds, c2s)
        elseif pg.subfamily == "d"
            if pg.n % 2 != 0
                push!(symels, Symel("i", i()))
            end
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n * 2, true)
            σds = generate_σd(pg.n)
            c2s = generate_C2(pg.n)
            symels = vcat(symels, cns, sns, σds, c2s)
        elseif isnothing(pg.subfamily)
            cns = generate_Cn(pg.n)
            c2s = generate_C2(pg.n)
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
            Cnv_irr(pg.n)
        elseif pg.subfamily == "h"
            Cnh_irr(pg.n)
        else
            Cn_irrmat(pg.n)
        end
    elseif pg.family == "D"
        if pg.subfamily == "d"
            Dnd_irr(pg.n)
        elseif pg.subfamily == "h"
            Dnh_irr(pg.n)
        else
            Dn_irr(pg.n)
        end
    elseif pg.family == "S"
        Sn_irr(pg.n)
    else
        if pg.family == "T"
            if pg.subfamily == "h"
                println("Th")
            elseif pg.subfamily == "d"
                println("Td")
            else
                println("T")
            end
        elseif pg.family == "O"
            if pg.subfamily == "h"
                println("Oh")
            else
                println("O")
            end
        elseif pg.family == "I"
            if pg.subfamily == "h"
                println("Ih")
            else
                println("I")
            end
        else
            throw(ArgumentError("Unrecognized Point Group"))
        end
    end
end
