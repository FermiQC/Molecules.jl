mutable struct Irrep
    name
    characters
end

function c_class_string(classes, pre, r, s)
    "Pushes class string to 'classes' for rotations, but ignores superscript if s is one"
    if s == 1
        push!(classes, pre*"_$(r)")
    else
        push!(classes, pre*"_$(r)^$(s)")
    end
end

function Cn_irr(n)
    irreps = [Irrep("A",ones(n))]
    if n % 2 == 0
        b = ones(n)
        for i = 1:n
            if i % 2 == 0
                b[i] *= -1
            end
        end
        push!(irreps, Irrep("B", b))
    end
    if 2 < n < 5
        θ = 2*pi / n
        v = zeros(n)
        for j = 0:n-1
            v[j+1] += 2*cos(j*θ)
        end
        push!(irreps, Irrep("E",v))
        
    elseif n > 2
        θ = 2*pi / n
        l = round(Int, (n - length(irreps)) / 2)
        for i = 1:l
            v = zeros(n)
            for j = 0:n-1
                v[j+1] += 2*cos(i*j*θ)
            end
            push!(irreps, Irrep("E$i",v))
        end
    end
    return irreps
end

function Cn_irrmat(n)
    names = ["A"]
    classes = ["E"]
    for c = 1:n-1
        r,s = reduce(n,c)
        c_class_string(classes, "C", r, s)
    end
    chars = ones(n)'
    if n % 2 == 0
        push!(names,"B")
        bi = ones(n)
        for i=1:n
            if i % 2 == 0
                bi[i] *= -1
            end
        end
        chars = vcat(chars, bi')
    end
    if 2 < n < 5
        push!(names, "E")
        θ = 2*pi / n
        v = zeros(n)
        for j = 0:n-1
            v[j+1] += 2*cos(j*θ)
        end
        chars = vcat(chars, v')
    elseif n ≥ 5
        θ = 2*π / n
        l = round(Int, ((n-length(names))/2))
        for i = 1:l
            push!(names, "E$i")
            v = zeros(n)
            for j = 0:n-1
                v[j+1] += 2*cos(i*j*θ)
            end
            chars = vcat(chars, v')
        end
    end
    return names, classes, chars
end

function Cnv_irr(n)
    names, classes, chars = Cn_irrmat(n)
    classes = ["E"]
    if n % 2 == 0
        for c = 1:(n>>1)-1
            r,s = reduce(n,c)
            c_class_string(classes, "2C", r, s)
        end
        push!(classes, "C_2")
        if n>>1 == 1
            push!(classes, "σv(xz)")
            push!(classes, "σd(yz)")
        else
            push!(classes, "$(n>>1)σv")
            push!(classes, "$(n>>1)σd")
        end
    else
        for c = 1:n>>1
            r,s = reduce(n,c)
            c_class_string(classes, "2C", r, s)
        end
        push!(classes, "$(n)σv")
    end
    names[1] = "A1"
    insert!(names, 2, "A2")
    chars = vcat(chars[1,:]', chars[1,:]', chars[2:end,:])
    for i = 2:n
        if i == n-i+2 || i > n-i+2
            break
        end
        #deleteat!(classes, n-i+2)
        chars = chars[:,1:end.!=n-i+2]
    end
    if n % 2 == 0
        nirr = round(Int,(n/2)+3)
        names[3] = "B1"
        insert!(names, 4, "B2")
        chars = vcat(chars[1:3,:], chars[3,:]', chars[4:end,:])
        σv = zeros(nirr)
        σd = zeros(nirr)
        σv[1:4] = [1,-1,1,-1]
        σd[1:4] = [1,-1,-1,1]
        chars = hcat(chars, σv, σd)
    else
        nirr = round(Int, (n-1)/2+2)
        σv = zeros(nirr)
        σv[1:2] = [1, -1]
        chars = hcat(chars, σv)
    end
    return names, classes, chars
end

function Cnh_irr(n)
    names, classes, cnchars = Cn_irrmat(n)
    if n % 2 == 0 
        push!(classes, "i")
        for i = 1:n-1
            if i == n>>1
                push!(classes, "σh")
            else
                r,s = reduce(n, (i+n>>1)%n)
                if s % 2 == 0
                    s += r
                end
                c_class_string(classes, "S", r, s)
            end
        end
    else
        push!(classes, "σh")
        for i = 1:n-1
            r,s = reduce(n,i)
            if i % 2 == 0
                c_class_string(classes, "S", r, s+n)
            else
                c_class_string(classes, "S", r, s)
            end
        end
    end
    if n % 2 == 0
        newnames = []
        for i = 1:length(names)
            push!(newnames, names[i]*"u")
            names[i] = names[i]*"g"
        end
        append!(names, newnames)
        cncharsi = -1 * cnchars
        top = hcat(cnchars, cnchars)
        bot = hcat(cnchars, cncharsi)
        chars = vcat(top, bot)
    else
        newnames = []
        for i = 1:length(names)
            push!(newnames, names[i]*"''")
            names[i] = names[i]*"'"
        end
        append!(names, newnames)
        cncharsi = -1 * cnchars
        top = hcat(cnchars, cnchars)
        bot = hcat(cnchars, cncharsi)
        chars = vcat(top, bot)
    end
    return names, classes, chars
end

function Sn_irr(n)
    if n % 4 == 0
        names, classes, chars = Cn_irrmat(n)
        for i = 1:n
            if i % 2 == 0
                classes[i] = "S"*classes[i][2:end]
            end
        end
    elseif n % 2 == 0
        ni = round(Int, n/2)
        names, classes, cnchars = Cn_irrmat(ni)
        classes = ["E"]
        for i = 1:n>>1-1
            r,s = reduce(n>>1, i)
            c_class_string(classes, "C", r, s)
        end
        push!(classes, "i")
        for i = 1:n>>1-1
            r,s = reduce(n, (i<<1+n>>1)%n)
            c_class_string(classes, "S", r, s)
        end
        newnames = []
        for i = 1:length(names)
            push!(newnames, names[i]*"u")
            names[i] = names[i]*"g"
        end
        append!(names, newnames)
        cncharsi = -1 * cnchars
        top = hcat(cnchars, cnchars)
        bot = hcat(cnchars, cncharsi)
        chars = vcat(top, bot)
    else
        throw(ArgumentError("Odd number n for S group"))
    end
    return names, classes, chars
end

function Dn_irr(n)
    if n == 2
        names = ["A", "B1", "B2", "B3"]
        classes = ["E", "C_2(z)", "C_2(y)", "C_2(x)"]
        chars = [1 1 1 1; 1 1 -1 -1; 1 -1 1 -1; 1 -1 -1 1]
        return names, classes, chars
    end
    names, garbage_classes, chars = Cn_irrmat(n)
    garbage_names, classes, garbage_chars = Cnv_irr(n)
    if n % 2 == 0
        classes[end-1] = classes[end-1][1]*"C_2'"
        classes[end] = classes[end][1]*"C_2''"
    else
        classes[end] = classes[end][1]*"C_2"
    end
    names[1] = "A1"
    insert!(names, 2, "A2")
    chars = vcat(chars[1,:]', chars[1,:]', chars[2:end,:])
    for i = 2:n
        if i == n-i+2 || i > n-i+2
            break
        end
        #deleteat!(classes, n-i+2)
        chars = chars[:,1:end.!=n-i+2]
    end
    if n % 2 == 0
        nirr = round(Int,(n/2)+3)
        names[3] = "B1"
        insert!(names, 4, "B2")
        chars = vcat(chars[1:3,:], chars[3,:]', chars[4:end,:])
        C2p = zeros(nirr)
        C2pp = zeros(nirr)
        C2p[1:4] = [1,-1,1,-1]
        C2pp[1:4] = [1,-1,-1,1]
        chars = hcat(chars, C2p, C2pp)
    else
        nirr = round(Int, (n-1)/2+2)
        C2p = zeros(nirr)
        C2p[1:2] = [1, -1]
        chars = hcat(chars, C2p)
    end
    return names, classes, chars
end

function Dnh_irr(n)
    names, classes, dnchars = Dn_irr(n)
    if n % 2 == 0
        push!(classes, "i")
        if n == 2
            push!(classes, "σ(xy)")
            push!(classes, "σ(xz)")
            push!(classes, "σ(yz)")
        else
            for i = 1:n>>1-1
                a = i+n>>1
                if a > n>>1
                    a = n - a
                end
                r,s = reduce(n, a)
                c_class_string(classes, "2S", r, s)
            end
            push!(classes, "σh")
            push!(classes, "$(n>>1)σv")
            push!(classes, "$(n>>1)σd")
        end
        newnames = []
        for i = 1:length(names)
            push!(newnames, names[i]*"u")
            names[i] = names[i]*"g"
        end
        append!(names, newnames)
        dncharsi = -1 * dnchars
        top = hcat(dnchars, dnchars)
        bot = hcat(dnchars, dncharsi)
        chars = vcat(top, bot)
    else
        push!(classes, "σh")
        for i = 1:n>>1
            if i % 2 == 0
                r,s = reduce(n, n-i)
            else
                r,s = reduce(n ,i)
            end
            c_class_string(classes, "2S", r, s)
        end
        push!(classes, "$(n)σv")
        newnames = []
        for i = 1:length(names)
            push!(newnames, names[i]*"''")
            names[i] = names[i]*"'"
        end
        append!(names, newnames)
        dncharsi = -1 * dnchars
        top = hcat(dnchars, dnchars)
        bot = hcat(dnchars, dncharsi)
        chars = vcat(top, bot)
    end
    return names, classes, chars
end

function Dnd_irr(n)
    if n % 2 == 0
        n2 = 2*n
        names, classes, chars = Sn_irr(n2)
        #classes = collect(1:2*n2)
        classes = classes[1:n+1]
        for i = 2:n
            classes[i] = "2"*classes[i]
        end
        push!(classes, "$(n)C_2'")
        push!(classes, "$(n)σd")
        names[1] = "A1"
        insert!(names, 2, "A2")
        chars = vcat(chars[1,:]', chars[1,:]', chars[2:end,:])
        for i = 2:n2
            if i == n2-i+2 || i > n2-i+2
                break
            end
            chars = chars[:,1:end.!=n2-i+2]
        end
        nirr = n+3
        names[3] = "B1"
        insert!(names, 4, "B2")
        chars = vcat(chars[1:3,:], chars[3,:]', chars[4:end,:])
        C2p = zeros(nirr)
        σd = zeros(nirr)
        C2p[1:4] = [1,-1,1,-1]
        σd[1:4] = [1,-1,-1,1]
        chars = hcat(chars, C2p, σd)
    else
        names, classes, dnchars = Dn_irr(n)
        push!(classes, "i")
        for i = 1:n>>1
            r,s = reduce(2*n, 2*i+n)
            if s > n
                s = 2*n-s
            end
            c_class_string(classes, "2S", r, s)
        end
        push!(classes, "$(n)σd")
        newnames = []
        for i = 1:length(names)
            push!(newnames, names[i]*"u")
            names[i] = names[i]*"g"
        end
        append!(names, newnames)
        dncharsi = -1 * dnchars
        top = hcat(dnchars, dnchars)
        bot = hcat(dnchars, dncharsi)
        chars = vcat(top, bot)
    end
    return names, classes, chars
end