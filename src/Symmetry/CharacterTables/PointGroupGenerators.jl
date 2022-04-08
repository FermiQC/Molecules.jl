function generate_Cn(n)
    symels = []
    axis = [0 0 1]'
    cn_r = Cn(axis, n)
    for i = 1:n-1
        a, b = reduce(n, i)
        push!(symels, Symel("C_$a^$b", cn_r^i)) # Cns
    end
    return symels
end

function generate_Sn(n)
    generate_Sn(n, false)
end

function generate_Sn(n, S2n)
    symels = []
    axis = [0 0 1]'
    σh = σ(axis)
    cn_r = Cn(axis, n)
    if S2n # Generating improper rotations for S2n PG
        for i = 1:n-1
            if i % 2 == 0
                continue
            else
                a, b = reduce(n, i)
                if a == 2
                    continue
                else
                    push!(symels, Symel("S_$a^$b", cn_r^i * σh))
                end
            end
        end
        return symels
    end
    for i = 1:n-1
        a, b = reduce(n, i)
        if b % 2 == 0
            b += n
        end
        if a == 2
            continue
        else
            push!(symels, Symel("S_$a^$b", cn_r^i * σh)) # Sns
        end
    end
    return symels    
end

function generate_σv(n)
    symels = []
    x_axis = [1 0 0]' # Orient C2 and σv along x-axis
    rot_mat = Cn([0 0 1]', n)
    for i = 1:n
        axis = (rot_mat ^ i) * x_axis
        push!(symels, Symel("sigmav_$(i % n + 1)", σ(axis)))
    end
    return symels
end

function generate_σd(n)
    symels = []
    x_axis = [1 0 0]' # Orient C2 and σv along x-axis
    z_axis = [0 0 1]'
    rot_mat = Cn(z_axis, 2*n)
    base_axis = Cn(z_axis, 4*n)*x_axis # Rotate x-axis by Cn/2 to produce an axis for σd's
    for i = 1:n
        axis = (rot_mat ^ i) * base_axis
        push!(symels, Symel("sigmad_$(i % n + 1)", σ(axis))) # Take moduli so that 360° rotated axis is 0
    end
    return symels
end

function generate_C2(n)
    if n % 2 == 0
        nn = 2*n
    else
        nn = n
    end
    symels = []
    x_axis = [1 0 0]' # Orient C2 and σv along x-axis
    rot_mat = Cn([0 0 1]', n)
    for i = 1:n
        axis = (rot_mat ^ i) * x_axis
        push!(symels, Symel("C_2($(i % n + 1))", Cn(axis, 2)))
    end
    return symels
end

function generate_T()
    """
        Assume a tetrahedron contained in a cube, then we can easily generate
        the vectors for the rotation elements.
    """
    # Generate C3's
    symels = []
    C3_1v = normalize!([1.0 1.0 1.0])
    C3_2v = normalize!([-1.0 1.0 -1.0])
    C3_3v = normalize!([-1.0 -1.0 1.0])
    C3_4v = normalize!([1.0 -1.0 -1.0])
    C3list = (C3_1v, C3_2v, C3_3v, C3_4v)
    namelist = ("α", "β", "γ", "δ")
    for i = 1:4
        C3 = Cn(C3list[i], 3)
        C3s = C3^2
        push!(symels, Symel("C_3($(namelist[i]))",C3))
        push!(symels, Symel("C_3^2($(namelist[i]))",C3s))
    end
    # Generate C2's
    C2_x = [1.0 0.0 0.0]
    C2_y = [0.0 1.0 0.0]
    C2_z = [0.0 0.0 1.0]
    C2list = (C2_x, C2_y, C2_z)
    namelist = ("x", "y", "z")
    for i = 1:3
        C2 = Cn(C2list[i], 2)
        push!(symels, Symel("C_2($(namelist[i]))", C2))
    end
    return symels
end

function generate_Td()
    symels = generate_T()
    # σd's
    σd_1v = normalize!([1.0 1.0 0.0])
    σd_2v = normalize!([1.0 -1.0 0.0])
    σd_3v = normalize!([1.0 0.0 1.0])
    σd_4v = normalize!([1.0 0.0 -1.0])
    σd_5v = normalize!([0.0 1.0 1.0])
    σd_6v = normalize!([0.0 1.0 -1.0])
    σs = (σd_1v,σd_2v,σd_3v,σd_4v,σd_5v,σd_6v)
    namelist = ("xyp","xym","xzp","xzm","yzp","yzm")
    for i = 1:6
        σd = σ(σs[i])
        push!(symels, Symel("σd($(namelist[i]))", σd))
    end
    # S4's
    S4_1v = [1.0 0.0 0.0]
    S4_2v = [0.0 1.0 0.0]
    S4_3v = [0.0 0.0 1.0]
    S4vlist = (S4_1v, S4_2v, S4_3v)
    namelist = ("x","y","z")
    for i = 1:3
        S4 = Sn(S4vlist[i], 4)
        S43 = S4^3
        push!(symels, Symel("S_4($(namelist[i]))", S4))
        push!(symels, Symel("S_4^3($(namelist[i]))", S43))
    end
    return symels
end

function generate_Th()
    symels = generate_T()
    # i
    push!(symels, Symel("i", i()))
    # S6
    S6_1v = normalize!([1.0 1.0 1.0])
    S6_2v = normalize!([-1.0 1.0 -1.0])
    S6_3v = normalize!([-1.0 -1.0 1.0])
    S6_4v = normalize!([1.0 -1.0 -1.0])
    S6list = (S6_1v, S6_2v, S6_3v, S6_4v)
    namelist = ("α", "β", "γ", "δ")
    for i = 1:4
        S6 = Sn(S6list[i], 6)
        S65 = S6^5
        push!(symels, Symel("S_6($(namelist[i]))", S6))
        push!(symels, Symel("S_6^5($(namelist[i]))", S65))
    end
    # 3σh
    σh_xv = [1.0 0.0 0.0]
    σh_yv = [0.0 1.0 0.0]
    σh_zv = [0.0 0.0 1.0]
    σlist = (σh_xv,σh_yv,σh_zv)
    namelist = ("x", "y", "z")
    for i = 1:3
        σh = σ(σlist[i])
        push!(symels, Symel("σh($(namelist[i]))", σh))
    end
    return symels
end

function generate_O()
    symels = []
    # C4
    C4_xv = [1.0 0.0 0.0]
    C4_yv = [0.0 1.0 0.0]
    C4_zv = [0.0 0.0 1.0]
    C4list = (C4_xv, C4_yv, C4_zv)
    namelist = ("x", "y", "z")
    for i = 1:3
        C4 = Cn(C4list[i], 4)
        C42 = C4^2
        C43 = C4^3
        push!(symels, Symel("C_4($(namelist[i]))", C4))
        push!(symels, Symel("C_2($(namelist[i]))", C42))
        push!(symels, Symel("C_4^3($(namelist[i]))", C43))
    end
    # C3
    C3_1v = normalize!([1.0 1.0 1.0])
    C3_2v = normalize!([1.0 -1.0 1.0])
    C3_3v = normalize!([1.0 1.0 -1.0])
    C3_4v = normalize!([1.0 -1.0 -1.0])
    C3list = (C3_1v, C3_2v, C3_3v, C3_4v)
    namelist = ("α", "β", "γ", "δ")
    for i = 1:4
        C3 = Cn(C3list[i], 3)
        C32 = C3^2
        push!(symels, Symel("C_3($(namelist[i]))", C3))
        push!(symels, Symel("C_3^2($(namelist[i]))", C32))
    end
    # C2
    C2_1v = normalize!([1.0 0.0 1.0])
    C2_2v = normalize!([1.0 0.0 -1.0])
    C2_3v = normalize!([1.0 1.0 0.0])
    C2_4v = normalize!([1.0 -1.0 0.0])
    C2_5v = normalize!([0.0 -1.0 1.0])
    C2_6v = normalize!([0.0 -1.0 -1.0])
    C2list = (C2_1v, C2_2v, C2_3v, C2_4v, C2_5v, C2_6v)
    namelist = ("xzp", "xzm", "xyp", "xym", "yzm", "yzp")
    for i = 1:6
        C2 = Cn(C2list[i],2)
        push!(symels, Symel("C_2($(namelist[i]))", C2))
    end
    return symels
end

function generate_Oh()
    symels = generate_O()
    symels = vcat(symels, Symel("i",i()))
    # S4 and σh
    S4_xv = [1.0 0.0 0.0]
    S4_yv = [0.0 1.0 0.0]
    S4_zv = [0.0 0.0 1.0]
    S4list = (S4_xv, S4_yv, S4_zv)
    namelist = ("x", "y", "z")
    for i = 1:3
        S4 = Sn(S4list[i], 4)
        σh = σ(S4list[i])
        S43 = S4^3
        push!(symels, Symel("S_4($(namelist[i]))", S4))
        push!(symels, Symel("σh($(namelist[i]))", σh))
        push!(symels, Symel("S_4^3($(namelist[i]))", S43))
    end
    # S6
    S6_1v = normalize!([1.0 1.0 1.0])
    S6_2v = normalize!([1.0 -1.0 1.0])
    S6_3v = normalize!([1.0 1.0 -1.0])
    S6_4v = normalize!([1.0 -1.0 -1.0])
    S6list = (S6_1v, S6_2v, S6_3v, S6_4v)
    namelist = ("α", "β", "γ", "δ")
    for i = 1:4
        S6 = Sn(S6list[i], 6)
        S65 = S6^5
        push!(symels, Symel("S_6($(namelist[i]))", S6))
        push!(symels, Symel("S_6^5($(namelist[i]))", S65))
    end
    # C2
    σd_1v = normalize!([1.0 0.0 1.0])
    σd_2v = normalize!([1.0 0.0 -1.0])
    σd_3v = normalize!([1.0 1.0 0.0])
    σd_4v = normalize!([1.0 -1.0 0.0])
    σd_5v = normalize!([0.0 -1.0 1.0])
    σd_6v = normalize!([0.0 -1.0 -1.0])
    σdlist = (σd_1v, σd_2v, σd_3v,σd_4v,σd_5v,σd_6v)
    namelist = ("xzp", "xzm", "xyp", "xym", "yzm", "yzp")
    for i = 1:6
        σd = σ(σdlist[i])
        push!(symels, Symel("σd($(namelist[i]))", σd))
    end
    return symels
end

function generate_I()
    symels = []
    faces, vertices, edgecenters = generate_I_vectors()
    
    # C5 (face vectors)
    for i = 1:6
        C5 = Cn(faces[i],5)
        C52 = C5^2
        C53 = C5^3
        C54 = C5^4
        push!(symels, Symel("C_5($i)", C5))
        push!(symels, Symel("C_5^2($i)", C52))
        push!(symels, Symel("C_5^3($i)", C53))
        push!(symels, Symel("C_5^4($i)", C54))
    end
    
    # C3 (vertex vectors)
    for i = 1:10
        C3 = Cn(vertices[i],3)
        C32 = C3^2
        push!(symels, Symel("C_3($i)", C3))
        push!(symels, Symel("C_3^2($i)", C32))
    end

    # C2 (edge vectors)
    for i = 1:15
        C2 = Cn(edgecenters[i],2)
        push!(symels, Symel("C_2($i)", C2))
    end
    
    return symels
end

function generate_Ih()
    symels = generate_I()
    faces, vertices, edgecenters = generate_I_vectors()
    push!(symels, Symel("i", i()))
    # S10 (face vectors)
    for i = 1:6
        S10 = Sn(faces[i],10)
        S103 = S10^3
        S107 = S10^7
        S109 = S10^9
        push!(symels, Symel("S_10($i)", S10))
        push!(symels, Symel("S_10^3($i)", S103))
        push!(symels, Symel("S_10^7($i)", S107))
        push!(symels, Symel("S_10^9($i)", S109))
    end
    
    # S6 (vertex vectors)
    for i = 1:10
        S6 = Sn(vertices[i],6)
        S65 = S6^5
        push!(symels, Symel("S_6($i)", S6))
        push!(symels, Symel("S_6^5($i)", S65))
    end

    # σ (edge vectors)
    for i = 1:15
        σi = σ(edgecenters[i])
        push!(symels, Symel("σ($i)", σi))
    end
    
    return symels
end

function generate_I_vectors()
    ϕ = (1+(5^0.5))/2
    ϕi = 1/ϕ
    faces = [[1.0 ϕ 0.0],[1.0 -ϕ 0.0],[-1.0 ϕ 0.0],[-1.0 -ϕ 0.0],
             [0.0 1.0 ϕ],[0.0 1.0 -ϕ],[0.0 -1.0 ϕ],[0.0 -1.0 -ϕ],
             [ϕ 0.0 1.0],[-ϕ 0.0 1.0],[ϕ 0.0 -1.0],[-ϕ 0.0 -1.0]]
    vertices = [[1.0 1.0 1.0],[1.0 1.0 -1.0],[1.0 -1.0 1.0],[-1.0 1.0 1.0],
                [1.0 -1.0 -1.0],[-1.0 1.0 -1.0],[-1.0 -1.0 1.0],[-1.0 -1.0 -1.0],
                [0.0 ϕ ϕi],[0.0 ϕ -ϕi],[0.0 -ϕ ϕi],[0.0 -ϕ -ϕi],
                [ϕi 0.0 ϕ],[-ϕi 0.0 ϕ],[ϕi 0.0 -ϕ],[-ϕi 0.0 -ϕ],
                [ϕ ϕi 0.0],[ϕ -ϕi 0.0],[-ϕ ϕi 0.0],[-ϕ -ϕi 0.0]]
    l = length(vertices)
    edglen = 2*ϕi
    edgecenters = []
    for i = 1:l
        for j = i+1:l
            if isapprox(distance(vertices[i], vertices[j]), edglen)
                v = normalize(vertices[i]+vertices[j])
                for k in edgecenters
                    if isapprox(abs(dot(k,v)), 1.0)
                        @goto here
                    end
                end
                push!(edgecenters, v)
                @label here
            end
        end
    end
    for i in faces
        normalize!(i)
    end
    face_vectors = []
    for i in faces
        v = normalize(i)
        for j in face_vectors
            if isapprox(abs(dot(v,j)), 1.0)
                @goto here2
            end
        end
        push!(face_vectors, v)
        @label here2
    end
    vertex_vectors = []
    for i in vertices
        v = normalize(i)
        for j in vertex_vectors
            if isapprox(abs(dot(v,j)), 1.0)
                @goto here3
            end
        end
        push!(vertex_vectors, v)
        @label here3
    end
    return (face_vectors, vertex_vectors, edgecenters)
end

function distance(a,b)
    return (sum((a-b).^2))^0.5
end

function reduce(n, i)
    g = gcd(n, i)
    return div(n, g), div(i, g)
end

function gcd(A, B)
    # A quick implementation of the Euclid algorithm for finding the greatest common divisor
    a = max(A,B)
    b = min(A,B)
    if a == 0
        return b
    elseif b == 0
        return a
    else
        r = a % b
        gcd(b, r)
    end
end
