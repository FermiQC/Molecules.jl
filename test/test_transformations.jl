import Molecules: translate, Cn, σ, Sn, i, transform, isequivalent

@testset "Transformations" begin
    @testset "Translate" begin
        atoms = [
            Atom(1, 1.008, [1.0, 2.0, 3.0])
        ]

        natoms = translate(atoms, [1.0, 2.0, 3.0])

        @test natoms[1].xyz == [0.0, 0.0, 0.0]
        end

    @testset "Operators" begin
        test = [[-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0],
                [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0],
                [-0.4999999999999998 -0.8660254037844387 0.0; 0.8660254037844387 -0.4999999999999998 0.0; 0.0 0.0 1.0], 
                [1.0 -0.0 -0.0; -0.0 1.0 -0.0; -0.0 -0.0 -1.0], 
                [-0.4999999999999998 -0.8660254037844387 0.0; 0.8660254037844387 -0.4999999999999998 0.0; 0.0 0.0 -1.0],
                Int8[-1 0 0; 0 -1 0; 0 0 -1]]
        C2 = Cn([0.0, 0.0, 1.0], 2)
        C2_deg = Molecules.rotation_matrixd([0.0, 0.0, 1.0], 180.0)
        C3 = Cn([0.0, 0.0, 1.0], 3)
        σxy = σ([0.0, 0.0, 1.0])
        S3 = Sn([0.0, 0.0, 1.0], 3)
        inv = i()
        out = [C2, C2_deg, C3, σxy, S3, inv]
        for i in 1:size(test,1)
            @test isapprox(out[i], test[i], rtol=1E-5)
        end

    end

    @testset "Transform" begin
        H2 = [Atom(1, 1.008, [-0.5, 0.0, 0.0]), Atom(1, 1.008, [0.5, 0.0, 0.0])]
        H2O = [Atom(1, 1.008, [-0.5, 0.0, 0.0]), Atom(1, 1.008, [0.5, 0.0, 0.0]), Atom(1, 1.008, [0.0, 1.0, 0.0])]
        ops = [Cn([0.0, 1.0, 0.0], 2), Cn([0.0, 1.0, 1.0], 2), Cn([1.0, 1.0, 0.0], 2), σ([1.0, 0.0, 0.0]), σ([1.0, 1.0, 0.0]), Sn([1.0, 0.0, 0.0], 5), Sn([1.0, 1.0, 0.0], 3), i()]
        H2ret = [true, true, false, true, false, true, false, true]
        H2Oret = [true, false, false, true, false, false, false, false]
        for i = 1:size(ops)[1]
            H2b = transform(H2, ops[i])
            H2Ob = transform(H2O, ops[i])
            @test isequivalent(H2, H2b) == H2ret[i]
            @test isequivalent(H2O, H2Ob) == H2Oret[i]
        end
        
    end
end
