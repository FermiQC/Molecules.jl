@testset "Character Tables" begin
    for i in allofem
        n, cl, ct = CharacterTables.pg_to_chartab(i)
        @test n == eval(Meta.parse(i*"cn"))
        @test isapprox(ct,eval(Meta.parse(i*"ct")))
    end
end