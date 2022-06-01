@testset "Character Tables" begin
    for i in allofem
        ct = CharacterTables.pg_to_chartab(i)
        @test ct.irreps == eval(Meta.parse(i*"cn"))
        @test isapprox(ct.characters ,eval(Meta.parse(i*"ct")))
    end
end