using OpticLobe
using Test

@testset "OpticLobe.jl" begin
    @test sum(W) == 54483992
    @test length(intrinsictypes) == 230
end
