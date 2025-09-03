using OpticLobe
using Test
using Preferences

@testset "OpticLobe.jl" begin
    # Test synapse counts based on version
    synapse_version = Preferences.load_preference(OpticLobe, "default_synapses", "Princeton")
    total_synapses = sum(W)

    # these numbers are after deleting T1 outgoing synapses
    if synapse_version == "Princeton"
        @test total_synapses == 76891137
    elseif synapse_version == "Buhmann"  
        @test total_synapses == 54483992
    end
    
    @test length(intrinsictypes) == 230
    
    # Test cell type consistency
    @test Set(alltypes) == Set(vcat(intrinsictypes, boundarytypes, othertypes))
    @test length(intersect(intrinsictypes, boundarytypes)) == 0  # No overlap
    
    println("Tests passed with $(synapse_version) synapses ($(total_synapses) total)")
end
