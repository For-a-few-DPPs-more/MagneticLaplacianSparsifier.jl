using MagneticLaplacianSparsifier
using Graphs, MetaGraphs
using Random
using Test

@testset "MagneticLaplacianSparsifier.jl" begin
    # Write your tests here.
end


@testset "RandomSpanningForests.jl" begin
    g = complete_graph(10)
    meta_g = MetaGraph(g, :angle, 0.0)

    rng = Random.default_rng()
    for e in edges(meta_g)
        θ = 2 * π * rand(rng)
        set_prop!(meta_g, e, :angle, θ)
    end

    @testset "crsf has only one cycle per component" begin
        q = 0.0
        mtsf = multi_type_spanning_forest(rng, meta_g, q).graph
        cc = [induced_subgraph(mtsf, cc)[1] for cc in connected_components(mtsf)]
        @test all(length(cycle_basis(c)) == 1 for c in cc)
    end
    # Write your tests here.
end