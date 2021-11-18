using MagneticLaplacianSparsifier
using Graphs, MetaGraphs
using Random
using LinearAlgebra
using Test

using MagneticLaplacianSparsifier: getRNG

@testset "MagneticLaplacianSparsifier.jl" begin
    rng = getRNG()

    n = 10
    p = 0.5
    eta = 0.3
    compGraph = gen_graph_mun(rng, n, p, eta)
    B = magnetic_incidence(compGraph)

    q = 0.0
    crsf = multi_type_spanning_forest(rng, compGraph, q)
    sparseB = magnetic_incidence(crsf)

    ind_e = mtsf_edge_indices(crsf, compGraph)
    B_sampled = B[:, ind_e]
    @test all(norm(sparseB - B_sampled) < 1e-10)
end

@testset "Magnetic vertex-edge incidence matrix" begin
    g = complete_graph(10)
    @testset "angle=0 recover oriented=$oriented incidence" for oriented in [true]
        graph = MetaGraph(g, :angle, 0.0)
        B = magnetic_incidence_matrix(graph; oriented=oriented)
        B_theo = Graphs.LinAlg.incidence_matrix(graph, eltype(B); oriented=oriented)
        @test B == B_theo
    end
end

@testset "Random spanning forests" begin
    g = complete_graph(10)
    meta_g = MetaGraph(g, :angle, 0.0)

    rng = Random.default_rng()
    for e in edges(meta_g)
        θ = 2 * π * rand(rng)
        set_prop!(meta_g, e, :angle, θ)
    end

    @testset "crsf has only one cycle per component" begin
        # crsf = cycle rooted spanning forest
        q = 0
        csrf = multi_type_spanning_forest(rng, meta_g, q).graph
        cc = [induced_subgraph(csrf, cc)[1] for cc in connected_components(csrf)]
        @test all(length(cycle_basis(c)) == 1 for c in cc)
    end
    # Write your tests here.
end
