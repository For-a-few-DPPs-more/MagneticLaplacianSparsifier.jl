@testset verbose = true "Magnetic vertex-edge incidence matrix" begin
    g = complete_graph(10)
    @testset "angle=0 recover oriented=$oriented incidence" for oriented in [true]
        graph = MetaGraph(g, :angle, 0.0)
        B = magnetic_incidence_matrix(graph; oriented=oriented)
        B_theo = Graphs.LinAlg.incidence_matrix(graph, eltype(B); oriented=oriented)
        @test B == B_theo
    end

    n = 100
    p = 0.5
    eta = 0.0
    rng = getRNG()

    @testset "Basic MUN: check least eigenvalue is zero if ranking trivial" begin
        meta_g = gen_graph_mun_basic(n, p, eta)
        B = magnetic_incidence(meta_g)
        L = B * B'
        @test abs(minimum(eigvals(L))) < 1e-12
    end
    # @testset "ReFactored MUN: check least eigenvalue is zero if ranking trivial" begin
    #     meta_g = gen_graph_mun(rng, n, p, eta)
    #     B = magnetic_incidence(meta_g)
    #     L = B * B'
    #     @test abs(minimum(eigvals(L))) < 1e-12
    # end

    @testset "Basic ERO: check least eigenvalue is zero if ranking trivial" begin
        meta_g = gen_graph_ero_basic(n, p, eta)
        B = magnetic_incidence(meta_g)
        L = B * B'
        @test abs(minimum(eigvals(L))) < 1e-12
    end
    # @testset "ReFactored ERO: check least eigenvalue is zero if ranking trivial" begin
    #     meta_g = gen_graph_ero(rng, n, p, eta)
    #     B = magnetic_incidence(meta_g)
    #     L = B * B'
    #     @test abs(minimum(eigvals(L))) < 1e-12
    # end
end
