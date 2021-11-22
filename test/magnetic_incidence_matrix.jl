@testset verbose = true "Magnetic vertex-edge incidence matrix B" begin
    @testset "angle=0 recover oriented=$oriented incidence" for oriented in [true, false]
        g = complete_graph(10)
        graph = MetaGraph(g, :angle, 0.0)
        B = magnetic_incidence_matrix(graph; oriented=oriented)
        B_theo = Graphs.LinAlg.incidence_matrix(graph, eltype(B); oriented=oriented)
        @test B == B_theo
    end

    @testset "λ_min(L = B B') ≈ 0 when trivial ranking $model" for model in [:mun, :ero]
        n = 100
        p = 0.5
        eta = 0.0
        rng = getRNG()

        meta_g = MagneticLaplacianSparsifier.erdos_renyi(rng, n, p, eta, model)
        B = magnetic_incidence(meta_g)
        L = B * B'
        @test abs(minimum(eigvals(L))) < 1e-12
    end
end
