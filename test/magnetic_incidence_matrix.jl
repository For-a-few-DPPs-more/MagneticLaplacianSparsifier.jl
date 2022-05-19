@testset verbose = true "Magnetic vertex-edge incidence matrix B" begin
    @testset "angle=0 recover oriented=$oriented incidence" for oriented in [true, false]
        g = complete_graph(10)
        graph = MetaGraph(g)
        for e in edges(graph)
            set_prop!(graph, e, :angle, 0.0)
        end
        B = magnetic_incidence_matrix(graph; oriented=oriented)
        B_theo = Graphs.LinAlg.incidence_matrix(graph, eltype(B); oriented=oriented)
        B_theo = transpose(B_theo)
        @test B == B_theo
    end

    @testset "correct node order for incidence of spanning subgraph" begin
        n = 50
        p = 0.8
        g = Graphs.erdos_renyi(n, p)
        meta_g = MetaGraph(g)
        for e in edges(meta_g)
            set_prop!(meta_g, e, :angle, 0.0)
        end

        m = ne(meta_g)
        B = magnetic_incidence(meta_g)

        q = 0.5
        rng = Random.default_rng()
        mtsf = multi_type_spanning_forest(rng, meta_g, q)

        ind_e = mtsf_edge_indices(mtsf, meta_g)
        B_mtsf = magnetic_incidence(mtsf)

        @test norm(abs.(Matrix(B[ind_e, :])) - abs.(Matrix(B_mtsf))) < 1e-12
    end

    @testset "λ_min(L = B' B) ≈ 0 when trivial ranking $model" for model in [:mun, :ero]
        n = 100
        p = 0.5
        eta = 0.0
        rng = getRNG()

        meta_g = MagneticLaplacianSparsifier.erdos_renyi(rng, n, p, eta, model)
        B = magnetic_incidence(meta_g)
        L = B' * B
        @test abs(minimum(eigvals(L))) < 1e-12
    end
end
