@testset verbose = true "Laplacian sparsifier" begin
    @testset "magnetic incidence sparsifier" begin
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
        B_sampled = B[ind_e, :]
        @test norm(sparseB - B_sampled) < 1e-10
    end

    @testset "average sparsifier mtsf is good" begin
        n = 20
        p = 0.5
        eta = 0.3
        q = 0

        rng = getRNG()
        meta_g = gen_graph_mun(rng, n, p, eta)
        B = magnetic_incidence(meta_g; oriented=true)
        L = B' * B
        lev = leverage_score(B, q)
        avgL, _, _, _ = average_sparsifier(rng, meta_g, lev, q, 10)
        relative_error = norm(avgL - L) / norm(L)
        @test relative_error < 0.4
    end

    @testset "average sparsifier mtsf is good" begin
        n = 20
        p = 0.5
        eta = 0.1
        q = 1

        rng = getRNG()
        meta_g = gen_graph_mun(rng, n, p, eta)
        B = magnetic_incidence(meta_g; oriented=true)
        L = B' * B
        lev = leverage_score(B, q)
        avgL, _, _, _ = average_sparsifier(rng, meta_g, lev, q, 10)
        relative_error = norm(avgL - L) / norm(L)
        @test relative_error < 0.4
    end
end
