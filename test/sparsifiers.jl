@testset verbose = true "Laplacian sparsifier" begin
    @testset "magnetic incidence sparsifier" begin
        rng = getRNG()

        n = 10
        p = 0.5
        eta = 0.3
        compGraph = gen_graph_mun(rng, n, p, eta)
        B = sp_magnetic_incidence(compGraph)

        q = 0.0
        crsf = multi_type_spanning_forest(rng, compGraph, q)
        sparseB = sp_magnetic_incidence(crsf)

        ind_e = mtsf_edge_indices(crsf, compGraph)
        B_sampled = B[ind_e, :]
        @test norm(sparseB - B_sampled) < 1e-10
    end

    @testset "average sparsifier crsf is good" begin
        n = 20
        p = 0.5
        eta = 0.3
        q = 0

        rng = getRNG()
        meta_g = gen_graph_mun(rng, n, p, eta)
        B = sp_magnetic_incidence(meta_g; oriented=true)
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
        B = sp_magnetic_incidence(meta_g; oriented=true)
        L = B' * B
        lev = leverage_score(B, q)
        avgL, _, _, _ = average_sparsifier(rng, meta_g, lev, q, 10)
        relative_error = norm(avgL - L) / norm(L)
        @test relative_error < 0.4
    end

    @testset "solver AX = B col by col" begin
        rng = getRNG()

        n = 30
        p = 0.5
        eta = 0.1
        compGraph = gen_graph_mun(rng, n, p, eta)
        B = sp_magnetic_incidence(compGraph)

        L = B' * B
        C = sparse((1.0 + 0im) .* (bitrand(n, 3)))

        X = linear_solve_matrix_system(L, C)

        X0 = Matrix(L) \ Matrix(C)

        accuracy = norm(Matrix(X) - X0)
        print("accuracy = ", accuracy)
        @test accuracy < 1e-8
    end

    @testset "sparse cholesky precond works" begin
        rng = getRNG()

        n = 100
        p = 0.5
        eta = 0.05
        q = 0

        meta_g = gen_graph_mun(rng, n, p, eta)
        B = sp_magnetic_incidence(meta_g; oriented=true)

        Lap = B' * B
        lev = leverage_score(B, q)
        avgL, _, _, _ = average_sparsifier(rng, meta_g, lev, q, 3)
        avgL = ((avgL + avgL') / 2)

        sp_pL, sp_R = sp_pcond_Lap(avgL, q, Lap)
        pL, R = pcond_Lap(Matrix(avgL), q, Matrix(Lap))
        rel_abs_diff = abs(cond_nb_pp(sp_pL) - cond_nb_pp(pL)) / cond_nb_pp(pL)
        print("relative abs difference on cond nb = ", rel_abs_diff)

        @test rel_abs_diff < 1e-2
    end
end
