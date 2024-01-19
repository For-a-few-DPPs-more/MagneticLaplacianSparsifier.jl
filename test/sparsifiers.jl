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

    @testset "sparsifier works by using vector of edge weights" begin
        rng = Random.default_rng()

        n = 100
        p = 0.5
        eta = 0.05
        q = 0
        meta_g = gen_graph_mun(rng, n, p, eta)
        B = sp_magnetic_incidence(meta_g)
        ls = leverage_score(B, q)

        weighted = false

        n = nv(meta_g)
        m = ne(meta_g)

        temp = spzeros(m)

        L = spzeros(n, n)

        weights = zeros(nb_samples, 1)
        w_tot = 0

        for i_sample in 1:nb_samples
            mtsf = multi_type_spanning_forest(rng, meta_g, q)

            D = props(mtsf)
            w = D[:weight]
            weights[i_sample] = w

            w_tot += w
            sparseB = sp_magnetic_incidence(mtsf; oriented=true)
            ind_e = mtsf_edge_indices(mtsf, meta_g)

            diag_elements = ones(length(ind_e))

            if ls === nothing
                nb_e = length(ind_e)
                W = I / (nb_e / m)
                diag_elements = diag_elements / (nb_e / m)
            else
                W = spdiagm(1 ./ ls[ind_e])
                diag_elements = diag_elements ./ ls[ind_e]
            end

            if weighted
                e_weights = get_edges_prop(meta_g, :e_weight, true, 1.0)
                W *= spdiagm(e_weights[ind_e])
                diag_elements = diag_elements .* e_weights[ind_e]
            end

            temp[ind_e] = temp[ind_e] + w * diag_elements

            L = L + w * sparseB' * W * sparseB
        end

        L = L / w_tot

        L_av = (1 / w_tot) * B' * spdiagm(temp) * B

        @test norm(L - L_av) < 1e-10
    end
end
