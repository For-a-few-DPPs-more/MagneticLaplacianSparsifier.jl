@testset verbose = true "Leverage scores approximation (magnetic case)" begin
    n = 20
    p = 0.5
    eta = 0.3

    rng = getRNG()
    meta_g = gen_graph_mun(rng, n, p, eta)
    B = sp_magnetic_incidence(meta_g; oriented=true)
    L = B' * B

    nb_samples = 100000
    @testset " for q = 0" begin
        rng = getRNG()
        q = 0
        emp_lev = emp_leverage_score(rng, meta_g, q, nb_samples)
        lev = leverage_score(B, q)

        relative_error = (norm(emp_lev - lev) / norm(lev))
        print("relative_error: ", relative_error, "\n")
        @test relative_error < 0.05
    end
    @testset " for q = 1" begin
        rng = getRNG()
        q = 1
        emp_lev = emp_leverage_score(rng, meta_g, q, nb_samples)
        lev = leverage_score(B, q)

        relative_error = (norm(emp_lev - lev) / norm(lev))
        print("relative_error: ", relative_error, "\n")
        @test relative_error < 0.05
    end
end

@testset verbose = true "Leverage scores approximation (UST case)" begin
    n = 20
    p = 0.5
    eta = 0.0

    rng = getRNG()
    meta_g = gen_graph_mun(rng, n, p, eta)
    B_ust = magnetic_incidence_matrix(meta_g; oriented=true, phases=false)

    nb_samples = 100000
    @testset " ust " begin
        rng = getRNG()
        q = 0
        emp_lev = emp_leverage_score(
            rng, meta_g, q, nb_samples; weighted=false, absorbing_node=true, ust=true
        )

        lev = leverage_score(B_ust, q)

        relative_error = (norm(emp_lev - lev) / norm(lev))
        print("relative_error: ", relative_error, "\n")
        @test relative_error < 0.05
    end
end

@testset "leverage scores approximation with weights and q = 0.1" begin
    weighted = true
    rng = Random.default_rng()

    k = 5
    x1 = randn(rng, (k, 2))
    x2 = randn(rng, (k, 2)) .+ 3
    x = [x1; x2]

    n = size(x, 1)
    meta_g = MetaGraph(n)

    bw = 2.0
    for i in 1:n
        for j in (i + 1):n
            w = exp(-norm(x[i, :] - x[j, :])^2 / bw^2)
            if w > 1e-6
                e = (i, j)
                add_edge!(meta_g, i, j, :angle, 0.0)
                set_prop!(meta_g, Edge(e), :e_weight, w)
            end
        end
    end

    e_weights = edge_weights(meta_g)

    W = diagm(e_weights)
    B = sp_magnetic_incidence(meta_g)

    Lap = B' * W * B
    nb_samples = Int(1e5)

    @testset "LS approximation with weights and q = 0.1" begin
        q = 0.1
        lev = leverage_score(B, q; e_weights)
        emp_lev = emp_leverage_score(rng, meta_g, q, nb_samples; weighted)

        relative_error = (norm(emp_lev - lev) / norm(lev))
        print("relative_error: ", relative_error, "\n")
        @test relative_error < 0.07
    end

    @testset "LS approximation with weights and q = 1." begin
        q = 1.0
        lev = leverage_score(B, q; e_weights)
        emp_lev = emp_leverage_score(rng, meta_g, q, nb_samples; weighted)

        relative_error = (norm(emp_lev - lev) / norm(lev))
        print("relative_error: ", relative_error, "\n")
        @test relative_error < 0.07
    end

    @testset "JL approximation of LS with weights and q = 1." begin
        # requiring much too large number of cols for the JL sketching
        # for checking consistency
        cst = 10000
        q = 1.0
        lev = leverage_score(B, q; e_weights)
        lev_JL = JL_lev_score_estimates(B, q; e_weights, cst)

        relative_error = (norm(lev_JL - lev) / norm(lev))
        print("relative_error: ", relative_error, "\n")
        @test relative_error < 0.01
    end
end
