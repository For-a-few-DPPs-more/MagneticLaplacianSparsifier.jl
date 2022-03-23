@testset verbose = true "syncrank" begin
    @testset "cumulate angles along tree" begin
        n = 8
        mtsf = MetaGraph(Graph(n))
        t = 1
        roots = [1; 6]
        add_edge!(mtsf, 1, 2, :angle, t)
        add_edge!(mtsf, 2, 3, :angle, t)
        add_edge!(mtsf, 3, 4, :angle, t)
        add_edge!(mtsf, 3, 5, :angle, t)
        add_edge!(mtsf, 6, 7, :angle, t)
        add_edge!(mtsf, 6, 8, :angle, t)

        set_prop!(mtsf, :roots, roots)

        vectors, _ = cumulate_angles(mtsf)
        v1 = vectors[1]
        v2 = vectors[2]
        angles1 = [0; 1; 2; 3; 3; 0; 0; 0]
        u1 = zeros(ComplexF64, n, 1)
        u1[1:5] = exp.(-im * angles1[1:5]) / sqrt(5)
        u2 = zeros(ComplexF64, n, 1)
        angles2 = [0; 0; 0; 0; 0; 0; 1; 1]
        u2[6:8] = exp.(-im * angles2[6:8]) / sqrt(3)

        @test norm(v1 - u1) < 1e-14
        @test norm(v2 - u2) < 1e-14
    end
    @testset "MC estimator on toy graph" begin
        n = 15
        p = 0.9
        eta = 0.1
        rng = Random.default_rng()
        meta_g = gen_graph_mun(rng, n, p, eta)

        B = magnetic_incidence(meta_g)
        L = B * B'

        q = 2
        H_target = q * inv(L + q * I)
        n_rep = 10000

        H = zeros(ComplexF64, n, n)
        for l in 1:n_rep
            mtsf = multi_type_spanning_forest(rng, meta_g, q)
            vectors, _ = cumulate_angles(mtsf)
            for v in vectors
                H += v * v'
            end
        end
        H /= n_rep

        @test norm(H_target - H) < 0.03
    end
    @testset "projector well approximated" begin
        n = 50
        p = 0.9
        eta = 0.2

        q = 5

        # planted ranking
        planted_score = randperm(rng, n)
        planted_ranking = ranking_from_score(planted_score)

        # graph model
        type = "ERO"

        if type == "MUN"
            meta_g = gen_graph_mun(rng, n, p, eta; planted_score)
        elseif type == "ERO"
            meta_g = gen_graph_ero(rng, n, p, eta; planted_score)
        end

        mtsf = multi_type_spanning_forest(rng, meta_g, q; weighted=false)

        vectors, pc_nodes = cumulate_angles(mtsf)

        P = zeros(n, n)
        for u in vectors
            P += u * u'
        end

        B_mtsf = magnetic_incidence(mtsf)
        L_mtsf = B_mtsf * B_mtsf'

        F = eigen(L_mtsf)
        val = F.values

        V = F.vectors

        P_mtsf = zeros(n, n)
        for l in 1:n
            if val[l] < 1e-13
                v = V[:, l]
                P_mtsf += v * v'
            end
        end

        @test norm(P_mtsf - P) < 1e-10
    end
    @testset "ranking from score" begin
        score = [2 3 1]
        ranking = ranking_from_score(score)
        @test norm(vec(ranking) - vec([2 1 3])) < 1e-14
    end

    @testset "syncrank retrieves noiseless ranking" begin
        rng = Random.default_rng()
        # graph parameters
        n = 10
        p = 0.7
        eta = 0.0
        # planted ranking
        planted_score = randperm(rng, n)
        planted_ranking = ranking_from_score(planted_score)
        # graph model
        meta_g = gen_graph_mun(rng, n, p, eta; planted_score)
        # unnormalized Laplacian
        B = magnetic_incidence(meta_g)
        L = B * B'
        ranking = syncrank(L, meta_g)
        @test corkendall(ranking, planted_ranking) > 0.999
    end

    @testset "least magnetic eigenvector noiseless case" begin
        rng = Random.default_rng()

        # graph parameters
        n = 10
        p = 0.7
        eta = 0.0
        # planted ranking
        planted_score = randperm(rng, n)

        # graph model
        type = "MUN"

        if type == "MUN"
            meta_g = gen_graph_mun(rng, n, p, eta; planted_score)
        elseif type == "ERO"
            meta_g = gen_graph_ero(rng, n, p, eta; planted_score)
        end

        B = magnetic_incidence(meta_g)
        L = B * B'

        temp = planted_score * π / (n - 1)
        y = exp.(im * temp)

        @test abs(norm(L * y)) < 1e-10
    end
end
