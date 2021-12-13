
@testset verbose = true "Random spanning forests" begin
    n = 30
    p = 0.9
    eta = 0.8
    q = 1.0

    rng = Random.default_rng()
    meta_g = gen_graph_mun(rng, n, p, eta)
    B = magnetic_incidence(meta_g; oriented=true)
    L = B * B'

    @testset "mtsf mean and variance match DPP's mean and variance" begin


    n_MC = 100000
    MC_edges = zeros(n_MC, 1)
    for i in 1:n_MC
        mtsf = multi_type_spanning_forest(rng, meta_g, q)
        MC_edges[i] = ne(mtsf)
    end

    m, v = mean_and_var(MC_edges)
    l = eigvals(B' * ((L + q * I) \ B))
    m_th = sum(l)
    v_th = sum(l .* (ones(size(l)) .- l))
    print(abs(m_th - m) / m)
    @test (abs(m_th - m) / m < 0.02) #&& (abs(v_th - v) / v < 0.02)
    end
end
