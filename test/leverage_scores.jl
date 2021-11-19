@testset verbose = true "Leverage scores approximation" begin
    n = 20
    p = 0.5
    eta = 0.3

    rng = getRNG()
    meta_g = gen_graph_mun(rng, n, p, eta)
    B = magnetic_incidence(meta_g; oriented=true)
    L = B * B'

    t = 100000
    @testset " for q = 0" begin
        q = 0
        emp_lev = emp_leverage_score(rng, meta_g, q, t)
        lev = leverage_score(B, q)

        relative_error = (norm(emp_lev - lev) / norm(lev))
        print("relative_error: ", relative_error, "\n")
        @test relative_error < 0.05
    end
    @testset " for q = 1" begin
        q = 1
        emp_lev = emp_leverage_score(rng, meta_g, q, t)
        lev = leverage_score(B, q)

        relative_error = (norm(emp_lev - lev) / norm(lev))
        print("relative_error: ", relative_error, "\n")
        @test relative_error < 0.05
    end
end
