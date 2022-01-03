@testset verbose = true "Random spanning forests" begin
    n_v = 10
    g = complete_graph(n_v)
    meta_g = MetaGraph(g, :angle, 0.0)

    rng = getRNG()

    for e in edges(meta_g)
        θ = 2 * π * rand(rng)
        set_prop!(meta_g, e, :angle, θ)
    end

    @testset "crsf has only one cycle per component" begin
        # crsf = cycle rooted spanning forest
        q = 0
        csrf = multi_type_spanning_forest(rng, meta_g, q).graph
        cc = [induced_subgraph(csrf, cc)[1] for cc in connected_components(csrf)]
        @test all(length(cycle_basis(c)) == 1 for c in cc)
    end
    @testset "mtsf has correct number of roots, cycles, branches" begin
        # mtsf = cycle rooted spanning forest
        q = 1
        mtsf = multi_type_spanning_forest(rng, meta_g, q)
        # get the roots
        roots = get_prop(mtsf, :roots)
        # get the nodes in the cycle(s)
        cycles = get_prop(mtsf, :cycle_nodes)
        # get the branches in the (reverse) order there were sampled
        branches = get_prop(mtsf, :branches)

        println("roots ", roots)
        flt_branches = collect(Iterators.flatten(branches))
        println("flt_branches ", flt_branches)
        println("cycles ", cycles)
        nb = isempty(cycles) ? 0 : sum(length(c) for c in cycles)
        nb += length(roots)
        nb += length(flt_branches)
        @test nb == nv(mtsf)
    end
    @testset "cholesky bound on number of non-zero entries" begin
        # upper bound on number of non-zero entries (Prop 2 in paper)
        # n + sum_{cycle c}(n_c - 3)
        q = 0
        crsf = multi_type_spanning_forest(rng, meta_g, q)

        ind = optimal_perm(crsf)

        B = sp_magnetic_incidence(crsf; oriented=true)
        L = B * B'

        # find cholesky with custom permutation
        lower_factor = cholesky(L; perm=ind).L

        chol_matrix = sparse(lower_factor)
        nb_off_diag_entries = nnz(chol_matrix) - n_v

        cycles = get_prop(crsf, :cycle_nodes)
        bound = isempty(cycles) ? 0 : sum(length(c) - 3 for c in cycles)
        bound += nv(crsf)

        @test nb_off_diag_entries <= bound
    end

    @testset "random walk steps have a correct distribution" begin
        rng = Random.default_rng()

        k = 5
        x1 = randn(rng, (k, 2))
        x2 = randn(rng, (k, 2)) .+ 3
        x = [x1; x2]

        n = size(x, 1)
        g = MetaGraph(n)

        bw = 1.0
        for i in 1:n
            for j in (i + 1):n
                w = exp(-norm(x[i, :] - x[j, :])^2 / bw^2)
                if w > 1e-6
                    e = (i, j)
                    add_edge!(g, i, j, :angle, 0.0)
                    set_prop!(g, Edge(e), :e_weight, w)
                end
            end
        end

        n0 = 1
        nb_list = neighbors(g, n0)

        @testset "without absording node" begin
            count = zeros(size(nb_list))
            nb_MC = 2 * 1e6

            weighted = true

            for i in 1:nb_MC
                nb = rand_step(rng, g, n0, weighted)
                id_nb = findfirst(nb_list .== nb)
                count[id_nb] += 1
            end

            empirical_p = count / nb_MC

            p = Vector{Float64}(undef, length(nb_list))
            it = 0
            for v in nb_list
                it += 1
                e = (n0, v)
                w = get_edge_prop(g, Edge(e), :e_weight)
                p[it] = abs(w)
            end
            p = p / sum(p)

            err = norm(empirical_p - p)
            println("error: ", err)

            @test err < 1e-3
        end
        @testset "with absording node" begin
            count = zeros(length(nb_list) + 1, 1)
            nb_MC = 2 * 1e6

            weighted = true
            q = 0.1

            for i in 1:nb_MC
                isroot = step_to_root(rng, g, n0, q, weighted)
                if isroot
                    count[end] += 1
                else
                    nb = rand_step(rng, g, n0, weighted)
                    id_nb = findfirst(nb_list .== nb)
                    count[id_nb] += 1
                end
            end

            empirical_p = count / nb_MC

            p = Vector{Float64}(undef, length(nb_list) + 1)
            it = 0
            for v in nb_list
                it += 1
                e = (n0, v)
                w = get_edge_prop(g, Edge(e), :e_weight)
                p[it] = abs(w)
            end
            p[end] = q
            p = p / sum(p)

            err = norm(empirical_p - p)
            println("error: ", err)

            @test err < 1e-3
        end
    end
end
