using MagneticLaplacianSparsifier
using Graphs, MetaGraphs
using Random
using LinearAlgebra
using SparseArrays

using Test

using MagneticLaplacianSparsifier: getRNG

@testset "MagneticLaplacianSparsifier.jl" begin
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
    B_sampled = B[:, ind_e]
    @test all(norm(sparseB - B_sampled) < 1e-10)
end

@testset "Magnetic vertex-edge incidence matrix" begin
    g = complete_graph(10)
    @testset "angle=0 recover oriented=$oriented incidence" for oriented in [true]
        graph = MetaGraph(g, :angle, 0.0)
        B = magnetic_incidence_matrix(graph; oriented=oriented)
        B_theo = Graphs.LinAlg.incidence_matrix(graph, eltype(B); oriented=oriented)
        @test B == B_theo
    end
end

@testset "Random spanning forests" begin
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

        flt_branches = collect(Iterators.flatten(branches))
        flt_cycles = collect(Iterators.flatten(cycles))
        nb = 0
        for i in 1:length(flt_cycles)
            nb += length(flt_cycles[i])
        end
        nb += length(roots)
        nb += length(flt_branches)
        @test nb == n_v
    end
    @testset "cholesky crsf is as sparse as possible" begin
        q = 0
        # sample a CRSF
        crsf = multi_type_spanning_forest(rng, meta_g, q)

        ind = optimal_perm(crsf)

        # construct Laplacian
        B = sp_magnetic_incidence(crsf; oriented=true)
        L = B * B'

        # find cholesky with custom permutation
        lower_factor = cholesky(L; perm=ind).L

        chol_matrix = sparse(lower_factor)
        # number of diag entries
        nb_off_diag_entries = nnz(chol_matrix) - n_v

        # upper bound on number of non-zero entries (Prop 2 in paper)
        # n + sum_{cycle c}(n_c - 3)
        cycles = get_prop(crsf, :cycle_nodes)
        flt_cycles = collect(Iterators.flatten(cycles))

        bound = 0
        for i in 1:length(flt_cycles)
            bound += length(flt_cycles[i]) - 3
        end
        bound += n_v

        @test nb_off_diag_entries <= bound
    end
end
