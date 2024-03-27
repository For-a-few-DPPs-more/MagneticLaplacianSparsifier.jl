function getRNG(seed=nothing)
    return isnothing(seed) ? Random.MersenneTwister() : Random.MersenneTwister(seed)
end
getRNG(seed::Random.AbstractRNG) = seed

consecutive_pairs(path) = partition(path, 2, 1)

function add_edges_from!(g::AbstractMetaGraph, edges)
    for e in edges
        edge = isa(e, Edge) ? e : Edge(e)
        add_edge!(g, edge)
    end
end

function add_edges_from_pairs!(edge_list,walk, ne_mtsf)
    for pair in consecutive_pairs(walk)
        edge_list[ne_mtsf + 1, 1] = pair[1]
        edge_list[ne_mtsf + 1, 2] = pair[2]
        ne_mtsf +=1
    end
end



function rem_edges_from!(g::AbstractMetaGraph, edges)
    for e in edges
        edge = isa(e, Edge) ? e : Edge(e)
        rem_edge!(g, edge)
    end
end

function set_edges_prop_from!(
    g1::AbstractMetaGraph,
    prop::Symbol,
    g2::AbstractMetaGraph,
    oriented::Bool=true,
    default::Real=1.0,
)
    for e in edges(g1)
        value = get_edge_prop(g2, e, prop, oriented, default)
        set_prop!(g1, e, prop, value)
    end
end

function get_edge_prop(
    g::AbstractMetaGraph, e::Edge, prop::Symbol, oriented::Bool=true, default::Real=1.0
)
    haskey(g.eprops, e) && haskey(g.eprops[e], prop) && return g.eprops[e][prop]
    oriented && return -get_edge_prop(g, reverse(e), prop, false, default)
    return default
end

function get_edges_prop(
    g::AbstractMetaGraph, prop::Symbol, oriented::Bool=true, default::Real=1.0
)
    return [get_edge_prop(g, e, prop, oriented, default) for e in edges(g)]
end

function cond_numbers(
    meta_g::AbstractMetaGraph,
    q::Real,
    max_batchsize::Integer,
    n_rep::Integer,
    rng::Random.AbstractRNG;
    q_system::Real=q,
    splg_methods::Vector{String}=nothing,
    weighted::Bool=false,
)::AbstractDict
    m = ne(meta_g)
    n = nv(meta_g)
    batch = n

    # magnetic Laplacian
    B = sp_magnetic_incidence(meta_g)

    e_weights = ones(m)
    W = I # weight matrix
    if weighted
        e_weights = get_edges_prop(meta_g, :e_weight, true, 1.0)
        W *= spdiagm(e_weights)
    end

    Lap = B' * W * B
    Lap = 0.5 * (Lap + Lap')

    # magnetic Laplacian condition number
    #cdL, l_min, _ = cond_nb_pp(Lap + q_system * I)

    cdL = arpack_cond(Lap + q_system * I)
    println("cond nb of Laplacian: ", cdL)

    # print least eigenvalues
    # println("least eigenvalue of Laplacian: ", real(l_min))

    if splg_methods === nothing
        splg_methods = [
            "DPP(K) unif",
            "DPP(K) JL-LS",
            "DPP(K) LS",
            "iid unif",
            "iid JL-LS",
            "iid LS",
            "ST unif",
            "ST JL-LS",
            "ST LS",
        ]
    end

    lev, lev_JL, lev_ust_JL, lev_ust = [0 for _ in 1:4]
    time_lev, time_lev_JL, time_lev_ust_JL, time_lev_ust = [0 for _ in 1:4]
    # magnetic leverage scores
    if ("iid LS" in splg_methods) || ("DPP(K) LS" in splg_methods)
        lev, time_lev = @timed leverage_score(B, q; e_weights)
    end
    # JL-estimates of magnetic leverage scores
    if ("DPP(K) JL-LS" in splg_methods) || ("iid JL-LS" in splg_methods)
        cst = 40
        lev_JL, time_lev_JL = @timed JL_lev_score_estimates(B, q; e_weights, cst)
    end
    # JL-estimates of combinatorial leverage scores
    if ("ST JL-LS" in splg_methods)
        cst = 40
        B_ust = magnetic_incidence_matrix(meta_g; oriented=true, phases=false)
        lev_ust_JL, time_lev_ust_JL = @timed JL_lev_score_estimates(
            B_ust, q; e_weights, cst
        )
    end
    # combinatorial leverage scores (non-magnetic)
    if ("ST LS" in splg_methods)
        B_ust = magnetic_incidence_matrix(meta_g; oriented=true, phases=false)
        lev_ust, time_lev_ust = @timed leverage_score(B_ust, 0; e_weights)
    end

    # initialization

    D_all = Dict()

    for method in splg_methods
        print("method: ", method)
        # initialization

        # metrics
        cnd, cnd_std = [zeros(max_batchsize) for _ in 1:2]

        # properties
        sparsity_L, sparsity_L_std = [zeros(max_batchsize) for _ in 1:2]
        timing, timing_std = [zeros(max_batchsize) for _ in 1:2]
        pc_edges, pc_edges_std = [zeros(max_batchsize) for _ in 1:2]
        cycles, cycles_std = [zeros(max_batchsize) for _ in 1:2]
        roots, roots_std = [zeros(max_batchsize) for _ in 1:2]
        weight, weight_std = [zeros(max_batchsize) for _ in 1:2]

        connected = ones(max_batchsize)

        for i in 1:max_batchsize

            # temporary arrays
            cnd_tp = zeros(n_rep)
            #
            sparsity_L_tp, time_tp, pc_edges_tp, roots_tp, cycles_tp, weight_tp  = [
                zeros(n_rep) for _ in 1:6
            ]
            connected_tp = ones(n_rep)

            for j in 1:n_rep
                sp_L = spzeros(n, n)
                time = 0
                n_cls = 0
                n_rts = 0
                isconnected = 1
                weights = zeros(i)

                if method == "DPP(K) unif"
                    # DPP(K) uniform weighting
                    vec = @timed average_sparsifier(rng, meta_g, nothing, q, i; weighted)
                    out = vec[1]
                    sp_L = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    weights = out[4]
                    isconnected = out[5]
                    time = vec[2]
                elseif method == "DPP(K) LS"
                    # DPP(K) leverage score weighting
                    vec = @timed average_sparsifier(rng, meta_g, lev, q, i; weighted)
                    out = vec[1]
                    sp_L = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    weights = out[4]
                    isconnected = out[5]
                    time = vec[2]
                elseif method == "DPP(K) JL-LS"
                    # DPP(K) leverage score weighting computed with JL sketching
                    vec = @timed average_sparsifier(rng, meta_g, lev_JL, q, i; weighted)
                    out = vec[1]
                    sp_L = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    weights = out[4]
                    isconnected = out[5]
                    time = vec[2]
                elseif method == "iid unif"
                    # iid uniform with uniform weighting
                    vec = @timed average_sparsifier_iid(
                        rng, meta_g, nothing, batch, i; weighted
                    )
                    out = vec[1]
                    sp_L = out[1]
                    isconnected = out[2]
                    time = vec[2]

                    # for avoiding laplacian matrix to have negative eigvals when
                    # graphs has many connected components
                    q_system = q_system + 1e-12
                elseif method == "iid LS"
                    # iid leverage score with leverage score weighting
                    vec = @timed average_sparsifier_iid(
                        rng, meta_g, lev, batch, i; weighted
                    )
                    out = vec[1]
                    sp_L = out[1]
                    isconnected = out[2]
                    time = vec[2]

                    # for avoiding laplacian matrix to have negative eigvals when
                    # graphs has many connected components
                    q_system = q_system + 1e-12
                elseif method == "iid JL-LS"
                    # iid leverage score with leverage score weighting
                    vec = @timed average_sparsifier_iid(
                        rng, meta_g, lev_JL, batch, i; weighted
                    )
                    out = vec[1]
                    sp_L = out[1]
                    isconnected = out[2]
                    time = vec[2]

                    # for avoiding laplacian matrix to have negative eigvals when
                    # graphs has many connected components
                    q_system = q_system + 1e-12
                elseif method == "ST unif"
                    # ST uniform weighting
                    absorbing_node = true
                    ust = true
                    q_ust = 0
                    vec = @timed average_sparsifier(
                        rng, meta_g, nothing, q_ust, i; weighted, absorbing_node, ust
                    )
                    out = vec[1]
                    sp_L = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    isconnected = out[5]
                    time = vec[2]
                elseif method == "ST LS"
                    # ST LS weighting
                    absorbing_node = true
                    ust = true
                    q_ust = 0
                    vec = @timed average_sparsifier(
                        rng, meta_g, lev_ust, q_ust, i; weighted, absorbing_node, ust
                    )
                    out = vec[1]
                    sp_L = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    isconnected = out[5]
                    time = vec[2]
                elseif method == "ST JL-LS"
                    # ST LS weighting with JL estimate
                    absorbing_node = true
                    ust = true
                    q_ust = 0
                    vec = @timed average_sparsifier(
                        rng, meta_g, lev_ust_JL, q_ust, i; weighted, absorbing_node, ust
                    )
                    out = vec[1]
                    sp_L = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    isconnected = out[5]
                    time = vec[2]
                end
                sp_L = 0.5 * (sp_L + sp_L')

                # by default q_system = q

                # sparse implementation but not fast
                #pcd_L, R = sp_pcond_Lap(sp_L, q_system, Lap) # not fast

                # pcd_L, R = pcond_Lap(sp_L, q_system, Lap)
                # pcd_L = 0.5 * (pcd_L + pcd_L')
                # cnd_tp[j], _, _, = cond_nb_pp(pcd_L) # rather fast

                cnd_tp[j] = arpack_rel_cond(Lap + q_system * I,sp_L  + q_system * I)
                sparsity_L_tp[j] = 0

                pc_edges_tp[j] = nb_of_edges(sp_L) / m
                time_tp[j] = time
                roots_tp[j] = n_rts
                cycles_tp[j] = n_cls
                weight_tp[j] = mean(weights)

                connected_tp[j] = isconnected
            end
            # computing mean and std's
            cnd[i], cnd_std[i] = mean(cnd_tp), std(cnd_tp)
            #
            sparsity_L[i], sparsity_L_std[i] = mean(sparsity_L_tp), std(sparsity_L_tp)
            pc_edges[i], pc_edges_std[i] = mean(pc_edges_tp), std(pc_edges_tp)
            timing[i], timing_std[i] = mean(time_tp), std(time_tp)
            roots[i], roots_std[i] = mean(roots_tp), std(roots_tp)
            cycles[i], cycles_std[i] = mean(cycles_tp), std(cycles_tp)
            weight[i], weight_std[i] = mean(weight_tp), std(weight_tp)
            connected[i] = mean(connected_tp)
        end
        D = Dict(
            "cdL" => cdL,
            #
            "cnd" => cnd,
            #
            "cnd_std" => cnd_std,
            #
            "sparsity_L" => sparsity_L,
            #
            "sparsity_L_std" => sparsity_L_std,
            #
            "pc_edges" => pc_edges,
            #
            "pc_edges_std" => pc_edges_std,
            #
            "timing" => timing,
            #
            "timing_std" => timing_std,
            #
            "roots" => roots,
            #
            "roots_std" => roots_std,
            #
            "cycles" => cycles,
            #
            "cycles_std" => cycles_std,
            #
            "n" => n,
            #
            "m" => m,
            #
            "time_lev" => time_lev,
            #
            "time_lev_JL" => time_lev_JL,
            #
            "time_lev_ust_JL" => time_lev_ust_JL,
            #
            "time_lev_ust" => time_lev_ust,
            #
            "connected" => connected,
            #
            "weight" => weight,
            #
            "weight_std" => weight_std,
        )
        push!(D_all, method => D)
    end

    return D_all
end

function benchmark_syncrank(
    meta_g::AbstractMetaGraph,
    planted_ranking_score::Array,
    max_batchsize::Integer,
    n_rep::Integer,
    rng::Random.AbstractRNG;
    splg_methods::Vector{String}=nothing,
    hkpv::Bool=false
)::AbstractDict
    n = nv(meta_g)
    m = ne(meta_g)
    #
    # compute ground truth
    planted_ranking = ranking_from_score(planted_ranking_score)

    #  include edge weights in meta_g: w_{uv} = 1/sqrt{d(u)d(v)}
    normalize_meta_g!(meta_g)

    # parameters for benchmarking
    batch = Int(floor(n))

    # technical parameters
    weighted = true # weighted graph is used
    singular = true # all eigenpairs are computed for more stability
    k = 10 # number of upsets in top k

    # incidence matrix
    B = sp_magnetic_incidence(meta_g)

    #######################################
    # syncrank with full magnetic Laplacian
    W = I # weight matrix
    e_weights = ones(m)
    if weighted
        e_weights = get_edges_prop(meta_g, :e_weight, true, 1.0)
        W *= spdiagm(e_weights)
    end
    L = B' * W * B

    # condL, _, _ = cond_nb_pp(L)

    # least eigenvector full Laplacian
    v, l_min = power_method_least_eigenvalue(L)

    println("least eigval of Laplacian= ", l_min)

    # recovered ranking full Laplacian
    ranking_full = syncrank(L, meta_g; singular)
    tau_full = corkendall(planted_ranking, ranking_full)
    spear_full = corspearman(planted_ranking, ranking_full)
    #######################################

    # start benchmarking
    if splg_methods === nothing
        splg_methods = [
            "DPP(K) unif",
            "DPP(K) JL-LS",
            "DPP(K) LS",
            "iid unif",
            "iid JL-LS",
            "iid LS",
            "ST unif",
            "ST JL-LS",
            "ST LS",
        ]
    end

    q = 0

    lev, lev_JL, lev_ust_JL, lev_ust = [0 for _ in 1:4]
    # magnetic leverage scores
    if ("iid LS" in splg_methods) || ("DPP(K) LS" in splg_methods)
        lev = leverage_score(B, q; e_weights)
    end
    # JL-estimates of magnetic leverage scores
    if ("DPP(K) JL-LS" in splg_methods) || ("iid JL-LS" in splg_methods)
        cst = 40
        lev_JL = JL_lev_score_estimates(B, q; e_weights, cst)
    end
    # JL-estimates of combinatorial leverage scores
    if ("ST JL-LS" in splg_methods)
        cst = 40
        B_ust = magnetic_incidence_matrix(meta_g; oriented=true, phases=false)
        lev_ust_JL = JL_lev_score_estimates(B_ust, q; e_weights, cst)
    end
    # combinatorial leverage scores (non-magnetic)
    if ("ST LS" in splg_methods)
        B_ust = magnetic_incidence_matrix(meta_g; oriented=true, phases=false)
        lev_ust = leverage_score(B_ust, 0; e_weights)
    end

    D_all = Dict()

    for method in splg_methods

        # initialization
        println("method: ", method)
        # metrics mean and stds
        err, err_std = [zeros(max_batchsize) for _ in 1:2]
        tau, tau_std = [zeros(max_batchsize) for _ in 1:2]
        spear, spear_std = [zeros(max_batchsize) for _ in 1:2]
        upsets_in_top, upsets_in_top_std = [zeros(max_batchsize) for _ in 1:2]

        # graph properties mean and stds
        pc_edges, pc_edges_std = [zeros(max_batchsize) for _ in 1:2]
        roots, roots_std = [zeros(max_batchsize) for _ in 1:2]
        cycles, cycles_std = [zeros(max_batchsize) for _ in 1:2]
        weight, weight_std = [zeros(max_batchsize) for _ in 1:2]

        connected = ones(max_batchsize)

        for i in 1:max_batchsize

            # metrics
            err_tp, tau_tp, spear_tp, upsets_in_top_tp = [zeros(n_rep) for _ in 1:4]
            # graph properties
            pc_edges_tp, roots_tp, cycles_tp, weight_tp = [zeros(n_rep) for _ in 1:4]

            connected_tp = ones(n_rep)
            # batch size (nb of spanning subgraphs)

            for j in 1:n_rep
                sp_L = spzeros(n, n)
                weights = zeros(i)

                n_cls = 0
                n_rts = 0

                if method == "DPP(K) unif"
                    # DPP(K) uniform weighting
                    vec = @timed average_sparsifier(rng, meta_g, nothing, q, i; weighted,hkpv)
                    out = vec[1]
                    sp_L = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    weights = out[4]
                    isconnected = out[5]
                    time = vec[2]

                elseif method == "DPP(K) JL-LS"
                    # DPP(K) leverage score weighting (JL approx)
                    vec = @timed average_sparsifier(rng, meta_g, lev_JL, q, i; weighted,hkpv)
                    out = vec[1]
                    sp_L = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    weights = out[4]
                    isconnected = out[5]
                    time = vec[2]

                elseif method == "DPP(K) LS"
                    # DPP(K) leverage score weighting
                    vec = @timed average_sparsifier(rng, meta_g, lev, q, i; weighted,hkpv)
                    out = vec[1]
                    sp_L = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    weights = out[4]
                    isconnected = out[5]
                    time = vec[2]

                elseif method == "iid unif"
                    # iid uniform with uniform weighting
                    vec = @timed average_sparsifier_iid(
                        rng, meta_g, nothing, batch, i; weighted
                    )
                    out = vec[1]
                    sp_L = out[1]
                    isconnected = out[2]
                    time = vec[2]
                elseif method == "iid JL-LS"
                    # iid leverage score with leverage score weighting (JL approx)
                    vec = @timed average_sparsifier_iid(
                        rng, meta_g, lev_JL, batch, i; weighted
                    )
                    out = vec[1]
                    sp_L = out[1]
                    isconnected = out[2]
                    time = vec[2]

                elseif method == "iid LS"
                    # iid leverage score with leverage score weighting
                    vec = @timed average_sparsifier_iid(
                        rng, meta_g, lev, batch, i; weighted
                    )
                    out = vec[1]
                    sp_L = out[1]
                    isconnected = out[2]
                    time = vec[2]

                elseif method == "ST unif"
                    # ST uniform weighting
                    absorbing_node = true
                    ust = true
                    q_ust = 0
                    vec = @timed average_sparsifier(
                        rng, meta_g, nothing, q_ust, i; weighted, absorbing_node, ust,hkpv
                    )
                    out = vec[1]
                    sp_L = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    isconnected = out[5]
                    time = vec[2]
                elseif method == "ST JL-LS"
                    # ST LS weighting (JL approx)
                    absorbing_node = true
                    ust = true
                    q_ust = 0
                    vec = @timed average_sparsifier(
                        rng, meta_g, lev_ust_JL, q_ust, i; weighted, absorbing_node, ust,hkpv
                    )
                    out = vec[1]
                    sp_L = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    isconnected = out[5]
                    time = vec[2]

                elseif method == "ST LS"
                    # ST LS weighting
                    absorbing_node = true
                    ust = true
                    q_ust = 0
                    vec = @timed average_sparsifier(
                        rng, meta_g, lev_ust, q_ust, i; weighted, absorbing_node, ust,hkpv
                    )
                    out = vec[1]
                    sp_L = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    isconnected = out[5]
                    time = vec[2]
                end

                # least eigenvector
                v_least, eig_least = power_method_least_eigenvalue(sp_L)
                #println("least eigenvalue of sparsifier: ", eig_least)

                sp_L = 0.5 * (sp_L + sp_L')

                ranking = syncrank(sp_L, meta_g; singular)

                # metrics
                err_tp[j] = eigenvec_dist(v, v_least)
                tau_tp[j] = corkendall(planted_ranking, ranking)
                spear_tp[j] = corspearman(planted_ranking, ranking)
                upsets_in_top_tp[j] = nb_upsets_in_top(meta_g, ranking, k)

                # graph properties
                pc_edges_tp[j] = nb_of_edges(sp_L) / m
                roots_tp[j] = n_rts
                cycles_tp[j] = n_cls
                weight_tp[j] = mean(weights)

                connected_tp[j] = isconnected
            end
            # metrics
            err[i], err_std[i] = mean(err_tp), std(err_tp)
            tau[i], tau_std[i] = mean(tau_tp), std(tau_tp)
            spear[i], spear_std[i] = mean(spear_tp), std(spear_tp)
            upsets_in_top[i], upsets_in_top_std[i] = mean(upsets_in_top_tp),
            std(upsets_in_top_tp)

            # graph properties
            pc_edges[i], pc_edges_std[i] = mean(pc_edges_tp), std(pc_edges_tp)
            roots[i], roots_std[i] = mean(roots_tp), std(roots_tp)
            cycles[i], cycles_std[i] = mean(cycles_tp), std(cycles_tp)
            weight[i], weight_std[i] = mean(weight_tp), std(weight_tp)
            connected[i] = mean(connected_tp)

            # cond number
            # cond_nb[i], cond_nb_std[i] = mean(cond_tp), std(cond_tp)
        end

        D = Dict(
            "err" => err,
            #
            "err_std" => err_std,
            #
            "tau" => tau,
            #
            "tau_std" => tau_std,
            #
            "spear" => spear,
            #
            "spear_std" => spear_std,
            #
            "pc_edges" => pc_edges,
            #
            "pc_edges_std" => pc_edges_std,
            #
            "tau_full" => tau_full,
            #
            "spear_full" => spear_full,
            #
            "upsets_in_top" => upsets_in_top,
            #
            "upsets_in_top_std" => upsets_in_top_std,
            #
            "roots" => roots,
            #
            "roots_std" => roots_std,
            #
            "cycles" => cycles,
            #
            "cycles_std" => cycles_std,
            #
            "weight" => weight,
            #
            "weight_std" => weight_std,
            #
            # "cond_nb" => cond_nb,
            # #
            # "cond_nb_std" => cond_nb_std,
            #
            # "condL" => condL,
            #
            "n" => n,
            #
            "m" => m,
            #
            "connected" => connected,
        )
        push!(D_all, method => D)
    end

    return D_all
end

# function eigenvalue_approx(
#     meta_g::Integer,
#     max_batchsize::Integer,
#     n_rep::Integer,
#     rng::Random.AbstractRNG;
#     splg_methods::Vector{String}=nothing,
# )::AbstractDict
#     n = nv(meta_g)
#     m = ne(meta_g)

#     #  include edge weights in meta_g: w_{uv} = 1/sqrt{d(u)d(v)}
#     normalize_meta_g!(meta_g)

#     # technical parameters
#     weighted = true # weighted graph is used

#     # incidence matrix
#     B = sp_magnetic_incidence(meta_g)

#     #######################################
#     # syncrank with full magnetic Laplacian
#     W = I # weight matrix
#     e_weights = ones(m)
#     if weighted
#         e_weights = get_edges_prop(meta_g, :e_weight, true, 1.0)
#         W *= diagm(e_weights)
#     end
#     L = B' * W * B

#     # leverage scores
#     q = 0
#     lev = leverage_score(B, q; e_weights)

#     # least eigenvalue full Laplacian
#     _, lambda_0 = power_method_least_eigenvalue(L)

#     # start benchmarking

#     rangebatch = 1:max_batchsize

#     if splg_methods === nothing
#         splg_methods = ["DPP(K) unif", "DPP(K) LS"]
#     end
#     D_all = Dict()

#     for method in splg_methods
#         print("method: ", method)
#         # initialization
#         lambda_sp, lambda_sp_std = zeros(max_batchsize), zeros(max_batchsize)
#         pc_edges, pc_edges_std = zeros(max_batchsize), zeros(max_batchsize)
#         cycles, cycles_std = zeros(max_batchsize), zeros(max_batchsize)
#         weight, weight_std = zeros(max_batchsize), zeros(max_batchsize)

#         for i in 1:max_batchsize

#             # initialization
#             lambda_sp_tp, cycles_tp, weight_tp, pc_edges_tp = [zeros(n_rep) for _ in 1:4]

#             t = rangebatch[i]

#             for j in 1:n_rep
#                 sp_L = spzeros(n, n)
#                 n_cls = 0
#                 n_rts = 0
#                 weights = zeros(t)
#                 if method == "DPP(K) unif"
#                     # DPP(K) uniform weighting
#                     sp_L, n_cls, n_rts, weights = average_sparsifier(
#                         rng, meta_g, nothing, q, t; weighted
#                     )

#                 elseif method == "DPP(K) LS"
#                     # DPP(K) leverage score weighting
#                     sp_L, n_cls, n_rts, weights = average_sparsifier(
#                         rng, meta_g, lev, q, t; weighted
#                     )
#                 end

#                 l = eigvals(sp_L)
#                 lambda_sp_tp[j] = real(l[1])
#                 cycles_tp[j] = n_cls
#                 weight_tp[j] = mean(weights)

#                 pc_edges_tp[j] = nb_of_edges(sp_L) / m
#             end
#             lambda_sp[i], lambda_sp_std[i] = mean(lambda_sp_tp), std(lambda_sp_tp)
#             pc_edges[i], pc_edges_std[i] = mean(pc_edges_tp), std(pc_edges_tp)
#             cycles[i], cycles_std[i] = mean(cycles_tp), std(cycles_tp)
#             weight[i], weight_std[i] = mean(weight_tp), std(weight_tp)
#         end

#         D = Dict(
#             "lambda_sp" => lambda_sp,
#             #
#             "lambda_sp_std" => lambda_sp_std,
#             #
#             "pc_edges" => pc_edges,
#             #
#             "pc_edges_std" => pc_edges_std,
#             #
#             "cycles" => cycles,
#             #
#             "cycles_std" => cycles_std,
#             #
#             "weight" => weight,
#             #
#             "weight_std" => weight_std,
#             #
#             "lambda" => lambda_0,
#             #
#             "n" => n,
#             #
#             "m" => m,
#         )
#         push!(D_all, method => D)
#     end

#     return D_all
# end

function plot_comparison_sync(
    metric::String,
    D_all::AbstractDict,
    y_limits;
    legendposition::Symbol=:topright,
    splg_methods::Vector{String}=nothing,
    check_connected=false,
)
    if splg_methods === nothing
        splg_methods = [
            "DPP(K) unif",
            "DPP(K) JL-LS",
            "DPP(K) LS",
            "iid unif",
            "iid JL-LS",
            "iid LS",
            "ST unif",
            "ST JL-LS",
            "ST LS",
        ]
    end
    D_color_method = Dict(
        "DPP(K) unif" => :blue,
        "DPP(K) JL-LS" => :skyblue,
        "DPP(K) LS" => :purple,
        "iid unif" => :cyan,
        "iid JL-LS" => :red,
        "iid LS" => :teal,
        "ST unif" => :gold4,
        "ST JL-LS" => :teal,
        "ST LS" => :maroon,
    )

    D_symbol_method = Dict(
        "DPP(K) unif" => :xcross,
        "DPP(K) JL-LS" => :plus,
        "DPP(K) LS" => :square,
        "iid unif" => :rtriangle,
        "iid JL-LS" => :utriangle,
        "iid LS" => :ltriangle,
        "ST unif" => :circle,
        "ST JL-LS" => :octagon,
        "ST LS" => :pentagon,
    )

    metric_std = metric * "_std"
    D = Dict()

    plt = Plots.plot()
    plt_ic = Plots.plot()

    for method in splg_methods
        if method == "DPP(K) unif"
            D = D_all[method]
            n = D["n"]
            m = D["m"]
            x = D["pc_edges"] * m / n
            y = D[metric]
            y_er = D[metric_std]

            Plots.plot!(
                plt,
                x,
                y;
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linewidth=2,
                markerstrokewidth=2,
                xtickfont=font(13),
                ytickfont=font(13),
                guidefont=font(13),
                legendfont=font(13),
                framestyle=:box,
                size=(500, 500),
                margins=0.1 * 2cm,
            )
            if check_connected
                connected = D["connected"]

                Plots.plot!(
                    plt_ic,
                    x,
                    connected;
                    legend=false,
                    labels=method,
                    markerstrokecolor=D_color_method[method],
                    markercolor=D_color_method[method],
                    linecolor=D_color_method[method],
                    markershape=D_symbol_method[method],
                    markersize=5,
                    xtickfontsize=13,
                    ytickfontsize=13,
                    xguidefontsize=13,
                    yguidefontsize=13,
                    legendfontsize=10,
                    linewidth=2,
                    framestyle=:box,
                    markerstrokewidth=2,
                    ylims=(-0.2, 1.2),
                    size=(500, 100),
                    ylabel="connectivity",
                    yticks=0:1:1,
                )
            end

        elseif method == "DPP(K) JL-LS"
            D = D_all[method]
            n = D["n"]
            m = D["m"]
            x = D["pc_edges"] * m / n
            y = D[metric]
            y_er = D[metric_std]

            Plots.plot!(
                plt,
                x,
                y;
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linewidth=2,
                markerstrokewidth=2,
                legend=legendposition,
            )
            if check_connected
                connected = D["connected"]

                Plots.plot!(
                    plt_ic,
                    x,
                    connected;
                    legend=false,
                    labels=method,
                    markerstrokecolor=D_color_method[method],
                    markercolor=D_color_method[method],
                    linecolor=D_color_method[method],
                    markershape=D_symbol_method[method],
                    markersize=5,
                    xtickfontsize=13,
                    ytickfontsize=13,
                    xguidefontsize=13,
                    yguidefontsize=13,
                    legendfontsize=10,
                    linewidth=2,
                    framestyle=:box,
                    markerstrokewidth=2,
                    ylims=(-0.2, 1.2),
                    size=(500, 100),
                    ylabel="connectivity",
                    yticks=0:1:1,
                )
            end

        elseif method == "DPP(K) LS"
            D = D_all[method]
            n = D["n"]
            m = D["m"]
            x = D["pc_edges"] * m / n
            y = D[metric]
            y_er = D[metric_std]

            Plots.plot!(
                plt,
                x,
                y;
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linewidth=2,
                markerstrokewidth=2,
                size=(500, 500),
                legend=legendposition,
            )

            if check_connected
                connected = D["connected"]

                Plots.plot!(
                    plt_ic,
                    x,
                    connected;
                    legend=false,
                    labels=method,
                    markerstrokecolor=D_color_method[method],
                    markercolor=D_color_method[method],
                    linecolor=D_color_method[method],
                    markershape=D_symbol_method[method],
                    markersize=5,
                    xtickfontsize=13,
                    ytickfontsize=13,
                    xguidefontsize=13,
                    yguidefontsize=13,
                    legendfontsize=10,
                    linewidth=2,
                    framestyle=:box,
                    markerstrokewidth=2,
                    ylims=(-0.2, 1.2),
                    size=(500, 100),
                    ylabel="connectivity",
                    yticks=0:1:1,
                )
            end

        elseif method == "iid unif"
            D = D_all[method]
            n = D["n"]
            m = D["m"]
            x = D["pc_edges"] * m / n
            #x_er = D["pc_edges_std"]
            n = D["n"]
            m = D["m"]

            y = D[metric]
            y_er = D[metric_std]

            Plots.plot!(
                plt,
                x,
                y;
                #xerror=x_er,
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linestyle=:dash,
                linewidth=2,
                markerstrokewidth=2,
                legend=legendposition,
            )
            if check_connected
                connected = D["connected"]

                Plots.plot!(
                    plt_ic,
                    x,
                    connected;
                    legend=false,
                    labels=method,
                    markerstrokecolor=D_color_method[method],
                    markercolor=D_color_method[method],
                    linecolor=D_color_method[method],
                    markershape=D_symbol_method[method],
                    markersize=5,
                    xtickfontsize=13,
                    ytickfontsize=13,
                    xguidefontsize=13,
                    yguidefontsize=13,
                    legendfontsize=10,
                    linewidth=2,
                    framestyle=:box,
                    markerstrokewidth=2,
                    ylims=(-0.2, 1.2),
                    size=(500, 100),
                    ylabel="connectivity",
                    yticks=0:1:1,
                    linestyle=:dash,
                )
            end

        elseif method == "iid JL-LS"
            D = D_all[method]
            n = D["n"]
            m = D["m"]
            x = D["pc_edges"] * m / n
            #x_er = D["pc_edges_std"]

            y = D[metric]
            y_er = D[metric_std]

            Plots.plot!(
                plt,
                x,
                y;
                #xerror=x_er,
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linestyle=:dash,
                linewidth=2,
                markerstrokewidth=2,
                legend=legendposition,
            )
            if check_connected
                connected = D["connected"]

                Plots.plot!(
                    plt_ic,
                    x,
                    connected;
                    legend=false,
                    labels=method,
                    markerstrokecolor=D_color_method[method],
                    markercolor=D_color_method[method],
                    linecolor=D_color_method[method],
                    markershape=D_symbol_method[method],
                    markersize=5,
                    xtickfontsize=13,
                    ytickfontsize=13,
                    xguidefontsize=13,
                    yguidefontsize=13,
                    legendfontsize=10,
                    linewidth=2,
                    framestyle=:box,
                    markerstrokewidth=2,
                    ylims=(-0.2, 1.2),
                    size=(500, 100),
                    ylabel="connectivity",
                    yticks=0:1:1,
                    linestyle=:dash,
                )
            end
        elseif method == "iid LS"
            D = D_all[method]
            n = D["n"]
            m = D["m"]
            x = D["pc_edges"] * m / n
            #x_er = D["pc_edges_std"]

            y = D[metric]
            y_er = D[metric_std]

            Plots.plot!(
                plt,
                x,
                y;
                #xerror=x_er,
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linestyle=:dash,
                linewidth=2,
                markerstrokewidth=2,
                legend=legendposition,
            )
            if check_connected
                connected = D["connected"]

                Plots.plot!(
                    plt_ic,
                    x,
                    connected;
                    legend=false,
                    labels=method,
                    markerstrokecolor=D_color_method[method],
                    markercolor=D_color_method[method],
                    linecolor=D_color_method[method],
                    markershape=D_symbol_method[method],
                    markersize=5,
                    xtickfontsize=13,
                    ytickfontsize=13,
                    xguidefontsize=13,
                    yguidefontsize=13,
                    legendfontsize=10,
                    linewidth=2,
                    framestyle=:box,
                    markerstrokewidth=2,
                    ylims=(-0.2, 1.2),
                    size=(500, 100),
                    ylabel="connectivity",
                    yticks=0:1:1,
                    linestyle=:dash,
                )
            end

        elseif method == "ST unif"
            D = D_all[method]
            n = D["n"]
            m = D["m"]

            x = D["pc_edges"] * m / n
            y = D[metric]
            y_er = D[metric_std]

            Plots.plot!(
                plt,
                x,
                y;
                #xerror=x_er,
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linewidth=2,
                markerstrokewidth=2,
                framestyle=:box,
                margins=0.1 * 2Plots.cm,
                legend=legendposition,
            )
            if check_connected
                connected = D["connected"]

                Plots.plot!(
                    plt_ic,
                    x,
                    connected;
                    legend=false,
                    labels=method,
                    markerstrokecolor=D_color_method[method],
                    markercolor=D_color_method[method],
                    linecolor=D_color_method[method],
                    markershape=D_symbol_method[method],
                    markersize=5,
                    xtickfontsize=13,
                    ytickfontsize=13,
                    xguidefontsize=13,
                    yguidefontsize=13,
                    legendfontsize=10,
                    linewidth=2,
                    framestyle=:box,
                    markerstrokewidth=2,
                    ylims=(-0.2, 1.2),
                    size=(500, 100),
                    ylabel="connectivity",
                    yticks=0:1:1,
                    linestyle=:dash,
                )
            end
        elseif method == "ST JL-LS"
            D = D_all[method]
            n = D["n"]
            m = D["m"]

            x = D["pc_edges"] * m / n
            y = D[metric]
            y_er = D[metric_std]

            Plots.plot!(
                plt,
                x,
                y;
                #xerror=x_er,
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linewidth=2,
                markerstrokewidth=2,
                framestyle=:box,
                margins=0.1 * 2Plots.cm,
                legend=legendposition,
            )
            if check_connected
                connected = D["connected"]

                Plots.plot!(
                    plt_ic,
                    x,
                    connected;
                    legend=false,
                    labels=method,
                    markerstrokecolor=D_color_method[method],
                    markercolor=D_color_method[method],
                    linecolor=D_color_method[method],
                    markershape=D_symbol_method[method],
                    markersize=5,
                    xtickfontsize=13,
                    ytickfontsize=13,
                    xguidefontsize=13,
                    yguidefontsize=13,
                    legendfontsize=10,
                    linewidth=2,
                    framestyle=:box,
                    markerstrokewidth=2,
                    ylims=(-0.2, 1.2),
                    size=(500, 100),
                    ylabel="connectivity",
                    yticks=0:1:1,
                    linestyle=:dash,
                )
            end

        elseif method == "ST LS"
            D = D_all[method]
            n = D["n"]
            m = D["m"]

            x = D["pc_edges"] * m / n
            y = D[metric]
            y_er = D[metric_std]

            Plots.plot!(
                plt,
                x,
                y;
                #xerror=x_er,
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linewidth=2,
                markerstrokewidth=2,
                framestyle=:box,
                margins=0.1 * 2Plots.cm,
                legend=legendposition,
            )
            if check_connected
                connected = D["connected"]

                Plots.plot!(
                    plt_ic,
                    x,
                    connected;
                    legend=false,
                    labels=method,
                    markerstrokecolor=D_color_method[method],
                    markercolor=D_color_method[method],
                    linecolor=D_color_method[method],
                    markershape=D_symbol_method[method],
                    markersize=5,
                    xtickfontsize=13,
                    ytickfontsize=13,
                    xguidefontsize=13,
                    yguidefontsize=13,
                    legendfontsize=10,
                    linewidth=2,
                    framestyle=:box,
                    markerstrokewidth=2,
                    ylims=(-0.2, 1.2),
                    size=(500, 100),
                    ylabel="connectivity",
                    yticks=0:1:1,
                    linestyle=:dash,
                )
            end
        end
    end
    xlabel!(plt, "number of edges over number of nodes")
    ylims!(plt, y_limits)


    if metric === "err"
        ylabel!(plt, "Eigenvector distance")
        yaxis!(plt, :log)

    elseif metric === "tau"
        n = D["n"]
        m = D["m"]
        x = D["pc_edges"] * m / n
        y = D["tau_full"] * ones(size(x))
        Plots.plot!(plt, x, y; labels="full",linewidth=2)
        ylabel!(plt, "Kendall's tau ")
        yaxis!(plt, :linear)

    elseif metric === "upsets_in_top"
        ylabel!(plt, "number of upsets in top 10 ")

    elseif metric === "spear"
        n = D["n"]
        m = D["m"]
        x = D["pc_edges"] * m / n
        y = D["spear_full"] * ones(size(x))
        Plots.plot!(plt, x, y; labels="full")
        ylabel!(plt, "Spearman")

    elseif metric === "least_eig"
        ylabel!("least eigenvalue")
        n = D["n"]
        m = D["m"]
        x = D["pc_edges"] * m / n
        y = D["exact_least_eig"] * ones(size(x))
        Plots.plot!(plt, x, y; labels="exact")

    elseif metric === "top_eig"
        ylabel!("top eigenvalue")
        n = D["n"]
        m = D["m"]
        x = D["pc_edges"] * m / n
        y = D["exact_top_eig"] * ones(size(x))
        Plots.plot!(plt, x, y; labels="exact")
    end

    yaxis!(plt_ic, :identity)
    xlabel!(plt_ic, "number of edges over number of nodes")

    if check_connected
        display(plot(plt, plt_ic; layout=Plots.grid(2, 1; heights=[0.8, 0.2])))
    else
        xlabel!(plt, "number of edges over number of nodes")
        display(plt)
    end

    return nothing
end

######

function plot_comparison_cond(
    D_all::AbstractDict, y_limits; legendposition::Symbol=:bottomright, splg_methods=nothing
)
    if splg_methods === nothing
        splg_methods = [
            "DPP(K) unif",
            "DPP(K) JL-LS",
            "DPP(K) LS",
            "iid unif",
            "iid JL-LS",
            "iid LS",
            "ST unif",
            "ST JL-LS",
            "ST LS",
        ]
    end

    D_color_method = Dict(
        "DPP(K) unif" => :blue,
        "DPP(K) JL-LS" => :skyblue,
        "DPP(K) LS" => :purple,
        "iid unif" => :cyan,
        "iid JL-LS" => :red,
        "iid LS" => :teal,
        "ST unif" => :gold4,
        "ST JL-LS" => :teal,
        "ST LS" => :maroon,
    )

    D_symbol_method = Dict(
        "DPP(K) unif" => :xcross,
        "DPP(K) JL-LS" => :plus,
        "DPP(K) LS" => :square,
        "iid unif" => :rtriangle,
        "iid JL-LS" => :utriangle,
        "iid LS" => :ltriangle,
        "ST unif" => :circle,
        "ST JL-LS" => :octagon,
        "ST LS" => :pentagon,
    )

    D = Dict()

    plt_ic = Plots.plot()
    plt = Plots.plot()

    for method in splg_methods
        if method == "DPP(K) unif"
            D = D_all[method]
            n = D["n"]
            m = D["m"]

            x = D["pc_edges"] * m / n
            y = D["cnd"]
            y_er = D["cnd_std"]

            Plots.plot!(
                plt,
                x,
                y;
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linewidth=2,
                markerstrokewidth=2,
                size=(500, 500),
            )

            connected = D["connected"]

            Plots.plot!(
                plt_ic,
                x,
                connected;
                legend=false,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                xtickfontsize=13,
                ytickfontsize=13,
                xguidefontsize=13,
                yguidefontsize=13,
                legendfontsize=10,
                linewidth=2,
                framestyle=:box,
                markerstrokewidth=2,
                ylims=(-0.2, 1.2),
                size=(500, 100),
                ylabel="connectivity",
                yticks=0:1:1,
            )

        elseif method == "DPP(K) JL-LS"
            D = D_all[method]
            n = D["n"]
            m = D["m"]

            x = D["pc_edges"] * m / n
            y = D["cnd"]
            y_er = D["cnd_std"]

            Plots.plot!(
                plt,
                x,
                y;
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linewidth=2,
                markerstrokewidth=2,
                size=(500, 500),
            )

            connected = D["connected"]

            Plots.plot!(
                plt_ic,
                x,
                connected;
                legend=false,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                xtickfontsize=13,
                ytickfontsize=13,
                xguidefontsize=13,
                yguidefontsize=13,
                legendfontsize=10,
                linewidth=2,
                framestyle=:box,
                markerstrokewidth=2,
                ylims=(-0.2, 1.2),
                ylabel="connectivity",
                yticks=0:1:1,
            )

        elseif method == "DPP(K) LS"
            D = D_all[method]

            x = D["pc_edges"] * m / n
            y = D["cnd"]
            y_er = D["cnd_std"]

            Plots.plot!(
                plt,
                x,
                y;
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linewidth=2,
                markerstrokewidth=2,
            )

            connected = D["connected"]

            Plots.plot!(
                plt_ic,
                x,
                connected;
                legend=false,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                xtickfontsize=13,
                ytickfontsize=13,
                xguidefontsize=13,
                yguidefontsize=13,
                legendfontsize=10,
                linewidth=2,
                framestyle=:box,
                markerstrokewidth=2,
                ylims=(-0.2, 1.2),
                ylabel="connectivity",
                yticks=0:1:1,
            )

        elseif method == "iid unif"
            D = D_all[method]
            n = D["n"]
            m = D["m"]

            x = D["pc_edges"] * m / n
            y = D["cnd"]
            y_er = D["cnd_std"]

            Plots.plot!(
                plt,
                x,
                y;
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linestyle=:dash,
                linewidth=2,
                markerstrokewidth=2,
            )

            connected = D["connected"]

            Plots.plot!(
                plt_ic,
                x,
                connected;
                legend=false,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                xtickfontsize=13,
                ytickfontsize=13,
                xguidefontsize=13,
                yguidefontsize=13,
                legendfontsize=10,
                linewidth=2,
                framestyle=:box,
                markerstrokewidth=2,
                ylims=(-0.2, 1.2),
                ylabel="connectivity",
                yticks=0:1:1,
                linestyle=:dash,
            )

        elseif method == "iid JL-LS"
            D = D_all[method]
            n = D["n"]
            m = D["m"]

            x = D["pc_edges"] * m / n
            y = D["cnd"]
            y_er = D["cnd_std"]

            Plots.plot!(
                plt,
                x,
                y;
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linestyle=:dash,
                linewidth=2,
                markerstrokewidth=2,
            )
            connected = D["connected"]

            Plots.plot!(
                plt_ic,
                x,
                connected;
                legend=false,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                xtickfontsize=13,
                ytickfontsize=13,
                xguidefontsize=13,
                yguidefontsize=13,
                legendfontsize=10,
                linewidth=2,
                framestyle=:box,
                markerstrokewidth=2,
                ylims=(-0.2, 1.2),
                ylabel="connectivity",
                yticks=0:1:1,
                linestyle=:dash,
            )

        elseif method == "iid LS"
            D = D_all[method]
            n = D["n"]
            m = D["m"]

            x = D["pc_edges"] * m / n
            y = D["cnd"]
            y_er = D["cnd_std"]

            Plots.plot!(
                plt,
                x,
                y;
                #yerror=y_er, # too large
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linestyle=:dash,
                linewidth=2,
                markerstrokewidth=2,
            )
            connected = D["connected"]

            Plots.plot!(
                plt_ic,
                x,
                connected;
                legend=false,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                xtickfontsize=13,
                ytickfontsize=13,
                xguidefontsize=13,
                yguidefontsize=13,
                legendfontsize=10,
                linewidth=2,
                framestyle=:box,
                markerstrokewidth=2,
                ylims=(-0.2, 1.2),
                ylabel="connectivity",
                yticks=0:1:1,
                linestyle=:dash,
            )

        elseif method == "ST unif"
            D = D_all[method]
            n = D["n"]
            m = D["m"]

            x = D["pc_edges"] * m / n
            y = D["cnd"]
            y_er = D["cnd_std"]

            Plots.plot!(
                plt,
                x,
                y;
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linewidth=2,
                markerstrokewidth=2,
            )
            connected = D["connected"]

            Plots.plot!(
                plt_ic,
                x,
                connected;
                legend=false,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                xtickfontsize=13,
                ytickfontsize=13,
                xguidefontsize=13,
                yguidefontsize=13,
                legendfontsize=10,
                linewidth=2,
                framestyle=:box,
                markerstrokewidth=2,
                ylims=(-0.2, 1.2),
                ylabel="connectivity",
                yticks=0:1:1,
            )

        elseif method == "ST JL-LS"
            D = D_all[method]
            n = D["n"]
            m = D["m"]

            x = D["pc_edges"] * m / n
            y = D["cnd"]
            y_er = D["cnd_std"]

            connected = D["connected"]

            Plots.plot!(
                plt,
                x,
                y;
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linewidth=2,
                markerstrokewidth=2,
            )
            Plots.plot!(
                plt_ic,
                x,
                connected;
                legend=false,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                xtickfont=font(13),
                ytickfont=font(13),
                guidefont=font(13),
                legendfont=font(13),
                linewidth=2,
                markerstrokewidth=2,
                ylims=(-0.2, 1.2),
                ylabel="connectivity",
                yticks=0:1:1,
            )

        elseif method == "ST LS"
            D = D_all[method]
            n = D["n"]
            m = D["m"]

            x = D["pc_edges"] * m / n
            y = D["cnd"]
            y_er = D["cnd_std"]

            Plots.plot!(
                plt,
                x,
                y;
                yerror=y_er,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=D_symbol_method[method],
                markersize=5,
                linewidth=2,
                markerstrokewidth=2,
            )
            Plots.plot!(
                plt_ic,
                x,
                connected;
                legend=false,
                labels=method,
                markerstrokecolor=D_color_method[method],
                markercolor=D_color_method[method],
                linecolor=D_color_method[method],
                markershape=:octagon,
                markersize=5,
                xtickfontsize=13,
                ytickfontsize=13,
                xguidefontsize=13,
                yguidefontsize=13,
                legendfontsize=10,
                linewidth=2,
                framestyle=:box,
                markerstrokewidth=2,
                ylims=(-0.2, 1.2),
                ylabel="connectivity",
                yticks=0:1:1,
            )
        end
    end
    x = D["pc_edges"] * m / n
    cdL = D["cdL"]
    n = D["n"]
    m = D["m"]

    y = cdL * ones(size(x))
    Plots.plot!(
        plt,
        x,
        y;
        labels="no precond.",
        ylabel="condition number",
        xtickfontsize=13,
        ytickfontsize=13,
        xguidefontsize=13,
        yguidefontsize=13,
        legendfontsize=10,
        linewidth=2,
        linecolor=:black,
        framestyle=:box,
        margins=0.1 * 2Plots.cm,
        markerstrokewidth=2,
        legend=legendposition,
    )
    ylims!(plt, y_limits)

    yaxis!(plt, :log)

    yaxis!(plt_ic, :identity)
    xlabel!(plt_ic, "number of edges over number of nodes")

    display(plot(plt, plt_ic; layout=Plots.grid(2, 1; heights=[0.8, 0.2])))
    return nothing
end

##########

"""
    flat_square_2d_grid(n, a, b)

constructs a regular grid in interval [a,b]^2

# Arguments
- `n:Integer`: total number of points sqrtn^2 with sqrtn = floor(sqrt(n))^2)
- `a:Float`: start point of interval [a,b].
- `b:Float`: end point of interval [a,b].

# Output
- `X::Array{Float64,2}`:  nx2 array with coordinates of n grid nodes with [a,b]
position (i,j) -> row = j + sqrtn (i-1) for i,j = 1, ..., sqrtn

"""
function flat_square_2d_grid(n::Integer, a::Real, b::Real)::Array{Float64,2}
    sqrtn = Int64(floor(sqrt(n)))
    X = zeros(sqrtn * sqrtn, 2)
    counter = 0
    for i in 1:sqrtn
        for j in 1:sqrtn
            counter += 1
            X[counter, 1] = a + (b - a) * (i - 1) / (sqrtn - 1)
            X[counter, 2] = a + (b - a) * (j - 1) / (sqrtn - 1)
        end
    end

    return X
end
