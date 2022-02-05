function getRNG(seed=nothing)
    return isnothing(seed) ? Random.MersenneTwister() : Random.MersenneTwister(seed)
end
getRNG(seed::Random.AbstractRNG) = seed

consecutive_pairs(path) = partition(path, 2, 1)

function add_edges_from!(g, edges)
    for e in edges
        edge = isa(e, Edge) ? e : Edge(e)
        add_edge!(g, edge)
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

function pcond_Lap(avgL, q, Lap)
    avgL = (avgL + avgL') / 2
    R = cholesky(avgL + q * I).L
    #pd_Lap = (avgL + q * I) \ (Lap + q * I)
    pd_Lap = R \ ((Lap + q * I) / R')
    return pd_Lap, R
end

function cond_numbers(meta_g, q, n_tot, n_rep, rng)
    m = ne(meta_g)
    n = nv(meta_g)
    batch = n

    cnd_nb_dpp_unif = zeros(n_tot, 1)
    cnd_nb_dpp_lev = zeros(n_tot, 1)
    cnd_nb_iid_unif = zeros(n_tot, 1)
    cnd_nb_iid_lev = zeros(n_tot, 1)

    sp_L_dpp_unif = zeros(n_tot, 1)
    sp_L_dpp_lev = zeros(n_tot, 1)
    sp_L_iid_unif = zeros(n_tot, 1)
    sp_L_iid_lev = zeros(n_tot, 1)

    percent_edges_dpp = zeros(n_tot, 1)
    percent_edges_iid = zeros(n_tot, 1)

    cnd_nb_dpp_unif_std = zeros(n_tot, 1)
    cnd_nb_dpp_lev_std = zeros(n_tot, 1)
    cnd_nb_iid_unif_std = zeros(n_tot, 1)
    cnd_nb_iid_lev_std = zeros(n_tot, 1)

    sp_L_dpp_unif_std = zeros(n_tot, 1)
    sp_L_dpp_lev_std = zeros(n_tot, 1)
    sp_L_iid_unif_std = zeros(n_tot, 1)
    sp_L_iid_lev_std = zeros(n_tot, 1)

    percent_edges_dpp_std = zeros(n_tot, 1)

    B = magnetic_incidence(meta_g)
    Lap = B * B'
    lev = leverage_score(B, q)

    for i in 1:n_tot
        cnd_nb_dpp_unif_tp = zeros(n_rep, 1)
        cnd_nb_dpp_lev_tp = zeros(n_rep, 1)
        cnd_nb_iid_unif_tp = zeros(n_rep, 1)
        cnd_nb_iid_lev_tp = zeros(n_rep, 1)

        sp_L_dpp_unif_tp = zeros(n_rep, 1)
        sp_L_dpp_lev_tp = zeros(n_rep, 1)
        sp_L_iid_unif_tp = zeros(n_rep, 1)
        sp_L_iid_lev_tp = zeros(n_rep, 1)

        percent_edges_dpp_tp = zeros(n_rep, 1)

        for j in 1:n_rep
            # DPP uniform weighting
            L_av = average_sparsifier(rng, meta_g, nothing, q, i)
            pcd_L, R = pcond_Lap(L_av, q, Lap)
            sp_L_dpp_unif_tp[j] = nnz(sparse(R))
            cnd_nb_dpp_unif_tp[j] = cond(pcd_L)
            percent_edges_dpp_tp[j] = nb_of_edges(L_av) / m

            # DPP leverage score weighting
            L_av = average_sparsifier(rng, meta_g, lev, q, i)
            pcd_L, R = pcond_Lap(L_av, q, Lap)
            sp_L_dpp_lev_tp[j] = nnz(sparse(R))
            cnd_nb_dpp_lev_tp[j] = cond(pcd_L)

            # iid uniform weighting
            L_av = average_sparsifier_iid(rng, meta_g, nothing, batch, i)
            pcd_L, R = pcond_Lap(L_av, q, Lap)
            sp_L_iid_unif_tp[j] = nnz(sparse(R))
            cnd_nb_iid_unif_tp[j] = cond(pcd_L)

            # iid leverage score weighting
            L_av = average_sparsifier_iid(rng, meta_g, lev, batch, i)
            pcd_L, R = pcond_Lap(L_av, q, Lap)
            sp_L_iid_lev_tp[j] = nnz(sparse(R))
            cnd_nb_iid_lev_tp[j] = cond(pcd_L)

            percent_edges_iid[i] = nb_of_edges(L_av) / m
        end

        cnd_nb_dpp_unif[i] = mean(cnd_nb_dpp_unif_tp)
        cnd_nb_dpp_lev[i] = mean(cnd_nb_dpp_lev_tp)
        cnd_nb_iid_unif[i] = mean(cnd_nb_iid_unif_tp)
        cnd_nb_iid_lev[i] = mean(cnd_nb_iid_lev_tp)

        sp_L_dpp_unif[i] = mean(sp_L_dpp_unif_tp)
        sp_L_dpp_lev[i] = mean(sp_L_dpp_lev_tp)
        sp_L_iid_unif[i] = mean(sp_L_iid_unif_tp)
        sp_L_iid_lev[i] = mean(sp_L_iid_lev_tp)

        percent_edges_dpp[i] = mean(percent_edges_dpp_tp)

        cnd_nb_dpp_unif_std[i] = std(cnd_nb_dpp_unif_tp)
        cnd_nb_dpp_lev_std[i] = std(cnd_nb_dpp_lev_tp)
        cnd_nb_iid_unif_std[i] = std(cnd_nb_iid_unif_tp)
        cnd_nb_iid_lev_std[i] = std(cnd_nb_iid_lev_tp)

        sp_L_dpp_unif_std[i] = std(sp_L_dpp_unif_tp)
        sp_L_dpp_lev_std[i] = std(sp_L_dpp_lev_tp)
        sp_L_iid_unif_std[i] = std(sp_L_iid_unif_tp)
        sp_L_iid_lev_std[i] = std(sp_L_iid_lev_tp)

        percent_edges_dpp_std[i] = std(percent_edges_dpp_tp)
    end

    return Dict(
        "cnd_nb_dpp_unif" => cnd_nb_dpp_unif,
        "cnd_nb_dpp_lev" => cnd_nb_dpp_lev,
        "cnd_nb_iid_unif" => cnd_nb_iid_unif,
        "cnd_nb_iid_lev" => cnd_nb_iid_lev,
        "sp_L_dpp_unif" => sp_L_dpp_unif,
        "sp_L_dpp_lev" => sp_L_dpp_lev,
        "sp_L_iid_unif" => sp_L_iid_unif,
        "sp_L_iid_lev" => sp_L_iid_lev,
        "cnd_nb_dpp_unif_std" => cnd_nb_dpp_unif_std,
        "cnd_nb_dpp_lev_std" => cnd_nb_dpp_lev_std,
        "cnd_nb_iid_unif_std" => cnd_nb_iid_unif_std,
        "cnd_nb_iid_lev_std" => cnd_nb_iid_lev_std,
        "sp_L_dpp_unif_std" => sp_L_dpp_unif_std,
        "sp_L_dpp_lev_std" => sp_L_dpp_lev_std,
        "sp_L_iid_unif_std" => sp_L_iid_unif_std,
        "sp_L_iid_lev_std" => sp_L_iid_lev_std,
        "percent_edges_dpp" => percent_edges_dpp,
        "percent_edges_dpp_std" => percent_edges_dpp_std,
        "percent_edges_iid" => percent_edges_iid,
    )
end

function benchmark_syncrank(meta_g, planted_ranking, n_batch, n_rep, rng)
    n = nv(meta_g)
    m = ne(meta_g)

    # include edge weights in meta_g
    normalize_meta_g!(meta_g)

    # parameters for benchmarking
    batch = Int(floor(n))

    # technical parameters
    weighted = true
    singular = true

    # incidence matrix
    B = magnetic_incidence(meta_g)

    # leverage score
    q = 0
    lev = leverage_score(B, q)

    # full magnetic Laplacian
    W = I # weight matrix
    if weighted
        e_weights = get_edges_prop(meta_g, :e_weight, true, 1.0)
        W *= diagm(e_weights)
    end
    L = B * W * B'

    # normalization of L
    normalize_Lap!(L)

    # least eigenvector full Laplacian
    v = least_eigenvector(L; singular)

    # recovered ranking full Laplacian
    ranking_full = syncrank(L, meta_g; singular)
    tau_full = corkendall(planted_ranking, ranking_full)

    rangebatch = 1:n_batch

    # initialization
    err_dpp_unif = zeros(size(rangebatch))
    err_dpp_lev = zeros(size(rangebatch))
    err_iid_unif = zeros(size(rangebatch))
    err_iid_lev = zeros(size(rangebatch))

    err_dpp_unif_std = zeros(size(rangebatch))
    err_dpp_lev_std = zeros(size(rangebatch))
    err_iid_unif_std = zeros(size(rangebatch))
    err_iid_lev_std = zeros(size(rangebatch))

    tau_dpp_unif = zeros(size(rangebatch))
    tau_dpp_lev = zeros(size(rangebatch))
    tau_iid_unif = zeros(size(rangebatch))
    tau_iid_lev = zeros(size(rangebatch))

    tau_dpp_unif_std = zeros(size(rangebatch))
    tau_dpp_lev_std = zeros(size(rangebatch))
    tau_iid_unif_std = zeros(size(rangebatch))
    tau_iid_lev_std = zeros(size(rangebatch))

    percent_edges_dpp = zeros(size(rangebatch))
    percent_edges_dpp_std = zeros(size(rangebatch))

    percent_edges_iid_unif = zeros(size(rangebatch))
    percent_edges_iid_lev = zeros(size(rangebatch))

    percent_edges_iid_unif_std = zeros(size(rangebatch))
    percent_edges_iid_lev_std = zeros(size(rangebatch))

    for i in 1:length(rangebatch)
        err_dpp_unif_tp = zeros(n_rep, 1)
        err_dpp_lev_tp = zeros(n_rep, 1)
        err_iid_unif_tp = zeros(n_rep, 1)
        err_iid_lev_tp = zeros(n_rep, 1)

        tau_dpp_unif_tp = zeros(n_rep, 1)
        tau_dpp_lev_tp = zeros(n_rep, 1)
        tau_iid_unif_tp = zeros(n_rep, 1)
        tau_iid_lev_tp = zeros(n_rep, 1)

        percent_edges_dpp_tp = zeros(n_rep, 1)
        percent_edges_iid_unif_tp = zeros(n_rep, 1)
        percent_edges_iid_lev_tp = zeros(n_rep, 1)

        t = rangebatch[i]

        for j in 1:n_rep

            # DPP uniform weighting
            L_av = average_sparsifier(rng, meta_g, nothing, q, t; weighted)
            v_av = least_eigenvector(L_av; singular)
            err_dpp_unif_tp[j] = eigenvec_dist(v, v_av)

            ranking = syncrank(L_av, meta_g; singular)
            tau_dpp_unif_tp[j] = corkendall(planted_ranking, ranking)

            percent_edges_dpp_tp[j] = nb_of_edges(L_av) / m

            # DPP leverage score weighting
            L_av = average_sparsifier(rng, meta_g, lev, q, t; weighted)
            v_av = least_eigenvector(L_av; singular)
            err_dpp_lev_tp[j] = eigenvec_dist(v, v_av)

            ranking = syncrank(L_av, meta_g; singular)
            tau_dpp_lev_tp[j] = corkendall(planted_ranking, ranking)

            # iid uniform with uniform weighting
            L_av = average_sparsifier_iid(rng, meta_g, nothing, batch, t; weighted)
            v_av = least_eigenvector(L_av; singular)
            err_iid_unif_tp[j] = eigenvec_dist(v, v_av)

            ranking = syncrank(L_av, meta_g; singular)
            tau_iid_unif_tp[j] = corkendall(planted_ranking, ranking)

            percent_edges_iid_unif_tp[j] = nb_of_edges(L_av) / m

            # iid leverage score with leverage score weighting
            L_av = average_sparsifier_iid(rng, meta_g, lev, batch, t; weighted)
            v_av = least_eigenvector(L_av; singular)
            err_iid_lev_tp[j] = eigenvec_dist(v, v_av)

            ranking = syncrank(L_av, meta_g; singular)
            tau_iid_lev_tp[j] = corkendall(planted_ranking, ranking)

            percent_edges_iid_lev_tp[j] = nb_of_edges(L_av) / m
        end
        err_dpp_unif[i] = mean(err_dpp_unif_tp)
        err_dpp_lev[i] = mean(err_dpp_lev_tp)
        err_iid_unif[i] = mean(err_iid_unif_tp)
        err_iid_lev[i] = mean(err_iid_lev_tp)

        err_dpp_unif_std[i] = std(err_dpp_unif_tp)
        err_dpp_lev_std[i] = std(err_dpp_lev_tp)
        err_iid_unif_std[i] = std(err_iid_unif_tp)
        err_iid_lev_std[i] = std(err_iid_lev_tp)

        tau_dpp_unif[i] = mean(tau_dpp_unif_tp)
        tau_dpp_lev[i] = mean(tau_dpp_lev_tp)
        tau_iid_unif[i] = mean(tau_iid_unif_tp)
        tau_iid_lev[i] = mean(tau_iid_lev_tp)

        tau_dpp_unif_std[i] = std(tau_dpp_unif_tp)
        tau_dpp_lev_std[i] = std(tau_dpp_lev_tp)
        tau_iid_unif_std[i] = std(tau_iid_unif_tp)
        tau_iid_lev_std[i] = std(tau_iid_lev_tp)

        percent_edges_dpp[i] = mean(percent_edges_dpp_tp)
        percent_edges_dpp_std[i] = std(percent_edges_dpp_tp)

        percent_edges_iid_unif[i] = mean(percent_edges_iid_unif_tp)
        percent_edges_iid_unif_std[i] = std(percent_edges_iid_unif_tp)

        percent_edges_iid_lev[i] = mean(percent_edges_iid_lev_tp)
        percent_edges_iid_lev_std[i] = std(percent_edges_iid_lev_tp)
    end

    return Dict(
        "err_dpp_unif" => err_dpp_unif,
        "err_dpp_lev" => err_dpp_lev,
        "err_iid_unif" => err_iid_unif,
        "err_iid_lev" => err_iid_lev,
        #
        "err_dpp_unif_std" => err_dpp_unif_std,
        "err_dpp_lev_std" => err_dpp_lev_std,
        "err_iid_unif_std" => err_iid_unif_std,
        "err_iid_lev_std" => err_iid_lev_std,
        #
        "tau_dpp_unif" => tau_dpp_unif,
        "tau_dpp_lev" => tau_dpp_lev,
        "tau_iid_unif" => tau_iid_unif,
        "tau_iid_lev" => tau_iid_lev,
        #
        "tau_dpp_unif_std" => tau_dpp_unif_std,
        "tau_dpp_lev_std" => tau_dpp_lev_std,
        "tau_iid_unif_std" => tau_iid_unif_std,
        "tau_iid_lev_std" => tau_iid_lev_std,
        #
        "percent_edges_dpp" => percent_edges_dpp,
        "percent_edges_iid_lev" => percent_edges_iid_lev,
        "percent_edges_iid_unif" => percent_edges_iid_unif,
        #
        "percent_edges_dpp_std" => percent_edges_dpp_std,
        "percent_edges_iid_unif_std" => percent_edges_iid_unif_std,
        "percent_edges_iid_lev_std" => percent_edges_iid_lev_std,
        #
        "tau_full" => tau_full,
    )
end
