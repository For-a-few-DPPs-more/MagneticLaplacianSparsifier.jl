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
    pd_Lap = (avgL + q * I) \ (Lap + q * I)
    #pd_Lap = R \ ((Lap + q * I) / R')
    return pd_Lap, R
end

function cond_numbers(meta_g, q, n_tot, n_rep, rng)
    m = ne(meta_g)
    n = nv(meta_g)
    batch = n

    cnd_number = zeros(n_tot, 1)
    cnd_number_no_lev = zeros(n_tot, 1)
    cnd_number_iid_unif = zeros(n_tot, 1)
    cnd_number_iid_lev = zeros(n_tot, 1)

    sp_L = zeros(n_tot, 1)
    sp_L_nl = zeros(n_tot, 1)
    sp_L_iid_unif = zeros(n_tot, 1)
    sp_L_iid_lev = zeros(n_tot, 1)

    percent_edges = zeros(n_tot, 1)
    percent_edges_iid = zeros(n_tot, 1)

    cnd_number_std = zeros(n_tot, 1)
    cnd_number_no_lev_std = zeros(n_tot, 1)
    cnd_number_iid_unif_std = zeros(n_tot, 1)
    cnd_number_iid_lev_std = zeros(n_tot, 1)

    sp_L_std = zeros(n_tot, 1)
    sp_L_nl_std = zeros(n_tot, 1)
    sp_L_iid_unif_std = zeros(n_tot, 1)
    sp_L_iid_lev_std = zeros(n_tot, 1)

    percent_edges_std = zeros(n_tot, 1)

    B = magnetic_incidence(meta_g)
    Lap = B * B'
    lev = leverage_score(B, q)

    for i in 1:n_tot
        cnd_number_tp = zeros(n_rep, 1)
        cnd_number_no_lev_tp = zeros(n_rep, 1)
        cnd_number_iid_unif_tp = zeros(n_rep, 1)
        cnd_number_iid_lev_tp = zeros(n_rep, 1)

        sp_L_tp = zeros(n_rep, 1)
        sp_L_nl_tp = zeros(n_rep, 1)
        sp_L_iid_unif_tp = zeros(n_rep, 1)
        sp_L_iid_lev_tp = zeros(n_rep, 1)

        percent_edges_tp = zeros(n_rep, 1)

        for j in 1:n_rep
            avgL = average_sparsifier(rng, meta_g, lev, q, i)
            avgL_no_lev = average_sparsifier(rng, meta_g, nothing, q, i)
            avgL_iid_unif = average_sparsifier_iid(rng, meta_g, nothing, batch, i)
            avgL_iid_lev = average_sparsifier_iid(rng, meta_g, lev, batch, i)

            precond_L, R = pcond_Lap(avgL, q, Lap)
            precond_L_nl, R_nl = pcond_Lap(avgL_no_lev, q, Lap)

            precond_L_iid_unif, R_iid_unif = pcond_Lap(avgL_iid_unif, q, Lap)
            precond_L_iid_lev, R_iid_lev = pcond_Lap(avgL_iid_lev, q, Lap)

            sp_L_tp[j] = nnz(sparse(R))
            sp_L_nl_tp[j] = nnz(sparse(R_nl))
            sp_L_iid_unif_tp[j] = nnz(sparse(R_iid_unif))
            sp_L_iid_lev_tp[j] = nnz(sparse(R_iid_lev))

            cnd_number_tp[j] = cond(precond_L)
            cnd_number_no_lev_tp[j] = cond(precond_L_nl)
            cnd_number_iid_unif_tp[j] = cond(precond_L_iid_unif)
            cnd_number_iid_lev_tp[j] = cond(precond_L_iid_lev)

            percent_edges_tp[j] = nb_of_edges(avgL) / m
            percent_edges_iid[i] = nb_of_edges(avgL_iid_unif) / m
        end

        cnd_number[i] = mean(cnd_number_tp)
        cnd_number_no_lev[i] = mean(cnd_number_no_lev_tp)
        cnd_number_iid_unif[i] = mean(cnd_number_iid_unif_tp)
        cnd_number_iid_lev[i] = mean(cnd_number_iid_lev_tp)

        sp_L[i] = mean(sp_L_tp)
        sp_L_nl[i] = mean(sp_L_nl_tp)
        sp_L_iid_unif[i] = mean(sp_L_iid_unif_tp)
        sp_L_iid_lev[i] = mean(sp_L_iid_lev_tp)

        percent_edges[i] = mean(percent_edges_tp)

        cnd_number_std[i] = std(cnd_number_tp)
        cnd_number_no_lev_std[i] = std(cnd_number_no_lev_tp)
        cnd_number_iid_unif_std[i] = std(cnd_number_iid_unif_tp)
        cnd_number_iid_lev_std[i] = std(cnd_number_iid_lev_tp)

        sp_L_std[i] = std(sp_L_tp)
        sp_L_nl_std[i] = std(sp_L_nl_tp)
        sp_L_iid_unif_std[i] = std(sp_L_iid_unif_tp)
        sp_L_iid_lev_std[i] = std(sp_L_iid_lev_tp)

        percent_edges_std[i] = std(percent_edges_tp)
    end

    return Dict(
        "cnd_number" => cnd_number,
        "cnd_number_no_lev" => cnd_number_no_lev,
        "cnd_number_iid_unif" => cnd_number_iid_unif,
        "cnd_number_iid_lev" => cnd_number_iid_lev,
        "sp_L" => sp_L,
        "sp_L_nl" => sp_L_nl,
        "sp_L_iid_unif" => sp_L_iid_unif,
        "sp_L_iid_lev" => sp_L_iid_lev,
        "percent_edges" => percent_edges,
        "cnd_number_std" => cnd_number_std,
        "cnd_number_no_lev_std" => cnd_number_no_lev_std,
        "cnd_number_iid_unif_std" => cnd_number_iid_unif_std,
        "cnd_number_iid_lev_std" => cnd_number_iid_lev_std,
        "sp_L_std" => sp_L_std,
        "sp_L_nl_std" => sp_L_nl_std,
        "sp_L_iid_unif_std" => sp_L_iid_unif_std,
        "sp_L_iid_lev_std" => sp_L_iid_lev_std,
        "percent_edges_std" => percent_edges_std,
        "percent_edges_iid" => percent_edges_iid,
    )
end
