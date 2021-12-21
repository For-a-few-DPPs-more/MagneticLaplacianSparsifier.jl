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

function nb_upsets(meta_g, ranking_score)
    oriented = true
    upsets = 0
    for e in edges(meta_g)
        a = get_edge_prop(meta_g, e, :angle, oriented)
        d = ranking_score[src(e)] - ranking_score[dst(e)]

        if sign(a) * sign(d) < 0
            upsets += 1
        end
    end
    return upsets
end

function syncrank(L, meta_g)
    n = nv(meta_g)
    lam, v = eigs(L; nev=1, which=:SM)
    ranking_score = angular_score(v)
    p = sortperm(vec(ranking_score); rev=true)

    upsets = ne(meta_g)

    score = zeros(n, 1)
    for shift in 1:n
        shifted_order = circshift(p, shift)
        upsets_score = nb_upsets(meta_g, shifted_order)
        if upsets_score <= upsets
            upsets = upsets_score
            score = shifted_order
        end
    end
    return score
end

function cond_numbers(meta_g, q, n_tot, n_rep, rng)
    m = ne(meta_g)
    cnd_number = zeros(n_tot, 1)
    cnd_number_no_lev = zeros(n_tot, 1)
    sp_L = zeros(n_tot, 1)
    sp_L_nl = zeros(n_tot, 1)
    percent_edges = zeros(n_tot, 1)

    cnd_number_std = zeros(n_tot, 1)
    cnd_number_no_lev_std = zeros(n_tot, 1)
    sp_L_std = zeros(n_tot, 1)
    sp_L_nl_std = zeros(n_tot, 1)
    percent_edges_std = zeros(n_tot, 1)

    B = magnetic_incidence(meta_g)
    Lap = B * B'
    lev = leverage_score(B, q)

    for i in 1:n_tot
        cnd_number_tp = zeros(n_rep, 1)
        cnd_number_no_lev_tp = zeros(n_rep, 1)
        sp_L_tp = zeros(n_rep, 1)
        sp_L_nl_tp = zeros(n_rep, 1)
        percent_edges_tp = zeros(n_rep, 1)

        for j in 1:n_rep
            avgL = average_sparsifier(rng, meta_g, lev, q, i)
            avgL_no_lev = average_sparsifier(rng, meta_g, nothing, q, i)
            avgL = (avgL + avgL') / 2
            avgL_no_lev = (avgL_no_lev + avgL_no_lev') / 2

            R = cholesky(avgL + q * I).L
            R_nl = cholesky(avgL_no_lev + q * I).L

            sp_L_tp[j] = nnz(sparse(R))
            sp_L_nl_tp[j] = nnz(sparse(R_nl))

            precond_L = R \ ((Lap + q * I) / R')
            precond_L_nl = R_nl \ ((Lap + q * I) / R_nl')

            cnd_number_tp[j] = cond(precond_L)
            cnd_number_no_lev_tp[j] = cond(precond_L_nl)
            percent_edges_tp[j] = nb_of_edges(avgL) / m
        end

        cnd_number[i] = mean(cnd_number_tp)
        cnd_number_no_lev[i] = mean(cnd_number_no_lev_tp)
        sp_L[i] = mean(sp_L_tp)
        sp_L_nl[i] = mean(sp_L_nl_tp)
        percent_edges[i] = mean(percent_edges_tp)

        cnd_number_std[i] = std(cnd_number_tp)
        cnd_number_no_lev_std[i] = std(cnd_number_no_lev_tp)
        sp_L_std[i] = std(sp_L_tp)
        sp_L_nl_std[i] = std(sp_L_nl_tp)
        percent_edges_std[i] = std(percent_edges_tp)
    end

    return cnd_number,
    cnd_number_no_lev,
    sp_L,
    sp_L_nl,
    percent_edges,
    cnd_number_std,
    cnd_number_no_lev_std,
    sp_L_std,
    sp_L_nl_std,
    percent_edges_std
end
