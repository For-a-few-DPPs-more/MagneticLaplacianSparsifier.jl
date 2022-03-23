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
    pd_Lap = R \ ((Lap + q * I) / R')
    return pd_Lap, R
end

function cond_numbers(meta_g, q, n_tot, n_rep, rng)
    m = ne(meta_g)
    n = nv(meta_g)
    batch = n

    weighted = false

    B = magnetic_incidence(meta_g)
    Lap = B * B'
    lev = leverage_score(B, q)

    methods = ["DPP unif", "DPP LS", "iid unif", "iid LS"]
    D_all = Dict()

    for method in methods

        # initialization
        cnd = zeros(n_tot, 1)
        sp_L = zeros(n_tot, 1)

        percent_edges = zeros(n_tot, 1)
        percent_edges_std = zeros(n_tot, 1)

        cnd_std = zeros(n_tot, 1)
        sp_L_std = zeros(n_tot, 1)

        for i in 1:n_tot

            # temporary arrays
            cnd_tp = zeros(n_rep, 1)
            sp_L_tp = zeros(n_rep, 1)
            percent_edges_tp = zeros(n_rep, 1)

            for j in 1:n_rep
                L_av = zeros(n, n)

                if method == "DPP unif"
                    # DPP uniform weighting
                    L_av = average_sparsifier(rng, meta_g, nothing, q, i; weighted)

                elseif method == "DPP LS"
                    # DPP leverage score weighting
                    L_av = average_sparsifier(rng, meta_g, lev, q, i; weighted)

                elseif method == "iid unif"
                    # iid uniform with uniform weighting
                    L_av = average_sparsifier_iid(rng, meta_g, nothing, batch, i; weighted)

                elseif method == "iid LS"
                    # iid leverage score with leverage score weighting
                    L_av = average_sparsifier_iid(rng, meta_g, lev, batch, i; weighted)
                end

                pcd_L, R = pcond_Lap(L_av, q, Lap)
                sp_L_tp[j] = nnz(sparse(R))
                cnd_tp[j] = cond(pcd_L)
                percent_edges_tp[j] = nb_of_edges(L_av) / m
            end

            cnd[i] = mean(cnd_tp)
            sp_L[i] = mean(sp_L_tp)
            percent_edges[i] = mean(percent_edges_tp)

            cnd_std[i] = std(cnd_tp)
            sp_L_std[i] = std(sp_L_tp)
            percent_edges_std[i] = std(percent_edges_tp)
        end
        D = Dict(
            "cnd" => cnd,
            #
            "cnd_std" => cnd_std,
            #
            "sp_L" => sp_L,
            #
            "sp_L_std" => sp_L_std,
            #
            "percent_edges" => percent_edges,
            #
            "percent_edges_std" => percent_edges_std,
        )
        push!(D_all, method => D)
    end

    return D_all
end

function benchmark_syncrank(meta_g, planted_ranking_score, n_batch, n_rep, rng)
    n = nv(meta_g)
    m = ne(meta_g)

    # compute ground truth
    planted_ranking = ranking_from_score(planted_ranking_score)

    #  include edge weights in meta_g: w_{uv} = 1/sqrt{d(u)d(v)}
    normalize_meta_g!(meta_g)

    # parameters for benchmarking
    batch = Int(floor(n))

    # technical parameters
    weighted = true # weighted graph is used
    extra_normalization = false # normalized Laplacian is used
    singular = true # all eigenpairs are computed for more stability

    # incidence matrix
    B = magnetic_incidence(meta_g)

    #######################################
    # syncrank with full magnetic Laplacian
    W = I # weight matrix
    if weighted
        e_weights = get_edges_prop(meta_g, :e_weight, true, 1.0)
        W *= diagm(e_weights)
    end
    L = B * W * B'

    # leverage score
    q = 0
    lev = leverage_score(B, q; W)

    N = 0
    if extra_normalization
        normalize_Lap!(L)
        Deg = Diagonal(L)
        N = inv(sqrt(Deg))
    end

    # least eigenvector full Laplacian
    v = least_eigenvector(L; singular)

    # recovered ranking full Laplacian
    ranking_full = syncrank(L, meta_g; singular)
    tau_full = corkendall(planted_ranking, ranking_full)
    #######################################

    # start benchmarking

    rangebatch = 1:n_batch

    methods = ["DPP unif", "DPP LS", "iid unif", "iid LS"]
    D_all = Dict()

    for method in methods

        # initialization
        err = zeros(size(rangebatch))
        err_std = zeros(size(rangebatch))

        tau = zeros(size(rangebatch))
        tau_std = zeros(size(rangebatch))

        percent_edges = zeros(size(rangebatch))
        percent_edges_std = zeros(size(rangebatch))

        for i in 1:length(rangebatch)
            err_tp = zeros(n_rep, 1)
            tau_tp = zeros(n_rep, 1)
            percent_edges_tp = zeros(n_rep, 1)

            t = rangebatch[i]

            for j in 1:n_rep
                L_av = zeros(n, n)

                if method == "DPP unif"
                    # DPP uniform weighting
                    L_av = average_sparsifier(rng, meta_g, nothing, q, t; weighted)

                elseif method == "DPP LS"
                    # DPP leverage score weighting
                    L_av = average_sparsifier(rng, meta_g, lev, q, t; weighted)

                elseif method == "iid unif"
                    # iid uniform with uniform weighting
                    L_av = average_sparsifier_iid(rng, meta_g, nothing, batch, t; weighted)

                elseif method == "iid LS"
                    # iid leverage score with leverage score weighting
                    L_av = average_sparsifier_iid(rng, meta_g, lev, batch, t; weighted)
                end

                if extra_normalization
                    # normalization step with N = inv(sqrt(Deg)) and Deg is the full degree matrix
                    L_av = N * L_av * N
                    #normalize_Lap!(L_av)
                end

                v_av = least_eigenvector(L_av; singular)
                err_tp[j] = eigenvec_dist(v, v_av)

                ranking = syncrank(L_av, meta_g; singular)
                tau_tp[j] = corkendall(planted_ranking, ranking)

                percent_edges_tp[j] = nb_of_edges(L_av) / m
            end
            err[i] = mean(err_tp)
            err_std[i] = std(err_tp)

            tau[i] = mean(tau_tp)
            tau_std[i] = std(tau_tp)

            percent_edges[i] = mean(percent_edges_tp)
            percent_edges_std[i] = std(percent_edges_tp)
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
            "percent_edges" => percent_edges,
            #
            "percent_edges_std" => percent_edges_std,
            #
            "tau_full" => tau_full,
        )
        push!(D_all, method => D)
    end

    return D_all
end
