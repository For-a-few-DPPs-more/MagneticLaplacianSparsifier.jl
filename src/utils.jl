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

function cond_numbers(meta_g, q, n_tot, n_rep, rng; q_system=q,methods=nothing)
    m = ne(meta_g)
    n = nv(meta_g)
    batch = n

    weighted = false

    B = magnetic_incidence(meta_g)
    Lap = B * B'
    lev = leverage_score(B, q)

    cdL = cond(Lap + q_system * I)

    B_ust = magnetic_incidence_matrix(meta_g; oriented=true, phases=false)
    lev_ust = leverage_score(Matrix(B_ust), 0)

    if methods===nothing
        methods = ["DPP unif", "DPP LS", "iid unif", "iid LS", "UST unif", "UST LS"]
    end

    D_all = Dict()

    for method in methods
        println(method)
        # initialization
        cnd = zeros(n_tot, 1)
        sp_L = zeros(n_tot, 1)
        timing = zeros(n_tot, 1)
        percent_edges = zeros(n_tot, 1)
        cycles = zeros(n_tot, 1)
        roots = zeros(n_tot, 1)

        cnd_std = zeros(n_tot, 1)
        sp_L_std = zeros(n_tot, 1)
        timing_std = zeros(n_tot, 1)
        percent_edges_std = zeros(n_tot, 1)
        roots_std = zeros(n_tot, 1)
        cycles_std = zeros(n_tot, 1)

        for i in 1:n_tot

            # temporary arrays
            cnd_tp = zeros(n_rep, 1)
            sp_L_tp = zeros(n_rep, 1)
            time_tp = zeros(n_rep, 1)
            percent_edges_tp = zeros(n_rep, 1)
            roots_tp = zeros(n_rep, 1)
            cycles_tp = zeros(n_rep, 1)

            for j in 1:n_rep
                L_av = zeros(n, n)
                time = 0
                n_cls = 0
                n_rts = 0
                if method == "DPP unif"
                    # DPP uniform weighting
                    # L_av = average_sparsifier(rng, meta_g, nothing, q, i; weighted)
                    vec = @timed average_sparsifier(rng, meta_g, nothing, q, i; weighted)
                    out = vec[1]
                    L_av = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    time = vec[2]
                elseif method == "DPP LS"
                    # DPP leverage score weighting
                    # L_av = average_sparsifier(rng, meta_g, lev, q, i; weighted)
                    vec = @timed average_sparsifier(rng, meta_g, lev, q, i; weighted)
                    out = vec[1]
                    L_av = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    time = vec[2]
                elseif method == "iid unif"
                    # iid uniform with uniform weighting
                    # L_av = average_sparsifier_iid(rng, meta_g, nothing, batch, i; weighted)
                    vec = @timed average_sparsifier_iid(
                        rng, meta_g, nothing, batch, i; weighted
                    )
                    L_av = vec[1]
                    time = vec[2]
                elseif method == "iid LS"
                    # iid leverage score with leverage score weighting
                    # L_av = average_sparsifier_iid(
                    #     rng, meta_g, lev, batch, i; weighted
                    # )
                    vec = @timed average_sparsifier_iid(
                        rng, meta_g, lev, batch, i; weighted
                    )
                    L_av = vec[1]
                    time = vec[2]
                elseif method == "UST unif"
                    # UST uniform weighting
                    absorbing_node = true
                    ust = true
                    q_ust = 0
                    # L_av = average_sparsifier(
                    #     rng, meta_g, nothing, q_ust, i; weighted, absorbing_node, ust
                    # )
                    vec = @timed average_sparsifier(
                        rng, meta_g, nothing, q_ust, i; weighted, absorbing_node, ust
                    )
                    out = vec[1]
                    L_av = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    time = vec[2]
                elseif method == "UST LS"
                    # UST LS weighting
                    absorbing_node = true
                    ust = true
                    q_ust = 0
                    # L_av = average_sparsifier(
                    #     rng, meta_g, lev_ust, q_ust, i; weighted, absorbing_node, ust
                    # )
                    vec = @timed average_sparsifier(
                        rng, meta_g, lev_ust, q_ust, i; weighted, absorbing_node, ust
                    )
                    out = vec[1]
                    L_av = out[1]
                    n_cls = out[2]
                    n_rts = out[3]
                    time = vec[2]
                end
                # by default q_system = q
                pcd_L, R = pcond_Lap(L_av, q_system, Lap)
                sp_L_tp[j] = nnz(sparse(R))
                cnd_tp[j] = cond(pcd_L)
                percent_edges_tp[j] = nb_of_edges(L_av) / m
                time_tp[j] = time
                roots_tp[j] = n_rts
                cycles_tp[j] = n_cls
            end

            cnd[i] = mean(cnd_tp)
            sp_L[i] = mean(sp_L_tp)
            percent_edges[i] = mean(percent_edges_tp)
            timing[i] = mean(time_tp)
            roots[i] = mean(roots_tp)
            cycles[i] = mean(cycles_tp)

            cnd_std[i] = std(cnd_tp)
            sp_L_std[i] = std(sp_L_tp)
            percent_edges_std[i] = std(percent_edges_tp)
            timing_std[i] = std(time_tp)
            roots_std[i] = std(roots_tp)
            cycles_std[i] = std(cycles_tp)
        end
        D = Dict(
            "cdL" => cdL,
            #
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
            #
            "timing" => timing,
            #
            "timing_std" => timing_std,
            #
            "roots" => roots,
            #
            "roots_std" => roots_std,
            #
            "cycles" => cycles_std,
            #
            "cycles_std" => cycles_std,
        )
        push!(D_all, method => D)
    end

    return D_all
end


function benchmark_syncrank(meta_g, planted_ranking_score, n_batch, n_rep, rng;methods=nothing)
    n = nv(meta_g)
    m = ne(meta_g)
    println("start")
    #
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
    k = 10 # number of upsets in top k

    # incidence matrix
    B = magnetic_incidence(meta_g)
    B_ust = magnetic_incidence_matrix(meta_g; oriented=true, phases=false)

    #######################################
    # syncrank with full magnetic Laplacian
    W = I # weight matrix
    if weighted
        e_weights = get_edges_prop(meta_g, :e_weight, true, 1.0)
        W *= diagm(e_weights)
    end
    L = B * W * B'

    # leverage scores
    println("computing leverage scores")
    q = 0
    lev = leverage_score(B, q; W)
    lev_ust = leverage_score(Matrix(B_ust), q; W)

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
    spear_full = corspearman(planted_ranking, ranking_full)
    #######################################

    # start benchmarking

    rangebatch = 1:n_batch

    if methods===nothing
        methods = ["DPP unif", "DPP LS", "iid unif", "iid LS", "UST unif", "UST LS"]
    end
    D_all = Dict()

    for method in methods
        println(method)

        # initialization
        err = zeros(size(rangebatch))
        err_std = zeros(size(rangebatch))

        tau = zeros(size(rangebatch))
        tau_std = zeros(size(rangebatch))

        spear = zeros(size(rangebatch))
        spear_std = zeros(size(rangebatch))

        percent_edges = zeros(size(rangebatch))
        percent_edges_std = zeros(size(rangebatch))

        upsets_in_top = zeros(size(rangebatch))
        upsets_in_top_std = zeros(size(rangebatch))

        roots = zeros(size(rangebatch))
        roots_std = zeros(size(rangebatch))

        cycles = zeros(size(rangebatch))
        cycles_std = zeros(size(rangebatch))

        for i in 1:length(rangebatch)
            err_tp = zeros(n_rep, 1)
            tau_tp = zeros(n_rep, 1)
            roots_tp = zeros(n_rep, 1)
            cycles_tp = zeros(n_rep, 1)
            spear_tp = zeros(n_rep, 1)
            upsets_in_top_tp = zeros(n_rep, 1)

            percent_edges_tp = zeros(n_rep, 1)

            t = rangebatch[i]

            for j in 1:n_rep
                L_av = zeros(n, n)
                n_cles = 0
                n_rts = 0
                if method == "DPP unif"
                    # DPP uniform weighting
                    L_av, n_cles, n_rts = average_sparsifier(
                        rng, meta_g, nothing, q, t; weighted
                    )

                elseif method == "DPP LS"
                    # DPP leverage score weighting
                    L_av, n_cles, n_rts = average_sparsifier(
                        rng, meta_g, lev, q, t; weighted
                    )

                elseif method == "iid unif"
                    # iid uniform with uniform weighting
                    L_av = average_sparsifier_iid(rng, meta_g, nothing, batch, t; weighted)

                elseif method == "iid LS"
                    # iid leverage score with leverage score weighting
                    L_av = average_sparsifier_iid(rng, meta_g, lev, batch, t; weighted)

                elseif method == "UST unif"
                    # UST uniform weighting
                    absorbing_node = true
                    ust = true
                    q_ust = 0
                    L_av, n_cles, n_rts = average_sparsifier(
                        rng, meta_g, nothing, q_ust, t; weighted, absorbing_node, ust
                    )
                elseif method == "UST LS"
                    # UST LS weighting
                    absorbing_node = true
                    ust = true
                    q_ust = 0
                    L_av, n_cles, n_rts = average_sparsifier(
                        rng, meta_g, lev_ust, q_ust, t; weighted, absorbing_node, ust
                    )
                end

                if extra_normalization
                    # normalization step with N = inv(sqrt(Deg)) and Deg is the full degree matrix
                    L_av = N * L_av * N
                    #normalize_Lap!(L_av)
                end

                v_av = least_eigenvector(L_av; singular)
                err_tp[j] = eigenvec_dist(v, v_av)
                roots_tp[j] = n_rts
                cycles_tp[j] = n_cles

                ranking = syncrank(L_av, meta_g; singular)
                tau_tp[j] = corkendall(planted_ranking, ranking)
                spear_tp[j] = corspearman(planted_ranking, ranking)
                upsets_in_top_tp[j] = nb_upsets_in_top(meta_g, ranking, k)

                percent_edges_tp[j] = nb_of_edges(L_av) / m
            end
            err[i] = mean(err_tp)
            err_std[i] = std(err_tp)

            tau[i] = mean(tau_tp)
            tau_std[i] = std(tau_tp)

            spear[i] = mean(spear_tp)
            spear_std[i] = std(spear_tp)

            upsets_in_top[i] = mean(upsets_in_top_tp)
            upsets_in_top_std[i] = std(upsets_in_top_tp)

            percent_edges[i] = mean(percent_edges_tp)
            percent_edges_std[i] = std(percent_edges_tp)

            roots[i] = mean(roots_tp)
            roots_std[i] = std(roots_tp)

            cycles[i] = mean(cycles_tp)
            cycles_std[i] = std(cycles_tp)
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
            "percent_edges" => percent_edges,
            #
            "percent_edges_std" => percent_edges_std,
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
        )
        push!(D_all, method => D)
    end

    return D_all
end

function plot_comparison_sync(
    metric::String, D_all, y_limits; legendposition::Symbol=:bottomright
)
    metric_std = metric * "_std"

    method = "DPP unif"
    D = D_all[method]
    x = D["percent_edges"]
    y = D[metric]
    y_er = D[metric_std]

    Plots.plot(
        x,
        y;
        yerror=y_er,
        labels=method,
        markerstrokecolor=:auto,
        markershape=:xcross,
        markersize=5,
        linewidth=2,
        markerstrokewidth=2,
        xtickfont=font(13),
        ytickfont=font(13),
        guidefont=font(13),
        legendfont=font(13),
        framestyle=:box,
        margins=0.1 * 2cm,
    )

    method = "DPP LS"
    D = D_all[method]
    x = D["percent_edges"]
    y = D[metric]
    y_er = D[metric_std]

    Plots.plot!(
        x,
        y;
        yerror=y_er,
        labels=method,
        markerstrokecolor=:auto,
        markershape=:circle,
        markersize=5,
        linewidth=2,
        markerstrokewidth=2,
    )

    method = "iid unif"
    D = D_all[method]
    x = D["percent_edges"]
    x_er = D["percent_edges_std"]

    y = D[metric]
    y_er = D[metric_std]

    Plots.plot!(
        x,
        y;
        xerror=x_er,
        yerror=y_er,
        labels=method,
        markerstrokecolor=:auto,
        markersize=5,
        linestyle=:dash,
        markershape=:rtriangle,
        linewidth=2,
        markerstrokewidth=2,
    )

    method = "iid LS"
    D = D_all[method]
    x = D["percent_edges"]
    x_er = D["percent_edges_std"]

    y = D[metric]
    y_er = D[metric_std]

    Plots.plot!(
        x,
        y;
        xerror=x_er,
        yerror=y_er,
        labels=method,
        markerstrokecolor=:auto,
        markersize=5,
        linestyle=:dash,
        markershape=:utriangle,
        linewidth=2,
        markerstrokewidth=2,
    )

    method = "UST unif"
    D = D_all[method]

    x = D["percent_edges"]
    y = D[metric]
    y_er = D[metric_std]

    Plots.plot!(
        x,
        y;
        xerror=x_er,
        yerror=y_er,
        labels=method,
        markerstrokecolor=:auto,
        markersize=5,
        markershape=:dtriangle,
        linewidth=2,
        markerstrokewidth=2,
        framestyle=:box,
        margins=0.1 * 2Plots.cm,
    )

    method = "UST LS"
    D = D_all[method]

    x = D["percent_edges"]
    y = D[metric]
    y_er = D[metric_std]

    Plots.plot!(
        x,
        y;
        xerror=x_er,
        yerror=y_er,
        labels=method,
        markerstrokecolor=:auto,
        markersize=5,
        markershape=:octagon,
        linewidth=2,
        markerstrokewidth=2,
        framestyle=:box,
        margins=0.1 * 2Plots.cm,
        legend=legendposition,
    )

    xlabel!("percentage of edges")
    ylims!(y_limits)

    if metric === "err"
        ylabel!("Distance between eigenvectors")
        yaxis!(:log)

    elseif metric === "tau"
        x = D["percent_edges"]
        y = D["tau_full"] * ones(size(x))
        Plots.plot!(x, y; labels="full")
        ylabel!("Kendall's tau distance ")

    elseif metric === "upsets_in_top"
        ylabel!("number of upsets in top 10 ")

    elseif metric === "spear"
        x = D["percent_edges"]
        y = D["spear_full"] * ones(size(x))
        Plots.plot!(x, y; labels="full")
        ylabel!("Spearman")
    end
end

function plot_comparison_cond(D_all, y_limits; legendposition::Symbol=:bottomright)
    method = "DPP unif"
    D = D_all[method]

    cdL = D["cdL"]

    x = D["percent_edges"]
    y = D["cnd"]
    y_er = D["cnd_std"]

    plot(
        x,
        y;
        yerror=y_er,
        xlabel="fraction of edges",
        yaxis=:log,
        labels=method,
        markerstrokecolor=:auto,
        markershape=:xcross,
        markersize=5,
        xtickfont=font(13),
        ytickfont=font(13),
        guidefont=font(13),
        legendfont=font(13),
        linewidth=2,
        markerstrokewidth=2,
    )

    method = "DPP LS"
    D = D_all[method]

    x = D["percent_edges"]
    y = D["cnd"]
    y_er = D["cnd_std"]

    plot!(
        x,
        y;
        yerror=y_er,
        yaxis=:log,
        labels=method,
        markerstrokecolor=:auto,
        markershape=:circle,
        markersize=5,
        linewidth=2,
        markerstrokewidth=2,
    )

    method = "iid unif"
    D = D_all[method]

    x = D["percent_edges"]
    y = D["cnd"]
    y_er = D["cnd_std"]

    plot!(
        x,
        y;
        yerror=y_er,
        yaxis=:log,
        labels=method,
        markerstrokecolor=:auto,
        markershape=:rtriangle,
        markersize=5,
        linestyle=:dash,
        linewidth=2,
        markerstrokewidth=2,
    )

    method = "iid LS"
    D = D_all[method]

    x = D["percent_edges"]
    y = D["cnd"]
    y_er = D["cnd_std"]

    plot!(
        x,
        y;
        yerror=y_er,
        yaxis=:log,
        labels=method,
        markerstrokecolor=:auto,
        markershape=:utriangle,
        markersize=5,
        linestyle=:dash,
        linewidth=2,
        markerstrokewidth=2,
    )

    method = "UST unif"
    D = D_all[method]

    x = D["percent_edges"]
    y = D["cnd"]
    y_er = D["cnd_std"]

    plot!(
        x,
        y;
        yerror=y_er,
        yaxis=:log,
        labels=method,
        markerstrokecolor=:auto,
        markershape=:dtriangle,
        markersize=5,
        linewidth=2,
        markerstrokewidth=2,
    )

    method = "UST LS"
    D = D_all[method]

    x = D["percent_edges"]
    y = D["cnd"]
    y_er = D["cnd_std"]

    plot!(
        x,
        y;
        yerror=y_er,
        yaxis=:log,
        labels=method,
        markerstrokecolor=:auto,
        markershape=:octagon,
        markersize=5,
        linewidth=2,
        markerstrokewidth=2,
    )

    x = D["percent_edges"]
    y = cdL * ones(size(x))
    plot!(
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
        framestyle=:box,
        margins=0.1 * 2Plots.cm,
        markerstrokewidth=2,
        legend=legendposition,
    )

    ylims!(y_limits)
    # foldername = "figures/"
    # type = "precond"
    # name = type*"n"*string(n)*"p"*string(p)*"eta"*string(eta)*"q"*string(q)*".pdf"
    # savefig(foldername*name)
end

function plot_nb_cycles(D_all, method; legendposition=:topleft)
    D = D_all[method]

    x = D["percent_edges"]
    y = D["cycles"]
    y_err = D["cycles_std"]

    n_batch = length(x)

    plot(
        x,
        y;
        ribbon=y_err,
        labels="average number of unicycles",
        xlabel="percentage of edges",
        markersize=5,
        markershape=:circle,
        markerstrokecolor=:auto,
        linewidth=2,
        markerstrokewidth=2,
        xtickfont=font(13),
        ytickfont=font(13),
        guidefont=font(13),
        legendfont=font(13),
        framestyle=:box,
        margins=0.1 * 2cm,
    )
    # baseline
    plot!(
        x,
        1:n_batch;
        linewidth=2,
        labels="minimum number of unicycles",
        legend=legendposition,
    )
    # end
end

function plot_nb_roots(D_all, method; legendposition=:topleft)
    D = D_all[method]

    x = D["percent_edges"]
    y = D["roots"]
    y_err = D["roots_std"]

    n_batch = length(x)

    plot(
        x,
        y;
        ribbon=y_err,
        xlabel="percentage of edges",
        labels="average number of roots",
        markersize=5,
        markershape=:circle,
        markerstrokecolor=:auto,
        linewidth=2,
        markerstrokewidth=2,
        xtickfont=font(13),
        ytickfont=font(13),
        guidefont=font(13),
        legendfont=font(13),
        framestyle=:box,
        margins=0.1 * 2cm,
    )
    # baseline
    plot!(
        x, 1:n_batch;
        linewidth=2,
        labels="minimum number of roots",
        legend=legendposition
    )
    # end
end
