# Following M. Cucuringu to build comparison graph
# SYNC-RANK: ROBUST RANKING, CONSTRAINED RANKING AND RANK AGGREGATION VIA EIGENVECTOR AND SDP SYNCHRONIZATION
# Multiplicative Uniform Noise model
function gen_graph_mun(rng, n, p, eta; planted_score=nothing)
    return erdos_renyi(rng, n, p, eta, :mun; planted_score)
end

# Following M. Cucuringu to build comparison graph
# SYNC-RANK: ROBUST RANKING, CONSTRAINED RANKING AND RANK AGGREGATION VIA EIGENVECTOR AND SDP SYNCHRONIZATION
# Erdos-Renyi Outliers model
function gen_graph_ero(rng, n, p, eta; planted_score=nothing)
    return erdos_renyi(rng, n, p, eta, :ero; planted_score)
end

function erdos_renyi(
    rng::Random.AbstractRNG,
    nv::Integer,
    p::Real,
    η::Real,
    model::Symbol;
    planted_score=nothing,
)
    m = div(nv * (nv - 1), 2)
    ne = rand(rng, Binomial(m, p)) # nb of edges is a sum of m Bernoulli's of parameter p.
    return erdos_renyi(rng, nv, ne, η, model; planted_score)
end

function erdos_renyi(
    rng::Random.AbstractRNG,
    n_v::Integer,
    n_e::Integer,
    η::Real,
    model::Symbol;
    planted_score=nothing,
)::AbstractMetaGraph
    g = MetaGraph(n_v)
    # convention: ranking score of node i is r_i = score[i]
    if planted_score === nothing
        planted_score = collect(1:n_v)
    end
    while ne(g) < n_e
        u = rand(rng, 1:n_v) # discrete uniform distribution
        v = rand(rng, 1:n_v)
        if u < v
            h_u = planted_score[u]
            h_v = planted_score[v]
            θ = (h_u - h_v) * π / (n_v - 1)
            if (model === :ero) && (rand(rng) < η) # Erdos-Renyi Outliers
                θ = rand(rng, (-n_v + 1):(n_v - 1)) * π / (n_v - 1)
            elseif model === :mun # Multiplicative Uniform Noise
                θ *= 1.0 + η * 2 * (rand(rng) - 0.5)
            end
            add_edge!(g, u, v, :angle, θ)
        end
    end
    return g
end

function ero_mun(
    rng::Random.AbstractRNG,
    n_v::Integer,
    p::Real,
    η::Real,
    noise::Real;
    planted_score=nothing,
)::AbstractMetaGraph
    g = MetaGraph(n_v)

    m = div(n_v * (n_v - 1), 2)
    n_e = rand(rng, Binomial(m, p))
    # convention: ranking score of node i is r_i = score[i]
    if planted_score === nothing
        planted_score = collect(1:n_v)
    end
    while ne(g) < n_e
        u = rand(rng, 1:n_v)
        v = rand(rng, 1:n_v)
        if u < v
            h_u = planted_score[u]
            h_v = planted_score[v]
            θ = (h_u - h_v) * π / (n_v - 1)
            if (rand(rng) < η)
                θ = rand(rng, (-n_v + 1):(n_v - 1)) * π / (n_v - 1)
            else
                θ *= 1.0 + noise * 2 * (rand(rng) - 0.5)
            end
            add_edge!(g, u, v, :angle, θ)
        end
    end
    return g
end

function ero_mun_sbm(
    rng::Random.AbstractRNG,
    n_v::Integer,
    p_in::Real,
    p_out::Real,
    η::Real,
    noise::Real;
    planted_score=nothing,
)::AbstractMetaGraph
    g = MetaGraph(n_v)

    # convention: ranking score of node i is r_i = score[i]
    if planted_score === nothing
        planted_score = collect(1:n_v)
    end

    n_half = div(n_v, 2)

    for u in 1:n_v
        for v in 1:n_v
            if u < v
                in_com_1 = (u <= n_half) && (v <= n_half)
                in_com_2 = (u > n_half) && (v > n_half)
                in_edge_exits = (in_com_1 || in_com_2) && rand(rng) < p_in
                out_edge_exits = (!in_com_1 && !in_com_2) && rand(rng) < p_out
                if in_edge_exits || out_edge_exits
                    h_u = planted_score[u]
                    h_v = planted_score[v]
                    θ = (h_u - h_v) * π / (n_v - 1)
                    if (rand(rng) < η)
                        θ = rand(rng, (-n_v + 1):(n_v - 1)) * π / (n_v - 1)
                    else
                        θ *= 1.0 + noise * 2 * (rand(rng) - 0.5)
                    end
                    add_edge!(g, u, v, :angle, θ)
                end
            end
        end
    end
    return g
end

function gen_graph_mun_basic(n, p, eta)
    # Following M. Cucuringu to build comparison graph
    # SYNC-RANK: ROBUST RANKING, CONSTRAINED RANKING AND RANK AGGREGATION VIA EIGENVECTOR AND SDP SYNCHRONIZATION
    # Multiplicative Uniform Noise model
    g = Graph(n)
    meta_g = MetaGraph(g, :angle, 0.0)

    # ranking
    r = collect(1:n)

    for i in r
        for j in r
            if i < j
                if rand(1)[1] < p
                    e = [i j]
                    edges = consecutive_pairs(e)
                    add_edges_from!(meta_g, edges)

                    err = eta * rand(1)[1]
                    a = pi * (i - j) * (1 + err) / (n - 1)

                    for e in edges
                        set_prop!(meta_g, Edge(e), :angle, a)
                    end
                end
            end
        end
    end

    return meta_g
end

function gen_graph_ero_basic(n, p, eta)
    # Following M. Cucuringu to build comparison graph
    # SYNC-RANK: ROBUST RANKING, CONSTRAINED RANKING AND RANK AGGREGATION VIA EIGENVECTOR AND SDP SYNCHRONIZATION
    # Erdos-Renyi Outliers model
    g = Graph(n)
    meta_g = MetaGraph(g, :angle, 0.0)

    # ranking
    r = collect(1:n)

    for i in r
        for j in r
            if i < j
                if rand(1)[1] < p
                    e = [i j]
                    edges = consecutive_pairs(e)
                    add_edges_from!(meta_g, edges)

                    if rand(1)[1] < 1 - eta
                        err = 0
                    else
                        err = (n - 1) * rand(1)[1]
                    end
                    a = pi * (i - j + err) / (n - 1)
                    for e in edges
                        set_prop!(meta_g, Edge(e), :angle, a)
                    end
                end
            end
        end
    end

    return meta_g
end
