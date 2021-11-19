
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

# Following M. Cucuringu to build comparison graph
# SYNC-RANK: ROBUST RANKING, CONSTRAINED RANKING AND RANK AGGREGATION VIA EIGENVECTOR AND SDP SYNCHRONIZATION
# Multiplicative Uniform Noise model
gen_graph_mun(rng, n, p, eta) = erdos_renyi(rng, n, p, eta, :mun)

# Following M. Cucuringu to build comparison graph
# SYNC-RANK: ROBUST RANKING, CONSTRAINED RANKING AND RANK AGGREGATION VIA EIGENVECTOR AND SDP SYNCHRONIZATION
# Erdos-Renyi Outliers model
gen_graph_ero(rng, n, p, eta) = erdos_renyi(rng, n, p, eta, :ero)

function erdos_renyi(rng::Random.AbstractRNG, nv::Integer, p::Real, η::Real, model::Symbol)
    m = div(nv * (nv - 1), 2)
    ne = rand(rng, Binomial(m, p))
    return erdos_renyi(rng, nv, ne, η, model)
end

function erdos_renyi(
    rng::Random.AbstractRNG, n_v::Integer, n_e::Integer, η::Real, model::Symbol
)
    g = MetaGraph(n_v)
    while ne(g) < n_e
        u = rand(rng, 1:n_v)
        v = rand(rng, 1:n_v)
        if u != v
            θ = π * (u - v) / (n_v - 1)
            if (model === :ero) && (rand(rng) < η) # Erdos-Renyi Outliers model
                θ += π * rand(rng)
            elseif model === :mun # Multiplicative Uniform Noise model
                θ *= (1.0 + η * rand(rng))
            end
            add_edge!(g, u, v, :angle, θ)
        end
    end
    return g
end
