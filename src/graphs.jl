# Following M. Cucuringu to build comparison graph
# SYNC-RANK: ROBUST RANKING, CONSTRAINED RANKING AND RANK AGGREGATION VIA EIGENVECTOR AND SDP SYNCHRONIZATION
# Multiplicative Uniform Noise model
function gen_graph_mun(rng, n, p, eta; planted_ranking=nothing)
    return erdos_renyi(rng, n, p, eta, :mun; planted_ranking)
end

# Following M. Cucuringu to build comparison graph
# SYNC-RANK: ROBUST RANKING, CONSTRAINED RANKING AND RANK AGGREGATION VIA EIGENVECTOR AND SDP SYNCHRONIZATION
# Erdos-Renyi Outliers model
function gen_graph_ero(rng, n, p, eta; planted_ranking=nothing)
    return erdos_renyi(rng, n, p, eta, :ero; planted_ranking)
end

function erdos_renyi(
    rng::Random.AbstractRNG,
    nv::Integer,
    p::Real,
    η::Real,
    model::Symbol;
    planted_ranking=nothing,
)
    m = div(nv * (nv - 1), 2)
    ne = rand(rng, Binomial(m, p)) # nb of edges is a sum of m Bernoulli's of parameter p.
    return erdos_renyi(rng, nv, ne, η, model; planted_ranking)
end

function erdos_renyi(
    rng::Random.AbstractRNG,
    n_v::Integer,
    n_e::Integer,
    η::Real,
    model::Symbol;
    planted_ranking=nothing,
)::AbstractMetaGraph
    g = MetaGraph(n_v)
    # convention: ranking score of node i is r_i = score[i]
    if planted_ranking === nothing
        planted_ranking = collect(1:n_v)
    end
    while ne(g) < n_e
        u = rand(rng, 1:n_v) # discrete uniform distribution
        v = rand(rng, 1:n_v)
        if u < v
            h_u = planted_ranking[u]
            h_v = planted_ranking[v]
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
