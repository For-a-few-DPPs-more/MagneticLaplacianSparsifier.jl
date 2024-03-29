# Following M. Cucuringu to build comparison graph
# SYNC-RANK: ROBUST RANKING, CONSTRAINED RANKING AND RANK AGGREGATION VIA EIGENVECTOR AND SDP SYNCHRONIZATION
# Multiplicative Uniform Noise model
function gen_graph_mun(rng, n, p, eta; planted_score=nothing, scaling=1)
    return erdos_renyi(rng, n, p, eta, :mun; planted_score, scaling)
end

# Following M. Cucuringu to build comparison graph
# SYNC-RANK: ROBUST RANKING, CONSTRAINED RANKING AND RANK AGGREGATION VIA EIGENVECTOR AND SDP SYNCHRONIZATION
# Erdos-Renyi Outliers model
function gen_graph_ero(rng, n, p, eta; planted_score=nothing, scaling=1)
    return erdos_renyi(rng, n, p, eta, :ero; planted_score, scaling)
end

function erdos_renyi(
    rng::Random.AbstractRNG,
    nv::Integer,
    p::Real,
    η::Real,
    model::Symbol;
    planted_score::Union{Array,Nothing}=nothing,
    scaling::Real=1.0,
)
    m = div(nv * (nv - 1), 2)
    ne = rand(rng, Binomial(m, p)) # nb of edges is a sum of m Bernoulli's of parameter p.
    return erdos_renyi(rng, nv, ne, η, model; planted_score, scaling)
end

function erdos_renyi(
    rng::Random.AbstractRNG,
    n_v::Integer,
    n_e::Integer,
    η::Real,
    model::Symbol;
    planted_score::Union{Array,Nothing}=nothing,
    scaling::Real=1.0,
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
            θ *= scaling
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
    planted_score::Union{Array,Nothing}=nothing,
    scaling::Real=1.0,
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
            θ *= scaling
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
    planted_score::Union{Array,Nothing}=nothing,
    scaling::Real=1.0,
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
                    θ *= scaling
                    add_edge!(g, u, v, :angle, θ)
                end
            end
        end
    end
    return g
end

function gen_graph_cliques(
    rng::Random.AbstractRNG,
    n::Integer,
    noise::Real,
    η::Real;
    planted_score::Union{Array,Nothing}=nothing,
)
    n_v = 2 * n
    g = Graph(n_v)
    meta_g = MetaGraph(g)

    # convention: ranking score of node i is r_i = score[i]
    if planted_score === nothing
        planted_score = collect(1:n_v)
    end

    # first clique
    for u in 1:n
        for v in 1:n
            if u < v
                h_u = planted_score[u]
                h_v = planted_score[v]
                if (rand(rng) < η)
                    θ = rand(rng, (-n_v + 1):(n_v - 1)) * π / (n_v - 1)
                    add_edge!(meta_g, u, v, :angle, θ)
                else
                    θ = (h_u - h_v) * π / (n_v - 1)
                    θ *= 1.0 + noise * 2 * (rand(rng) - 0.5)
                    add_edge!(meta_g, u, v, :angle, θ)
                end
            end
        end
    end
    # bottleneck
    u = n
    v = n + 1

    h_u = planted_score[u]
    h_v = planted_score[v]
    θ = (h_u - h_v) * π / (n_v - 1)
    θ *= 1.0 + noise * 2 * (rand(rng) - 0.5)
    add_edge!(meta_g, u, v, :angle, θ)

    # second clique
    for u in (n + 1):(2 * n)
        for v in (n + 1):(2 * n)
            if u < v
                h_u = planted_score[u]
                h_v = planted_score[v]
                if (rand(rng) < η)
                    θ = rand(rng, (-n_v + 1):(n_v - 1)) * π / (n_v - 1)
                    add_edge!(meta_g, u, v, :angle, θ)
                else
                    θ = (h_u - h_v) * π / (n_v - 1)
                    θ *= 1.0 + noise * 2 * (rand(rng) - 0.5)
                    add_edge!(meta_g, u, v, :angle, θ)
                end
            end
        end
    end
    return meta_g
end

function gen_graph_planted_triangles(
    rng::Random.AbstractRNG,
    n_v::Integer,
    noise::Real,
    p_edge::Real,
    p_triangles::Real;
    planted_score::Union{Array,Nothing}=nothing,
)
    g = Graph(n_v)
    meta_g = MetaGraph(g)

    # convention: ranking score of node i is r_i = score[i]
    if planted_score === nothing
        planted_score = collect(1:n_v)
    end

    for u in 1:n_v
        for v in 1:n_v
            if u < v
                h_u = planted_score[u]
                h_v = planted_score[v]
                if (rand(rng) < p_edge)
                    θ = (h_u - h_v) * π / (n_v - 1)
                    θ *= 1.0 + noise * 2 * (rand(rng) - 0.5)
                    add_edge!(meta_g, u, v, :angle, θ)
                end
            end
        end
    end
    nb_t = 0
    for u in 1:n_v
        for v in 1:n_v
            for w in 1:n_v
                if u < v < w
                    if (rand(rng) < p_triangles)
                        add_edge!(meta_g, u, v, :angle, π / 6)
                        add_edge!(meta_g, v, w, :angle, π / 6)
                        add_edge!(meta_g, w, u, :angle, π / 6)
                        nb_t += 1
                    end
                end
            end
        end
    end
    println(nb_t)
    return meta_g
end

function ero_located(
    rng::Random.AbstractRNG,
    n_v::Integer,
    p::Real,
    η::Real;
    planted_score::Union{Array,Nothing}=nothing,
)
    g = MetaGraph(n_v)

    m = div(n_v * (n_v - 1), 2)
    n_e = rand(rng, Binomial(m, p))

    noisy_edges = []
    err_edges = []

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
                push!(noisy_edges, [u v])

                erroneous_θ = rand(rng, (-n_v + 1):(n_v - 1)) * π / (n_v - 1)
                diff = erroneous_θ - θ

                push!(err_edges, diff)
                θ = rand(rng, (-n_v + 1):(n_v - 1)) * π / (n_v - 1)
            end
            add_edge!(g, u, v, :angle, θ)
        end
    end
    return g, noisy_edges, err_edges
end


function turn_into_connection_graph(rng,meta_g,eta,type,planted_score)

    n = nv(meta_g)

    model = :mun
    if type=="MUN"
        model = :mun
    elseif type=="ERO"
        model = :ero
    end

    meta_g = MetaGraph(meta_g)
    for e in edges(meta_g)
        u = src(e)
        v = dst(e)
        h_u = planted_score[u]
        h_v = planted_score[v]
        θ = (h_u - h_v) * π / (n - 1)
        if (model === :ero) && (rand(rng) < eta) # Erdos-Renyi Outliers
            θ = rand(rng, (-n + 1):(n - 1)) * π / (n - 1)
        elseif model === :mun # Multiplicative Uniform Noise
            θ *= 1.0 + eta * 2 * (rand(rng) - 0.5)
        end
        set_prop!(meta_g, e, :angle, θ)
    end

    return meta_g
end


function main_component(g)
    c = connected_components(g)
    _, i = findmax(length.(c))
    return g[c[i]]
end
