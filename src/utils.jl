getRNG(seed=nothing) = isnothing(seed) ? Random.default_rng() : Random.default_rng(seed)
getRNG(seed::Random.AbstractRNG) = seed

consecutive_pairs(path) = partition(path, 2, 1)

function add_edges_from!(g, edges)
    for e in edges
        edge = isa(e, Edge) ? e : Edge(e)
        add_edge!(g, edge)
    end
end

function add_angles_from!(mtsf, edges, mg)
    for e in edges
        edge = isa(e, Edge) ? e : Edge(e)
        θ = get_edge_property(mg, edge, :angle)
        set_prop!(mtsf, edge, :angle, θ)
    end
end
