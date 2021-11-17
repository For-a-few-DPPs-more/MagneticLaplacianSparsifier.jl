getRNG(seed=nothing) = isnothing(seed) ? Random.default_rng() : Random.default_rng(seed)
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
    g::AbstractMetaGraph, e::Edge, property::Symbol, oriented::Bool=true, default::Real=1.0
)
    if haskey(g.eprops, e) && haskey(g.eprops[e], property)
        return g.eprops[e][property]
    end
    _e = reverse(e)
    if oriented && haskey(g.eprops, _e) && haskey(g.eprops[_e], property)
        return -g.eprops[_e][property]
    end
    return default
end

function get_edges_prop(
    g::AbstractMetaGraph, prop::Symbol, oriented::Bool=true, default::Real=1.0
)
    return [get_edge_prop(g, e, prop, oriented, default) for e in edges(g)]
end
