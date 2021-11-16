# todo write the docstrings

function multi_type_spanning_forest(g::MetaGraph{T,U}, q::Real)::MetaGraph{T,U} where {T,U}
    rng = getRNG()
    return multi_type_spanning_forest(rng, g, q)
end

function multi_type_spanning_forest(
    rng::Random.AbstractRNG, g::MetaGraph{T,U}, q::Real
)::MetaGraph{T,U} where {T,U}
    # Initialize the multi type spanning forest
    mtsf = MetaGraph{T,U}(nv(g))
    nv_mtsf = 0
    weight = 1.0
    roots = T[]

    # Initialize the random walk
    walk = T[]
    unvisited = Set{T}(vertices(g))

    # Start the random walk
    rng = getRNG(rng)
    n0 = rand(rng, unvisited)
    push!(walk, n0)
    setdiff!(unvisited, n0)

    while nv_mtsf < nv(g)
        n0_is_root = rand(rng) < q / (q + degree(g, n0))
        if n0_is_root
            push!(roots, n0)
            add_edges_from!(mtsf, consecutive_pairs(walk))
            nv_mtsf += length(walk) - 1
            setdiff!(unvisited, walk)
            n0 = restart_walk_from_unvisited_node!(rng, walk, unvisited)
            continue
        end

        n1 = rand(rng, neighbors(g, n0))
        push!(walk, n1)

        if n1 in unvisited
            setdiff!(unvisited, n1)
            n0 = n1  # continue the walk

        elseif degree(mtsf, n1) > 0  # n1 in mtsf
            add_edges_from!(mtsf, consecutive_pairs(walk))
            nv_mtsf += length(walk) - 1
            setdiff!(unvisited, walk)
            n0 = restart_walk_from_unvisited_node!(rng, walk, unvisited)

        else  # if n1 in walk: identify unique cycle/loop in walk with knot n1
            idx_n1 = findfirst(x -> x == n1, walk)
            cycle_nodes = @view walk[idx_n1:end]
            keep, alpha = keep_cycle(rng, g, consecutive_pairs(cycle_nodes))

            if keep  # cycle
                weight *= max(alpha, 1)
                add_edges_from!(mtsf, consecutive_pairs(walk))
                nv_mtsf += length(walk) - 1 # since walk contains twice the knot
                setdiff!(unvisited, walk)
                n0 = restart_walk_from_unvisited_node!(rng, walk, unvisited)

            else  # pop cycle but keep loopy node
                union!(unvisited, cycle_nodes)
                setdiff!(unvisited, n1)  # remove n1 which was part of cycle_nodes
                resize!(walk, idx_n1)
                n0 = n1  # continue the walk
            end
        end
    end
    set_prop!(mtsf, :weight, weight)
    set_prop!(mtsf, :roots, roots)
    return mtsf
end

function restart_walk_from_unvisited_node!(rng, walk, unvisited)
    empty!(walk)
    n0 = -1
    if !isempty(unvisited)
        rng = getRNG(rng)
        n0 = rand(rng, unvisited)
        setdiff!(unvisited, n0)
        push!(walk, n0)
    end
    return n0
end

function keep_cycle(rng, graph::AbstractMetaGraph, cycle_edges)
    alpha = 1.0 - cos(curvature(graph, cycle_edges))
    keep = rand(rng) < min(alpha, 1)
    return keep, alpha
end

function curvature(
    g::AbstractMetaGraph, edges, oriented::Bool=true, default_angle::Real=0.0
)
    curvature_ = 0.0
    for e in edges
        angle = get_edge_property(g, Edge(e), :angle, oriented, default_angle)
        curvature_ += angle
    end
    return curvature_
end

function get_edge_property(
    g::AbstractMetaGraph,
    e::Edge,
    property::Symbol=:angle,
    oriented::Bool=true,
    default::Real=1.0,
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
