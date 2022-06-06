"""
    multi_type_spanning_forest(g, q)

samples an mtsf with default random number generator.

# Arguments
- `g::MetaGraph{T,U}` connection graph.
- `q::Real` positive parameter playing the role of regularizer.
# Output
- `mtsf::MetaGraph{T,U}` connection graph associated with mtsf.

"""
function multi_type_spanning_forest(g::MetaGraph{T,U}, q::Real)::MetaGraph{T,U} where {T,U}
    return multi_type_spanning_forest(getRNG(), g, q)
end

"""
    multi_type_spanning_forest(
        rng,
        g,
        q;
        weighted,
        absorbing_node,
        ust,
    )

samples an mtsf with a given random number generator, and possibly weighted links, an absorbing node. For sampling a spanning tree, put absorbing_node = true and ust = true. This algorithm uses cycle-popping with capped cycle weights

# Arguments
- `rng::Random.AbstractRNG` random number generator.
- `g::MetaGraph{T,U}` connection graph.
- `q::Real` positive parameter playing the role of regularizer.
- `weighted::Bool=false` (optional) uses the edge weights in the random walk.
- `absorbing_node::Bool=false` (optional) absorbing node sampled uniformly.
- `ust::Bool=false` (optional) allows to a sample spanning tree.
# Output
- `mtsf::MetaGraph{T,U}` connection graph associated with mtsf.

"""
function multi_type_spanning_forest(
    rng::Random.AbstractRNG,
    g::MetaGraph{T,U},
    q::Real;
    weighted::Bool=false,
    absorbing_node::Bool=false,
    ust::Bool=false,
)::MetaGraph{T,U} where {T,U}
    # Initialize the multi type spanning forest
    mtsf = MetaGraph{T,U}(nv(g))
    nv_mtsf = 0
    weight = 1.0
    roots = T[]
    nodes_in_cycles = Vector{T}[]
    reverse_order_branches = Vector{T}[]

    # Initialize the random walk
    walk = T[]
    unvisited = Set{T}(vertices(g))

    #fix a root if necessary
    if absorbing_node
        ab_node = rand(rng, unvisited)
        push!(roots, ab_node)
        setdiff!(unvisited, ab_node)
        nv_mtsf += 1
    end

    # Start the random walk
    n0 = rand(rng, unvisited)
    push!(walk, n0)
    setdiff!(unvisited, n0)

    while nv_mtsf < nv(g)
        n0_is_root = false
        if q > 1e-10
            #n0_is_root = rand(rng) < q / (q + degree(g, n0))
            n0_is_root = step_to_root(rng, g, n0, q, weighted)
        end

        if n0_is_root
            push!(roots, n0)
            add_edges_from!(mtsf, consecutive_pairs(walk))
            nv_mtsf += length(walk)
            setdiff!(unvisited, walk)

            # record branch w/o  root
            push!(reverse_order_branches, walk[1:(end - 1)])

            n0 = restart_walk_from_unvisited_node!(rng, walk, unvisited)
            continue
        end

        #n1 = rand(rng, neighbors(g, n0))
        n1 = rand_step(rng, g, n0, weighted)

        push!(walk, n1)

        if n1 in unvisited
            setdiff!(unvisited, n1)
            n0 = n1  # continue the walk

        elseif (degree(mtsf, n1) > 0 || n1 in roots) || (n1 in roots && ust)
            add_edges_from!(mtsf, consecutive_pairs(walk))
            nv_mtsf += length(walk) - 1
            setdiff!(unvisited, walk)

            push!(reverse_order_branches, walk[1:(end - 1)])

            n0 = restart_walk_from_unvisited_node!(rng, walk, unvisited)

        else  # if n1 in walk: identify unique cycle/loop in walk with knot n1
            idx_n1 = findfirst(x -> x == n1, walk)
            cycle_nodes = @view walk[idx_n1:end]
            keep = false
            alpha = 0
            if !ust
                keep, alpha = keep_cycle(rng, g, consecutive_pairs(cycle_nodes))
            end

            if keep # cycle
                weight *= max(alpha, 1)
                add_edges_from!(mtsf, consecutive_pairs(walk))
                nv_mtsf += length(walk) - 1 # since walk contains twice the knot
                setdiff!(unvisited, walk)

                # record cycle nodes: remove starting node so that it appears only once
                push!(nodes_in_cycles, cycle_nodes[2:end])
                # record branch without the knot
                push!(reverse_order_branches, walk[1:(idx_n1 - 1)])

                n0 = restart_walk_from_unvisited_node!(rng, walk, unvisited)

            else  # pop cycle but keep loopy node
                union!(unvisited, cycle_nodes)
                setdiff!(unvisited, n1)  # remove n1 which was part of cycle_nodes
                resize!(walk, idx_n1)

                n0 = n1
            end
        end
    end
    set_edges_prop_from!(mtsf, :angle, g, true)
    set_prop!(mtsf, :weight, weight)
    set_prop!(mtsf, :roots, roots)
    set_prop!(mtsf, :cycle_nodes, nodes_in_cycles)
    set_prop!(mtsf, :branches, reverse(reverse_order_branches))

    return mtsf
end

function restart_walk_from_unvisited_node!(rng::Random.AbstractRNG, walk, unvisited)
    empty!(walk)
    n0 = -1
    if !isempty(unvisited)
        n0 = rand(rng, unvisited)
        setdiff!(unvisited, n0)
        push!(walk, n0)
    end
    return n0
end

function keep_cycle(rng::Random.AbstractRNG, graph::AbstractMetaGraph, cycle_edges)
    alpha = 1.0 - cos(curvature(graph, cycle_edges))
    keep = rand(rng) < min(alpha, 1)
    return keep, alpha
end

function curvature(
    g::AbstractMetaGraph, edges, oriented::Bool=true, default_angle::Real=0.0
)
    curvature = 0.0
    for e in edges
        angle = get_edge_prop(g, Edge(e), :angle, oriented, default_angle)
        curvature += angle
    end
    return curvature
end

function rand_step(rng::Random.AbstractRNG, g, n0, weighted::Bool=false)
    if weighted
        nb_list = neighbors(g, n0)
        p = Vector{Float64}(undef, length(nb_list))
        it = 0
        for v in nb_list
            it += 1
            e = (n0, v)
            w = get_edge_prop(g, Edge(e), :e_weight)
            p[it] = abs(w)
        end
        p = p / sum(p)

        ind_nb = rand(rng, Categorical(p), 1)[1]
        n1 = nb_list[ind_nb]
    else
        n1 = rand(rng, neighbors(g, n0))
    end

    return n1
end

function step_to_root(rng::Random.AbstractRNG, g, n0, q::Real, weighted::Bool=false)
    if weighted
        nb_list = neighbors(g, n0)
        sum_of_w = 0
        for v in nb_list
            e = (n0, v)
            w = get_edge_prop(g, Edge(e), :e_weight)
            sum_of_w += abs(w)
        end
        n0_is_root = rand(rng) < q / (q + sum_of_w)
    else
        n0_is_root = rand(rng) < q / (q + degree(g, n0))
    end

    return n0_is_root
end
