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
    multi_type_spanning_forest(rng,g,q;weighted,absorbing_node,ust)

samples an mtsf with a given random number generator, and possibly weighted links, an absorbing node. For sampling a spanning tree, put absorbing_node = true and ust = true. This algorithm uses cycle-popping with capped cycle weights.

The output meta graph has the following properties:
- :weight self-normalized sampling weight of the mtsf.
- :roots array containing the mtsf roots.
- :cycle_nodes array containing the nodes in each of the cycles.
- :branches array containing the branches rooted at the cycles or roots.

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

    # fix an absorbing root if necessary
    if absorbing_node
        ab_node = rand(rng, unvisited)
        push!(roots, ab_node)
        setdiff!(unvisited, ab_node)
        nv_mtsf += 1
    end

    # Start the random walk
    n0 = rand(rng, unvisited)
    # add n0 to walk
    push!(walk, n0)
    # mark n0 as visited
    setdiff!(unvisited, n0)

    while nv_mtsf < nv(g)
        n0_is_root = false
        # check if n0 is a root
        if q > 1e-10
            # if q not too small
            n0_is_root = step_to_root(rng, g, n0, q, weighted)
            # explicitly n0_is_root = rand(rng) < q / (q + degree(g, n0))
        end

        if n0_is_root
            # if n0 is indeed a root record the rooted branch in mtsf
            push!(roots, n0)
            add_edges_from!(mtsf, consecutive_pairs(walk))
            nv_mtsf += length(walk)
            # mark nodes in walk as visited
            setdiff!(unvisited, walk)

            # record branch w/o  root
            push!(reverse_order_branches, walk[1:(end - 1)])
            # restart with a new n0 uniformly among the unvisited nodes
            n0 = restart_walk_from_unvisited_node!(rng, walk, unvisited)
            continue
        end

        # walk to a nb with proba propto edge weight
        n1 = rand_step(rng, g, n0, weighted)
        # add this nb to walk
        push!(walk, n1)

        if n1 in unvisited
            # if n1 unvisited, mark it as visited
            setdiff!(unvisited, n1)
            n0 = n1  # and continue the walk

        elseif (degree(mtsf, n1) > 0 || n1 in roots) || (n1 in roots && ust)
            # otherwise if n1 is already visited and is in the mtsf or is a root
            # record the walk as a branch
            add_edges_from!(mtsf, consecutive_pairs(walk))
            nv_mtsf += length(walk) - 1
            setdiff!(unvisited, walk)
            # mark nodes in walk as visited

            push!(reverse_order_branches, walk[1:(end - 1)])
            # restart with a new n0 uniformly among the unvisited nodes
            n0 = restart_walk_from_unvisited_node!(rng, walk, unvisited)

        else  # else n1 is in walk.
            # We identify unique cycle in walk with knot n1
            idx_n1 = findfirst(x -> x == n1, walk)
            cycle_nodes = @view walk[idx_n1:end]

            keep = false # by default, cycle is popped
            alpha = 0
            if !ust # if not spanning tree, toss a coin to pop or not
                keep, alpha = keep_cycle(rng, g, consecutive_pairs(cycle_nodes))
            end

            if keep # cycle is kept
                # compute cycle weight in view of reweighted MC
                weight *= max(alpha, 1)
                add_edges_from!(mtsf, consecutive_pairs(walk))
                nv_mtsf += length(walk) - 1 # since walk contains twice the knot
                setdiff!(unvisited, walk)

                # record cycle nodes: remove starting node so that it appears only once
                push!(nodes_in_cycles, cycle_nodes[2:end])
                # record branch without the knot
                push!(reverse_order_branches, walk[1:(idx_n1 - 1)])
                # restart with a new n0 uniformly among the unvisited nodes
                n0 = restart_walk_from_unvisited_node!(rng, walk, unvisited)

            else  # pop cycle but keep loopy node
                union!(unvisited, cycle_nodes)
                # mark cycle nodes as unvisited
                setdiff!(unvisited, n1)
                # remove n1 which was part of cycle_nodes
                resize!(walk, idx_n1)
                # remove cycle nodes from the walk
                n0 = n1
            end
        end
    end
    # store outputs
    set_edges_prop_from!(mtsf, :angle, g, true)
    set_prop!(mtsf, :weight, weight)
    set_prop!(mtsf, :roots, roots)
    set_prop!(mtsf, :cycle_nodes, nodes_in_cycles)
    set_prop!(mtsf, :branches, reverse(reverse_order_branches))

    return mtsf
end

"""
    restart_walk_from_unvisited_node!(rng, walk, unvisited)

returns an unvisited node.

# Arguments
- `rng::Random.AbstractRNG` random number generator.
- `walk` array containing the current path of the walk.
- `unvisited` array indicating which node is visited/unvisited/in the path.
# Output
- `n0` an unvisited node.

"""
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

"""
    keep_cycle(rng, graph, cycle_edges)

keeps or pops a cycle.

# Arguments
- `rng::Random.AbstractRNG` random number generator.
- `graph` connection graph.
- `cycle_edges` array of nodes in the cycle.

# Output
- `keep` boolean determining if the cycle is kept (true) or not (false).
- `alpha` cycle weight.

"""
function keep_cycle(rng::Random.AbstractRNG, graph::AbstractMetaGraph, cycle_edges)
    alpha = 1.0 - cos(curvature(graph, cycle_edges))
    keep = rand(rng) < min(alpha, 1)
    return keep, alpha
end

function optim_keep_cycle(rng::Random.AbstractRNG, Laplacian, cycle_edges)
    alpha = 1.0 - real(optim_curvature(Laplacian, cycle_edges))
    keep = rand(rng) < min(alpha, 1)
    return keep, alpha
end

function optim_curvature(Laplacian, edges)
    curvature = 1.0
    for e in edges
        phase = -Laplacian[e[1], e[2]]
        curvature *= phase
    end
    return curvature
end

"""
    curvature(g, edges, oriented, default_angle)

sum of the angles along a cycle.
# Arguments
- `g::AbstractMetaGraph` connection graph.
- `edges` cycle edges.
- `oriented` boolean indicating if edges are oriented.
- `default_angle::Real=0.0` angle by default.

# Output
- `curvature` sum of angles along cycle edges.

"""
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

function crsf(rng::Random.AbstractRNG, g::MetaGraph{T,U}) where {T,U}
    n = nv(g)
    # Initialize outputs
    node_list = zeros(Bool, n)
    edge_list = zeros(Int, n, 2)
    nv_subgraph = 0
    ne_subgraph = 0
    weight = 1.0

    # Initialize the random walk
    stps_in_walk = 0
    # walk = zeros(Int, n) # largest possible cycle contains n nodes
    walk = T[]
    visited = zeros(Bool, n)

    # Start the random walk
    # v0 = minimum(ordering[.!visited])
    v0 = findfirst(x -> x == 0, visited)

    # println("v0: ", v0)
    stps_in_walk += 1
    # add v0 to walk
    # walk[stps_in_walk] = v0
    push!(walk, v0)

    # mark v0 as visited
    visited[v0] = 1
    # println("walk: ", walk)

    while nv_subgraph < n
        # walk to a nb with proba propto edge weight
        v1 = rand(rng, neighbors(g, v0))
        # println("v1: ", v1)
        stps_in_walk += 1
        # add this nb to walk
        # walk[stps_in_walk] = v1
        push!(walk, v1)

        # println("walk: ", walk)

        if !visited[v1] # v1 unvisited,
            # println("not visited")
            # mark it as visited
            visited[v1] = 1
            v0 = v1  # and continue the walk

        elseif node_list[v1] == 1  # v1 is in the subgraph
            # println("in subgraph")
            # mark branch in subgraph

            # branch = @view walk[1:stps_in_walk]
            branch = walk
            node_list[branch] .= 1
            # println("store branch in node list: ", branch)

            # increment nb vertices in subgraph
            nv_subgraph += stps_in_walk - 1 # since last node is in subgraph

            # mark nodes in branch as visited
            visited[branch] .= 1

            # add edges to elist
            for pair in consecutive_pairs(branch)
                edge_list[ne_subgraph + 1, 1] = pair[1]
                edge_list[ne_subgraph + 1, 2] = pair[2]
                ne_subgraph += 1
            end

            if nv_subgraph < n
                # restart with a new v0
                # v0 = minimum(ordering[.!visited])
                v0 = findfirst(x -> x == 0, visited)

                # println("restart from :", v0)
                visited[v0] = 1
                # empty walk
                # walk[1:stps_in_walk] .= 0
                empty!(walk)
                stps_in_walk = 1
                # walk[stps_in_walk] = v0
                push!(walk, v0)
                # println("walk :", walk)
            end

        else  # else v1 is in walk.
            # We identify unique cycle in walk with knot v1
            # println("v1 makes a cycle")
            cr_branch = @view walk[1:stps_in_walk]
            # println("cycle-rooted branch: ", cr_branch)
            idx_v1 = findfirst(x -> x == v1, cr_branch)
            cycle_nodes = @view cr_branch[idx_v1:stps_in_walk]
            # println("nodes in cycle", cycle_nodes)
            keep, alpha = keep_cycle(rng, g, consecutive_pairs(cycle_nodes))

            if keep # cycle is kept
                # compute cycle weight in view of reweighted MC
                weight *= max(alpha, 1)
                # println("keep")
                for pair in consecutive_pairs(cr_branch)
                    edge_list[ne_subgraph + 1, 1] = pair[1]
                    edge_list[ne_subgraph + 1, 2] = pair[2]
                    ne_subgraph += 1
                end

                # mark cr branch as in subgraph
                node_list[cr_branch] .= 1
                nv_subgraph += length(cr_branch) - 1 # since walk contains twice the knot
                # mark cr branch as visited
                visited[cr_branch] .= 1

                if nv_subgraph < n
                    # restart with a new v0
                    # walk[1:stps_in_walk] .= 0
                    empty!(walk)
                    v0 = findfirst(x -> x == 0, visited)
                    # v0 = minimum(ordering[.!visited])
                    # println("restart from: ", v0)
                    visited[v0] = 1
                    stps_in_walk = 1
                    # walk[stps_in_walk] = v0
                    push!(walk, v0)
                    # println("walk: ", walk)
                end

            else  # pop cycle but keep loopy node
                # println("pop")
                visited[cycle_nodes[1:(end - 1)]] .= 0
                # mark cycle nodes as visited
                visited[v1] = 1
                # v1 was part of cycle_nodes
                resize!(walk, idx_v1)
                # walk[(idx_v1 + 1):stps_in_walk] .= 0
                # println("walk: ", walk)
                stps_in_walk = idx_v1
                # remove cycle nodes from the walk
                v0 = v1
                # println("restart from knot: ", v0)
            end
        end
    end

    return node_list, edge_list, weight
end

function simple_multi_type_spanning_forest(
    rng::Random.AbstractRNG,
    g::MetaGraph{T,U},
    q::Real;
    weighted::Bool=false,
    absorbing_node::Bool=false,
    ust::Bool=false,
) where {T,U} #::MetaGraph{T,U} where {T,U}
    # Initialize the multi type spanning forest
    mtsf = MetaGraph{T,U}(nv(g))
    nv_mtsf = 0
    weight = 1.0
    roots = T[]

    # Initialize the random walk
    walk = T[]
    unvisited = Set{T}(vertices(g))

    # fix an absorbing root if necessary
    if absorbing_node
        ab_node = rand(rng, unvisited)
        push!(roots, ab_node)
        setdiff!(unvisited, ab_node)
        nv_mtsf += 1
    end

    # Start the random walk
    n0 = rand(rng, unvisited)
    # add n0 to walk
    push!(walk, n0)
    # mark n0 as visited
    setdiff!(unvisited, n0)

    while nv_mtsf < nv(g)
        n0_is_root = false
        # check if n0 is a root
        if q > 1e-10
            # if q not too small
            n0_is_root = step_to_root(rng, g, n0, q, weighted)
            # explicitly n0_is_root = rand(rng) < q / (q + degree(g, n0))
        end

        if n0_is_root
            # if n0 is indeed a root record the rooted branch in mtsf
            push!(roots, n0)
            add_edges_from!(mtsf, consecutive_pairs(walk))
            nv_mtsf += length(walk)
            # mark nodes in walk as visited
            setdiff!(unvisited, walk)

            # restart with a new n0 uniformly among the unvisited nodes
            n0 = restart_walk_from_unvisited_node!(rng, walk, unvisited)
            continue
        end

        # walk to a nb with proba propto edge weight
        n1 = rand_step(rng, g, n0, weighted)
        # add this nb to walk
        push!(walk, n1)

        if n1 in unvisited
            # if n1 unvisited, mark it as visited
            setdiff!(unvisited, n1)
            n0 = n1  # and continue the walk

        elseif (degree(mtsf, n1) > 0 || n1 in roots) || (n1 in roots && ust)
            # otherwise if n1 is already visited and is in the mtsf or is a root
            # record the walk as a branch
            add_edges_from!(mtsf, consecutive_pairs(walk))
            nv_mtsf += length(walk) - 1
            setdiff!(unvisited, walk)
            # mark nodes in walk as visited
            # restart with a new n0 uniformly among the unvisited nodes
            n0 = restart_walk_from_unvisited_node!(rng, walk, unvisited)

        else  # else n1 is in walk.
            # We identify unique cycle in walk with knot n1
            idx_n1 = findfirst(x -> x == n1, walk)
            cycle_nodes = @view walk[idx_n1:end]

            keep = false # by default, cycle is popped
            alpha = 0
            if !ust # if not spanning tree, toss a coin to pop or not
                keep, alpha = keep_cycle(rng, g, consecutive_pairs(cycle_nodes))
            end

            if keep # cycle is kept
                # compute cycle weight in view of reweighted MC
                weight *= max(alpha, 1)
                add_edges_from!(mtsf, consecutive_pairs(walk))
                nv_mtsf += length(walk) - 1 # since walk contains twice the knot
                setdiff!(unvisited, walk)

                # restart with a new n0 uniformly among the unvisited nodes
                n0 = restart_walk_from_unvisited_node!(rng, walk, unvisited)

            else  # pop cycle but keep loopy node
                union!(unvisited, cycle_nodes)
                # mark cycle nodes as unvisited
                setdiff!(unvisited, n1)
                # remove n1 which was part of cycle_nodes
                resize!(walk, idx_n1)
                # remove cycle nodes from the walk
                n0 = n1
            end
        end
    end
    # store outputs
    # set_edges_prop_from!(mtsf, :angle, g, true)
    # set_prop!(mtsf, :weight, weight)

    return mtsf, weight
end
