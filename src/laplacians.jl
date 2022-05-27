
function sp_magnetic_incidence(graph; oriented::Bool=true)
    return magnetic_incidence_matrix(graph; oriented=oriented)
end

function magnetic_incidence(graph; oriented::Bool=true)::Matrix{Complex{Float64}}
    return Array(sp_magnetic_incidence(graph; oriented=oriented))
end

function magnetic_incidence_matrix(
    graph::AbstractGraph; oriented::Bool=true, phases::Bool=true
)::SparseMatrixCSC{Complex{Float64}}
    n_v, n_e = nv(graph), ne(graph)

    I = vcat(src.(edges(graph)), dst.(edges(graph)))
    J = vcat(1:n_e, 1:n_e)

    θ = get_edges_prop(graph, :angle, true, 0.0)
    w = @. exp(-im * 0.5 * θ)
    if !phases
        w = abs.(w)
    end
    V = vcat(w, oriented ? -conj.(w) : conj.(w))
    # here different sign wrt paper but unimportant
    B = transpose(sparse(I, J, V, n_v, n_e))
    return B
end

# todo naming mtsf, csrf
function mtsf_edge_indices(mtsf, graph)
    return [i for (i, e) in enumerate(edges(graph)) if has_edge(mtsf, src(e), dst(e))]
end

function average_sparsifier(
    rng::Random.AbstractRNG,
    meta_g::AbstractMetaGraph,
    ls::Union{Array,Nothing},
    q::Real,
    nb_samples::Integer;
    weighted::Bool=false,
    absorbing_node::Bool=false,
    ust::Bool=false,
)
    n = nv(meta_g)
    m = ne(meta_g)
    L = zeros(n, n)
    nb_cycles = zeros(nb_samples, 1)
    nb_roots = zeros(nb_samples, 1)
    weights = zeros(nb_samples, 1)
    w_tot = 0

    for i_sample in 1:nb_samples
        mtsf = multi_type_spanning_forest(rng, meta_g, q; weighted, absorbing_node, ust)
        # check nb roots and cycles
        cycles = get_prop(mtsf, :cycle_nodes)
        nb_cycles[i_sample] = length(cycles)
        roots = get_prop(mtsf, :roots)
        nb_roots[i_sample] = length(roots)

        D = props(mtsf)
        w = D[:weight]
        weights[i_sample] = w

        w_tot += w
        sparseB = magnetic_incidence(mtsf; oriented=true)
        ind_e = mtsf_edge_indices(mtsf, meta_g)
        if ls === nothing
            nb_e = length(ind_e)
            W = I / (nb_e / m)
        else
            W = diagm(1 ./ ls[ind_e])
        end
        if weighted
            e_weights = get_edges_prop(meta_g, :e_weight, true, 1.0)
            W *= diagm(e_weights[ind_e])
        end

        L = L + w * sparseB' * W * sparseB
    end
    nb_sampled_cycles = sum(nb_cycles)
    nb_sampled_roots = sum(nb_roots)

    L = L / w_tot

    return L, nb_sampled_cycles, nb_sampled_roots, weights
end

function sample_subgraph_iid(
    rng::Random.AbstractRNG,
    meta_g::AbstractMetaGraph,
    ls::Union{Array,Nothing},
    batch::Integer,
)
    n = nv(meta_g)
    m = ne(meta_g)
    subgraph = MetaGraph(n)
    if ls === nothing
        ind_rd = rand(rng, 1:m, (batch, 1))
    else
        p = vec(ls / sum(ls))
        ind_rd = rand(rng, Categorical(p), (batch, 1))
    end
    all_edges = collect(edges(meta_g))
    subset_edges = all_edges[ind_rd]

    for e in subset_edges
        add_edge!(subgraph, e)
        angle = get_edge_prop(meta_g, e, :angle, true)
        set_prop!(subgraph, e, :angle, angle)
    end

    return subgraph
end

function average_sparsifier_iid(
    rng::Random.AbstractRNG,
    meta_g::AbstractMetaGraph,
    ls::Union{Array,Nothing},
    batch::Integer,
    nb_samples::Integer;
    weighted::Bool=false,
)
    n = nv(meta_g)
    m = ne(meta_g)
    L = zeros(n, n)
    w_tot = 0

    for _ in 1:nb_samples
        subgraph = sample_subgraph_iid(rng::Random.AbstractRNG, meta_g, ls, batch)
        w = 1
        w_tot += w
        sparseB = magnetic_incidence(subgraph; oriented=true)
        ind_e = mtsf_edge_indices(subgraph, meta_g)
        if ls === nothing
            nb_e = length(ind_e)
            W = I / (nb_e / m)
        else
            W = diagm(1 ./ ls[ind_e])
        end

        if weighted
            e_weights = edge_weights(meta_g)
            W *= diagm(e_weights[ind_e])
        end

        L = L + w * sparseB' * W * sparseB
    end
    L = L / w_tot

    return L
end

function leverage_score(B::Array, q::Real; W=I::Array)
    if q > 1e-13
        levScores = real(diag(W * B * ((B' * W * B + q * I) \ B')))
    else
        levScores = real(diag(W * B * pinv(B' * W * B + q * I) * B'))
    end

    return levScores
end

function emp_leverage_score(
    rng::Random.AbstractRNG,
    meta_g::AbstractMetaGraph,
    q::Real,
    t::Integer;
    weighted::Bool=false,
    absorbing_node::Bool=false,
    ust::Bool=false,
)
    m = ne(meta_g)
    emp_lev = zeros(m, 1)
    for _ in 1:t
        mtsf = multi_type_spanning_forest(rng, meta_g, q; weighted, absorbing_node, ust)
        ind_e = mtsf_edge_indices(mtsf, meta_g)
        emp_lev[ind_e] = emp_lev[ind_e] .+ 1
    end
    emp_lev /= t

    return emp_lev
end

function edge_weights(g)
    e_weights = get_edges_prop(g, :e_weight, true, 1.0)
    return e_weights
end

nb_of_edges(L::AbstractMatrix) = (nnz(sparse(L)) - size(L, 1)) / 2

function optimal_perm(mtsf::AbstractMetaGraph)
    n_v = nv(mtsf)
    #get the roots
    roots = get_prop(mtsf, :roots)
    # get the branches in the (reverse) order there were sampled
    branches = get_prop(mtsf, :branches)
    flt_branches = collect(Iterators.flatten(branches))

    # indices array putting the nodes in the right order
    non_cycle_nodes = [roots; flt_branches]
    ind_perm = [non_cycle_nodes; setdiff(1:n_v, non_cycle_nodes)]
    return ind_perm
end
