
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
    w = @. exp(im * 0.5 * θ)
    if !phases
        w = abs.(w)
    end
    V = vcat(oriented ? -w : w, conj.(w))
    # here different sign wrt paper but unimportant
    return sparse(I, J, V, n_v, n_e)
end

# todo naming mtsf, csrf
function mtsf_edge_indices(mtsf, graph)
    return [i for (i, e) in enumerate(edges(graph)) if has_edge(mtsf, src(e), dst(e))]
end

function average_sparsifier(
    rng,
    meta_g,
    ls,
    q,
    nb_samples;
    weighted::Bool=false,
    absorbing_node::Bool=false,
    ust::Bool=false,
)
    n = nv(meta_g)
    m = ne(meta_g)
    L = zeros(n, n)
    nb_cycles = zeros(nb_samples, 1)
    nb_roots = zeros(nb_samples, 1)
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

        L = L + w * sparseB * W * sparseB'
    end
    nb_sampled_cycles = sum(nb_cycles)
    nb_sampled_roots = sum(nb_roots)

    L = L / w_tot

    return L, nb_sampled_cycles, nb_sampled_roots
end

# function average_sparsifier!(
#     L,
#     rng,
#     meta_g,
#     ls,
#     q,
#     nb_samples;
#     weighted::Bool=false,
#     absorbing_node::Bool=false,
#     ust::Bool=false,
# )
#     n = nv(meta_g)
#     m = ne(meta_g)
#     w_tot = 0

#     for _ in 1:nb_samples
#         mtsf = multi_type_spanning_forest(rng, meta_g, q; weighted, absorbing_node, ust)
#         D = props(mtsf)
#         w = D[:weight]
#         w_tot += w
#         sparseB = magnetic_incidence(mtsf; oriented=true)
#         ind_e = mtsf_edge_indices(mtsf, meta_g)
#         if ls === nothing
#             nb_e = length(ind_e)
#             W = I / (nb_e / m)
#         else
#             W = diagm(1 ./ ls[ind_e])
#         end
#         if weighted
#             e_weights = get_edges_prop(meta_g, :e_weight, true, 1.0)
#             W *= diagm(e_weights[ind_e])
#         end

#         L = L + w * sparseB * W * sparseB'
#     end
#     return L = L / w_tot
# end

function sample_subgraph_iid(rng, meta_g, ls, batch)
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

function average_sparsifier_iid(rng, meta_g, ls, batch, nb_samples; weighted::Bool=false)
    n = nv(meta_g)
    m = ne(meta_g)
    L = zeros(n, n)
    w_tot = 0

    for _ in 1:nb_samples
        subgraph = sample_subgraph_iid(rng, meta_g, ls, batch)
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

        L = L + w * sparseB * W * sparseB'
    end
    L = L / w_tot

    return L
end

# function average_sparsifier_iid!(
#     L, rng, meta_g, ls, batch, nb_samples; weighted::Bool=false
# )
#     n = nv(meta_g)
#     m = ne(meta_g)
#     L = zeros(n, n)
#     w_tot = 0

#     for _ in 1:nb_samples
#         subgraph = sample_subgraph_iid(rng, meta_g, ls, batch)
#         w = 1
#         w_tot += w
#         sparseB = magnetic_incidence(subgraph; oriented=true)
#         ind_e = mtsf_edge_indices(subgraph, meta_g)
#         if ls === nothing
#             nb_e = length(ind_e)
#             W = I / (nb_e / m)
#         else
#             W = diagm(1 ./ ls[ind_e])
#         end

#         if weighted
#             e_weights = edge_weights(meta_g)
#             W *= diagm(e_weights[ind_e])
#         end

#         L = L + w * sparseB * W * sparseB'
#     end
#     return L = L / w_tot
# end

function leverage_score(B, q; W=I)
    if q > 1e-13
        levScores = real(diag(W * B' * ((B * W * B' + q * I) \ B)))
    else
        levScores = real(diag(W * B' * pinv(B * W * B' + q * I) * B))
    end

    return levScores
end

function emp_leverage_score(
    rng, meta_g, q, t; weighted::Bool=false, absorbing_node::Bool=false, ust::Bool=false
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

function optimal_perm(mtsf)
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
# Refactoring Proposition

# function average_sparsifier(rng, meta_g, ls, useLS, q, t)
#     return average_sparsifier_bis(rng, meta_g, q, t, useLS ? ls : nothing)
# end

# function average_sparsifier_bis(
#     rng,
#     graph::AbstractMetaGraph,
#     q::Real=0.0,
#     nb_samples::Integer=1,
#     leverage_scores::Union{Dict{Edge,Float64},Nothing}=nothing,
# )
#     L = spzeros(ComplexF64, nv(graph), nv(graph))
#     weight_total = 0.0
#     for _ in 1:nb_samples
#         mtsf = multi_type_spanning_forest(rng, graph, q)
#         w_csrf = get_prop(mtsf, :weight)
#         weight_total += w_csrf
#         B = magneticIncidence(mtsf; oriented=true)
#         # L += B * W * B'
#         if isnothing(leverage_scores)
#             w = w_csrf * ne(graph) / ne(csrf)
#             mul!(L, B, B', w, true)
#         else
#             W = diagm([w_csrf / leverage_scores[e] for e in edges(csrf)])
#             mul!(L, mul(B, W), B', true, true)
#         end
#     end
#     L /= weight_total
#     return L
# end

# function leverage_score(B, q)
#     return leverage_scores_from_incidence(B, q)
# end

# function emp_leverage_score(rng, meta_g, q, t)
#     ls_dict = leverage_scores_empirical(rng, meta_g; q=q, nb_samples=t)
#     return collect(values(ls_dict))
# end

# nb_of_edges(L::AbstractMatrix) = (nnz(sparse(L)) - size(L, 1)) / 2

# function leverage_scores_from_incidence(B::AbstractMatrix, q::Real=0.0)
#     return real(diag(B' * ((B * B' + q * I) \ B)))
# end

# function leverage_scores(g::AbstractGraph; q::Real=0.0)::Dict{Edge,Float64}
#     B = magnetic_incidence_matrix(g; oriented=true)
#     vals = leverage_scores_from_incidence(B, q)
#     leverage_scores = Dict(edges(g) .=> vals)
#     return leverage_scores
# end

# function leverage_scores_empirical(
#     rng::Random.AbstractRNG, graph::AbstractGraph; q::Real=0.0, nb_samples::Integer=1
# )::Dict{Edge,Float64}
#     leverage_scores = Dict{Edge,Float64}(edges(graph) .=> 0.0)
#     for _ in 1:nb_samples
#         mtsf = multi_type_spanning_forest(rng, graph, q)
#         for e in edges(mtsf)
#             leverage_scores[e] += 1.0
#         end
#     end
#     map!(x->x/nb_samples, values(leverage_scores))
#     return leverage_scores
# end
