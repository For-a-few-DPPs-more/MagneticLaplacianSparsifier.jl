
function sp_magnetic_incidence(graph; oriented::Bool=false)
    return magnetic_incidence_matrix(graph; oriented=oriented)
end

function magnetic_incidence(graph; oriented::Bool=false)::Matrix{Complex{Float64}}
    return Array(sp_magnetic_incidence(graph; oriented=oriented))
end

function magnetic_incidence_matrix(
    graph::AbstractGraph; oriented::Bool=true
)::SparseMatrixCSC{Complex{Float64}}
    #B = spzeros(Complex, nv(graph), ne(graph))
    #for (idx_e, e) in enumerate(edges(graph))
    #    u = src(e)
    #    v = dst(e)
    #    angle = get_edge_prop(graph, e, :angle)
    #    print(angle)
    #    B[u, idx_e] = exp(0.5 * angle * im)
    #    B[v, idx_e] = -exp(-0.5 * angle * im)
    #end
    #return B
    n_v, n_e = nv(graph), ne(graph)

    I = vcat(src.(edges(graph)), dst.(edges(graph)))
    J = vcat(1:n_e, 1:n_e)

    θ = get_edges_prop(graph, :angle, true, 0.0)
    w = @. exp(im * 0.5 * θ)
    V = vcat(oriented ? -conj.(w) : w, w)
    return sparse(I, J, V, n_v, n_e)
end

# todo naming mtsf, csrf
function mtsf_edge_indices(crsf, graph)
    return [i for (i, e) in enumerate(edges(graph)) if has_edge(crsf, src(e), dst(e))]
end

function average_sparsifier(rng, compGraph, ls, useLS, q, t)
    n = nv(compGraph)
    m = ne(compGraph)
    sparseL = zeros(n, n)
    w_tot = 0

    for i in 1:t
        crsf = multi_type_spanning_forest(rng, compGraph, q)
        D = props(crsf)
        w = D[:weight]
        w_tot += w
        sparseB = magnetic_incidence(crsf; oriented=true)
        ind_e = mtsf_edge_indices(crsf, compGraph)
        if useLS
            W = diagm(1 ./ ls[ind_e])
        else
            nb_e = length(ind_e)
            W = I / (nb_e / m)
        end
        sparseL = sparseL + w * sparseB * W * sparseB'
    end
    sparseL = sparseL / w_tot

    return sparseL
end

function leverage_score(B, q)
    levScores = real(diag(B' * ((B * B' + q * I) \ B)))
    return levScores
end

function emp_leverage_score(rng, compGraph, q, t)
    m = ne(compGraph)
    empLev = zeros(m, 1)

    for i in 1:t
        crsf = multi_type_spanning_forest(rng, compGraph, q)
        ind_e = mtsf_edge_indices(crsf, compGraph)
        empLev[ind_e] = empLev[ind_e] .+ 1
    end
    empLev = empLev / t

    return empLev
end

function nb_of_edges(L)
    n = size(L)[1]
    nb_e = ((nnz(sparse(L)) - n) / 2)
    return nb_e
end

# Refactoring Proposition

# function average_sparsifier(rng, compGraph, ls, useLS, q, t)
#     return average_sparsifier_bis(rng, compGraph, q, t, useLS ? ls : nothing)
# end

# function leverage_score(B, q)
#     return leverage_scores_from_incidence(B, q)
# end

# function emp_leverage_score(rng, compGraph, q, t)
#     ls_dict = leverage_scores_empirical(rng, compGraph; q=q, nb_samples=t)
#     return collect(values(ls_dict))
# end

# nb_of_edges(L::AbstractMatrix) = (nnz(sparse(L)) - size(L, 1)) / 2

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
#         crsf = multi_type_spanning_forest(rng, graph, q)
#         w_csrf = get_prop(crsf, :weight)
#         weight_total += w_csrf
#         B = magneticIncidence(crsf; oriented=true)
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
#             leverage_scores[e] += 1 / nb_samples
#         end
#     end
#     return leverage_scores
# end
