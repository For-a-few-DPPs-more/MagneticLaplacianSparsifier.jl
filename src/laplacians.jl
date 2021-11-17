
function spMagneticIncidence(graph; oriented::Bool=false)
    return magnetic_incidence_matrix(graph; oriented=oriented)
end

function magneticIncidence(graph; oriented::Bool=false)::Matrix{Complex{Float64}}
    return Array(spMagneticIncidence(graph; oriented=oriented))
end

function magnetic_incidence_matrix(
    graph::AbstractGraph; oriented::Bool=false
)::SparseMatrixCSC{Complex{Float64}}
    # B = spzeros(Complex, nv(graph), ne(graph))
    # for (idx_e, e) in enumerate(edges(graph))
    #     u = src(e)
    #     v = dst(e)
    #     angle = get_edge_prop(graph, e, :angle)
    #     B[u, idx_e] = exp(0.5 * angle * im)
    #     B[v, idx_e] = -exp(-0.5 * angle * im)
    # end

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
