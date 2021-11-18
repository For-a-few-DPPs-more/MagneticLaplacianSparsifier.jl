
function spMagneticIncidence(graph; oriented::Bool=false)
    return magnetic_incidence_matrix(graph; oriented=oriented)
end

function magneticIncidence(graph; oriented::Bool=false)::Matrix{Complex{Float64}}
    return Array(spMagneticIncidence(graph; oriented=oriented))
end

function magnetic_incidence_matrix(
    graph::AbstractGraph; oriented::Bool=false
)::SparseMatrixCSC{Complex{Float64}}
    B = spzeros(Complex, nv(graph), ne(graph))
    for (idx_e, e) in enumerate(edges(graph))
        u = src(e)
        v = dst(e)
        angle = get_edge_prop(graph, e, :angle)
        B[u, idx_e] = exp(0.5 * angle * im)
        B[v, idx_e] = -exp(-0.5 * angle * im)
    end
    return B
    #n_v, n_e = nv(graph), ne(graph)
    #
    #I = vcat(src.(edges(graph)), dst.(edges(graph)))
    #J = vcat(1:n_e, 1:n_e)
    #
    #θ = get_edges_prop(graph, :angle, true, 0.0)
    #w = @. exp(im * 0.5 * θ)
    #V = vcat(oriented ? -conj.(w) : w, w)
    #return sparse(I, J, V, n_v, n_e)
end

# todo naming mtsf, csrf
function mtsf_edge_indices(crsf, graph)
    return [i for (i, e) in enumerate(edges(graph)) if has_edge(crsf, src(e), dst(e))]
end

function averageSparsifier(rng, compGraph, ls, useLS, q, t)
    n = nv(compGraph)
    m = ne(compGraph)
    sparseL = zeros(n, n)
    w_tot = 0

    for i in 1:t
        crsf = multi_type_spanning_forest(rng, compGraph, q)
        D = props(crsf)
        w = D[:weight]
        w_tot += w
        sparseB = magneticIncidence(crsf)
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

function leverageScore(B, q)
    levScores = real(diag(B' * ((B * B' + q * I) \ B)))
    return levScores
end

function empLeverageScore(rng, compGraph, q, t)
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
    nb_e = ((nnz(sparse(L))-n)/2);
    return nb_e
end
