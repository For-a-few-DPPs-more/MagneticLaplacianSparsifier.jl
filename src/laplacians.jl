
function spMagneticIncidence(compGraph)
    n = nv(compGraph)
    m = ne(compGraph)

    B = spzeros(Complex, n, m)
    edges_compGraph = edges(compGraph)

    ind_e = 0
    for e in edges_compGraph
        ind_e += 1
        u = src(e)
        v = dst(e)
        angle = get_edge_prop(compGraph, e, :angle)
        B[u, ind_e] = exp(-0.5 * angle * im)
        B[v, ind_e] = -exp(0.5 * angle * im)
    end

    return B
end

function magneticIncidence(compGraph)::Array{Complex,2}
    n = nv(compGraph)
    m = ne(compGraph)

    B = zeros(Complex, n, m)
    edges_compGraph = edges(compGraph)

    ind_e = 0
    for e in edges_compGraph
        ind_e += 1
        u = src(e)
        v = dst(e)
        angle = get_edge_prop(compGraph, e, :angle)
        B[u, ind_e] = exp(-0.5 * angle * im)
        B[v, ind_e] = -exp(0.5 * angle * im)
    end

    return B
end

function mtsf_edge_indices(crsf, compGraph)
    edges_compGraph = edges(compGraph)
    edges_crsf = edges(crsf)

    ind_e = []
    it = 0
    for e in edges_compGraph
        it += 1
        if e in edges_crsf
            push!(ind_e, it)
        end
    end

    return ind_e
end
