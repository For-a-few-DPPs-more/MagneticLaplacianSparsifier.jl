
function magneticIncidence(compGraph)
    n = nv(compGraph)
    m = ne(compGraph)

    B = spzeros(Complex,n,m)
    edges_compGraph = edges(compGraph)

    ind_e = 0
    for e in edges_compGraph
        ind_e += 1
        u = src(e)
        v = dst(e)
        angle = get_edge_property(compGraph, e, :angle)
        B[u, ind_e] = exp(-0.5 * angle*im)
        B[v, ind_e] = -exp(0.5 * angle*im)
    end

    return B
end