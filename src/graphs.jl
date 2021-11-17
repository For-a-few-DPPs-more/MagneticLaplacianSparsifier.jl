
function generateGraphMUN(n, p, eta)
    # Following M. Cucuringu to build comparison graph
    # SYNC-RANK: ROBUST RANKING, CONSTRAINED RANKING AND RANK AGGREGATION VIA EIGENVECTOR AND SDP SYNCHRONIZATION
    # Multiplicative Uniform Noise model
    g = Graph(n)
    meta_g = MetaGraph(g, :angle, 0.0)

    # ranking
    r = collect(1:n)

    for i in r
        for j in r
            if i < j
                if rand(1)[1] < p
                    e = [i j]
                    edges = consecutive_pairs(e)
                    add_edges_from!(meta_g, edges)

                    err = eta * rand(1)[1]
                    a = pi * (i - j) * (1 + err) / (n - 1)
                    for e in edges
                        set_prop!(meta_g, Edge(e), :angle, a)
                    end
                end
            end
        end
    end

    return meta_g
end

function generateGraphERO(n, p, eta)
    # Following M. Cucuringu to build comparison graph
    # SYNC-RANK: ROBUST RANKING, CONSTRAINED RANKING AND RANK AGGREGATION VIA EIGENVECTOR AND SDP SYNCHRONIZATION
    # Erdos-Renyi Outliers model
    g = Graph(n)
    meta_g = MetaGraph(g, :angle, 0.0)

    # ranking
    r = collect(1:n)

    for i in r
        for j in r
            if i < j
                if rand(1)[1] < p
                    e = [i j]
                    edges = consecutive_pairs(e)
                    add_edges_from!(meta_g, edges)

                    if rand(1)[1] < 1 - eta
                        err = 0
                    else
                        err = (n - 1) * rand(1)[1]
                    end
                    a = pi * (i - j + err) / (n - 1)
                    for e in edges
                        set_prop!(meta_g, Edge(e), :angle, a)
                    end
                end
            end
        end
    end

    return meta_g
end
