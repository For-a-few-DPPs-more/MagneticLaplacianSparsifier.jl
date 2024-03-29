angular_score(v) = mod.(angle.(v), 2 * pi)
modulus_entries(v) = abs.(v)

function nb_upsets(meta_g, ranking)
    oriented = true
    upsets = 0
    for e in edges(meta_g)
        a = get_edge_prop(meta_g, e, :angle, oriented)
        d = ranking[src(e)] - ranking[dst(e)]
        # upset if score assignment contradicts the pairwise comparison
        if sign(a) * sign(d) > 0
            upsets += 1
        end
    end
    return upsets
end

function nb_upsets_in_top(meta_g, ranking, k)
    oriented = true
    upsets = 0
    for e in edges(meta_g)
        a = get_edge_prop(meta_g, e, :angle, oriented)
        r_scr = ranking[src(e)]
        r_dst = ranking[dst(e)]
        if (r_scr <= k) && (r_dst <= k)
            d = r_scr - r_dst
            # upset if score assignment contradicts the pairwise comparison
            if sign(a) * sign(d) > 0
                upsets += 1
            end
        end
    end
    return upsets
end

function syncrank(spL, meta_g; singular=true)
    # least eigenvector
    v, _ = power_method_least_eigenvalue(spL)# v = least_eigenvector(L; singular)
    # entry uv of L is -exp(i * theta(u,v))
    # with theta(u,v) approx h(u) - h(v) where h(u) is the score of node u
    score = angular_score(v)
    # since v(u) approx exp( i h(u) )
    p = ranking_from_score(score)
    ranking = best_shift(p, meta_g)
    return ranking
end

function ranking_from_score(score)
    # score is such that score[1] is the score of node 1
    # we find permutation putting ranking_score in descending order
    p = sortperm(vec(score); rev=true)
    # such that score[p[1]] > score[p[2]] > ...
    #
    # find inverse permutation containing score of each entry
    ranking = invperm(p)
    # such that ranking[1] is the ranking of node 1, etc ...
    return ranking
end

function best_shift(p, meta_g)
    n = nv(meta_g)
    # best circular shift to minimize upsets
    upsets = ne(meta_g)
    ranking = p
    for shift in 0:(n - 1)
        # shift ranking by 'shift'
        shifted_ranking = mod.(p .+ shift, n) .+ 1
        # compute # of upsets
        upsets_score = nb_upsets(meta_g, shifted_ranking)

        if upsets_score <= upsets
            # if # of upsets is lower, update score
            upsets = upsets_score
            ranking = shifted_ranking
        end
    end
    return ranking
end

# least eigenvector thanks to std julia eigensolver (not sparse)
function least_eigenvector(L; singular=true)
    # least eigenvector
    if singular
        # if Laplacian is singular
        F = eigen(L)
        V = F.vectors
        v = V[:, 1]
    else
        # otherwise simply compute the 6 smallest
        _, v = eigs(L; nev=6, which=:SM)
    end
    return v[:, 1]
end

# least eigenvector thanks to pow meth with sp matrices
function power_method_least_eigenvalue(sparseL; iters=20)
    sparseL = Hermitian(sparseL) # + 1e-10*I # regularization to use KrylovJL_CG
    n = size(sparseL)[1]
    est_eigenvec = rand(Complex{Float64}, n)
    est_eigenvec /= sqrt(est_eigenvec' * est_eigenvec)
    #Power iteration
    for _ in 1:iters
        prob = LinearProblem(sparseL, est_eigenvec)
        #est_eigenvec = solve(prob, KrylovJL_CG())# KrylovJL_CG (error when ill cond)
        est_eigenvec = solve(prob, KrylovJL_MINRES())
        est_eigenvec /= sqrt(est_eigenvec' * est_eigenvec)
    end
    est_eigenval = real(est_eigenvec' * sparseL * est_eigenvec)
    return est_eigenvec, est_eigenval
end

function power_method_top_eigenvalue(sparseL; iters=50)
    n = size(sparseL)[1]
    est_eigenvec = rand(Complex{Float64}, n)
    est_eigenvec /= sqrt(est_eigenvec' * est_eigenvec)
    #Power iteration
    for _ in 1:iters
        est_eigenvec = sparseL * est_eigenvec
        est_eigenvec /= sqrt(est_eigenvec' * est_eigenvec)
    end
    est_eigenval = real(est_eigenvec' * sparseL * est_eigenvec)
    return est_eigenvec, est_eigenval
end

function cond_nb_pp(spL)
    _, est_least_eigenval = power_method_least_eigenvalue(spL)
    _, est_top_eigenval = power_method_top_eigenvalue(spL)
    return est_top_eigenval / est_least_eigenval, est_least_eigenval, est_top_eigenval
end

function eigenvec_dist(u, v)
    # for normalized vectors
    scalar = v' * u
    dist = 1 - abs.(scalar[1])
    return dist
end

function normalize_Lap!(L)
    Deg = Diagonal(L)
    N = inv(sqrt(Deg))
    return L = N * L * N
end

function normalize_meta_g!(meta_g)
    for e in edges(meta_g)
        # compute # of neighbors
        deg_scr = length(neighbors(meta_g, src(e)))
        deg_dst = length(neighbors(meta_g, dst(e)))
        # standard weight normalization
        w = 1 / sqrt(deg_scr * deg_dst)
        set_prop!(meta_g, e, :e_weight, w)
    end
end

function cumulate_angles(mtsf)
    n = nv(mtsf)
    roots = get_prop(mtsf, :roots)
    vectors = []
    nb_nodes_in_trees = 0
    for cc in connected_components(mtsf)
        angles = zeros(Float64, n, 1)
        tree, nodes = induced_subgraph(mtsf, cc)
        for r in roots
            if r in cc
                s = 1
                # breadth first search
                parents = bfs_parents(tree, s)
                nb_nodes_in_trees += length(nodes)
                for i in 1:length(parents)
                    src_t = i
                    tgt_t = parents[i]
                    angle = 0
                    while src_t !== s
                        scr_node = nodes[src_t]
                        tgt_node = nodes[tgt_t]
                        theta = get_edge_prop(mtsf, Edge(scr_node, tgt_node), :angle)
                        angle += theta
                        src_t = tgt_t
                        tgt_t = parents[src_t]
                    end
                    angles[nodes[i]] = angle
                end
                v = zeros(ComplexF64, n, 1)
                v[nodes] = exp.(im * angles[nodes])
                v /= sqrt(length(nodes))
                push!(vectors, v)
            end
        end
    end
    percent_nodes_in_trees = nb_nodes_in_trees / n
    return vectors, percent_nodes_in_trees
end
