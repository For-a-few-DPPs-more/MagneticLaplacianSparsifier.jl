angular_score(v) = angle.(v)
modulus_entries(v) = abs.(v)

function nb_upsets(meta_g, ranking_score)
    oriented = true
    upsets = 0
    for e in edges(meta_g)
        a = get_edge_prop(meta_g, e, :angle, oriented)
        d = ranking_score[src(e)] - ranking_score[dst(e)]

        if sign(a) * sign(d) < 0
            upsets += 1
        end
    end
    return upsets
end

function syncrank(L, meta_g, singular=false)
    n = nv(meta_g)

    # permutation
    id_p = randperm(n)
    inv_id_p = invperm(id_p)

    # least eigenvector
    v = least_eigenvector(L[id_p, id_p]; singular)

    ranking_score = angular_score(v)

    # inverse perturbation
    ranking_score = ranking_score[inv_id_p]
    p = sortperm(vec(ranking_score); rev=true)

    # best cyclic shift to minimize upsets
    upsets = ne(meta_g)
    score = zeros(n, 1)
    for shift in 1:n
        shifted_order = circshift(p, shift)
        upsets_score = nb_upsets(meta_g, shifted_order)
        if upsets_score <= upsets
            upsets = upsets_score
            score = shifted_order
        end
    end

    return score
end

function least_eigenvector(L; singular=false)
    # least eigenvector
    if singular
        F = eigen(L)
        V = F.vectors
        v = V[:, 1]
    else
        _, v = eigs(L; nev=1, which=:SM)
    end
    return v
end

function eigenvec_dist(u, v)
    scalar = v' * u
    dist = 1 - abs.(scalar[1])
    return dist
end
