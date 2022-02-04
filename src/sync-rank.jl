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

function syncrank(L, meta_g; singular=true)
    n = nv(meta_g)

    ## permutation
    #id_p = randperm(n)
    #inv_id_p = invperm(id_p)
    #v = least_eigenvector(L[id_p, id_p]; singular)
    #ranking_score = angular_score(v)
    #ranking_score = ranking_score[inv_id_p]    # inverse perturbation
    ##

    # least eigenvector
    v = least_eigenvector(L; singular)
    ranking_score = angular_score(v)

    # find permutation putting ranking_score in descending order
    p = sortperm(vec(ranking_score); rev=true)
    # find inverse permutation containing score of each entry
    p = invperm(p)

    # best circular shift to minimize upsets
    upsets = ne(meta_g)
    score = zeros(n, 1)
    for shift in 0:n
        shifted_order = circshift(p, shift)
        upsets_score = nb_upsets(meta_g, shifted_order)
        if upsets_score <= upsets
            upsets = upsets_score
            score = shifted_order
        end
    end

    return score
end

function least_eigenvector(L; singular=true)
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

function normalize_Lap!(L)
    Deg = Diagonal(L)
    N = inv(sqrt(Deg))
    return L = N * L * N
end

function normalize_meta_g!(meta_g)
    for e in edges(meta_g)
        deg_scr = length(neighbors(meta_g, src(e)))
        deg_dst = length(neighbors(meta_g, dst(e)))
        w = 1 / sqrt(deg_scr * deg_dst)
        set_prop!(meta_g, e, :e_weight, w)
    end
end
