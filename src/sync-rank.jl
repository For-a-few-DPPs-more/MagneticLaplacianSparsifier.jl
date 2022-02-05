angular_score(v) = angle.(v)
modulus_entries(v) = abs.(v)

function nb_upsets(meta_g, ranking_score)
    oriented = true
    upsets = 0
    for e in edges(meta_g)
        a = get_edge_prop(meta_g, e, :angle, oriented)
        d = ranking_score[src(e)] - ranking_score[dst(e)]
        # upset if score assignment contradicts the pairwise comparison
        if sign(a) * sign(d) < 0
            upsets += 1
        end
    end
    return upsets
end

function syncrank(L, meta_g; singular=true)
    n = nv(meta_g)

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
        # shift p by 'shift'
        shifted_order = circshift(p, shift)
        # compute # of upsets
        upsets_score = nb_upsets(meta_g, shifted_order)
        if upsets_score <= upsets
            # if # of upsets is lower, update score
            upsets = upsets_score
            score = shifted_order
        end
    end

    return score
end

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
