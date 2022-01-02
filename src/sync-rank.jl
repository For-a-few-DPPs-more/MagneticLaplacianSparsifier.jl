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
    if singular
        F = eigen(L)
        V = F.vectors
        v = V[:, 1]
    else
        _, v = eigs(L; nev=1, which=:SM)
    end

    ranking_score = angular_score(v)
    p = sortperm(vec(ranking_score); rev=true)

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
