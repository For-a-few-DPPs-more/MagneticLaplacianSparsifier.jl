function angular_score(v)
    n = length(v)
    score = zeros(n, 1)
    for i in 1:n
        score[i] = angle(v[i])
    end
    return score
end

function modulus_entries(v)
    n = length(v)
    mod_ent = zeros(n, 1)
    for i in 1:n
        mod_ent[i] = abs(v[i])
    end
    return mod_ent
end
