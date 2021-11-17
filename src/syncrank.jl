function angular_score(v)
    n = length(v)
    score = zeros(n, 1)
    for i in 1:n
        score[i] = angle(v[i])
    end
    return score
end
