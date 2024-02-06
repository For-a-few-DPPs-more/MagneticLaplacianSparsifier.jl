# borrowed from Determinantal.jl
#
# https://github.com/dahtah/Determinantal.jl
# written by Simon Barthemé, Nicolas Tremblay, Guillaume Gautier
# Implementation of HKPV algorithm

function lvg(U)
    p = abs.(U[:, 1]) .^ 2
    if (size(U, 2) > 1)
        for j in 2:size(U, 2)
            for i in 1:size(U, 1)
                p[i] += abs(U[i, j])^2
            end
        end
    end
    return p
end

function wsample(w::AbstractVector, ttl)
    u = rand() * ttl
    i = 0
    s = 0
    while s < u
        i += 1
        s += w[i]
    end
    return i
end

function sample_pdpp(U::AbstractMatrix)
    return sample_pdpp(U, lvg(U))
end

function sample_pdpp(U::AbstractMatrix, lvg::AbstractVector)
    n = size(U, 1)
    m = size(U, 2)
    #Initial distribution
    #Um = Matrix(U)
    p = copy(lvg)
    F = zeros(ComplexF64, m, m)
    f = zeros(ComplexF64, m)
    v = zeros(ComplexF64, m)
    tmp = zeros(ComplexF64, n)
    inds = BitSet()
    ss = sum(lvg)
    @inbounds for i in 1:m
        itm = wsample(p, ss)
        push!(inds, itm)
        #v = vec(Matrix(U[itm,:]))
        #v = @view U[itm,:]
        #v = U[itm,:]
        copyto!(v, U[itm, :])
        if i == 1
            copyto!(f, v)
        else
            Fv = @view F[:, 1:(i - 1)]
            copyto!(f, v - Fv * (Fv' * v))
        end
        F[:, i] = f / sqrt(dot(v, f))
        mul!(tmp, U, @view F[:, i])
        ss = 0.0
        @inbounds for j in 1:n
            s = p[j] - abs(tmp[j])^2
            p[j] = (s > 0 ? s : 0)
            ss += p[j]
        end
        @inbounds for j in inds
            ss -= p[j]
            p[j] = 0
        end
    end
    return collect(inds)
end
