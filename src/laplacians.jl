
function sp_magnetic_incidence(graph; oriented::Bool=true)
    return magnetic_incidence_matrix(graph; oriented=oriented)
end

function magnetic_incidence(graph; oriented::Bool=true)::Matrix{Complex{Float64}}
    return Array(sp_magnetic_incidence(graph; oriented=oriented))
end

function magnetic_incidence_matrix(
    graph::AbstractGraph; oriented::Bool=true, phases::Bool=true
)::SparseMatrixCSC{Complex{Float64}}
    n_v, n_e = nv(graph), ne(graph)

    I = vcat(src.(edges(graph)), dst.(edges(graph)))
    J = vcat(1:n_e, 1:n_e)

    θ = get_edges_prop(graph, :angle, true, 0.0)
    w = @. exp(-im * 0.5 * θ)
    if !phases
        w = abs.(w)
    end
    V = vcat(w, oriented ? -conj.(w) : conj.(w))
    # here different sign wrt paper but unimportant
    B = transpose(sparse(I, J, V, n_v, n_e))
    return B
end

function mtsf_edge_indices(mtsf, graph)
    return [i for (i, e) in enumerate(edges(graph)) if has_edge(mtsf, src(e), dst(e))]
end

function average_sparsifier(
    rng::Random.AbstractRNG,
    meta_g::AbstractMetaGraph,
    ls::Union{Array,Nothing},
    q::Real,
    nb_samples::Integer;
    weighted::Bool=false,
    absorbing_node::Bool=false,
    ust::Bool=false,
    hkpv::Bool=false,
)
    m = ne(meta_g)
    edge_weights = get_edges_prop(meta_g, :e_weight, true, 1.0)

    # Initialization
    nb_cycles = zeros(nb_samples, 1)
    nb_roots = zeros(nb_samples, 1)
    subgraph_weights = zeros(nb_samples, 1)
    w_tot = 0

    # for storing weights of edges in the sparsifier
    sp_e_weight_diag_el = spzeros(m)

    sparseB = sp_magnetic_incidence(meta_g; oriented=true)

    if hkpv == true
        if q < 1e-10
            # QR
            V = Matrix(qr(Matrix(sparseB)).Q)
            # full = false
            # V = Matrix(svd(Matrix(sparseB);full).U)
        else
            # dilating to a DPP on the edges + nodes
            # so that it becomes projective
            BigB = [sparseB; sqrt(q) * sparse(I, n, n)]
            V = Matrix(qr(Matrix(BigB)).Q)
        end
    end

    for i_sample in 1:nb_samples
        #
        if hkpv == true
            # standard algorithm by Hough et al (HKPV)
            if q < 1e-10
                id = sample_pdpp(V)
            else
                ind_e = sample_pdpp(V)
                # restricting the extension to the edges only
                ind_e = ind_e[ind_e .< (m + 1)]
            end
            ind_e = vec(id)
            w = 1
            w_tot = 1
        else
            # cycle popping for weakly inconsistent cycles
            mtsf = multi_type_spanning_forest(rng, meta_g, q; weighted, absorbing_node, ust)
            # check nb roots and cycles
            cycles = get_prop(mtsf, :cycle_nodes)
            nb_cycles[i_sample] = length(cycles)
            roots = get_prop(mtsf, :roots)
            nb_roots[i_sample] = length(roots)

            D = props(mtsf)
            w = D[:weight]
            subgraph_weights[i_sample] = w
            w_tot += w
            ind_e = mtsf_edge_indices(mtsf, meta_g)
        end
        #

        diag_elements = ones(length(ind_e))

        if ls === nothing
            diag_elements /= (length(ind_e) / m)
        else
            diag_elements ./= ls[ind_e]
        end

        if weighted
            diag_elements = diag_elements .* edge_weights[ind_e]
        end

        sp_e_weight_diag_el[ind_e] += w * diag_elements
    end
    nb_sampled_cycles = sum(nb_cycles)
    nb_sampled_roots = sum(nb_roots)

    # checking connectivity

    ind_tot = vec(1:m) # vec with all edge indices
    ind_edges_sparsifier = ind_tot[sp_e_weight_diag_el .!= 0] # vec with edge indices in sparsifier
    subgraph = MetaGraph(nv(meta_g))
    all_edges = collect(edges(meta_g))
    subset_edges = all_edges[ind_edges_sparsifier]

    for e in subset_edges
        add_edge!(subgraph, e)
    end

    isconnected = is_connected(subgraph)

    L = (1 / w_tot) * sparseB' * spdiagm(sp_e_weight_diag_el) * sparseB

    return L, nb_sampled_cycles, nb_sampled_roots, subgraph_weights, isconnected
end

function average_sparsifier_iid(
    rng::Random.AbstractRNG,
    meta_g::AbstractMetaGraph,
    ls::Union{Array,Nothing},
    batch::Integer,
    nb_samples::Integer;
    weighted::Bool=false,
)
    m = ne(meta_g)
    e_weights = get_edges_prop(meta_g, :e_weight, true, 1.0)

    if ls !== nothing
        p = vec(ls / sum(ls))
    end

    w = 1
    w_tot = 0
    # for storing weights of edges in the sparsifier
    sp_e_weight_diag_el = spzeros(m)

    for _ in 1:nb_samples
        w_tot += w
        ind_e = zeros(batch)

        if ls === nothing
            ind_e = rand(rng, 1:m, batch)
        else
            ind_e = rand(rng, Categorical(p), batch)
        end

        diag_elements = ones(length(ind_e))

        if ls === nothing
            diag_elements /= (length(ind_e) / m)
        else
            diag_elements ./= ls[ind_e]
        end

        if weighted
            diag_elements = diag_elements .* e_weights[ind_e]
        end

        sp_e_weight_diag_el[ind_e] += w * diag_elements
    end

    # checking connectivity

    ind_tot = vec(1:m) # vec with all edge indices
    ind_edges_sparsifier = ind_tot[sp_e_weight_diag_el .!= 0] # vec with edge indices in sparsifier
    subgraph = MetaGraph(nv(meta_g))
    all_edges = collect(edges(meta_g))
    subset_edges = all_edges[ind_edges_sparsifier]

    for e in subset_edges
        add_edge!(subgraph, e)
    end

    isconnected = is_connected(subgraph)

    sparseB = sp_magnetic_incidence(meta_g; oriented=true)
    L = (1 / w_tot) * sparseB' * spdiagm(sp_e_weight_diag_el) * sparseB

    return L, isconnected
end

function leverage_score(
    spB::SparseMatrixCSC{ComplexF64,Int64}, q::Real; e_weights=ones(size(spB)[1])
)
    W = spdiagm(e_weights)
    spL_reg = spB' * W * spB + q * I

    # nb: linear systems with sparse rhs not yet available in julia
    lev_scores = real(diag(W * spB * (spL_reg \ Matrix(spB'))))

    return lev_scores
end

function JL_lev_score_estimates(spB, q; e_weights=ones(size(spB)[1]), cst=40)
    # Johnson-Lindenstrauss estimate of leverage scores
    # with Rademacher sketching
    n = size(spB)[2]
    m = size(spB)[1]
    if q < 1e-12
        # by courtesy of an anonymous reviewer of ACHA
        k = Int(ceil(cst * log(m) + 1)) # number of samples
        Q = (2 * bitrand((m, k)) .- 1) # sketching with Rademacher random variables
        spL = spB' * diagm(e_weights) * spB + q * sparse(I, n, n)
        M = sqrt.(e_weights) .* (spB * (spL \ (spB' * (sqrt.(e_weights) .* Q))))
        lev_score_estimates = (abs.(M) .^ 2) * ones(k)
        lev_score_estimates /= k # normalization of sketching matrix
    else
        k = Int(ceil(cst * log(m + n) + 1)) # number of samples
        Q = (2 * bitrand((m + n, k)) .- 1) / sqrt(k) # sketching with Rademacher random variables
        wB = diagm(sqrt.(e_weights)) * spB
        BigB = [sqrt(q) * sparse(I, n, n); wB]
        reg_spL = BigB' * BigB
        M = wB * (reg_spL \ (BigB' * Q))
        lev_score_estimates = (abs.(M) .^ 2) * ones(k)
    end
    println("k = ", k, " vs nb edges= ", m, "\n ")

    return lev_score_estimates
end

function emp_leverage_score(
    rng::Random.AbstractRNG,
    meta_g::AbstractMetaGraph,
    q::Real,
    t::Integer;
    weighted::Bool=false,
    absorbing_node::Bool=false,
    ust::Bool=false,
    hkpv::Bool=false,
)
    m = ne(meta_g)
    n = nv(meta_g)
    emp_lev = zeros(m, 1)

    if hkpv
        sparseB = sp_magnetic_incidence(meta_g; oriented=true)
        if q < 1e-10
            V = Matrix(qr(Matrix(sparseB)).Q)
        else
            # dilating to a DPP on the edges + nodes
            # so that it becomes projective
            BigB = [sparseB; sqrt(q) * sparse(I, n, n)]
            V = Matrix(qr(Matrix(BigB)).Q)
        end
        println("HKPV")
    end
    for _ in 1:t
        if hkpv
            if q < 1e-10
                ind_e = sample_pdpp(V)
            else
                ind_e = sample_pdpp(V)
                # restricting the dilation to the edges only
                ind_e = ind_e[ind_e .< (m + 1)]
            end
        else
            mtsf = multi_type_spanning_forest(rng, meta_g, q; weighted, absorbing_node, ust)
            ind_e = mtsf_edge_indices(mtsf, meta_g)
        end
        emp_lev[ind_e] = emp_lev[ind_e] .+ 1
    end
    emp_lev /= t

    return emp_lev
end

function edge_weights(g)
    e_weights = get_edges_prop(g, :e_weight, true, 1.0)
    return e_weights
end

nb_of_edges(L::AbstractMatrix) = (nnz(sparse(L)) - size(L, 1)) / 2

function optimal_perm(mtsf::AbstractMetaGraph)
    n_v = nv(mtsf)
    #get the roots
    roots = get_prop(mtsf, :roots)
    # get the branches in the (reverse) order there were sampled
    branches = get_prop(mtsf, :branches)
    flt_branches = collect(Iterators.flatten(branches))

    # indices array putting the nodes in the right order
    non_cycle_nodes = [roots; flt_branches]
    ind_perm = [non_cycle_nodes; setdiff(1:n_v, non_cycle_nodes)]
    return ind_perm
end

function arpack_rel_cond(L, spL)
    # warning ! "Arpack", version="0.5.3" NOT "0.5.4" (broken)
    λ_min, _ = eigs(L, spL; nev=2, which=:SM, tol=1e-1, maxiter=600)
    λ_max, _ = eigs(L, spL; nev=2, which=:LM, tol=1e-1, maxiter=600)

    cnd = real(λ_max[1]) / real(λ_min[1])

    return cnd[1]
end


function arpack_cond(L)
    # warning ! "Arpack", version="0.5.3" NOT "0.5.4" (broken)
    λ_min, _ = eigs(L; nev=1, which=:SM)
    λ_max, _ = eigs(L; nev=1, which=:LM)

    cnd = real(λ_max[1]) / real(λ_min[1])

    return cnd
end

function sp_pcond_Lap(spL::SparseMatrixCSC{ComplexF64,Int64}, q, L)
    temp = spL + q * I
    temp = 0.5 * (temp + temp')
    C = cholesky(temp)
    R = sparse(C.L)[invperm(C.p), :] # since sparse cholesky is pivoted

    T = linear_solve_matrix_system(R, L + q * I) # sparse matrix AX=B
    pL = linear_solve_matrix_system(R, T') # sparse matrix AX=B

    return pL, R
end

function pcond_Lap(spL, q, L)
    C = cholesky(spL + q * I)
    R = sparse(C.L)
    R = R[invperm(C.p), :]
    #R = sparse(C.L)[invperm(C.p), :]
    T = Matrix(R) \ Matrix(L + q * I)
    pL = Matrix(R) \ (T')

    return pL, R
end

function linear_solve_matrix_system(A, B)
    if size(A)[1] != size(A)[2]
        error("A is not square in AX = B")
    end
    if size(A)[1] != size(B)[1]
        error("nb of rows of A and B should be equal in AX = B")
    end
    if size(B)[2] == 1
        error("RHS is a col vector. Please use LinearSolve.jl directly")
    end

    n = size(A)[1]
    k = size(B)[2]
    X = spzeros(Complex{Float64}, n, k)

    for i in 1:k
        prob = LinearProblem(A, B[:, i])
        x = solve(prob, KrylovJL_GMRES())
        X[:, i] = x
    end
    return X
end
