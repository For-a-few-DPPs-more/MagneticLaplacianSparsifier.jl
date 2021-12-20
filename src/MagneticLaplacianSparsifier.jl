module MagneticLaplacianSparsifier

using Distributions
using Graphs
using MetaGraphs
using Random
using SparseArrays
using LinearAlgebra
using Arpack

using IterTools: partition

include("utils.jl")
include("mtsf.jl")
include("graphs.jl")
include("laplacians.jl")
include("sync-rank.jl")

# todo cleamtidy up exports
export multi_type_spanning_forest,
    add_edges_from!,
    consecutive_pairs,
    get_edge_prop,
    erdos_renyi,
    gen_graph_mun,
    gen_graph_ero,
    magnetic_incidence,
    mtsf_edge_indices,
    angular_score,
    magnetic_incidence_matrix,
    sp_magnetic_incidence,
    average_sparsifier,
    leverage_score,
    emp_leverage_score,
    modulus_entries,
    nb_of_edges,
    optimal_perm,
    nb_upsets,
    syncrank,
    cond_numbers
end
