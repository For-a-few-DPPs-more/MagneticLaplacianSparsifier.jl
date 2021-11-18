module MagneticLaplacianSparsifier

using Graphs
using MetaGraphs
using Random
using SparseArrays
using LinearAlgebra

using IterTools: partition

include("utils.jl")
include("mtsf.jl")
include("graphs.jl")
include("laplacians.jl")
include("sync-rank.jl")

export multi_type_spanning_forest,
    add_edges_from!,
    consecutive_pairs,
    get_edge_prop,
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
    nb_of_edges

end
