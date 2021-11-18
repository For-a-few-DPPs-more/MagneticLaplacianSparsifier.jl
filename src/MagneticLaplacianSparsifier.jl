module MagneticLaplacianSparsifier

using Distributions
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
    generateGraphMUN,
    generateGraphERO,
    magneticIncidence,
    mtsf_edge_indices,
    angular_score,
    magnetic_incidence_matrix,
    averageSparsifier,
    leverageScore,
    empLeverageScore,
    modulus_entries,
    nb_of_edges

end
