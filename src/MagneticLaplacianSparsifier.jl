module MagneticLaplacianSparsifier

using Graphs
using MetaGraphs
using Random
using SparseArrays

using IterTools: partition

include("utils.jl")
include("mtsf.jl")
include("graphs.jl")
include("laplacians.jl")
include("syncrank.jl")

export multi_type_spanning_forest,
    add_edges_from!,
    consecutive_pairs,
    get_edge_property,
    set_prop!,
    generateGraphMUN,
    generateGraphERO,
    magneticIncidence,
    mtsf_edge_indices,
    angular_score

end
