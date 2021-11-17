module MagneticLaplacianSparsifier

include("RandomSpanningForests.jl")
using .RandomSpanningForests # note the dot

using Graphs
using MetaGraphs
using Random
using SparseArrays
using Arpack

include("graphs.jl")
include("laplacians.jl")

export multi_type_spanning_forest, add_edges_from!, consecutive_pairs,get_edge_property, set_prop!, generateGraphMUN, generateGraphERO,magneticIncidence

end