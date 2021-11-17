module RandomSpanningForests

using Graphs
using MetaGraphs
using Random

# Write your package code here.
using IterTools: partition

export multi_type_spanning_forest,
    add_edges_from!, consecutive_pairs, set_prop!, get_edge_property
include("utils.jl")
include("mtsf.jl")

end
