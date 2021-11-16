module RandomSpanningForests

using Graphs
using MetaGraphs
using Random

# Write your package code here.
using IterTools: partition

export multi_type_spanning_forest
include("utils.jl")
include("mtsf.jl")

end
