module MagneticLaplacianSparsifier

include("RandomSpanningForests.jl")
using .RandomSpanningForests # note the dot

include("graphs.jl")

export multi_type_spanning_forest, add_edges_from!, consecutive_pairs, set_prop!, generateGraphMUN, generateGraphERO

end