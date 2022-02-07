using Test, MagneticLaplacianSparsifier

using Graphs, MetaGraphs
using Random
using LinearAlgebra
using SparseArrays
using StatsBase

using MagneticLaplacianSparsifier: getRNG

const testdir = dirname(@__FILE__)

const name_test_file = [
    "random_spanning_forests",
    "magnetic_incidence_matrix",
    "leverage_scores",
    "sparsifiers",
    "edge_statistics",
    "syncrank",
]

@testset verbose = true "MagneticLaplacianSparsifier.jl" begin
    for name in name_test_file
        test_path = joinpath(testdir, "$(name).jl")
        include(test_path)
    end
end
