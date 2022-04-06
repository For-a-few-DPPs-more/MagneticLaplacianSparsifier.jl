module MagneticLaplacianSparsifier

using Distributions
using Graphs
using MetaGraphs
using Random
using SparseArrays
using LinearAlgebra
using Arpack
using Statistics
using Measures
using StatsBase
using Plots

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
    cond_numbers,
    average_sparsifier_iid,
    pcond_Lap,
    rand_step,
    step_to_root,
    get_edges_prop,
    edge_weights,
    least_eigenvector,
    eigenvec_dist,
    normalize_Lap!,
    normalize_meta_g!,
    benchmark_syncrank,
    cumulate_angles,
    best_shift,
    ranking_from_score,
    timings_cond_numbers,
    number_of_upsets_in_top,
    gen_graph_mun_basic,
    gen_graph_ero_basic,
    ero_mun,
    ero_mun_sbm,
    plot_comparison_sync,
    plot_comparison_cond,
    plot_nb_cycles,
    plot_nb_roots,
    gen_graph_cliques,
    gen_graph_planted_triangles
end
