# MagneticLaplacianSparsifier

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://For-a-few-DPPs-more.github.io/MagneticLaplacianSparsifier.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://For-a-few-DPPs-more.github.io/MagneticLaplacianSparsifier.jl/dev)
[![Build Status](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl)

## Documentation
This is the code associated with the manuscript 
[Sparsification of the regularized magnetic Laplacian  thanks to multi-type spanning forests](https://arxiv.org/pdf/2208.14797.pdf)
by [Michaël Fanuel](https://mrfanuel.github.io/) and [Rémi Bardenet](https://rbardenet.github.io/).

This code uses multiple samples of Multi-Type Spanning Forests (MTSF) to sparsify the regularized magnetic Laplacian.

Here is an illustration of MTSF sampling with Wilson's algorithm

![](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/main/notebooks/figures/MTSF_grid.gif)

The output is an MTSF with one cycle-rooted tree and four rooted tree (roots are in red).

![](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/main/notebooks/figures/MTSF_137.png)


## Use of the code  
### Install Julia

If you do not have Julia installed, please visit [julialang.org](https://julialang.org/learning/getting-started/)
### Installation

[`MagneticLaplacianSparsifier.jl`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl) is not registered.
The way to use it is to type

```julia
julia> ]add https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl
```

### Jupyter notebooks for reproducing the paper figures

You can execute the Jupyter [`notebooks`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/master/notebooks) to generate the paper figures.

- [`plotting random spanning graphs`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/master/notebooks/plots.ipynb) See Figure 1.
- [`sparsify-and-eigensolve`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/master/notebooks/syncrank.ipynb) See Figure 3 and 12.
- [`sparsify-and-precondition`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/master/notebooks/preconditioning.ipynb) See Figure 4, 5, 6 and 7.
- [`timings`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/master/notebooks/timings.ipynb) See Figure 8 and 9.
- [`leverage score estimation`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/master/notebooks/demo_lev_scores.ipynb) See Figure 11.



### Usage

To use this functions of this package, simply type

```julia
julia> using MagneticLaplacianSparsifier
```
