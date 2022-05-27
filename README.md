# MagneticLaplacianSparsifier

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://For-a-few-DPPs-more.github.io/MagneticLaplacianSparsifier.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://For-a-few-DPPs-more.github.io/MagneticLaplacianSparsifier.jl/dev)
[![Build Status](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl)

### Documentation
This is the code associated with the manuscript 
[Sparsification of the regularized magnetic Laplacian  thanks to multi-type spanning forests]()
by Michaël Fanuel and Rémi Bardenet
### Install Julia

If you do not have Julia installed, please visit [julialang.org](https://julialang.org/learning/getting-started/)
### Installation

[`MagneticLaplacianSparsifier.jl`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl) is not registered.
The way to use it is to type

```julia
julia> ]add https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl
```

### Jupyter notebooks for reproducing the paper figures

You can execute the Jupyter [`notebooks`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/master/notebooks) to play with the code.

- [`sparsify-and-eigensolve`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/master/notebooks/syncrank.ipynb)
- [`sparsify-and-precondition`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/master/notebooks/preconditioning.ipynb)
- [`timings`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/master/notebooks/timings.ipynb)
- [`leverage score estimation`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/master/notebooks/demo_lev_scores.ipynb)
- [`plotting random spanning graphs`](https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/master/notebooks/plots.ipynb)


### Usage

To use this functions of this package, simply type

```julia
julia> using MagneticLaplacianSparsifier
```
