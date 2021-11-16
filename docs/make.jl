using MagneticLaplacianSparsifier
using Documenter

DocMeta.setdocmeta!(MagneticLaplacianSparsifier, :DocTestSetup, :(using MagneticLaplacianSparsifier); recursive=true)

makedocs(;
    modules=[MagneticLaplacianSparsifier],
    authors="Michaël Fanuel <mrfanuel@hotmail.fr>, Guillaume Gautier <guillaume.gga@gmail.com>, and contributors",
    repo="https://github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl/blob/{commit}{path}#{line}",
    sitename="MagneticLaplacianSparsifier.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://For-a-few-DPPs-more.github.io/MagneticLaplacianSparsifier.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/For-a-few-DPPs-more/MagneticLaplacianSparsifier.jl",
    devbranch="main",
)
