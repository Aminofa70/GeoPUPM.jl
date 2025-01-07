using GeoPUPM
using Documenter

DocMeta.setdocmeta!(GeoPUPM, :DocTestSetup, :(using GeoPUPM); recursive=true)

makedocs(;
    modules=[GeoPUPM],
    authors="Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    sitename="GeoPUPM.jl",
    format=Documenter.HTML(;
        canonical="https://Aminofa70.github.io/GeoPUPM.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Aminofa70/GeoPUPM.jl",
    devbranch="main",
)
