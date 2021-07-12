using BlackScholesGreeks
using Documenter

DocMeta.setdocmeta!(BlackScholesGreeks, :DocTestSetup, :(using BlackScholesGreeks); recursive=true)

makedocs(;
    modules=[BlackScholesGreeks],
    authors="Yefeng Wang",
    repo="https://github.com/Stochastic-Adventure/BlackScholesGreeks.jl/blob/{commit}{path}#{line}",
    sitename="BlackScholesGreeks.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Stochastic-Adventure.github.io/BlackScholesGreeks.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Stochastic-Adventure/BlackScholesGreeks.jl",
)
