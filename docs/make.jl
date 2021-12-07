using MotifPvalue
using Documenter

DocMeta.setdocmeta!(MotifPvalue, :DocTestSetup, :(using MotifPvalue); recursive=true)

makedocs(;
    modules=[MotifPvalue],
    authors="Shane Kuei Hsien Chu (skchu@wustl.edu)",
    repo="https://github.com/kchu25/MotifPvalue.jl/blob/{commit}{path}#{line}",
    sitename="MotifPvalue.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kchu25.github.io/MotifPvalue.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/MotifPvalue.jl",
    devbranch="main",
)
