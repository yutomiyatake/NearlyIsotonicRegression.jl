using NearlyIsotonicRegression
using Documenter

DocMeta.setdocmeta!(NearlyIsotonicRegression, :DocTestSetup, :(using NearlyIsotonicRegression); recursive=true)

makedocs(;
    modules=[NearlyIsotonicRegression],
    authors="yutomiyatake <yuto.miyatake.cmc@osaka-u.ac.jp> and contributors",
    repo="https://github.com/yutomiyatake/NearlyIsotonicRegression.jl/blob/{commit}{path}#{line}",
    sitename="NearlyIsotonicRegression.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yutomiyatake.github.io/NearlyIsotonicRegression.jl",
        repolink="https://github.com/yutomiyatake/NearlyIsotonicRegression.jl",
        edit_link="main",
        assets=String[],
    ),
    # pages=[
    #     "Home" => "index.md",
    # ],
    pages = ["Home" => "index.md", "Functions" => Any["func/iso.md", "func/neariso.md"],  "examples.md"]
)

deploydocs(;
    repo="github.com/yutomiyatake/NearlyIsotonicRegression.jl",
    devbranch="main",
)
