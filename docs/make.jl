using Documenter, ImagePhaseCongruency
using TestImages
using DemoCards

testimage("cameraman") # used to trigger artifact downloading

# generate
demopage, postprocess_cb, demo_assets = makedemos("examples") # this is the relative path to docs/
assets = []
isnothing(demo_assets) || (push!(assets, demo_assets))

format = Documenter.HTML(edit_link = "master",
                         prettyurls = get(ENV, "CI", nothing) == "true",
                         assets = assets)
makedocs(
    format=format,
    sitename = "ImagePhaseCongruency",
    pages = [
        "index.md",
        demopage,
        "functions.md"
    ]
)

postprocess_cb()

deploydocs(
    repo = "github.com/peterkovesi/ImagePhaseCongruency.jl.git",
)
