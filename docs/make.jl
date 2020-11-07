#push!(LOAD_PATH, "../src/")

using Documenter, ImagePhaseCongruency

makedocs(
    sitename = "ImagePhaseCongruency",
    pages = [
        "index.md",
        "examples.md",
        "functions.md"
    ]
)

deploydocs(
    repo = "github.com/peterkovesi/ImagePhaseCongruency.jl.git",
)
