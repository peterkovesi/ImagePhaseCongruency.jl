#push!(LOAD_PATH,"/Users/pk/JuliaProjects/ImagePhaseCongruency/src/")

using Documenter, ImagePhaseCongruency

makedocs(
    format = :html,
    sitename = "ImagePhaseCongruency",
    pages = [
        "index.md",
        "examples.md",
        "functions.md"
    ]
)
