using Documenter, alphaPredictableComponents

makedocs(
    modules = [alphaPredictableComponents],
    clean = false,
    sitename = "alphaPredictable\nComponents.jl",
    pages = [
        "Home" => "index.md",
        "Guide" => "man/Example1.md",
        "Library" => "lib/library.md"
    ],
)


deploydocs(
    repo = "github.com/Mattriks/alphaPredictableComponents.jl.git",
)

