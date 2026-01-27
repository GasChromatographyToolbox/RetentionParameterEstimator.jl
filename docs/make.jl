using Documenter
using RetentionParameterEstimator

makedocs(
    sitename = "RetentionParameterEstimator",
    modules = [RetentionParameterEstimator],
    pages = Any[
                "Home" => "index.md",
                "Docstrings" => "docstrings.md"
            ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/GasChromatographyToolbox/RetentionParameterEstimator.jl",
    devbranch = "main"
)
