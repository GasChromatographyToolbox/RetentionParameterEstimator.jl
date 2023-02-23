using Documenter
using RetentionParameterEstimator

makedocs(
    sitename = "RetentionParameterEstimator",
    format = Documenter.HTML(),
    modules = [RetentionParameterEstimator]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
