module RetentionParameterEstimator

using Reexport

@reexport using CSV
@reexport using DataFrames
using ForwardDiff
using GasChromatographySimulator
@reexport using Measurements
using Optimization
using OptimizationBBO
using OptimizationCMAEvolutionStrategy
using OptimizationOptimJL
using OptimizationOptimisers
@reexport using Plots
@reexport using PlutoUI
@reexport using StatsPlots
using UrlDownload
@reexport using Statistics

include("Load.jl")
include("Loss.jl")
include("Estimate_Start_Values.jl")
include("Optimization.jl")
include("Simulate_Test.jl")
include("Misc.jl")
include("Notebook.jl")


const θref = 30.0
const rT_nom = 0.69
const Tst = 273.15
const R = 8.31446261815324
const std_opt = GasChromatographySimulator.Options(alg=Tsit5(),abstol=1e-8, reltol=1e-5, ng=true, odesys=false)

end # module