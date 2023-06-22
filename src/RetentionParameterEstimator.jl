module RetentionParameterEstimator

#using Reexport
using GasChromatographySimulator
#using ForwardDiff
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using OptimizationOptimisers
using OptimizationCMAEvolutionStrategy
using Interpolations
using DataFrames
using CSV
using Statistics
using Measurements
using ForwardDiff
#using BenchmarkTools

include("Load.jl")
include("Loss.jl")
include("Estimate_Start_Values.jl")
include("Optimization.jl")
include("Simulate_Test.jl")
include("Misc.jl")

const θref = 30.0
const rT_nom = 0.69
const Tst = 273.15
const R = 8.31446261815324
const std_opt = GasChromatographySimulator.Options(abstol=1e-8, reltol=1e-5, ng=true, odesys=false)

end # module