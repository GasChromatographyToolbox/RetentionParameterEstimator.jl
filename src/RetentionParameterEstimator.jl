module RetentionDataEstimator

#using Reexport
using GasChromatographySimulator
#using ForwardDiff
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using OptimizationOptimisers
using OptimizationCMAEvolutionStrategy
using Dierckx
#using BenchmarkTools

include("Load.jl")
include("Loss.jl")
include("Estimate_Start_Values.jl")
include("Optimization.jl")
include("Simulate_Test.jl")

const θref = 30.0
const rT_nom = 0.6
const Tst = 273.15

end # module