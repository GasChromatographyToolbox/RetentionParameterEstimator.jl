#!/usr/bin/env julia
"""
Quick benchmark script for method_m1 only.

Usage:
    julia benchmark_m1.jl [data_file] [--parallel|--no-parallel] [maxiters] [maxtime]

Examples:
    julia benchmark_m1.jl
    julia -t 4 benchmark_m1.jl --parallel
    julia benchmark_m1.jl --no-parallel 5000 300.0
"""

# Activate the package environment
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using RetentionParameterEstimator
using Dates
using Printf

# Default values
default_data_file = joinpath(@__DIR__, "..", "data", "meas_df05_Rxi5SilMS.csv")
default_maxiters = 10000
default_maxtime = 600.0

# Parse command-line arguments
data_file = default_data_file
use_parallel = nothing
maxiters = default_maxiters
maxtime = default_maxtime

for arg in ARGS
    if arg == "--parallel"
        use_parallel = true
    elseif arg == "--no-parallel"
        use_parallel = false
    elseif endswith(arg, ".csv")
        data_file = arg
    elseif maxiters == default_maxiters
        try
            maxiters = parse(Int, arg)
        catch
            println("Warning: Could not parse maxiters: $arg, using default: $maxiters")
        end
    else
        try
            maxtime = parse(Float64, arg)
        catch
            println("Warning: Could not parse maxtime: $arg, using default: $maxtime")
        end
    end
end

# Auto-detect parallelization if not specified
if use_parallel === nothing
    use_parallel = Threads.nthreads() > 1
end

println("=" ^ 80)
println("Benchmark: method_m1 only")
println("=" ^ 80)
println("\nConfiguration:")
println("  Data file: $data_file")
println("  Parallelization: $(use_parallel ? "ENABLED ($(Threads.nthreads()) threads)" : "DISABLED")")
println("  Optimization limits: maxiters=$maxiters, maxtime=$maxtime")
println()

# Load data
println("Loading data...")
meas = RetentionParameterEstimator.load_chromatograms(data_file)
col_input = (L=meas[1].L, d=meas[1].d*1e3)  # Convert d from m to mm
println("  Loaded $(length(meas[4])) substances with $(length(meas[3].measurement)) measurements")
println()

# Warm-up run
println("Warming up (compiling code)...")
try
    RetentionParameterEstimator.method_m1(meas, col_input, se_col=false, parallel=use_parallel, maxiters=10, maxtime=1.0)
    println("  Warm-up complete")
catch e
    println("  Warning: Warm-up failed: $e")
end
println()

# Run benchmark
println("Running method_m1...")
println("-" ^ 80)
time_m1 = @elapsed begin
    res_m1, Telu_max_m1 = RetentionParameterEstimator.method_m1(meas, col_input, se_col=true, parallel=use_parallel, maxiters=maxiters, maxtime=maxtime)
end

println("-" ^ 80)
println("\nResults:")
println("  Execution time: $(@sprintf("%.2f", time_m1)) seconds")
println("  Substances: $(length(res_m1.Name))")
Telu_max_val = isa(Telu_max_m1, Vector) ? maximum(Telu_max_m1) : Telu_max_m1
println("  Maximum elution temperature: $(@sprintf("%.2f", Telu_max_val)) K")
println()

# Show some results
println("Sample results (first 3 substances):")
for i in 1:min(3, length(res_m1.Name))
    println("  $(res_m1.Name[i]):")
    println("    Tchar = $(@sprintf("%.2f", res_m1.Tchar[i])) K")
    println("    θchar = $(@sprintf("%.4f", res_m1.θchar[i]))")
    println("    ΔCp = $(@sprintf("%.2f", res_m1.ΔCp[i])) J/(mol·K)")
    println("    Loss = $(@sprintf("%.6f", res_m1.min[i]))")
end

println("\n" * "=" ^ 80)
println("Benchmark complete!")
println("=" ^ 80)
