#!/usr/bin/env julia
"""
Benchmark script for multistart optimization with method_m1 using single chromatograms.

For each chromatogram in the data file, creates a single-chromatogram measurement
and compares optimization with and without multistart.

Usage:
    julia benchmark_multistart_m1.jl [data_file] [multistart_n] [maxiters] [maxtime]

Examples:
    julia benchmark_multistart_m1.jl
    julia benchmark_multistart_m1.jl data/meas_df05_Rxi5SilMS.csv 10
    julia benchmark_multistart_m1.jl data/meas_df05_Rxi5SilMS.csv 10 5000 300.0
"""

# Activate the package environment
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using RetentionParameterEstimator
using GasChromatographySimulator
using Printf
using Dates

# Default values
default_data_file = joinpath(@__DIR__, "..", "data", "meas_df05_Rxi5SilMS.csv")
default_multistart_n = 10
default_maxiters = 10000
default_maxtime = 600.0

# Parse command-line arguments
data_file = default_data_file
multistart_n = default_multistart_n
maxiters = default_maxiters
maxtime = default_maxtime

for arg in ARGS
    if endswith(arg, ".csv")
        data_file = arg
    elseif multistart_n == default_multistart_n
        try
            multistart_n = parse(Int, arg)
        catch
            println("Warning: Could not parse multistart_n: $arg, using default: $multistart_n")
        end
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

println("=" ^ 80)
println("Benchmark: method_m1 with multistart (single chromatograms)")
println("=" ^ 80)
println("\nConfiguration:")
println("  Data file: $data_file")
println("  Multistart runs: $multistart_n")
println("  Optimization limits: maxiters=$maxiters, maxtime=$maxtime")
println()

# Load data
println("Loading data...")
meas_full = RetentionParameterEstimator.load_chromatograms(data_file)
col_input = (L=meas_full[1].L, d=meas_full[1].d*1e3)  # Convert d from m to mm
n_measurements = length(meas_full[2])
n_substances = length(meas_full[4])
println("  Loaded $n_substances substances with $n_measurements measurements")
println()

# Results storage
results = []

# Process each chromatogram separately
for meas_idx in 1:n_measurements
    println("-" ^ 80)
    println("Processing measurement $meas_idx of $n_measurements: $(meas_full[3].measurement[meas_idx])")
    println("-" ^ 80)
    
    # Create single-chromatogram measurement
    meas_single = (
        meas_full[1],  # Column info (same)
        [meas_full[2][meas_idx]],  # Single program
        meas_full[3][meas_idx:meas_idx, :],  # Single row
        meas_full[4],  # All substances
        meas_full[5],  # pout
        meas_full[6]   # time_unit
    )
    
    # Define column
    col = GasChromatographySimulator.Column(col_input.L, col_input.d*1e-3, meas_full[1].df, meas_full[1].sp, meas_full[1].gas)
    
    # Get initial estimates
    Tchar_est, θchar_est, ΔCp_est, Telu_max = RetentionParameterEstimator.estimate_start_parameter(
        meas_single[3], col, meas_single[2]; time_unit=meas_single[6]
    )
    
    # Benchmark without multistart
    println("  Running without multistart...")
    time_no_multistart = @elapsed begin
        res_no_multistart = RetentionParameterEstimator.estimate_parameters(
            meas_single[3], meas_single[4], col, meas_single[2], 
            Tchar_est, θchar_est, ΔCp_est;
            mode="Kcentric_single", pout=meas_single[5], time_unit=meas_single[6],
            method=RetentionParameterEstimator.NewtonTrustRegion(), 
            opt=RetentionParameterEstimator.std_opt,
            maxiters=maxiters, maxtime=maxtime, multistart_n=0
        )[1]
    end
    
    # Benchmark with multistart
    println("  Running with multistart (n=$multistart_n)...")
    time_multistart = @elapsed begin
        res_multistart = RetentionParameterEstimator.estimate_parameters(
            meas_single[3], meas_single[4], col, meas_single[2], 
            Tchar_est, θchar_est, ΔCp_est;
            mode="Kcentric_single", pout=meas_single[5], time_unit=meas_single[6],
            method=RetentionParameterEstimator.NewtonTrustRegion(), 
            opt=RetentionParameterEstimator.std_opt,
            maxiters=maxiters, maxtime=maxtime, multistart_n=multistart_n
        )[1]
    end
    
    # Calculate average loss improvement
    avg_loss_no_multistart = mean(res_no_multistart.min)
    avg_loss_multistart = mean(res_multistart.min)
    loss_improvement = (avg_loss_no_multistart - avg_loss_multistart) / avg_loss_no_multistart * 100
    
    # Store results
    push!(results, (
        measurement=meas_single[3].measurement[1],
        time_no_multistart=time_no_multistart,
        time_multistart=time_multistart,
        avg_loss_no_multistart=avg_loss_no_multistart,
        avg_loss_multistart=avg_loss_multistart,
        loss_improvement=loss_improvement,
        speedup=time_no_multistart / time_multistart
    ))
    
    println("  Results:")
    println("    Time (no multistart): $(@sprintf("%.2f", time_no_multistart)) s")
    println("    Time (multistart):    $(@sprintf("%.2f", time_multistart)) s")
    println("    Speedup:              $(@sprintf("%.2f", time_no_multistart / time_multistart))x")
    println("    Avg loss (no multistart): $(@sprintf("%.6f", avg_loss_no_multistart))")
    println("    Avg loss (multistart):    $(@sprintf("%.6f", avg_loss_multistart))")
    println("    Loss improvement:         $(@sprintf("%.2f", loss_improvement))%")
    println()
end

# Summary
println("=" ^ 80)
println("Summary")
println("=" ^ 80)
println("\nAverage across all measurements:")
avg_time_no_multistart = mean([r.time_no_multistart for r in results])
avg_time_multistart = mean([r.time_multistart for r in results])
avg_loss_improvement = mean([r.loss_improvement for r in results])
avg_speedup = mean([r.speedup for r in results])

println("  Average time (no multistart): $(@sprintf("%.2f", avg_time_no_multistart)) s")
println("  Average time (multistart):    $(@sprintf("%.2f", avg_time_multistart)) s")
println("  Average speedup:              $(@sprintf("%.2f", avg_speedup))x")
println("  Average loss improvement:    $(@sprintf("%.2f", avg_loss_improvement))%")

println("\nPer-measurement details:")
for r in results
    println("  $(r.measurement):")
    println("    Time: $(@sprintf("%.2f", r.time_no_multistart))s → $(@sprintf("%.2f", r.time_multistart))s ($(@sprintf("%.2f", r.speedup))x)")
    println("    Loss: $(@sprintf("%.6f", r.avg_loss_no_multistart)) → $(@sprintf("%.6f", r.avg_loss_multistart)) ($(@sprintf("%+.2f", r.loss_improvement))%)")
end

println("\n" * "=" ^ 80)
println("Benchmark complete!")
println("=" ^ 80)
