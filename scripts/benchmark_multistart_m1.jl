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
using Logging
using CSV
using DataFrames

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

# Suppress warnings during benchmarks, but allow info/warn for multistart logging
old_logger = global_logger()
# Use a logger that filters out most warnings but allows our multistart messages
# We'll use ConsoleLogger with a custom filter or just use the default logger
# For now, keep NullLogger to suppress warnings, but multistart verbose logging
# can be enabled by setting verbose=true in estimate_parameters (if we add that parameter)
global_logger(NullLogger())

println("=" ^ 80)
println("Benchmark: method_m1 with multistart (single chromatograms)")
println("=" ^ 80)
println("\nConfiguration:")
println("  Data file: $data_file")
println("  Multistart runs: $multistart_n")
println("  Optimization limits: maxiters=$maxiters, maxtime=$maxtime")
println("\nNote: Multistart can sometimes be faster if non-multistart hits time limits.")
println("      Multistart should never give worse loss than non-multistart (it tries the")
println("      original starting point first). If it does, it's likely due to:")
println("      - Non-determinism when hitting time/iteration limits")
println("      - Different convergence paths due to numerical precision")
println("      - The first optimization in multistart failing (caught silently)")
println()

# Load database for comparison (if available)
db_file = joinpath(@__DIR__, "..", "data", "database_Rxi5SilMS_beta125.csv")
global db = nothing
if isfile(db_file)
    try
        global db = DataFrame(CSV.File(db_file))
        println("  Database file found: $db_file")
        println("  Will compare results with database parameters")
        println("  Loaded $(nrow(db)) substances from database")
    catch e
        println("  Warning: Could not load database file: $db_file")
        println("  Error: $(typeof(e)) - $(e)")
        global db = nothing  # Explicitly set to nothing on error
    end
else
    println("  Database file not found: $db_file")
end
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
            maxiters=maxiters, maxtime=maxtime, multistart_n=0, coupled_perturbation=true
        )[1]
    end
    
    # Benchmark with multistart
    println("  Running with multistart (n=$multistart_n, coupled_perturbation=true)...")
    time_multistart = @elapsed begin
        res_multistart = RetentionParameterEstimator.estimate_parameters(
            meas_single[3], meas_single[4], col, meas_single[2], 
            Tchar_est, θchar_est, ΔCp_est;
            mode="Kcentric_single", pout=meas_single[5], time_unit=meas_single[6],
            method=RetentionParameterEstimator.NewtonTrustRegion(), 
            opt=RetentionParameterEstimator.std_opt,
            maxiters=maxiters, maxtime=maxtime, multistart_n=multistart_n, coupled_perturbation=true
        )[1]
    end
    
    # Calculate average loss improvement
    avg_loss_no_multistart = mean(res_no_multistart.min)
    avg_loss_multistart = mean(res_multistart.min)
    loss_improvement = if avg_loss_no_multistart > 0
        (avg_loss_no_multistart - avg_loss_multistart) / avg_loss_no_multistart * 100
    else
        0.0
    end
    
    # Compare with database if available (for both no multistart and multistart)
    db_comparison_no_multistart = nothing
    db_comparison_multistart = nothing
    
    # Helper function to compute database comparison
    function compute_db_comparison(res_df, label)
        if db !== nothing
            try
                # Extract values from Measurements if needed
                res_values = DataFrame(
                    Name=res_df.Name,
                    Tchar=[isa(x, Measurements.Measurement) ? Measurements.value(x) : x for x in res_df.Tchar],
                    θchar=[isa(x, Measurements.Measurement) ? Measurements.value(x) : x for x in res_df.θchar],
                    ΔCp=[isa(x, Measurements.Measurement) ? Measurements.value(x) : x for x in res_df.ΔCp]
                )
                diff = RetentionParameterEstimator.difference_estimation_to_alternative_data(res_values, db)
                
                # Check if any matches were found (not all NaN)
                valid_relΔTchar = collect(skipmissing(abs.(diff.relΔTchar)))
                valid_relΔθchar = collect(skipmissing(abs.(diff.relΔθchar)))
                valid_relΔΔCp = collect(skipmissing(abs.(diff.relΔΔCp)))
                
                if length(valid_relΔTchar) > 0 && length(valid_relΔθchar) > 0 && length(valid_relΔΔCp) > 0
                    return (
                        mean_abs_relΔTchar=mean(valid_relΔTchar) * 100,
                        mean_abs_relΔθchar=mean(valid_relΔθchar) * 100,
                        mean_abs_relΔΔCp=mean(valid_relΔΔCp) * 100,
                        n_matched=length(valid_relΔTchar)
                    )
                else
                    # No matches found - substances in results don't match database
                    return (n_matched=0,)
                end
            catch e
                # If comparison fails, show error for debugging
                println("    Warning: Database comparison failed for $label: $(typeof(e)) - $(e)")
                return nothing
            end
        else
            return nothing
        end
    end
    
    db_comparison_no_multistart = compute_db_comparison(res_no_multistart, "no multistart")
    db_comparison_multistart = compute_db_comparison(res_multistart, "multistart")
    
    # Store results
    push!(results, (
        measurement=meas_single[3].measurement[1],
        time_no_multistart=time_no_multistart,
        time_multistart=time_multistart,
        avg_loss_no_multistart=avg_loss_no_multistart,
        avg_loss_multistart=avg_loss_multistart,
        loss_improvement=loss_improvement,
        speedup=time_no_multistart / time_multistart,
        db_comparison_no_multistart=db_comparison_no_multistart,
        db_comparison_multistart=db_comparison_multistart
    ))
    
    # Format loss values appropriately
    format_loss(loss) = if loss < 1e-6
        @sprintf("%.2e", loss)
    elseif loss < 1.0
        @sprintf("%.8f", loss)
    else
        @sprintf("%.6f", loss)
    end
    
    println("  Results:")
    println("    Time (no multistart): $(@sprintf("%.2f", time_no_multistart)) s")
    println("    Time (multistart):    $(@sprintf("%.2f", time_multistart)) s")
    println("    Speedup:              $(@sprintf("%.2f", time_no_multistart / time_multistart))x")
    println("    Avg loss (no multistart): $(format_loss(avg_loss_no_multistart))")
    println("    Avg loss (multistart):    $(format_loss(avg_loss_multistart))")
    if loss_improvement != 0.0
        println("    Loss improvement:         $(@sprintf("%.2f", loss_improvement))%")
    end
    # Database comparison for no multistart
    if db_comparison_no_multistart !== nothing
        if haskey(db_comparison_no_multistart, :n_matched) && db_comparison_no_multistart.n_matched == 0
            println("    Database comparison (no multistart): No matching substances found (check CAS numbers)")
        elseif haskey(db_comparison_no_multistart, :mean_abs_relΔTchar)
            println("    Database comparison (no multistart, $(db_comparison_no_multistart.n_matched) substances):")
            println("      Mean |rel ΔTchar|: $(@sprintf("%.2f", db_comparison_no_multistart.mean_abs_relΔTchar))%")
            println("      Mean |rel Δθchar|: $(@sprintf("%.2f", db_comparison_no_multistart.mean_abs_relΔθchar))%")
            println("      Mean |rel ΔΔCp|:   $(@sprintf("%.2f", db_comparison_no_multistart.mean_abs_relΔΔCp))%")
        end
    end
    
    # Database comparison for multistart
    if db_comparison_multistart !== nothing
        if haskey(db_comparison_multistart, :n_matched) && db_comparison_multistart.n_matched == 0
            println("    Database comparison (multistart): No matching substances found (check CAS numbers)")
        elseif haskey(db_comparison_multistart, :mean_abs_relΔTchar)
            println("    Database comparison (multistart, $(db_comparison_multistart.n_matched) substances):")
            println("      Mean |rel ΔTchar|: $(@sprintf("%.2f", db_comparison_multistart.mean_abs_relΔTchar))%")
            println("      Mean |rel Δθchar|: $(@sprintf("%.2f", db_comparison_multistart.mean_abs_relΔθchar))%")
            println("      Mean |rel ΔΔCp|:   $(@sprintf("%.2f", db_comparison_multistart.mean_abs_relΔΔCp))%")
        end
    end
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

# Format loss values appropriately
format_loss(loss) = if loss < 1e-6
    @sprintf("%.2e", loss)
elseif loss < 1.0
    @sprintf("%.8f", loss)
else
    @sprintf("%.6f", loss)
end

println("\nPer-measurement details:")
for r in results
    println("  $(r.measurement):")
    println("    Time: $(@sprintf("%.2f", r.time_no_multistart))s → $(@sprintf("%.2f", r.time_multistart))s ($(@sprintf("%.2f", r.speedup))x)")
    println("    Loss: $(format_loss(r.avg_loss_no_multistart)) → $(format_loss(r.avg_loss_multistart))", 
            r.loss_improvement != 0.0 ? " ($(@sprintf("%+.2f", r.loss_improvement))%)" : "")
    if r.db_comparison_no_multistart !== nothing && haskey(r.db_comparison_no_multistart, :mean_abs_relΔTchar)
        println("    DB no multistart ($(r.db_comparison_no_multistart.n_matched) matched): |rel ΔTchar|=$(@sprintf("%.2f", r.db_comparison_no_multistart.mean_abs_relΔTchar))%, |rel Δθchar|=$(@sprintf("%.2f", r.db_comparison_no_multistart.mean_abs_relΔθchar))%, |rel ΔΔCp|=$(@sprintf("%.2f", r.db_comparison_no_multistart.mean_abs_relΔΔCp))%")
    end
    if r.db_comparison_multistart !== nothing && haskey(r.db_comparison_multistart, :mean_abs_relΔTchar)
        println("    DB multistart ($(r.db_comparison_multistart.n_matched) matched): |rel ΔTchar|=$(@sprintf("%.2f", r.db_comparison_multistart.mean_abs_relΔTchar))%, |rel Δθchar|=$(@sprintf("%.2f", r.db_comparison_multistart.mean_abs_relΔθchar))%, |rel ΔΔCp|=$(@sprintf("%.2f", r.db_comparison_multistart.mean_abs_relΔΔCp))%")
    end
end

# Restore logger
global_logger(old_logger)

println("\n" * "=" ^ 80)
println("Benchmark complete!")
println("=" ^ 80)
