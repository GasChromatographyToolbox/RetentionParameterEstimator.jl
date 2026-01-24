#!/usr/bin/env julia
"""
    benchmark_methods.jl

Script to benchmark and compare all optimization methods (method_m1, m2, m3, m4) performance and results.
Supports parallelization if Julia is started with multiple threads (e.g., `julia -t 4 benchmark_methods.jl`).

Usage:
    julia benchmark_methods.jl [data_file]
    julia -t 4 benchmark_methods.jl [data_file]  # With 4 threads for parallelization

If no data_file is provided, uses ../data/meas_df05_Rxi5SilMS.csv (relative to scripts/)
"""

# Activate the package environment (go up one level from scripts/ to package root)
using Pkg
Pkg.activate(dirname(@__DIR__))

using RetentionParameterEstimator
using Printf
using CSV
using DataFrames

function benchmark_and_compare(data_file::String; use_parallel::Bool=true)
    println("=" ^ 80)
    println("Benchmarking All Methods: method_m1, m2, m3, m4")
    println("=" ^ 80)
    println("\nData file: $data_file")
    
    # Check thread availability
    nthreads = Threads.nthreads()
    if use_parallel && nthreads > 1
        println("\nParallelization: ENABLED ($nthreads threads)")
        parallel_flag = true
    else
        if use_parallel && nthreads == 1
            println("\nParallelization: REQUESTED but only 1 thread available")
            println("  Start Julia with multiple threads: julia -t 4 benchmark_methods.jl")
        else
            println("\nParallelization: DISABLED")
        end
        parallel_flag = false
    end
    
    # Load data
    println("\nLoading data...")
    meas = RetentionParameterEstimator.load_chromatograms(data_file)
    println("  Loaded $(length(meas[4])) substances")
    println("  Loaded $(length(meas[3].measurement)) measurements")
    
    # Prepare col_input for method_m1 (needs L and d in mm)
    col_input = (L=meas[1].L, d=meas[1].d*1000)  # Convert d from m to mm
    
    # Run method_m1
    println("\n" * "-" ^ 80)
    println("Running method_m1 (requires known d)...")
    if parallel_flag
        println("  Using parallelization: $(parallel_flag)")
    end
    println("-" ^ 80)
    time_m1 = @elapsed begin
        res_m1, Telu_max_m1 = RetentionParameterEstimator.method_m1(meas, col_input, se_col=true, parallel=parallel_flag)
    end
    
    println("  Completed in $(@sprintf("%.2f", time_m1)) seconds")
    println("  Estimated parameters for $(length(res_m1.Name)) substances")
    println("  Note: method_m1 uses fixed d = $(@sprintf("%.6f", meas[1].d)) m")
    
    # Run method_m2
    println("\n" * "-" ^ 80)
    println("Running method_m2 (two-pass heuristic)...")
    if parallel_flag
        println("  Using parallelization: $(parallel_flag)")
    end
    println("-" ^ 80)
    time_m2 = @elapsed begin
        res_m2, Telu_max_m2 = RetentionParameterEstimator.method_m2(meas, se_col=true, parallel=parallel_flag)
    end
    
    println("  Completed in $(@sprintf("%.2f", time_m2)) seconds")
    println("  Estimated parameters for $(length(res_m2.Name)) substances")
    
    # Run method_m3
    println("\n" * "-" ^ 80)
    println("Running method_m3 (joint optimization - may be slow for many substances)...")
    if parallel_flag
        println("  Using parallelization: $(parallel_flag) (for standard error calculation)")
    end
    println("-" ^ 80)
    time_m3 = @elapsed begin
        res_m3, Telu_max_m3 = RetentionParameterEstimator.method_m3(meas, se_col=true, parallel=parallel_flag)
    end
    
    println("  Completed in $(@sprintf("%.2f", time_m3)) seconds")
    println("  Estimated parameters for $(length(res_m3.Name)) substances")
    
    # Run method_m4
    println("\n" * "-" ^ 80)
    println("Running method_m4 (alternating optimization)...")
    if parallel_flag
        println("  Using parallelization: $(parallel_flag) (for Block 2 and standard error calculation)")
    end
    println("-" ^ 80)
    time_m4 = @elapsed begin
        res_m4, Telu_max_m4 = RetentionParameterEstimator.method_m4(meas, se_col=true, parallel=parallel_flag)
    end
    
    println("  Completed in $(@sprintf("%.2f", time_m4)) seconds")
    println("  Estimated parameters for $(length(res_m4.Name)) substances")
    
    # Run method_m4_ (previous implementation)
    println("\n" * "-" ^ 80)
    println("Running method_m4_ (alternating optimization - previous implementation)...")
    if parallel_flag
        println("  Using parallelization: $(parallel_flag) (for Block 2 and standard error calculation)")
    end
    println("-" ^ 80)
    time_m4_ = @elapsed begin
        res_m4_, Telu_max_m4_ = RetentionParameterEstimator.method_m4_(meas, se_col=true, parallel=parallel_flag)
    end
    
    println("  Completed in $(@sprintf("%.2f", time_m4_)) seconds")
    println("  Estimated parameters for $(length(res_m4_.Name)) substances")
    
    # Compare results
    println("\n" * "=" ^ 80)
    println("Performance Comparison")
    println("=" ^ 80)
    
    times = [time_m1, time_m2, time_m3, time_m4, time_m4_]
    methods = ["method_m1", "method_m2", "method_m3", "method_m4", "method_m4_"]
    fastest_idx = argmin(times)
    
    println("\nExecution Times:")
    for (i, (method, time)) in enumerate(zip(methods, times))
        speedup = time / times[fastest_idx]
        if i == fastest_idx
            println("  $method: $(@sprintf("%8.2f", time)) seconds (fastest)")
        else
            println("  $method: $(@sprintf("%8.2f", time)) seconds ($(@sprintf("%.2fx", speedup)) slower)")
        end
    end
    
    # Parameter comparison (compare m2, m3, m4 against m1 as baseline, or m4 as reference)
    println("\n" * "=" ^ 80)
    println("Parameter Comparison (using method_m4 as reference)")
    println("=" ^ 80)
    
    # Compare each method against m4
    methods_to_compare = [
        ("m1", res_m1, Telu_max_m1, false),  # m1 doesn't optimize d
        ("m2", res_m2, Telu_max_m2, true),
        ("m3", res_m3, Telu_max_m3, true),
        ("m4", res_m4, Telu_max_m4, true),
        ("m4_", res_m4_, Telu_max_m4_, true)
    ]
    
    for (method_name, res, Telu_max, has_d) in methods_to_compare
        println("\n" * "-" ^ 80)
        println("method_$method_name vs method_m4:")
        println("-" ^ 80)
        
        max_tchar_diff = 0.0
        max_thetachar_diff = 0.0
        max_deltacp_diff = 0.0
        max_d_diff = 0.0
        max_min_diff = 0.0
        
        for i in 1:length(res_m4.Name)
            # Find corresponding substance in res
            idx = findfirst(res.Name .== res_m4.Name[i])
            if isnothing(idx)
                continue
            end
            
            tchar_diff = abs(res.Tchar[idx] - res_m4.Tchar[i])
            thetachar_diff = abs(res.θchar[idx] - res_m4.θchar[i])
            deltacp_diff = abs(res.ΔCp[idx] - res_m4.ΔCp[i])
            min_diff = abs(res.min[idx] - res_m4.min[i])
            
            max_tchar_diff = max(max_tchar_diff, tchar_diff)
            max_thetachar_diff = max(max_thetachar_diff, thetachar_diff)
            max_deltacp_diff = max(max_deltacp_diff, deltacp_diff)
            max_min_diff = max(max_min_diff, min_diff)
            
            if has_d && "d" in names(res)
                d_diff = abs(res.d[idx] - res_m4.d[i])
                max_d_diff = max(max_d_diff, d_diff)
            end
        end
        
        println("  Maximum Differences:")
        println("    Tchar:    $(@sprintf("%.2f", max_tchar_diff)) K")
        println("    θchar:    $(@sprintf("%.4f", max_thetachar_diff))")
        println("    ΔCp:      $(@sprintf("%.2f", max_deltacp_diff)) J/(mol·K)")
        if has_d && "d" in names(res)
            println("    d:        $(@sprintf("%.2e", max_d_diff)) m")
        else
            println("    d:        N/A (not optimized)")
        end
        println("    min:      $(@sprintf("%.2e", max_min_diff))")
    end
    
    # Column diameter consistency check
    println("\n" * "=" ^ 80)
    println("Column Diameter (d) Consistency Check")
    println("=" ^ 80)
    
    for (method_name, res, _, has_d) in methods_to_compare
        if has_d && "d" in names(res)
            d_unique = length(unique(res.d))
            d_mean = mean(res.d)
            d_std = std(res.d)
            
            println("\nmethod_$method_name:")
            println("  Unique d values: $d_unique")
            if d_unique == 1
                println("  ✓ Correctly enforces same d for all substances")
                println("  d = $(@sprintf("%.6f", res.d[1])) m")
            else
                println("  ✗ WARNING: d varies across substances (violates constraint)")
                println("  Mean d = $(@sprintf("%.6f", d_mean)) m, std = $(@sprintf("%.6e", d_std)) m")
            end
        else
            println("\nmethod_$method_name:")
            println("  d is fixed (not optimized)")
        end
    end
    
    # Telu_max comparison
    println("\n" * "=" ^ 80)
    println("Maximum Elution Temperature (Telu_max) Comparison")
    println("=" ^ 80)
    
    Telu_max_values = []
    for (method_name, _, Telu_max, _) in methods_to_compare
        Telu_max_val = isa(Telu_max, Vector) ? maximum(Telu_max) : Telu_max
        push!(Telu_max_values, Telu_max_val)
        println("  method_$method_name: $(@sprintf("%.2f", Telu_max_val)) K")
    end
    
    if all(x -> abs(x - Telu_max_values[1]) < 0.01, Telu_max_values)
        println("  ✓ All methods agree on Telu_max")
    else
        println("  ✗ WARNING: Different Telu_max values across methods!")
    end
    
    # Detailed parameter comparison table
    println("\n" * "=" ^ 80)
    println("Detailed Parameter Comparison (All Substances)")
    println("=" ^ 80)
    
    # Create a comparison DataFrame
    comparison_data = []
    for i in 1:length(res_m4.Name)
        substance_name = res_m4.Name[i]
        
        # Find corresponding indices in each result
        idx_m1 = findfirst(res_m1.Name .== substance_name)
        idx_m2 = findfirst(res_m2.Name .== substance_name)
        idx_m3 = findfirst(res_m3.Name .== substance_name)
        idx_m4 = i
        idx_m4_ = findfirst(res_m4_.Name .== substance_name)
        
        row = Dict(
            "Substance" => substance_name,
            "m1_Tchar" => isnothing(idx_m1) ? missing : res_m1.Tchar[idx_m1],
            "m1_θchar" => isnothing(idx_m1) ? missing : res_m1.θchar[idx_m1],
            "m1_ΔCp" => isnothing(idx_m1) ? missing : res_m1.ΔCp[idx_m1],
            "m1_min" => isnothing(idx_m1) ? missing : res_m1.min[idx_m1],
            "m2_Tchar" => isnothing(idx_m2) ? missing : res_m2.Tchar[idx_m2],
            "m2_θchar" => isnothing(idx_m2) ? missing : res_m2.θchar[idx_m2],
            "m2_ΔCp" => isnothing(idx_m2) ? missing : res_m2.ΔCp[idx_m2],
            "m2_d" => isnothing(idx_m2) ? missing : res_m2.d[idx_m2],
            "m2_min" => isnothing(idx_m2) ? missing : res_m2.min[idx_m2],
            "m3_Tchar" => isnothing(idx_m3) ? missing : res_m3.Tchar[idx_m3],
            "m3_θchar" => isnothing(idx_m3) ? missing : res_m3.θchar[idx_m3],
            "m3_ΔCp" => isnothing(idx_m3) ? missing : res_m3.ΔCp[idx_m3],
            "m3_d" => isnothing(idx_m3) ? missing : res_m3.d[idx_m3],
            "m3_min" => isnothing(idx_m3) ? missing : res_m3.min[idx_m3],
            "m4_Tchar" => res_m4.Tchar[idx_m4],
            "m4_θchar" => res_m4.θchar[idx_m4],
            "m4_ΔCp" => res_m4.ΔCp[idx_m4],
            "m4_d" => res_m4.d[idx_m4],
            "m4_min" => res_m4.min[idx_m4],
            "m4__Tchar" => isnothing(idx_m4_) ? missing : res_m4_.Tchar[idx_m4_],
            "m4__θchar" => isnothing(idx_m4_) ? missing : res_m4_.θchar[idx_m4_],
            "m4__ΔCp" => isnothing(idx_m4_) ? missing : res_m4_.ΔCp[idx_m4_],
            "m4__d" => isnothing(idx_m4_) ? missing : res_m4_.d[idx_m4_],
            "m4__min" => isnothing(idx_m4_) ? missing : res_m4_.min[idx_m4_],
        )
        push!(comparison_data, row)
    end
    
    comparison_df = DataFrame(comparison_data)
    
    # Print summary table
    println("\nSummary Table (first 5 substances):")
    println("=" ^ 80)
    for i in 1:min(5, nrow(comparison_df))
        substance_name = comparison_df.Substance[i]
        idx_m1 = findfirst(res_m1.Name .== substance_name)
        idx_m2 = findfirst(res_m2.Name .== substance_name)
        idx_m3 = findfirst(res_m3.Name .== substance_name)
        
        println("\n$(substance_name):")
        println("  Tchar (K):    m1=$(ismissing(comparison_df.m1_Tchar[i]) ? "N/A" : @sprintf("%.2f", comparison_df.m1_Tchar[i])), m2=$(ismissing(comparison_df.m2_Tchar[i]) ? "N/A" : @sprintf("%.2f", comparison_df.m2_Tchar[i])), m3=$(ismissing(comparison_df.m3_Tchar[i]) ? "N/A" : @sprintf("%.2f", comparison_df.m3_Tchar[i])), m4=$(@sprintf("%.2f", comparison_df.m4_Tchar[i])), m4_=$(ismissing(comparison_df.m4__Tchar[i]) ? "N/A" : @sprintf("%.2f", comparison_df.m4__Tchar[i]))")
        println("  θchar:        m1=$(ismissing(comparison_df.m1_θchar[i]) ? "N/A" : @sprintf("%.2f", comparison_df.m1_θchar[i])), m2=$(ismissing(comparison_df.m2_θchar[i]) ? "N/A" : @sprintf("%.2f", comparison_df.m2_θchar[i])), m3=$(ismissing(comparison_df.m3_θchar[i]) ? "N/A" : @sprintf("%.2f", comparison_df.m3_θchar[i])), m4=$(@sprintf("%.2f", comparison_df.m4_θchar[i])), m4_=$(ismissing(comparison_df.m4__θchar[i]) ? "N/A" : @sprintf("%.2f", comparison_df.m4__θchar[i]))")
        println("  ΔCp (J/mol·K): m1=$(ismissing(comparison_df.m1_ΔCp[i]) ? "N/A" : @sprintf("%.2f", comparison_df.m1_ΔCp[i])), m2=$(ismissing(comparison_df.m2_ΔCp[i]) ? "N/A" : @sprintf("%.2f", comparison_df.m2_ΔCp[i])), m3=$(ismissing(comparison_df.m3_ΔCp[i]) ? "N/A" : @sprintf("%.2f", comparison_df.m3_ΔCp[i])), m4=$(@sprintf("%.2f", comparison_df.m4_ΔCp[i])), m4_=$(ismissing(comparison_df.m4__ΔCp[i]) ? "N/A" : @sprintf("%.2f", comparison_df.m4__ΔCp[i]))")
        if !ismissing(comparison_df.m2_d[i])
            println("  d (m):        m2=$(@sprintf("%.6f", comparison_df.m2_d[i])), m3=$(@sprintf("%.6f", comparison_df.m3_d[i])), m4=$(@sprintf("%.6f", comparison_df.m4_d[i])), m4_=$(ismissing(comparison_df.m4__d[i]) ? "N/A" : @sprintf("%.6f", comparison_df.m4__d[i]))")
        end
        println("  min:          m1=$(ismissing(comparison_df.m1_min[i]) ? "N/A" : @sprintf("%.4e", comparison_df.m1_min[i])), m2=$(ismissing(comparison_df.m2_min[i]) ? "N/A" : @sprintf("%.4e", comparison_df.m2_min[i])), m3=$(ismissing(comparison_df.m3_min[i]) ? "N/A" : @sprintf("%.4e", comparison_df.m3_min[i])), m4=$(@sprintf("%.4e", comparison_df.m4_min[i])), m4_=$(ismissing(comparison_df.m4__min[i]) ? "N/A" : @sprintf("%.4e", comparison_df.m4__min[i]))")
    end
    
    if nrow(comparison_df) > 5
        println("\n... (showing first 5 of $(nrow(comparison_df)) substances)")
        println("Full comparison data available in returned DataFrame")
    end
    
    # Save full comparison to CSV for detailed analysis
    comparison_file = replace(data_file, r"\.csv$" => "_comparison.csv")
    try
        CSV.write(comparison_file, comparison_df)
        println("\nFull comparison saved to: $comparison_file")
    catch e
        println("\nWarning: Could not save comparison file: $e")
    end
    
    # Summary
    println("\n" * "=" ^ 80)
    println("Summary")
    println("=" ^ 80)
    println("\nPerformance Ranking (fastest to slowest):")
    sorted_indices = sortperm(times)
    for (rank, idx) in enumerate(sorted_indices)
        println("  $rank. $(methods[idx]): $(@sprintf("%.2f", times[idx]))s")
    end
    
    println("\nKey Observations:")
    println("  • method_m1: Fastest, but requires known d")
    println("  • method_m2: Two-pass heuristic (legacy method)")
    println("  • method_m3: Joint optimization (slow for many substances)")
    println("  • method_m4: Alternating optimization (uses estimate_parameters mode=\"d_only\")")
    println("  • method_m4_: Alternating optimization (previous implementation with direct optimize_d_only call)")
    
    if parallel_flag
        println("\n  ✓ Parallelization was used ($nthreads threads)")
        println("    Expected speedup: 2-4x for methods with many substances")
    else
        if length(meas[4]) > 5
            println("\n  💡 Tip: Start Julia with multiple threads for better performance:")
            println("    julia -t 4 benchmark_methods.jl")
        end
    end
    
    if length(meas[4]) > 10
        println("\n  ⚠️  Note: With $(length(meas[4])) substances, method_m3 may be impractical.")
        println("     Consider using method_m4 for better performance.")
    end
    
    return (m1=(res_m1, Telu_max_m1, time_m1), 
            m2=(res_m2, Telu_max_m2, time_m2),
            m3=(res_m3, Telu_max_m3, time_m3),
            m4=(res_m4, Telu_max_m4, time_m4),
            m4_=(res_m4_, Telu_max_m4_, time_m4_))
end

# Main execution
# When run directly from command line
if abspath(PROGRAM_FILE) == @__FILE__
    data_file = length(ARGS) > 0 ? ARGS[1] : joinpath(dirname(@__DIR__), "data", "meas_df05_Rxi5SilMS.csv")
    
    if !isfile(data_file)
        println("Error: File not found: $data_file")
        exit(1)
    end
    
    try
        benchmark_and_compare(data_file)
    catch e
        println("\nError during benchmarking:")
        println(e)
        rethrow(e)
    end
end

# When included, automatically run with default file (unless already executed above)
if !(abspath(PROGRAM_FILE) == @__FILE__)
    data_file = joinpath(dirname(@__DIR__), "data", "meas_df05_Rxi5SilMS.csv")
    if isfile(data_file)
        try
            benchmark_and_compare(data_file)
        catch e
            println("\nError during benchmarking:")
            println(e)
            rethrow(e)
        end
    else
        println("Note: Default data file not found: $data_file")
        println("Call benchmark_and_compare(\"path/to/data.csv\") to run manually")
    end
end


#=
================================================================================
Performance Comparison
================================================================================

Execution Times:
  method_m1:    26.45 seconds (fastest)
  method_m2:   380.82 seconds (14.40x slower)
  method_m3:   653.17 seconds (24.70x slower)
  method_m4:   102.53 seconds (3.88x slower)

================================================================================
Parameter Comparison (using method_m4 as reference)
================================================================================

--------------------------------------------------------------------------------
method_m1 vs method_m4:
--------------------------------------------------------------------------------
  Maximum Differences:
    Tchar:    0.78 K
    θchar:    0.5293
    ΔCp:      6.76 J/(mol·K)
    d:        N/A (not optimized)
    min:      1.07e-02

--------------------------------------------------------------------------------
method_m2 vs method_m4:
--------------------------------------------------------------------------------
  Maximum Differences:
    Tchar:    0.06 K
    θchar:    0.0446
    ΔCp:      0.56 J/(mol·K)
    d:        9.78e-08 m
    min:      7.19e-04

--------------------------------------------------------------------------------
method_m3 vs method_m4:
--------------------------------------------------------------------------------
  Maximum Differences:
    Tchar:    0.30 K
    θchar:    0.4905
    ΔCp:      9.29 J/(mol·K)
    d:        1.98e-07 m
    min:      3.10e-02

--------------------------------------------------------------------------------
method_m4 vs method_m4:
--------------------------------------------------------------------------------
  Maximum Differences:
    Tchar:    0.00 K
    θchar:    0.0000
    ΔCp:      0.00 J/(mol·K)
    d:        0.00e+00 m
    min:      0.00e+00

================================================================================
Column Diameter (d) Consistency Check
================================================================================

method_m1:
  d is fixed (not optimized)

method_m2:
  Unique d values: 1
  ✓ Correctly enforces same d for all substances
  d = 0.000249 m

method_m3:
  Unique d values: 1
  ✓ Correctly enforces same d for all substances
  d = 0.000249 m

method_m4:
  Unique d values: 1
  ✓ Correctly enforces same d for all substances
  d = 0.000249 m

================================================================================
Maximum Elution Temperature (Telu_max) Comparison
================================================================================
  method_m1: 511.88 K
  method_m2: 511.88 K
  method_m3: 511.88 K
  method_m4: 511.88 K
  ✓ All methods agree on Telu_max

================================================================================
Summary
================================================================================

Performance Ranking (fastest to slowest):
  1. method_m1: 26.45s
  2. method_m4: 102.53s
  3. method_m2: 380.82s
  4. method_m3: 653.17s

Key Observations:
  • method_m1: Fastest, but requires known d
  • method_m2: Two-pass heuristic (legacy method)
  • method_m3: Joint optimization (slow for many substances)
  • method_m4: Alternating optimization (recommended when d unknown)

  ⚠️  Note: With 12 substances, method_m3 may be impractical.
     Consider using method_m4 for better performance.
(m1 = (12×8 DataFrame
 Row │ Name                 min         Tchar    Tchar_std  θchar    θchar_std  ΔCp       ΔCp_std  
     │ String               Float64     Float64  Float64    Float64  Float64    Float64   Float64  
─────┼─────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 2-Octanone           0.0164255   400.155    9.5489   33.916     6.26267   85.1842   65.2932
   2 │ 4-Fluoroaniline      0.00904903  404.695    9.4731   35.8462    6.2739    89.1576   60.1928
   3 │ 2,6-Dimethylphenol   0.0095027   425.133    5.88949  37.4436    3.29483   89.9029   28.923
   4 │ Octan-1-ol           0.0186872   414.889    6.28717  34.7332    3.60564  100.086    37.8009
   5 │ 1,3-Dichlorobenzene  0.0144235   408.75    10.0252   37.4897    6.58484   70.7259   53.6445
   6 │ Cyclohexanone        0.00495986  384.657   16.1867   35.5913   13.8822    74.3097  124.73
   7 │ 4-Methylphenol       0.0132451   415.87     6.20243  35.5661    3.64663  112.797    38.2814
   8 │ Aniline              0.00394122  400.697   11.0122   36.1166    7.68737   80.7639   70.2848
   9 │ n-Butylbenzene       0.0170755   416.243    7.56401  36.5514    4.45674   74.95     38.6359
  10 │ Nonan-2-one          0.0176068   419.028    5.88357  34.9426    3.22974   88.4841   31.99
  11 │ 1-Phenyl-2-propanol  0.0331945   430.899    5.34944  38.133     2.88332   84.7351   23.9226
  12 │ Hexan-1-ol           0.00577887  374.77    16.1206   31.8053   14.4595    88.7813  167.834, [479.36199999999997, 483.93399999999997, 506.038, 495.87399999999997, 487.234, 457.25199999999995, 497.39799999999997, 478.6, 496.38399999999996, 500.19399999999996, 511.88199999999995, 446.32599999999996], 26.447403125), m2 = (12×10 DataFrame
 Row │ Name                 min         Tchar    Tchar_std  θchar    θchar_std  ΔCp       ΔCp_std  d          ⋯
     │ String               Float64     Float64  Float64    Float64  Float64    Float64   Float64  Float64    ⋯
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 2-Octanone           0.0136734   399.509   0.704344  33.4689    1.14913   80.0556  23.5274  0.00024886 ⋯
   2 │ 4-Fluoroaniline      0.00879331  404.019   0.743752  35.3709    1.21616   84.6064  21.2541  0.00024886
   3 │ 2,6-Dimethylphenol   0.00665955  424.426   0.780516  36.9597    1.18921   85.7352  17.8909  0.00024886
   4 │ Octan-1-ol           0.0201893   414.229   0.702162  34.2788    1.09176   95.2299  20.3788  0.00024886
   5 │ 1,3-Dichlorobenzene  0.0140492   408.044   0.79956   37.0021    1.27883   66.5719  19.5556  0.00024886 ⋯
   6 │ Cyclohexanone        0.00535502  383.988   0.790737  35.1059    1.45147   69.4594  29.0472  0.00024886
   7 │ 4-Methylphenol       0.0149483   415.196   0.716893  35.0895    1.13892  108.113   20.0241  0.00024886
   8 │ Aniline              0.00688216  400.016   0.760557  35.6361    1.26398   76.2267  21.9607  0.00024886
   9 │ n-Butylbenzene       0.0153246   415.556   0.754314  36.0912    1.1545    70.7256  18.5971  0.00024886 ⋯
  10 │ Nonan-2-one          0.0115542   418.363   0.721137  34.4943    1.09421   83.7053  19.9385  0.00024886
  11 │ 1-Phenyl-2-propanol  0.0232539   430.179   0.801979  37.6451    1.19745   80.693   17.0887  0.00024886
  12 │ Hexan-1-ol           0.00470332  374.17    0.697393  31.3682    1.31761   82.5814  38.0477  0.00024886
                                                                                               1 column omitted, [479.36199999999997, 483.93399999999997, 506.038, 495.87399999999997, 487.234, 457.25199999999995, 497.39799999999997, 478.6, 496.38399999999996, 500.19399999999996, 511.88199999999995, 446.32599999999996], 380.819940792), m3 = (12×10 DataFrame
 Row │ Name                 min        Tchar    Tchar_std  θchar    θchar_std  ΔCp      ΔCp_std   d           ⋯
     │ String               Float64    Float64  Float64    Float64  Float64    Float64  Float64   Float64     ⋯
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 2-Octanone           0.0356077  399.682   13.0123   33.7153    9.13163  84.5767  106.4     0.00024896  ⋯
   2 │ 4-Fluoroaniline      0.0356077  404.034    8.55628  35.3518    5.98093  84.1589   60.2494  0.00024896
   3 │ 2,6-Dimethylphenol   0.0356077  424.453    9.05762  36.9518    6.18103  85.4235   55.6338  0.00024896
   4 │ Octan-1-ol           0.0356077  414.192   11.4504   34.1493    7.85423  92.4858   88.4792  0.00024896
   5 │ 1,3-Dichlorobenzene  0.0356077  408.281    7.76492  37.3945    5.40039  72.4731   46.6477  0.00024896  ⋯
   6 │ Cyclohexanone        0.0356077  384.12    14.4678   35.3316   10.5875   73.6836  108.261   0.00024896
   7 │ 4-Methylphenol       0.0356077  414.929   10.562    34.5619    7.27786  98.4706   78.2301  0.00024896
   8 │ Aniline              0.0356077  400.138   12.1145   35.8058    8.59881  78.9026   83.3591  0.00024896
   9 │ n-Butylbenzene       0.0356077  415.743    7.82723  36.3996    5.32996  75.7836   50.2746  0.00024896  ⋯
  10 │ Nonan-2-one          0.0356077  418.527    8.81693  34.7124    5.98724  87.3465   65.1923  0.00024896
  11 │ 1-Phenyl-2-propanol  0.0356077  430.348   10.064    37.8542    6.85536  83.3474   58.0011  0.00024896
  12 │ Hexan-1-ol           0.0356077  374.163   51.7132   31.2302   37.4613   77.713   562.029   0.00024896
                                                                                               1 column omitted, [479.36199999999997, 483.93399999999997, 506.038, 495.87399999999997, 487.234, 457.25199999999995, 497.39799999999997, 478.6, 496.38399999999996, 500.19399999999996, 511.88199999999995, 446.32599999999996], 653.165687), m4 = (12×10 DataFrame
 Row │ Name                 min         Tchar    Tchar_std  θchar    θchar_std  ΔCp       ΔCp_std  d          ⋯
     │ String               Float64     Float64  Float64    Float64  Float64    Float64   Float64  Float64    ⋯
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 2-Octanone           0.0135588   399.453   0.7005    33.4293    1.14144   79.5811  23.4567  0.00024876 ⋯
   2 │ 4-Fluoroaniline      0.00894906  403.961   0.740625  35.3291    1.20997   84.1841  21.2223  0.00024876
   3 │ 2,6-Dimethylphenol   0.00655225  424.365   0.779054  36.9182    1.18654   85.3666  17.9091  0.00024876
   4 │ Octan-1-ol           0.0204682   414.173   0.700765  34.2403    1.08915   94.8092  20.397   0.00024876
   5 │ 1,3-Dichlorobenzene  0.0141249   407.984   0.797927  36.9614    1.27582   66.2222  19.5734  0.00024876 ⋯
   6 │ Cyclohexanone        0.00543339  383.931   0.786812  35.0643    1.44227   69.0311  28.9757  0.00024876
   7 │ 4-Methylphenol       0.0152366   415.141   0.718867  35.0524    1.14272  107.76    20.1491  0.00024876
   8 │ Aniline              0.00686713  399.958   0.759161  35.5954    1.26131   75.8352  21.9893  0.00024876
   9 │ n-Butylbenzene       0.0152159   415.494   0.752759  36.0466    1.15159   70.2713  18.617   0.00024876 ⋯
  10 │ Nonan-2-one          0.0111875   418.306   0.719821  34.4562    1.09181   83.2908  19.9588  0.00024876
  11 │ 1-Phenyl-2-propanol  0.0225349   430.117   0.803225  37.6037    1.19946   80.3431  17.1698  0.00024876
  12 │ Hexan-1-ol           0.00462178  374.118   0.695968  31.3304    1.3144    82.0204  38.0954  0.00024876
                                                                                              2 columns omitted, [479.36199999999997, 483.93399999999997, 506.038, 495.87399999999997, 487.234, 457.25199999999995, 497.39799999999997, 478.6, 496.38399999999996, 500.19399999999996, 511.88199999999995, 446.32599999999996], 102.526260333))
=#