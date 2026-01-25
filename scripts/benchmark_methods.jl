#!/usr/bin/env julia
"""
    benchmark_methods.jl

Script to benchmark and compare all optimization methods (method_m1, m2, m3, m4, m4_) performance and results.
Supports optional parallelization controlled via command-line argument.

Usage:
    julia benchmark_methods.jl [data_file] [--parallel|--no-parallel] [maxiters] [maxtime]
    julia -t 4 benchmark_methods.jl [data_file] --parallel  # With 4 threads for parallelization
    julia benchmark_methods.jl [data_file] --no-parallel    # Force no parallelization

Arguments:
    data_file: Path to CSV file with measurement data (default: ../data/meas_df05_Rxi5SilMS.csv)
    --parallel: Enable parallelization (default if multiple threads available)
    --no-parallel: Disable parallelization (for debugging)
    maxiters: Maximum iterations for optimization (default: 10000)
    maxtime: Maximum time in seconds for optimization (default: 600.0)

If no parallelization flag is provided, parallelization is enabled if multiple threads are available.
"""

# Activate the package environment (go up one level from scripts/ to package root)
using Pkg
Pkg.activate(dirname(@__DIR__))

using RetentionParameterEstimator
using Printf
using CSV
using DataFrames
using Dates
using LibGit2

"""
    get_git_info()

Get git repository information for tracking code state.
"""
function get_git_info()
    repo_path = dirname(@__DIR__)
    try
        repo = LibGit2.GitRepo(repo_path)
        head = LibGit2.head(repo)
        commit_oid = LibGit2.GitHash(head)
        commit_hash = string(commit_oid)[1:7]  # Short hash
        
        # Get branch name
        branch_ref = LibGit2.head(repo)
        branch_name = try
            LibGit2.shortname(branch_ref)
        catch
            "detached"
        end
        
        # Check for uncommitted changes
        has_changes = false
        try
            status_dict = LibGit2.status(repo, ".")
            has_changes = status_dict !== nothing && !isempty(status_dict)
        catch
            has_changes = false
        end
        
        # Get commit message
        commit_obj = LibGit2.get(LibGit2.GitCommit, repo, commit_oid)
        commit_msg = LibGit2.message(commit_obj)
        commit_msg_short = length(commit_msg) > 50 ? commit_msg[1:47] * "..." : commit_msg
        
        LibGit2.close(repo)
        
        return Dict(
            "commit_hash" => commit_hash,
            "commit_hash_full" => string(commit_oid),
            "branch" => branch_name,
            "has_uncommitted_changes" => has_changes,
            "commit_message" => strip(commit_msg_short)
        )
    catch e
        return Dict(
            "commit_hash" => "unknown",
            "commit_hash_full" => "unknown",
            "branch" => "unknown",
            "has_uncommitted_changes" => false,
            "commit_message" => "Could not read git info: $e"
        )
    end
end

function benchmark_and_compare(data_file::String; use_parallel::Union{Bool,Nothing}=nothing, maxiters::Int=10000, maxtime::Float64=600.0)
    println("=" ^ 80)
    println("Benchmarking All Methods: method_m1, m2, m3, m4, m4_")
    println("=" ^ 80)
    println("\nData file: $data_file")
    
    # Determine parallelization setting
    nthreads = Threads.nthreads()
    if use_parallel === nothing
        # Auto-detect: use parallelization if multiple threads available
        if nthreads > 1
            parallel_flag = true
            println("\nParallelization: AUTO-ENABLED ($nthreads threads available)")
        else
            parallel_flag = false
            println("\nParallelization: AUTO-DISABLED (only 1 thread available)")
            println("  Start Julia with multiple threads: julia -t 4 benchmark_methods.jl")
        end
    elseif use_parallel
        if nthreads > 1
            parallel_flag = true
            println("\nParallelization: ENABLED ($nthreads threads)")
        else
            parallel_flag = false
            println("\nParallelization: REQUESTED but only 1 thread available")
            println("  Start Julia with multiple threads: julia -t 4 benchmark_methods.jl")
            println("  Continuing without parallelization...")
        end
    else
        parallel_flag = false
        println("\nParallelization: DISABLED (requested)")
        println("  Threads available: $nthreads (not used)")
    end
    
    # Load data
    println("\nLoading data...")
    meas = RetentionParameterEstimator.load_chromatograms(data_file)
    println("  Loaded $(length(meas[4])) substances")
    println("  Loaded $(length(meas[3].measurement)) measurements")
    
    # Prepare col_input for method_m1 (needs L and d in mm)
    col_input = (L=meas[1].L, d=meas[1].d*1000)  # Convert d from m to mm
    
    # Warm-up run to trigger JIT compilation
    println("\n" * "-" ^ 80)
    println("Warming up (compiling code)...")
    println("-" ^ 80)
    try
        # Quick warm-up with minimal iterations/time to compile all code paths
        RetentionParameterEstimator.method_m1(meas, col_input, se_col=false, parallel=parallel_flag, maxiters=10, maxtime=1.0)
        RetentionParameterEstimator.method_m2(meas, se_col=false, parallel=parallel_flag, maxiters=10, maxtime=1.0)
        RetentionParameterEstimator.method_m3(meas, se_col=false, parallel=parallel_flag, maxiters=10, maxtime=1.0)
        RetentionParameterEstimator.method_m4(meas, se_col=false, parallel=parallel_flag, maxiters=10, maxtime=1.0)
        RetentionParameterEstimator.method_m4_(meas, se_col=false, parallel=parallel_flag, maxiters=10, maxtime=1.0)
        println("  Warm-up completed")
    catch e
        println("  Warning: Warm-up failed (this is usually fine): $e")
    end
    println()
    
    # Run method_m1
    println("\n" * "-" ^ 80)
    println("Running method_m1 (requires known d)...")
    if parallel_flag
        println("  Using parallelization: $(parallel_flag)")
    end
    println("-" ^ 80)
    time_m1 = @elapsed begin
        res_m1, Telu_max_m1 = RetentionParameterEstimator.method_m1(meas, col_input, se_col=true, parallel=parallel_flag, maxiters=maxiters, maxtime=maxtime)
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
        res_m2, Telu_max_m2 = RetentionParameterEstimator.method_m2(meas, se_col=true, parallel=parallel_flag, maxiters=maxiters, maxtime=maxtime)
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
        res_m3, Telu_max_m3 = RetentionParameterEstimator.method_m3(meas, se_col=true, parallel=parallel_flag, maxiters=maxiters, maxtime=maxtime)
    end
    
    println("  Completed in $(@sprintf("%.2f", time_m3)) seconds")
    println("  Estimated parameters for $(length(res_m3.Name)) substances")
    
    # Run method_m4
    println("\n" * "-" ^ 80)
    println("Running method_m4 (alternating optimization)...")
    if parallel_flag
        println("  Using parallelization: $(parallel_flag) (for Block 2 and standard error calculation)")
        println("  ⚠️  WARNING: method_m4 may have issues with parallelization - results may differ")
    end
    println("-" ^ 80)
    time_m4 = @elapsed begin
        res_m4, Telu_max_m4 = RetentionParameterEstimator.method_m4(meas, se_col=true, parallel=parallel_flag, maxiters=maxiters, maxtime=maxtime)
    end
    
    println("  Completed in $(@sprintf("%.2f", time_m4)) seconds")
    println("  Estimated parameters for $(length(res_m4.Name)) substances")
    
    # Run method_m4_ (previous implementation)
    println("\n" * "-" ^ 80)
    println("Running method_m4_ (alternating optimization - previous implementation)...")
    if parallel_flag
        println("  Using parallelization: $(parallel_flag) (for Block 2 and standard error calculation)")
        println("  ⚠️  WARNING: method_m4_ may have issues with parallelization - results may differ")
    end
    println("-" ^ 80)
    time_m4_ = @elapsed begin
        res_m4_, Telu_max_m4_ = RetentionParameterEstimator.method_m4_(meas, se_col=true, parallel=parallel_flag, maxiters=maxiters, maxtime=maxtime)
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
    
    # Get git info for tracking
    git_info = get_git_info()
    
    # Create timestamp for unique filename
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    data_basename = basename(replace(data_file, r"\.csv$" => ""))
    parallel_suffix = parallel_flag ? "parallel" : "no_parallel"
    
    # Save full comparison to CSV with timestamp
    comparison_file = joinpath(dirname(data_file), "$(data_basename)_comparison_$(parallel_suffix)_$(timestamp).csv")
    try
        CSV.write(comparison_file, comparison_df)
        println("\nFull comparison saved to: $comparison_file")
    catch e
        println("\nWarning: Could not save comparison file: $e")
    end
    
    # Save benchmark summary to text file
    output_file = joinpath(dirname(data_file), "$(data_basename)_benchmark_output_$(parallel_suffix)_$(timestamp).txt")
    try
        open(output_file, "w") do f
            println(f, "=" ^ 80)
            println(f, "Benchmark Summary - $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
            println(f, "=" ^ 80)
            println(f, "\nGit Repository Information:")
            println(f, "  Commit Hash: $(git_info["commit_hash"]) ($(git_info["commit_hash_full"]))")
            println(f, "  Branch: $(git_info["branch"])")
            println(f, "  Uncommitted Changes: $(git_info["has_uncommitted_changes"] ? "YES" : "NO")")
            println(f, "  Commit Message: $(git_info["commit_message"])")
            println(f, "\n" * "=" ^ 80)
            println(f, "\nBenchmark Configuration:")
            println(f, "  Data file: $data_file")
            println(f, "  Parallelization: $(parallel_flag ? "ENABLED ($nthreads threads)" : "DISABLED")")
            println(f, "  Optimization limits: maxiters=$maxiters, maxtime=$maxtime")
            println(f, "  Substances: $(length(meas[4]))")
            println(f, "  Measurements: $(length(meas[3].measurement))")
            println(f, "\n" * "=" ^ 80)
            println(f, "\nExecution Times:")
            sorted_indices = sortperm(times)
            for (rank, idx) in enumerate(sorted_indices)
                speedup = times[idx] / times[sorted_indices[1]]
                if rank == 1
                    println(f, "  $rank. $(methods[idx]): $(@sprintf("%.2f", times[idx])) seconds (fastest)")
                else
                    println(f, "  $rank. $(methods[idx]): $(@sprintf("%.2f", times[idx])) seconds ($(@sprintf("%.2fx", speedup)) slower)")
                end
            end
            println(f, "\n" * "=" ^ 80)
            println(f, "\nKey Results:")
            println(f, "  Column diameter (d):")
            for (method_name, res, _, has_d) in methods_to_compare
                if has_d && "d" in names(res)
                    d_unique = length(unique(res.d))
                    if d_unique == 1
                        println(f, "    method_$method_name: $(@sprintf("%.6f", res.d[1])) m (consistent)")
                    else
                        d_mean = mean(res.d)
                        println(f, "    method_$method_name: $(@sprintf("%.6f", d_mean)) m (varies, std=$(@sprintf("%.2e", std(res.d))) m)")
                    end
                end
            end
            println(f, "\n  Maximum Elution Temperature (Telu_max):")
            for (method_name, _, Telu_max, _) in methods_to_compare
                Telu_max_val = isa(Telu_max, Vector) ? maximum(Telu_max) : Telu_max
                println(f, "    method_$method_name: $(@sprintf("%.2f", Telu_max_val)) K")
            end
            println(f, "\n" * "=" ^ 80)
            println(f, "\nFull comparison CSV: $comparison_file")
        end
        println("\nBenchmark summary saved to: $output_file")
    catch e
        println("\nWarning: Could not save summary file: $e")
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
        println("    ⚠️  Note: method_m4 and method_m4_ may have issues with parallelization")
    else
        if length(meas[4]) > 5 && nthreads > 1
            println("\n  💡 Tip: Start Julia with multiple threads for better performance:")
            println("    julia -t 4 benchmark_methods.jl --parallel")
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
            m4_=(res_m4_, Telu_max_m4_, time_m4_)), parallel_flag
end

# Main execution
# When run directly from command line
if abspath(PROGRAM_FILE) == @__FILE__
    # Parse command-line arguments
    local data_file = nothing
    local use_parallel = nothing  # nil means auto-detect
    local maxiters = 10000
    local maxtime = 600.0
    
    for arg in ARGS
        if arg == "--parallel"
            use_parallel = true
        elseif arg == "--no-parallel"
            use_parallel = false
        elseif startswith(arg, "--")
            println("Warning: Unknown option: $arg")
        elseif data_file === nothing
            data_file = arg
        elseif maxiters == 10000  # First numeric arg after data_file
            try
                maxiters = parse(Int, arg)
            catch
                println("Warning: Could not parse maxiters: $arg, using default: $maxiters")
            end
        else  # Second numeric arg
            try
                maxtime = parse(Float64, arg)
            catch
                println("Warning: Could not parse maxtime: $arg, using default: $maxtime")
            end
        end
    end
    
    # Set default data file if not provided
    if data_file === nothing
        data_file = joinpath(dirname(@__DIR__), "data", "meas_df05_Rxi5SilMS.csv")
    end
    
    if !isfile(data_file)
        println("Error: File not found: $data_file")
        exit(1)
    end
    
    try
        # Create output file path
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        data_basename = basename(replace(data_file, r"\.csv$" => ""))
        
        # Determine parallelization for file naming
        nthreads = Threads.nthreads()
        if use_parallel === nothing
            parallel_flag_for_file = nthreads > 1
        else
            parallel_flag_for_file = use_parallel && nthreads > 1
        end
        parallel_suffix = parallel_flag_for_file ? "parallel" : "no_parallel"
        output_file = joinpath(dirname(data_file), "$(data_basename)_benchmark_output_$(parallel_suffix)_$(timestamp).txt")
        
        # Get git info
        git_info = get_git_info()
        
        # Run benchmark (output goes to console)
        benchmark_and_compare(data_file, use_parallel=use_parallel, maxiters=maxiters, maxtime=maxtime)
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
