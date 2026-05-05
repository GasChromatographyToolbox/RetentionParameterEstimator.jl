# test the multistart optimization
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using RetentionParameterEstimator
using CSV
using DataFrames


# load data from data/meas_df05_Rxi5SilMS.csv
data_file = joinpath(@__DIR__, "..", "data", "meas_df05_Rxi5SilMS.csv")
meas_full = RetentionParameterEstimator.load_chromatograms(data_file)

# load retention parameter data from data/database_Rxi5SilMS_beta125.csv
db_file = joinpath(@__DIR__, "..", "data", "database_Rxi5SilMS_beta125.csv")
db = DataFrame(CSV.File(db_file))

# different estimates for the retention parameters for each measurement separately:
meas_select = Array{Tuple{RetentionParameterEstimator.GasChromatographySimulator.Column, Vector{RetentionParameterEstimator.GasChromatographySimulator.Program}, DataFrame, Vector{String}, String, String}}(undef, length(meas_full[2]))
for i=1:length(meas_full[2])
    meas_select[i] = RetentionParameterEstimator.filter_selected_measurements(meas_full, ["meas$(i)"], meas_full[4])
end
# 1. estimate_start_parameter_single_ramp
est_1 = Array{NTuple{4, Vector{Float64}}}(undef, length(meas_full[2]))
for i=1:length(meas_full[2])
    est_1[i] = RetentionParameterEstimator.estimate_start_parameter_single_ramp(meas_select[i][3], meas_select[i][1], meas_select[i][2]; time_unit=meas_select[i][6])
end

# 1b. corrected initial estimate for Tchar: Tchar_est = Telu ./ (0.25*sqrt(rT) + 0.8)
est_1_corrected = Array{NTuple{4, Vector{Float64}}}(undef, length(meas_full[2]))
for i=1:length(meas_full[2])
    est_1_corrected[i] = RetentionParameterEstimator.estimate_start_parameter_single_measurement_corrected(
        meas_select[i][3],  # DataFrame with retention times
        meas_select[i][1],  # Column
        meas_select[i][2][1],  # Single program (first element of program array)
        time_unit=meas_select[i][6]
    )
end



# 2. [] estimate_start_parameter_single_ramp_weighted -> remove this function
# 3. [] estimate_start_parameter_mean_elu_temp -> is this function still needed if we use approximation for complex temperature programs by using a mean heating ramp between first and last temperature plateau?
est_3 = Array{NTuple{4, Vector{Float64}}}(undef, length(meas_full[2]))
for i=1:length(meas_full[2])
    est_3[i] = RetentionParameterEstimator.estimate_start_parameter_mean_elu_ramp(meas_select[i][3], meas_select[i][1], meas_select[i][2]; time_unit=meas_select[i][6])
end
# using only one chromatogram results in the same initial estimates as `est_1`
# separate tests to compare the different functions needed for
# a) single chromatograms with complex temperature programs
# b) multiple chromatograms  
# 4. [x] use the retention parameter data from the database -> already done with `opt_1_db` (`Tchar_db`, `θchar_db`, `ΔCp_db`)

# compare the different initial estimates for the retention parameters

# use the different initial estimates for optimization
# 1. without multistart
# 1.1. method_m1
#       basicly done with `opt_1`, `opt_1_corrected` and `opt_1_db`
#
# => make the optimization for all separate measurements using:
# A. the initial estimates `est_1`
# B. the corrected initial estimates `est_1_corrected`
# C. the database values as initial estimates

# Prepare database lookup
db_name_to_row = Dict(name => row for (row, name) in enumerate(db.Name))

# Create output directory
results_dir = joinpath(@__DIR__, "..", "data", "multistart_test")
mkpath(results_dir)

# Loop over all measurements
n_meas = length(meas_full[2])
println("Running optimization for $n_meas measurements...")

opt_1 = Array{Any}(undef, n_meas)
opt_1_corrected = Array{Any}(undef, n_meas)
opt_1_db = Array{Any}(undef, n_meas)
for i=1:n_meas
    println("\n=== Processing measurement $i ===")
    
    # A. Optimize using initial estimates `est_1`
    println("  Running opt_1 (est_1 initial estimates)...")
    opt_1[i] = RetentionParameterEstimator.estimate_parameters(
        meas_select[i][3], 
        meas_select[i][4],
        meas_select[i][1], 
        meas_select[i][2], 
        est_1[i][1], 
        est_1[i][2], 
        est_1[i][3]; 
        method=RetentionParameterEstimator.NewtonTrustRegion(), 
        opt=RetentionParameterEstimator.std_opt, 
        maxiters=10000, 
        maxtime=600.0, 
        mode="Kcentric_single", 
        metric="squared", 
        pout=meas_select[i][5], 
        time_unit=meas_select[i][6], 
        parallel=false, 
        multistart_n=0, 
        coupled_perturbation=true
    )
    
    # B. Optimize using corrected initial estimates `est_1_corrected`
    println("  Running opt_1_corrected (est_1_corrected initial estimates)...")
    opt_1_corrected[i] = RetentionParameterEstimator.estimate_parameters(
        meas_select[i][3], 
        meas_select[i][4],
        meas_select[i][1], 
        meas_select[i][2], 
        est_1_corrected[i][1], 
        est_1_corrected[i][2], 
        est_1_corrected[i][3]; 
        method=RetentionParameterEstimator.NewtonTrustRegion(), 
        opt=RetentionParameterEstimator.std_opt, 
        maxiters=10000, 
        maxtime=600.0, 
        mode="Kcentric_single", 
        metric="squared", 
        pout=meas_select[i][5], 
        time_unit=meas_select[i][6], 
        parallel=false, 
        multistart_n=0, 
        coupled_perturbation=true
    )
    
    # C. Prepare database values as initial estimates
    Tchar_db = Array{Float64}(undef, length(meas_select[i][4]))
    θchar_db = Array{Float64}(undef, length(meas_select[i][4]))
    ΔCp_db = Array{Float64}(undef, length(meas_select[i][4]))
    
    for (j, subst_name) in enumerate(meas_select[i][4])
        db_idx = get(db_name_to_row, subst_name, nothing)
        if isnothing(db_idx)
            # If substance not in database, use est_1_corrected as fallback
            println("    Warning: Substance $subst_name not found in database, using est_1_corrected")
            Tchar_db[j] = est_1_corrected[i][1][j]
            θchar_db[j] = est_1_corrected[i][2][j]
            ΔCp_db[j] = est_1_corrected[i][3][j]
        else
            # Database Tchar is in °C, convert to K
            Tchar_db[j] = db.Tchar[db_idx] + 273.15
            # Database θchar is in °C
            θchar_db[j] = db.thetachar[db_idx]
            # Database ΔCp is in J/mol/K
            ΔCp_db[j] = db.DeltaCp[db_idx]
        end
    end
    
    # Optimize using database values as initial estimates
    println("  Running opt_1_db (database initial estimates)...")
    opt_1_db[i] = RetentionParameterEstimator.estimate_parameters(
        meas_select[i][3], 
        meas_select[i][4],
        meas_select[i][1], 
        meas_select[i][2], 
        Tchar_db, 
        θchar_db, 
        ΔCp_db; 
        method=RetentionParameterEstimator.NewtonTrustRegion(), 
        opt=RetentionParameterEstimator.std_opt, 
        maxiters=10000, 
        maxtime=600.0, 
        mode="Kcentric_single", 
        metric="squared", 
        pout=meas_select[i][5], 
        time_unit=meas_select[i][6], 
        parallel=false, 
        multistart_n=0, 
        coupled_perturbation=true
    )
    
    # Save results to CSV files
    # estimate_parameters returns a tuple, the first element is the DataFrame
    CSV.write(joinpath(results_dir, "opt_1_meas$(i).csv"), opt_1[i][1])
    CSV.write(joinpath(results_dir, "opt_1_corrected_meas$(i).csv"), opt_1_corrected[i][1])
    CSV.write(joinpath(results_dir, "opt_1_db_meas$(i).csv"), opt_1_db[i][1])
    
    println("  Saved results for measurement $i")
end

println("\n✓ All optimization results saved to CSV files in $(results_dir)")
println("  Files: opt_1_meas1.csv through opt_1_meas$(n_meas).csv")
println("         opt_1_corrected_meas1.csv through opt_1_corrected_meas$(n_meas).csv")
println("         opt_1_db_meas1.csv through opt_1_db_meas$(n_meas).csv")
# 1.2. method_m2
# 2. with multistart
# 2.1. method_m1
# 2.1.1. with coupled perturbation
# 2.1.2. without coupled perturbation
# 2.2. method_m2
# 2.2.1. with coupled perturbation
# 2.2.2. without coupled perturbation
# note: also record the perturbation initial estimates for the retention parameters

# compare the results

# plot the results

# save the results

# save the plots