### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ df61062c-0eaa-448c-84a2-aa0d507a5f2a
begin
    import Pkg
    # activate a temporary environment
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="Plots", version="1"),
        Pkg.PackageSpec(name="PlutoUI", version="0.7"),
		Pkg.PackageSpec(name="RetentionParameterEstimator", rev="multistart"),
		Pkg.PackageSpec(name="CSV"),
		Pkg.PackageSpec(name="DataFrames"),
		Pkg.PackageSpec(name="LsqFit")
    ])
    using Plots, PlutoUI, RetentionParameterEstimator, CSV, DataFrames, LsqFit
end


# ╔═╡ 6746f24d-5d7c-4856-be20-296a74c604b0
TableOfContents(depth=5)

# ╔═╡ 26ac8579-9ed8-40b2-b554-46e45e46a547
md"""
# Test for Multistart
"""

# ╔═╡ 523227ff-db3e-49d2-874a-bb6e8795e394
md"""
## I. Measured example data
Use of the measured example data from the RetentionParameterEstimator package.
"""

# ╔═╡ 4c01b1bc-7ded-4861-bf74-7bf2bd2d4a02
md"""
### I.1 Load the measurement data
"""

# ╔═╡ 0f574073-c05b-4a93-a831-de8edb0ed2c3
begin
    data_file = joinpath(@__DIR__, "..", "data", "meas_df05_Rxi5SilMS.csv")
    meas_full = RetentionParameterEstimator.load_chromatograms(data_file)
end

# ╔═╡ c5f3d1fd-360f-47ac-a45f-ff5c5df18f68
a = RetentionParameterEstimator.time_unit_conversion_factor(meas_full[6])

# ╔═╡ 3eb27671-8398-4b43-b4c4-4ffe64b0e268
n_meas = length(meas_full[2])

# ╔═╡ 8f927f6e-565b-4465-b918-26c93241636b
md"""
### I.2 Load the retention data 
"""

# ╔═╡ 172f6380-5ff9-4375-8982-402ee1c1312e
begin
    db_file = joinpath(@__DIR__, "..", "data", "database_Rxi5SilMS_beta125.csv")
    db = DataFrame(CSV.File(db_file))
end

# ╔═╡ 2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6c7d
md"""
### I.3 Prepare separate measurements
Only single chromatograms are used
"""

# ╔═╡ 3b4c5d6e-7f8a-9b0c-1d2e-3f4a5b6c7d8e
begin
	# Different estimates for the retention parameters for each measurement separately:
	meas_select = Array{Tuple{RetentionParameterEstimator.GasChromatographySimulator.Column, Vector{RetentionParameterEstimator.GasChromatographySimulator.Program}, DataFrame, Vector{String}, String, String}}(undef, length(meas_full[2]))
	for i=1:length(meas_full[2])
	    meas_select[i] = RetentionParameterEstimator.filter_selected_measurements(meas_full, ["meas$(i)"], meas_full[4])
	end
	meas_select
end

# ╔═╡ 4c5d6e7f-8a9b-0c1d-2e3f-4a5b6c7d8e9f
md"""
### I.4 Initial parameter estimation
"""

# ╔═╡ 5d6e7f8a-9b0c-1d2e-3f4a-5b6c7d8e9f0a
md"""
#### I.4.1 `estimate_start_parameter_single_ramp`
"""

# ╔═╡ 6e7f8a9b-0c1d-2e3f-4a5b-6c7d8e9f0a1b
begin
	# 1. estimate_start_parameter_single_ramp
	est_1 = Array{NTuple{4, Vector{Float64}}}(undef, length(meas_full[2]))
	for i=1:length(meas_full[2])
	    est_1[i] = RetentionParameterEstimator.estimate_start_parameter_single_ramp(meas_select[i][3], meas_select[i][1], meas_select[i][2]; time_unit=meas_select[i][6])
	end
	# -> Tchar_est = Telu_max
	# -> so for single ramp the initial estimate for Tchar is the elution temperature of the solutes
	est_1
end

# ╔═╡ 892d0a54-1247-421c-a24c-b9b15304f9b7
# reorder the estimations in the same order as in the database

# ╔═╡ dc692369-c139-4792-b127-6b31da0665cd
index_estimates = [findfirst(name .== meas_full[4]) for name in db.Name]

# ╔═╡ 24d3b224-0a22-4622-a8a1-11a2a9913af6
meas_full[4][index_estimates] == db.Name

# ╔═╡ 7f8a9b0c-1d2e-3f4a-5b6c-7d8e9f0a1b2c
begin
	Plots.plot(ylabel="Telu in K")
	for i=1:length(meas_full[2])
	    Plots.plot!(est_1[i][4][index_estimates], label="meas$(i)")
	end
	Plots.plot!(title="estimate_start_parameter_single_ramp")
	Plots.scatter!(db.Tchar.+273.15, label="database")
end

# ╔═╡ 3dadde84-a664-4c00-845c-1a6ed97c24a6
md"""
- `Tchar_est = Telu_max`
- so for single ramp the initial estimate for Tchar is the elution temperature of the solutes
- for lower heating rates the elution temperature is lower than for higher heating rates
- elution temperatures of the fourth measurement best matches the database values
- have a look at the holdup time and dimensionless heating rate
"""

# ╔═╡ 8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2c3d
begin
	Plots.plot(ylabel="θchar in °C")
	for i=1:length(meas_full[2])
	    Plots.plot!(est_1[i][2][index_estimates], label="meas$(i)")
	end
	Plots.plot!(title="estimate_start_parameter_single_ramp")
	Plots.scatter!(db.thetachar, label="database")
	# -> thetachar values with offset
end

# ╔═╡ 668e5789-a707-44ed-a103-c758067a1b8c
md"""
- estimates for `thetachar` underestimated for all exept the last measurement (highest heating rate -> highest elution temperatures)
"""

# ╔═╡ 9b0c1d2e-3f4a-5b6c-7d8e-9f0a1b2c3d4e
begin
	Plots.plot(ylabel="Delta C_p in J/mol/K")
	for i=1:length(meas_full[2])
	    Plots.plot!(est_1[i][3][index_estimates], label="meas$(i)")
	end
	Plots.plot!(title="estimate_start_parameter_single_ramp")
	Plots.scatter!(db.DeltaCp, label="database")
end

# ╔═╡ dc353cbd-d55a-4fc8-b8f6-26827869dcf3
md"""
- estimates for initial `ΔCp` for most measurements in similar magnitude
- `ΔCp` variies much more in general
- for highest heating rates the initial estimate is much higher
"""

# ╔═╡ 0c1d2e3f-4a5b-6c7d-8e9f-0a1b2c3d4e5f
md"""
##### I.4.1.1 Heating rate analysis

- Hold-up time at reference temperature (150°C): `tMref`
- reference thermal constant: `θref = 30`°C (as rough approximation for the characteristic thermal constant θchar)
- dimensionless heating rate: `rT = RT * tMref/θref`
- nominel dimensionless heating rate: `rT_nom = 0.69` (dimensionless heating rate where `Tchar ≈ Telu` for retained substances eluting in a heating ramp)
"""

# ╔═╡ 1d2e3f4a-5b6c-7d8e-9f0a-1b2c3d4e5f6a
begin
	tMref = Array{Float64}(undef, length(meas_full[2]))
	RT = Array{Float64}(undef, length(meas_full[2])) 
	for i=1:length(meas_full[2])
	    tMref[i] = RetentionParameterEstimator.reference_holdup_time(meas_select[i][1], meas_select[i][2][1])/a
	    # single-ramp temperature programs with ramp between time_steps 2 and 3 are assumed
	    RT[i] = (meas_select[i][2][1].temp_steps[3] - meas_select[i][2][1].temp_steps[2])/meas_select[i][2][1].time_steps[3]*a 
	end 
	rT = RT.*tMref./RetentionParameterEstimator.θref
	# -> the fourth heating rate is the closest to the nominal heating rate rT_nom
	# -> elution temperatures of this measurement also best matches the database values

end

# ╔═╡ 2b4aed30-bd17-4908-82d1-a60b6b07e3f6
rT_diff = rT .- RetentionParameterEstimator.rT_nom

# ╔═╡ c5f60d98-d39c-4217-8f80-c65f5f4a7d88
md"""
The fourth heating rate is the closest to the nominal rate `rT_norm`. Elution temperatures for this measurement closely match the database values.

**Idea**
- look at relative change of elution temperature `Telu` to database `db.Tchar`
- plot this over the dimensionless heating rate `rT`
- with the relationship `Telu/db.Tchar = f(rT)`, we can create a correction factor for the initial guess, if only one measurement is available
"""

# ╔═╡ 2e3f4a5b-6c7d-8e9f-0a1b-2c3d4e5f6a7b
begin
	RetentionParameterEstimator.Plots.plot(ylabel="Telu/db.Tchar", xlabel="rT")
	mean_ratio_Telu_over_Tchar = Array{Float64}(undef, length(meas_full[2]))
	for i=1:length(meas_full[2])
	    mean_ratio_Telu_over_Tchar[i] = mean(est_1[i][4]./(db.Tchar.+273.15))
	end
	RetentionParameterEstimator.Plots.scatter!(rT, mean_ratio_Telu_over_Tchar, label="mean")
	# => non-linear relationship, looks more like root function
	model_sqrt(x, p) = p[1] .* sqrt.(x) .+ p[2]
	rT_range = 0.0:0.01:2.7
	RetentionParameterEstimator.Plots.plot!(rT_range, model_sqrt(rT_range, [0.25, 0.8]), label="model")
end

# ╔═╡ 36cd926c-0714-4e66-85e5-d4ffeb277a43
md"""
- non-linear relationship -> square root
- **is this true in general?** -> more data needed for tests (use simulated chromatograms and other measurements)
- for now make a fit to the data
"""

# ╔═╡ 5040fee7-69da-440c-8c86-d0170e026ee0
md"""
##### LsqFit
"""

# ╔═╡ ab9bd046-f289-4561-a28f-6aae57520b06
fit = LsqFit.curve_fit(model_sqrt, rT, mean_ratio_Telu_over_Tchar, [0.25, 0.8])

# ╔═╡ 15b8c131-c54c-4268-995c-ed956b7cf7d0
fit.param

# ╔═╡ 3f4a5b6c-7d8e-9f0a-1b2c-3d4e5f6a7b8c
md"""
#### I.4.2 Corrected initial estimate for Tchar
"""

# ╔═╡ 4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8c9d
begin
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
	est_1_corrected
end

# ╔═╡ 5b6c7d8e-9f0a-1b2c-3d4e-5f6a7b8c9d0e
begin
	RetentionParameterEstimator.Plots.plot(ylabel="Tchar in K")
	for i=1:length(meas_full[2])
	    RetentionParameterEstimator.Plots.plot!(est_1_corrected[i][1][index_estimates], label="meas$(i)")
	end
	RetentionParameterEstimator.Plots.plot!(title="estimate_start_parameter_single_measurement_corrected")
	RetentionParameterEstimator.Plots.scatter!(db.Tchar.+273.15, label="database")
	# -> estimates for Tchar much closer to the database values
end

# ╔═╡ 01e4eb87-d378-4493-9024-a42c9e4fb6d3
md"""
- similar estimates for `Tchar` for all heating rates
- lower heating rates still underestimating
- but less severe
"""

# ╔═╡ 6c7d8e9f-0a1b-2c3d-4e5f-6a7b8c9d0e1f
begin
	RetentionParameterEstimator.Plots.plot(ylabel="θchar in °C")
	for i=1:length(meas_full[2])
	    RetentionParameterEstimator.Plots.plot!(est_1_corrected[i][2][index_estimates], label="meas$(i)")
	end
	RetentionParameterEstimator.Plots.plot!(title="estimate_start_parameter_single_measurement_corrected")
	RetentionParameterEstimator.Plots.scatter!(db.thetachar, label="database")
	# -> thetachar values with offset
end

# ╔═╡ 9e1acb4c-57f3-442a-8743-420734cc2a22
md"""
- initial estimates for `θchar` with correction factor similar for all heating rates
- still relative big offset to database values
"""

# ╔═╡ 7d8e9f0a-1b2c-3d4e-5f6a-7b8c9d0e1f2a
begin
	RetentionParameterEstimator.Plots.plot(ylabel="Delta C_p in J/mol/K")
	for i=1:length(meas_full[2])
	    RetentionParameterEstimator.Plots.plot!(est_1_corrected[i][3][index_estimates], label="meas$(i)")
	end
	RetentionParameterEstimator.Plots.plot!(title="estimate_start_parameter_single_measurement_corrected")
	RetentionParameterEstimator.Plots.scatter!(db.DeltaCp, label="database")
end

# ╔═╡ 95427b94-976c-4af5-9446-40d18630342f
md"""
- initial estimate for `ΔCp` with correction factor gives similar values for all heating rates
- database values scatter with relativ big deviations around these estimates
"""

# ╔═╡ 2390a8b7-5dd3-49c7-8d49-89f7b4e0200e
md"""
#### I.4.3 `estimate_start_parameter_mean_elu_temp`
"""

# ╔═╡ a8fa25a7-f03a-43af-8214-9f45aa6af41e
begin
	est_3 = Array{NTuple{4, Vector{Float64}}}(undef, length(meas_full[2]))
	for i=1:length(meas_full[2])
	    est_3[i] = RetentionParameterEstimator.estimate_start_parameter_mean_elu_temp(meas_select[i][3], meas_select[i][1], meas_select[i][2]; time_unit=meas_select[i][6])
	end
	est_3
end

# ╔═╡ 4959b43d-8727-4ed6-8cab-8a5f776d5506
begin
	RetentionParameterEstimator.Plots.plot(ylabel="Tchar in K")
	for i=1:length(meas_full[2])
	    RetentionParameterEstimator.Plots.plot!(est_3[i][1][index_estimates], label="meas$(i)")
	end
	RetentionParameterEstimator.Plots.plot!(title="estimate_start_parameter_mean_elu_temp")
	RetentionParameterEstimator.Plots.scatter!(db.Tchar.+273.15, label="database")
	# -> estimates for Tchar much closer to the database values
end

# ╔═╡ 97ecbc60-95b8-4eec-93d2-f3fc1fa42c10
est_1 == est_3

# ╔═╡ a941154e-91a6-4f0b-9e82-a9076de068b4
md"""
Using only one chromatogram to estimate the initial retention parameters results in the same parameters using functions `estimate_start_parameter_single_ramp` and `estimate_start_parameter_mean_elu_temp`.

To compare these functions, chromatograms with complex programs should be used, or the combination of multiple chromatograms.
"""

# ╔═╡ 454b19ba-e9ed-4a0f-bc28-5fe1da017911
md"""
#### 1.4.4 Using all measurements for estimates
"""

# ╔═╡ a8027040-afce-499a-8d2c-03d87a651861
# 1. estimate_start_parameter_single_ramp using all measurements
est_1_full =  RetentionParameterEstimator.estimate_start_parameter_single_ramp(meas_full[3], meas_full[1], meas_full[2]; time_unit=meas_full[6])

# ╔═╡ c752a216-ce4e-4391-8bf7-3e401d4f5541
begin
	RetentionParameterEstimator.Plots.plot(ylabel="Tchar in K")
	RetentionParameterEstimator.Plots.plot!(est_1_full[1][index_estimates], label="all measurements")
	RetentionParameterEstimator.Plots.plot!(title="estimate_start_parameter_single_ramp")
	RetentionParameterEstimator.Plots.scatter!(db.Tchar.+273.15, label="database")
	# -> estimates for Tchar much closer to the database values
end

# ╔═╡ bb3196f9-37fa-4a24-ace3-88d8bc8c4f81
est_1_full_corrected = RetentionParameterEstimator.estimate_start_parameter_single_measurement_corrected(
	        meas_full[3],  # DataFrame with retention times
	        meas_full[1],  # Column
	        meas_full[2],  # Single program (first element of program array)
	        time_unit=meas_full[6]
	    )
# gives error
# - either only allow one measurement
# - or fix the problem

# ╔═╡ 7e7eef6e-5f5f-4924-85f8-11b60da63aff
est_3_full = RetentionParameterEstimator.estimate_start_parameter_mean_elu_temp(meas_full[3], meas_full[1], meas_full[2]; time_unit=meas_full[6])

# ╔═╡ a7a779a1-5273-4576-8e32-335c23a7eeed
est_3_full == est_1_full

# ╔═╡ bb4e78db-9835-4586-b62b-4c1b9a089fad
begin
	RetentionParameterEstimator.Plots.plot(ylabel="Tchar in K")
	RetentionParameterEstimator.Plots.plot!(est_3_full[1][index_estimates], label="all measurements")
	RetentionParameterEstimator.Plots.plot!(title="eestimate_start_parameter_mean_elu_temp")
	RetentionParameterEstimator.Plots.scatter!(db.Tchar.+273.15, label="database")
	# -> estimates for Tchar much closer to the database values
end

# ╔═╡ 0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4c5d
md"""
### []I.5 Load optimization results from CSV files
"""

# ╔═╡ 00ff8636-b2cf-4a3b-9e85-f2c3d904b242
results_dir = joinpath(@__DIR__, "..", "data", "multistart_test")

# ╔═╡ 2c3d4e5f-6a7b-8c9d-0e1f-2a3b4c5d6e7f
md"""
#### I.5.1 Comparison of initial estimates and optimization results
"""

# ╔═╡ c36c7e5e-c524-479f-8da8-625b81a4d566
md"""
##### I.5.1.1 Optimization results using `est_1`
"""

# ╔═╡ 98359aa4-6d1c-423c-b988-b3de09f42da7
begin
	opt_1 = Array{Any}(undef, n_meas)
	for i=1:n_meas
		opt_1[i] = try
			DataFrame(CSV.File(joinpath(results_dir, "opt_1_meas$(i).csv")))[index_estimates,:]
		catch
			nothing
		end
	end
	opt_1
end

# ╔═╡ 9c2be544-4723-4a78-b927-0a581f0de306
md"""
##### I.5.1.2 Optimization results using `est_1_corrected`
"""

# ╔═╡ ab5fbbaf-c193-4fc0-9e39-fce5a2267bc6
begin
	opt_1_corrected = Array{Any}(undef, n_meas)
	for i=1:n_meas
		opt_1_corrected[i] = try
			DataFrame(CSV.File(joinpath(results_dir, "opt_1_corrected_meas$(i).csv")))[index_estimates,:]
		catch
			nothing
		end
	end
	opt_1_corrected
end

# ╔═╡ 7ae3c0b8-b24e-4466-9a50-24d0cb4805bc
md"""
##### I.5.1.3 Optimization results using database `db`
"""

# ╔═╡ 96086869-fdd0-49f8-9761-5d29b7f43a6b
begin
	opt_1_db = Array{Any}(undef, n_meas)
	for i=1:n_meas
		opt_1_db[i] = try
			DataFrame(CSV.File(joinpath(results_dir, "opt_1_corrected_meas$(i).csv")))[index_estimates,:]
		catch
			nothing
		end
	end
	opt_1_db
end

# ╔═╡ 252290a2-f4c5-4d4b-9927-7af2db1f586a
md"""
### I.6 Comparison
"""

# ╔═╡ 4e5f6a7b-8c9d-0e1f-2a3b-4c5d6e7f8a9b
md"""
#### []I.6.1 `Tchar` comparison
"""

# ╔═╡ 5f6a7b8c-9d0e-1f2a-3b4c-5d6e7f8a9b0c
begin
	RetentionParameterEstimator.Plots.scatter(est_1[1][1][opt_1_indices], c=:blue, markershape=:diamond, label="est_1")
	RetentionParameterEstimator.Plots.scatter!(est_1_corrected[1][1][opt_1_corrected_indices], c=:red, markershape=:diamond, label="est_1_corrected")
	if !isnothing(opt_1_final)
		RetentionParameterEstimator.Plots.scatter!(opt_1_final.Tchar, c=:blue, label="opt_1")
	end
	if !isnothing(opt_1_corrected_final)
		RetentionParameterEstimator.Plots.scatter!(opt_1_corrected_final.Tchar, c=:red, label="opt_1_corrected")
	end
	if !isnothing(opt_1_db_final)
		RetentionParameterEstimator.Plots.scatter!(opt_1_db_final.Tchar, c=:green, label="opt_db")
	end
	RetentionParameterEstimator.Plots.scatter!(db_final.Tchar.+273.15, c=:black, markershape=:x, label="database")
	RetentionParameterEstimator.Plots.plot!(title="Tchar", ylabel="Tchar in K")
end

# ╔═╡ 6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0c1d
md"""
#### []I.6.2 `θchar` comparison
"""

# ╔═╡ 7b8c9d0e-1f2a-3b4c-5d6e-7f8a9b0c1d2e
begin
	RetentionParameterEstimator.Plots.scatter(est_1[1][2][opt_1_indices], c=:blue, markershape=:diamond, label="est_1")
	RetentionParameterEstimator.Plots.scatter!(est_1_corrected[1][2][opt_1_corrected_indices], c=:red, markershape=:diamond, label="est_1_corrected")
	if !isnothing(opt_1_final)
		RetentionParameterEstimator.Plots.scatter!(opt_1_final.θchar, c=:blue, label="opt_1")
	end
	if !isnothing(opt_1_corrected_final)
		RetentionParameterEstimator.Plots.scatter!(opt_1_corrected_final.θchar, c=:red, label="opt_1_corrected")
	end
	if !isnothing(opt_1_db_final)
		RetentionParameterEstimator.Plots.scatter!(opt_1_db_final.θchar, c=:green, label="opt_db")
	end
	RetentionParameterEstimator.Plots.scatter!(db_final.thetachar, c=:black, markershape=:x, label="database")
	RetentionParameterEstimator.Plots.plot!(title="θchar", ylabel="θchar in °C")
end

# ╔═╡ 8c9d0e1f-2a3b-4c5d-6e7f-8a9b0c1d2e3f
md"""
#### []I.6.3 `ΔCp` comparison
"""

# ╔═╡ 9d0e1f2a-3b4c-5d6e-7f8a-9b0c1d2e3f4a
begin
	RetentionParameterEstimator.Plots.scatter(est_1[1][3][opt_1_indices], c=:blue, markershape=:diamond, label="est_1")
	RetentionParameterEstimator.Plots.scatter!(est_1_corrected[1][3][opt_1_corrected_indices], c=:red, markershape=:diamond, label="est_1_corrected")
	if !isnothing(opt_1_final)
		RetentionParameterEstimator.Plots.scatter!(opt_1_final.ΔCp, c=:blue, label="opt_1")
	end
	if !isnothing(opt_1_corrected_final)
		RetentionParameterEstimator.Plots.scatter!(opt_1_corrected_final.ΔCp, c=:red, label="opt_1_corrected")
	end
	if !isnothing(opt_1_db_final)
		RetentionParameterEstimator.Plots.scatter!(opt_1_db_final.ΔCp, c=:green, label="opt_db")
	end
	RetentionParameterEstimator.Plots.scatter!(db_final.DeltaCp, c=:black, markershape=:x, label="database")
	RetentionParameterEstimator.Plots.plot!(title="ΔCp", ylabel="ΔCp in J/mol/K")
end

# ╔═╡ 8e9f0a1b-2c3d-4e5f-6a7b-8c9d0e1f2a3b
md"""
#### []I.6.4 Loss calculations
"""

# ╔═╡ 9f0a1b2c-3d4e-5f6a-7b8c-9d0e1f2a3b4c
begin
	# Calculate loss for each substance using initial estimates
	tR_meas = Array(meas_select[1][3][:,2:end]).*a
	ns = size(tR_meas, 2)  # number of substances

	loss_est_1 = Array{Float64}(undef, ns)
	loss_est_1_corrected = Array{Float64}(undef, ns)

	for j=1:ns
	    # Prepare data for this substance (filter missing values)
	    tRs_, prog_, subst_list_ = RetentionParameterEstimator.prepare_single_substance_data(
	        tR_meas[:,j], 
	        meas_select[1][2], 
	        meas_select[1][4][j]
	    )
	    
	    # Calculate loss with est_1 initial estimates
	    Tchar_j = est_1[1][1][j]
	    θchar_j = est_1[1][2][j]
	    ΔCp_j = est_1[1][3][j]
	    
	    loss_est_1[j] = RetentionParameterEstimator.loss(
	        tRs_, 
	        [Tchar_j],  # Single value as vector
	        [θchar_j], 
	        [ΔCp_j], 
	        subst_list_, 
	        meas_select[1][1].L, 
	        meas_select[1][1].d, 
	        prog_, 
	        meas_select[1][1].gas; 
	        opt=RetentionParameterEstimator.std_opt, 
	        metric="squared"
	    )
	    
	    # Calculate loss with est_1_corrected initial estimates
	    Tchar_j_corrected = est_1_corrected[1][1][j]
	    θchar_j_corrected = est_1_corrected[1][2][j]
	    ΔCp_j_corrected = est_1_corrected[1][3][j]
	    
	    loss_est_1_corrected[j] = RetentionParameterEstimator.loss(
	        tRs_, 
	        [Tchar_j_corrected], 
	        [θchar_j_corrected], 
	        [ΔCp_j_corrected], 
	        subst_list_, 
	        meas_select[1][1].L, 
	        meas_select[1][1].d, 
	        prog_, 
	        meas_select[1][1].gas; 
	        opt=RetentionParameterEstimator.std_opt, 
	        metric="squared"
	    )
	end

	# Create DataFrames with loss values for comparison
	loss_df_initial = DataFrame(
	    Name=meas_select[1][4],
	    loss_est_1=loss_est_1,
	    loss_est_1_corrected=loss_est_1_corrected
	)
	loss_df_initial
end

# ╔═╡ 0e1f2a3b-4c5d-6e7f-8a9b-0c1d2e3f4a5b
md"""
##### []I.6.4.1 Loss comparison
"""

# ╔═╡ 1f2a3b4c-5d6e-7f8a-9b0c-1d2e3f4a5b6c
begin
	if !isnothing(loss_df_full)
		# Reorder loss_df to match db_final order
		loss_df_name_to_row = Dict(name => row for (row, name) in enumerate(loss_df_full.Name))
		loss_df_indices = [get(loss_df_name_to_row, name, nothing) for name in db_final.Name]
		if any(isnothing, loss_df_indices)
			println("Warning: Some substances in db_final not found in loss_df")
			valid = findall(!isnothing, loss_df_indices)
			loss_df_indices = loss_df_indices[valid]
			db_final_loss = db_final[valid, :]
		else
			db_final_loss = db_final
		end
		loss_df_final = loss_df_full[loss_df_indices, :]
		
		println("\nLoss comparison (using db order):")
		println(loss_df_final)
		
		# Plot loss comparison
		RetentionParameterEstimator.Plots.scatter(loss_df_final.loss_est_1, c=:blue, markershape=:diamond, label="est_1")
		RetentionParameterEstimator.Plots.scatter!(loss_df_final.loss_est_1_corrected, c=:red, markershape=:diamond, label="est_1_corrected")
		if hasproperty(loss_df_final, :loss_est_db)
			RetentionParameterEstimator.Plots.scatter!(loss_df_final.loss_est_db, c=:green, markershape=:diamond, label="est_db")
		end
		if hasproperty(loss_df_final, :loss_opt_1)
			RetentionParameterEstimator.Plots.scatter!(loss_df_final.loss_opt_1, c=:blue, label="opt_1")
		end
		if hasproperty(loss_df_final, :loss_opt_1_corrected)
			RetentionParameterEstimator.Plots.scatter!(loss_df_final.loss_opt_1_corrected, c=:red, label="opt_1_corrected")
		end
		if hasproperty(loss_df_final, :loss_opt_db)
			RetentionParameterEstimator.Plots.scatter!(loss_df_final.loss_opt_db, c=:green, label="opt_db")
		end
		RetentionParameterEstimator.Plots.plot!(title="Loss comparison", ylabel="Loss", yscale=:log10)
	else
		# Plot only initial estimates if optimization results not loaded
		RetentionParameterEstimator.Plots.scatter(loss_df_initial.loss_est_1, c=:blue, markershape=:diamond, label="est_1")
		RetentionParameterEstimator.Plots.scatter!(loss_df_initial.loss_est_1_corrected, c=:red, markershape=:diamond, label="est_1_corrected")
		RetentionParameterEstimator.Plots.plot!(title="Loss comparison (initial estimates only)", ylabel="Loss", yscale=:log10)
	end
end

# ╔═╡ 8ba57a4c-c300-40b6-bea5-b8a950fd9d63
md"""
## []II. Simulated data
"""

# ╔═╡ Cell order:
# ╠═df61062c-0eaa-448c-84a2-aa0d507a5f2a
# ╠═6746f24d-5d7c-4856-be20-296a74c604b0
# ╠═26ac8579-9ed8-40b2-b554-46e45e46a547
# ╠═523227ff-db3e-49d2-874a-bb6e8795e394
# ╠═4c01b1bc-7ded-4861-bf74-7bf2bd2d4a02
# ╠═0f574073-c05b-4a93-a831-de8edb0ed2c3
# ╠═c5f3d1fd-360f-47ac-a45f-ff5c5df18f68
# ╠═3eb27671-8398-4b43-b4c4-4ffe64b0e268
# ╠═8f927f6e-565b-4465-b918-26c93241636b
# ╠═172f6380-5ff9-4375-8982-402ee1c1312e
# ╠═2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6c7d
# ╠═3b4c5d6e-7f8a-9b0c-1d2e-3f4a5b6c7d8e
# ╠═4c5d6e7f-8a9b-0c1d-2e3f-4a5b6c7d8e9f
# ╠═5d6e7f8a-9b0c-1d2e-3f4a-5b6c7d8e9f0a
# ╠═6e7f8a9b-0c1d-2e3f-4a5b-6c7d8e9f0a1b
# ╠═892d0a54-1247-421c-a24c-b9b15304f9b7
# ╠═dc692369-c139-4792-b127-6b31da0665cd
# ╠═24d3b224-0a22-4622-a8a1-11a2a9913af6
# ╟─7f8a9b0c-1d2e-3f4a-5b6c-7d8e9f0a1b2c
# ╟─3dadde84-a664-4c00-845c-1a6ed97c24a6
# ╟─8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2c3d
# ╟─668e5789-a707-44ed-a103-c758067a1b8c
# ╟─9b0c1d2e-3f4a-5b6c-7d8e-9f0a1b2c3d4e
# ╟─dc353cbd-d55a-4fc8-b8f6-26827869dcf3
# ╟─0c1d2e3f-4a5b-6c7d-8e9f-0a1b2c3d4e5f
# ╠═1d2e3f4a-5b6c-7d8e-9f0a-1b2c3d4e5f6a
# ╠═2b4aed30-bd17-4908-82d1-a60b6b07e3f6
# ╠═c5f60d98-d39c-4217-8f80-c65f5f4a7d88
# ╠═2e3f4a5b-6c7d-8e9f-0a1b-2c3d4e5f6a7b
# ╠═36cd926c-0714-4e66-85e5-d4ffeb277a43
# ╠═5040fee7-69da-440c-8c86-d0170e026ee0
# ╠═ab9bd046-f289-4561-a28f-6aae57520b06
# ╠═15b8c131-c54c-4268-995c-ed956b7cf7d0
# ╠═3f4a5b6c-7d8e-9f0a-1b2c-3d4e5f6a7b8c
# ╠═4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8c9d
# ╟─5b6c7d8e-9f0a-1b2c-3d4e-5f6a7b8c9d0e
# ╠═01e4eb87-d378-4493-9024-a42c9e4fb6d3
# ╟─6c7d8e9f-0a1b-2c3d-4e5f-6a7b8c9d0e1f
# ╟─9e1acb4c-57f3-442a-8743-420734cc2a22
# ╟─7d8e9f0a-1b2c-3d4e-5f6a-7b8c9d0e1f2a
# ╟─95427b94-976c-4af5-9446-40d18630342f
# ╠═2390a8b7-5dd3-49c7-8d49-89f7b4e0200e
# ╠═a8fa25a7-f03a-43af-8214-9f45aa6af41e
# ╠═4959b43d-8727-4ed6-8cab-8a5f776d5506
# ╠═97ecbc60-95b8-4eec-93d2-f3fc1fa42c10
# ╠═a941154e-91a6-4f0b-9e82-a9076de068b4
# ╠═454b19ba-e9ed-4a0f-bc28-5fe1da017911
# ╠═a8027040-afce-499a-8d2c-03d87a651861
# ╠═c752a216-ce4e-4391-8bf7-3e401d4f5541
# ╠═bb3196f9-37fa-4a24-ace3-88d8bc8c4f81
# ╠═7e7eef6e-5f5f-4924-85f8-11b60da63aff
# ╠═a7a779a1-5273-4576-8e32-335c23a7eeed
# ╠═bb4e78db-9835-4586-b62b-4c1b9a089fad
# ╠═0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4c5d
# ╠═00ff8636-b2cf-4a3b-9e85-f2c3d904b242
# ╠═2c3d4e5f-6a7b-8c9d-0e1f-2a3b4c5d6e7f
# ╠═c36c7e5e-c524-479f-8da8-625b81a4d566
# ╠═98359aa4-6d1c-423c-b988-b3de09f42da7
# ╠═9c2be544-4723-4a78-b927-0a581f0de306
# ╠═ab5fbbaf-c193-4fc0-9e39-fce5a2267bc6
# ╠═7ae3c0b8-b24e-4466-9a50-24d0cb4805bc
# ╠═96086869-fdd0-49f8-9761-5d29b7f43a6b
# ╠═252290a2-f4c5-4d4b-9927-7af2db1f586a
# ╠═4e5f6a7b-8c9d-0e1f-2a3b-4c5d6e7f8a9b
# ╠═5f6a7b8c-9d0e-1f2a-3b4c-5d6e7f8a9b0c
# ╠═6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0c1d
# ╠═7b8c9d0e-1f2a-3b4c-5d6e-7f8a9b0c1d2e
# ╠═8c9d0e1f-2a3b-4c5d-6e7f-8a9b0c1d2e3f
# ╠═9d0e1f2a-3b4c-5d6e-7f8a-9b0c1d2e3f4a
# ╠═8e9f0a1b-2c3d-4e5f-6a7b-8c9d0e1f2a3b
# ╠═9f0a1b2c-3d4e-5f6a-7b8c-9d0e1f2a3b4c
# ╠═0e1f2a3b-4c5d-6e7f-8a9b-0c1d2e3f4a5b
# ╠═1f2a3b4c-5d6e-7f8a-9b0c-1d2e3f4a5b6c
# ╠═8ba57a4c-c300-40b6-bea5-b8a950fd9d63
