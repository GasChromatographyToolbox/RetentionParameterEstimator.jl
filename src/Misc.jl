# additional misc. functions 
function compare!(df, db)
	ΔTchar = Array{Float64}(undef, length(df.Name))
	Δθchar = Array{Float64}(undef, length(df.Name))
	ΔΔCp = Array{Float64}(undef, length(df.Name))
	relΔTchar = Array{Float64}(undef, length(df.Name))
	relΔθchar = Array{Float64}(undef, length(df.Name))
	relΔΔCp = Array{Float64}(undef, length(df.Name))
	for i=1:length(df.Name)
		ii = findfirst(df.Name[i].==db.Name)
        if isnothing(ii)
            ΔTchar[i] = NaN
            Δθchar[i] = NaN
            ΔΔCp[i] = NaN
            relΔTchar[i] = NaN
            relΔθchar[i] = NaN
            relΔΔCp[i] = NaN
        else
		    ΔTchar[i] = df.Tchar[i] - (db.Tchar[ii] + 273.15)
		    Δθchar[i] = df.θchar[i] - db.thetachar[ii]
		    ΔΔCp[i] = df.ΔCp[i] - db.DeltaCp[ii]
		    relΔTchar[i] = ΔTchar[i]/(db.Tchar[ii] + 273.15)
		    relΔθchar[i] = Δθchar[i]/db.thetachar[ii]
		    relΔΔCp[i] = ΔΔCp[i]/db.DeltaCp[ii]
        end
	end
	df[!, :ΔTchar] = ΔTchar
	df[!, :Δθchar] = Δθchar
	df[!, :ΔΔCp] = ΔΔCp
	df[!, :relΔTchar] = relΔTchar
	df[!, :relΔθchar] = relΔθchar
	df[!, :relΔΔCp] = relΔΔCp
	return df
end

function simulate_measurements(compare, meas, result::DataFrame; opt=std_opt)
	
	L = compare[1].L
	sp = compare[1].sp
	φ₀ = compare[1].df/compare[1].d
	gas = compare[1].gas
	prog = compare[2]
	solute_names = result.Name
	a = time_unit_conversion_factor(compare[6])
	CAS = GasChromatographySimulator.CAS_identification(string.(result.Name)).CAS

	sub = Array{GasChromatographySimulator.Substance}(undef, length(solute_names))
	for i=1:length(solute_names)
		
		#ii = findfirst(solute_names[i].==result[j].Name)
		if "θchar" in names(result) 
			Cag = GasChromatographySimulator.diffusivity(CAS[i], gas)
			sub[i] = GasChromatographySimulator.Substance(result.Name[i], CAS[i], result.Tchar[i], result.θchar[i], result.ΔCp[i], φ₀, "", Cag, 0.0, 0.0)
		else
			Cag = GasChromatographySimulator.diffusivity(CAS[i], gas)
			sub[i] = GasChromatographySimulator.Substance(result.Name[i], CAS[i], result.Tchar[i]+273.15, result.thetachar[i], result.DeltaCp[i], φ₀, "", Cag, 0.0, 0.0)
		end
	end
	
	par = Array{GasChromatographySimulator.Parameters}(undef, length(prog))
	pl = Array{DataFrame}(undef, length(prog))
	loss = Array{Float64}(undef, length(prog))
	for i=1:length(compare[2])
		#ii = findfirst(solute_names[i].==solute_names_compare)
		if "d" in names(result)
			d = mean(result.d) # if d was estimate use the mean value
		else
			d = compare[1].d
		end
		λ = meas[1].L/d
		col = GasChromatographySimulator.Column(λ*d, d, φ₀*d, sp, gas)
		par[i] = GasChromatographySimulator.Parameters(col, compare[2][i], sub, opt)
		try
			pl[i] = GasChromatographySimulator.simulate(par[i])[1]
		catch
			pl[i] = DataFrame(Name=solute_names, tR=NaN.*ones(length(solute_names)), τR=NaN.*ones(length(solute_names)))
		end
		#CAS = GasChromatographySimulator.CAS_identification(string.(result.Name)).CAS
		ΔtR = Array{Float64}(undef, length(solute_names))
		relΔtR = Array{Float64}(undef, length(solute_names))
		for k=1:length(solute_names)
			kk = findfirst(pl[i].Name[k].==compare[4])
			tR_compare = Array(compare[3][i, 2:(length(compare[4])+1)]).*a
			ΔtR[k] = pl[i].tR[k] - tR_compare[kk]
			relΔtR[k] = (pl[i].tR[k] - tR_compare[kk])/tR_compare[kk]		
		end
		pl[i][!, :ΔtR] = ΔtR
		pl[i][!, :relΔtR] = relΔtR
		loss[i] = sum(ΔtR.^2)/length(ΔtR)
	end
	return pl, loss, par
end

"""
	separate_error_columns(res)

If the result dataframe `res` contains columns with `Measurements` typed values, these columns will be split in two. The first column contains the value and uses the original column name. The second column contains the uncertainty and the name of the column is the original name with an added "_uncertainty". If the column is not of the `Measurements` type, it will be copied as is for the new dataframe.
"""
function separate_error_columns(res)
	new_res = DataFrame()
	for i=1:size(res)[2]
		if typeof(res[!,i]) == Array{Measurement{Float64}, 1}
			new_res[!, names(res)[i]] = Measurements.value.(res[!,i])
			new_res[!, names(res)[i]*"_uncertainty"] = Measurements.uncertainty.(res[!,i])
		else
			new_res[!, names(res)[i]] = res[!,i]
		end
	end
	return new_res
end

"""
    substance_list(substance_names, tRs)

Generate a vector of substance names, excluding rows with missing values in `tRs`.

# Arguments
- `substance_names`: A vector of substance names.
- `tRs`: A 2D array where each row corresponds to a set of retention times and columns to different substances.

# Returns
- A vector of substance names corresponding to non-missing rows in `tRs`.
"""
function substance_list(substance_names, tRs)
	substance_list_2d = Array{String}(undef, size(tRs))
	for i=1:size(tRs)[1]
		substance_list_2d[i,:] = substance_names
	end
	index_notmissing = Not(findall(ismissing.(tRs[:,:])))
	substance_list = substance_list_2d[index_notmissing]
	return substance_list
end

"""
    prog_list(prog, tRs)

Generate a vector of program conditions, excluding rows with missing values in `tRs`.

# Arguments
- `prog`: A program condition or vector of programs.
  - If a single `Program`: repeated for all measurements and substances.
  - If a `Vector{Program}`: one program per measurement (length must match number of rows in `tRs`).
- `tRs`: A 2D array where each row corresponds to a set of retention times and columns to different substances.

# Returns
- A vector of program conditions corresponding to non-missing rows in `tRs`.

# Examples
```jldoctest
julia> using RetentionParameterEstimator, GasChromatographySimulator

julia> tRs = [100.0 150.0; 200.0 250.0]
2×2 Matrix{Float64}:
 100.0  150.0
 200.0  250.0

julia> prog_single = GasChromatographySimulator.Program([0.0, 5.0], [40.0, 340.0], [150000.0, 250000.0], "vacuum", "min")
...

julia> prog_vec = [prog_single, prog_single]
2-element Vector{GasChromatographySimulator.Program}:
 ...

julia> # Single program (repeated)
julia> prog_list(prog_single, tRs)
4-element Vector{GasChromatographySimulator.Program}:
 ...

julia> # Vector of programs (one per measurement)
julia> prog_list(prog_vec, tRs)
4-element Vector{GasChromatographySimulator.Program}:
 ...
```
"""
function prog_list(prog, tRs)
	prog_list_2d = Array{GasChromatographySimulator.Program}(undef, size(tRs))
	
	# Handle both single program and vector of programs
	if isa(prog, Vector)
		# Vector of programs: one per measurement (row)
		# Each column gets the same vector of programs
		for j=1:size(tRs)[2]
			prog_list_2d[:,j] = prog
		end
	else
		# Single program: repeat for all measurements
		for j=1:size(tRs)[2]
			prog_list_2d[:,j] = fill(prog, size(tRs)[1])
		end
	end
	
	index_notmissing = Not(findall(ismissing.(tRs[:,:])))
	prog_list = prog_list_2d[index_notmissing]
	return prog_list
end

"""
    time_unit_conversion_factor(time_unit)

Get the conversion factor for time units.

# Arguments
- `time_unit`: String indicating time unit, either `"min"` or `"s"`.

# Returns
- Conversion factor: `60.0` for `"min"`, `1.0` for `"s"` (or any other value).

# Examples
```jldoctest
julia> RetentionParameterEstimator.time_unit_conversion_factor("min")
60.0

julia> RetentionParameterEstimator.time_unit_conversion_factor("s")
1.0
```
"""
function time_unit_conversion_factor(time_unit)
    if time_unit == "min"
        a = 60.0
    else
        a = 1.0
    end
	return a
end

"""
    prepare_optimization_data(tRs, solute_names, prog, time_unit)

Prepare data for optimization by converting time units, extracting retention times,
filtering missing values, and creating vectorized program and substance lists.

# Arguments
- `tRs`: DataFrame with retention times (first column is measurement names, subsequent columns are retention times for each solute).
- `solute_names`: Vector of solute names (one per column in `tRs` after the first column).
- `prog`: Vector of programs (one per measurement/row in `tRs`).
- `time_unit`: String indicating time unit, either `"min"` or `"s"`.

# Returns
- `tRs_`: Vector of retention times with missing values filtered out.
- `prog_`: Vector of programs matching the filtered `tRs_` (expanded to 2D and filtered).
- `subst_list_`: Vector of substance names matching the filtered `tRs_`.

# Examples
```jldoctest
julia> using RetentionParameterEstimator, DataFrames

julia> tRs = DataFrame(measurement=["meas1", "meas2"], solute1=[100.0, 200.0], solute2=[150.0, 250.0])
2×3 DataFrame
 Row │ measurement  solute1   solute2  
     │ String       Float64   Float64  
─────┼─────────────────────────────────
   1 │ meas1           100.0    150.0
   2 │ meas2           200.0    250.0

julia> solute_names = ["solute1", "solute2"]
2-element Vector{String}:
 "solute1"
 "solute2"

julia> prog = [GasChromatographySimulator.Program([0.0, 5.0], [40.0, 340.0], [150000.0, 250000.0], "vacuum", "min") for _ in 1:2]
2-element Vector{GasChromatographySimulator.Program}:
 ...

julia> tRs_, prog_, subst_list_ = RetentionParameterEstimator.prepare_optimization_data(tRs, solute_names, prog, "min")
([100.0, 200.0, 150.0, 250.0], [...], ["solute1", "solute1", "solute2", "solute2"])
```
"""
function prepare_optimization_data(tRs, solute_names, prog, time_unit)
    a = time_unit_conversion_factor(time_unit)
    tR_meas = Array(tRs[:,2:end]).*a
    index_notmissing = Not(findall(ismissing.(tR_meas[:,:])))
    tRs_ = collect(skipmissing(tR_meas[:,:]))
    prog_ = prog_list(prog, tR_meas)
    subst_list_ = substance_list(solute_names, tR_meas)
    return tRs_, prog_, subst_list_
end

"""
    prepare_single_substance_data(tR_meas, prog, solute_name)

Prepare data for single-substance optimization by filtering missing values and creating 
vectorized program and substance lists.

# Arguments
- `tR_meas`: 1D array of retention times (already time-unit converted) for a single substance.
- `prog`: Vector of programs (one per measurement).
- `solute_name`: Name of the single substance.

# Returns
- `tRs_`: Vector of retention times with missing values filtered out.
- `prog_`: Vector of programs matching the filtered `tRs_`.
- `subst_list_`: Vector with repeated substance name matching the filtered `tRs_`.

# Examples
```jldoctest
julia> using RetentionParameterEstimator

julia> tR_meas = [100.0, missing, 200.0, 250.0]
4-element Vector{Union{Missing, Float64}}:
 100.0
  missing
 200.0
 250.0

julia> prog = [GasChromatographySimulator.Program([0.0, 5.0], [40.0, 340.0], [150000.0, 250000.0], "vacuum", "min") for _ in 1:4]
4-element Vector{GasChromatographySimulator.Program}:
 ...

julia> tRs_, prog_, subst_list_ = RetentionParameterEstimator.prepare_single_substance_data(tR_meas, prog, "solute1")
([100.0, 200.0, 250.0], [...], ["solute1", "solute1", "solute1"])
```
"""
function prepare_single_substance_data(tR_meas, prog, solute_name)
    index_notmissing = Not(findall(ismissing.(tR_meas[:])))
    tRs_ = collect(skipmissing(tR_meas[:]))
    prog_ = prog[index_notmissing]
    subst_list_ = fill(solute_name, length(tRs_))
    return tRs_, prog_, subst_list_
end