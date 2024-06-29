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
	if compare[6] == "min"
		a = 60.0
	else
		a = 1.0
	end
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