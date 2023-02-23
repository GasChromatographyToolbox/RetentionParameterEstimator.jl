### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 8a1861d4-4acf-11ed-39b0-7b1af92df0ee
using DrWatson

# ╔═╡ 59dd4f69-4a2a-4fd6-b679-6f807c297948
@quickactivate "RetentionParameterEstimator_Paper"

# ╔═╡ 09422105-a747-40ac-9666-591326850d8f
begin 
	using RetentionParameterEstimator, GasChromatographySimulator, Plots, OptimizationOptimJL, PlutoUI, CSV, DataFrames, Statistics, LaTeXStrings, StatsPlots, Measurements
	TableOfContents()
end

# ╔═╡ ee641f1d-4840-4250-9eec-0f0166f49ff2
using ForwardDiff

# ╔═╡ 5be3a351-f6e0-4cc6-800d-a55f8decf889
plotly()

# ╔═╡ eb5fc23c-2151-47fa-b56c-5771a4f8b9c5
html"""
<style>
  main {
	max-width: 800px;
  }
</style>
"""

# ╔═╡ 6d4ec54b-01b2-4208-9b9e-fcb70d236c3e
md"""
# Estimation of K-centric retention parameters
by temperature programmed GC **simulations** with different temperature programs
"""

# ╔═╡ 3e808e31-081e-48d2-908e-d66726e270e8
# if nothing is selected -> load the data for the paper from the github
# for alternativa data and isothermal data load these from the github only if the measurements are loaded from the github

# ╔═╡ ebc2a807-4413-4721-930a-6328ae72a1a9
md"""
## Select measurement data
Measured chromatograms used to **estimate** the parameters by optimization:

$(@bind file_meas FilePicker([MIME("text/csv")]))
"""

# ╔═╡ eb14e619-82d4-49ac-ab2a-28e56230dbc6
begin
	if isnothing(file_meas)
		md"""
		Selected chromatograms for **estimation**: _nothing_
		
		Please select a file of chromatograms for **estimation**!
		"""
	else
		meas = RetentionParameterEstimator.load_chromatograms(file_meas);
		md"""
		Selected chromatograms for **estimation**: $(file_meas["name"])
		"""
	end
end

# ╔═╡ d745c22b-1c96-4a96-83da-abb1df91ab87
begin
	if !isnothing(file_meas)
		md"""
		Select measurements:
		$(@bind selected_measurements confirm(MultiSelect(meas[3].measurement; default=meas[3].measurement)))

		
		Select solutes:
		$(@bind selected_solutes confirm(MultiSelect(meas[4]; default=meas[4])))
		"""
	end
end

# ╔═╡ 365ec346-65e4-4bb8-9b8f-824ed0e5e3f2
# reduce the modes to the two methods m1 and m2 of the paper + test

# ╔═╡ f3ffd4ce-a378-4033-88e9-bc1fb8cc4bbe
md"""
## Select mode

$(@bind select_mode confirm(Select(["check measurement for plausibility", "estimate retention parameters", "estimate retention parameters and diameter, 2x single", "estimate retention parameters and diameter, single", "estimate retention parameters and diameter, all", "estimate retention parameters and diameter, direct", "estimate retention parameters and diameter, all, mod", "estimate retention parameters and diameter, multistart"])))
"""

# ╔═╡ e98f4b1b-e577-40d0-a7d8-71c53d99ee1b
if select_mode == "estimate retention parameters"
	#md"""
	#Column parameters:

	#L in m: $(@bind L_input NumberField(0.0:0.01:1000.0; default=meas[1].L))

	#d in mm: $(@bind d_input NumberField(0.0:0.0001:1.0; default=meas[1].d*1000.0))

	@bind col_input confirm(
		PlutoUI.combine() do Child
		md"""
		## Column dimensions 
		L in m: $(Child("L", NumberField(0.0:0.01:1000.0; default=meas[1].L)))

		d in mm: $(Child("d", NumberField(0.0:0.0001:1.0; default=meas[1].d*1000.0))) 
		"""
		end
	)
end	

# ╔═╡ 882b9cf0-7320-408d-b7aa-3e20fffb0ff2
md"""
## Test confidence
"""

# ╔═╡ 1f34c58b-5f35-422b-a645-4aa9e96a51f6
"squared" # metric

# ╔═╡ 0d61fd05-c0c6-4764-9f96-3b867a456ad3
md"""
## Select verification data

Measured chromatograms used to **verify** the estimated parameters:

$(@bind file_comp FilePicker([MIME("text/csv")]))

Select folder of chromatograms:
$(@bind folder_chromatograms confirm(TextField(default="/Users/janleppert/Documents/GitHub/RetentionParameterEstimator_Paper/papers/data/chromatograms")))
"""

# ╔═╡ 762c877d-3f41-49de-a7ea-0eef1456ac11
begin
	if isnothing(file_comp)
		md"""
		Selected chromatograms for **verification**: _nothing_

		Please select a file of chromatograms for **verification**!
		"""
	else
		if file_meas == file_comp
			md"""
			Selected chromatograms for **verification**: $(file_comp["name"])
	
			_**Attention**_: The selected data for estimation and verification is the same!
			"""
		else
			comp = RetentionParameterEstimator.load_chromatograms(file_comp);
			md"""
			Selected chromatograms for **verification**: $(file_comp["name"])
			"""
		end
	end
end

# ╔═╡ 07a7e45a-d73e-4a83-9323-700d3e2a88cc
begin
	if !isnothing(file_comp)
		md"""
		Select comparison measurements:
		$(@bind selected_comparison confirm(MultiSelect(comp[3].measurement; default=comp[3].measurement)))
		"""
	end
end

# ╔═╡ af2aa109-c514-40dc-a5a9-c7b9974dd8f1
begin
	chroms = Array{DataFrame}(undef, length(selected_comparison))
	for i=1:length(selected_comparison)
		chroms[i] = DataFrame(CSV.File(joinpath(folder_chromatograms,string("Chrom_", selected_comparison[i],".csv")); header=4))
	end
	chroms
end

# ╔═╡ ae424251-f4f7-48aa-a72b-3be85a193705
md"""
## Verification

Using the estimated parameters (``T_{char}``, ``θ_{char}``, ``ΔC_p``) resp. (``T_{char}``, ``θ_{char}``, ``ΔC_p``, ``d``) to simulate other temperature programs. Compare the simulated retention times with measured retention times.
"""

# ╔═╡ f83b26e7-f16d-49ac-ab3b-230533ac9c82
md"""
## Alternative parameters

Select file with alternative parameters:

$(@bind file_db FilePicker([MIME("text/csv")]))
"""

# ╔═╡ 6187b9a9-83e4-49e0-be6f-8aea1a09da7c
# add test for correct data format

# ╔═╡ 903d12f1-6d6f-4a71-8095-7457cffcafc4
md"""
## Isothermal ``\ln{k}`` data

Select file with isothermal measured ``\ln{k}`` data.

$(@bind file_isolnk FilePicker([MIME("text/csv")]))
"""

# ╔═╡ 347c53d6-9c94-47de-814a-aef43fa5aae2
# add test for correct data format

# ╔═╡ 66287dfc-fe77-4fa8-8892-79d2d3de6cb3
begin
	if isnothing(file_isolnk)
		md"""
		Selected isothermal measured ``\ln{k}`` data: _nothing_

		Please select a file of isothermal measured ``\ln{k}`` data!
		"""
	else
		isolnk = DataFrame(CSV.File(file_isolnk["data"]))
		md"""
		Selected isothermal measured ``\ln{k}`` data: $(file_isolnk["name"])

		$(embed_display(isolnk))
		"""
	end
end


# ╔═╡ b4fc2f6e-4500-45da-852e-59fb084e365a
md"""
## Plot ``\ln{k}`` over ``T``

Select solute for ``\ln{k}`` plot:
$(@bind select_solute confirm(Select(meas[4]; default=meas[4][1])))
"""

# ╔═╡ 9178967d-26dc-43be-b6e4-f35bbd0b0b04
md"""
# End
"""

# ╔═╡ 46fab3fe-cf88-4cbf-b1cc-a232bb7520db
# filter function for selected measurements and selected solutes
function filter_selected_measurements(meas, selected_measurements, selected_solutes)
	index_measurements = Array{Int}(undef, length(selected_measurements))
	for i=1:length(selected_measurements)
		index_measurements[i] = findfirst(selected_measurements[i].==meas[3].measurement)
	end
	index_solutes = Array{Int}(undef, length(selected_solutes))
	for i=1:length(selected_solutes)
		index_solutes[i] = findfirst(selected_solutes[i].==names(meas[3]))
	end
	meas_select = (meas[1], meas[2][index_measurements], meas[3][index_measurements, [1; index_solutes]], selected_solutes, meas[5], meas[6])
	return meas_select
end

# ╔═╡ b2c254a2-a5d6-4f18-803a-75d048fc7cdf
meas_select = filter_selected_measurements(meas, selected_measurements, selected_solutes)

# ╔═╡ 31b176dc-b667-4604-b490-046d55e545c9
meas_select[3][!,2] # tRs

# ╔═╡ 1ed16860-6899-4a17-aa2e-dffd32a4cf82
meas_select[1].L # L

# ╔═╡ 2c418101-9184-49b7-8ad8-ae3be25d7cdb
meas_select[1].d # d

# ╔═╡ fcc0c0cf-aee7-42bb-9ad9-790325db31aa
meas_select[2] # prog

# ╔═╡ cdef487b-1a2b-470a-9f13-e08b88e079a2
meas_select[1].gas # gas

# ╔═╡ fb4ccf8c-18af-4b7d-8bff-50028afdcb6a
@isdefined col_input

# ╔═╡ 06d5dd5a-6883-4b62-8942-dc19e7b08e4d
function comparison(res, meas, comp)
	opt = GasChromatographySimulator.Options(ng=true, odesys=false)
	CAS_comp = GasChromatographySimulator.CAS_identification(comp[4]).CAS
	CAS_meas = GasChromatographySimulator.CAS_identification(meas[4]).CAS
	i_sub = findall(x->x in CAS_meas, CAS_comp) # indices of common elements of CAS_comp in CAS_meas (indeices relative to CAS_meas)
	sub = Array{GasChromatographySimulator.Substance}(undef, length(i_sub))
	for i=i_sub
		CAS = CAS_meas[i]
		ii = findfirst(meas[4][i].==res.Name)
		Cag = GasChromatographySimulator.diffusivity(CAS, comp[1].gas)
		if "θchar" in names(res) 
			sub[i] = GasChromatographySimulator.Substance(res.Name[ii], CAS, res.Tchar[ii], res.θchar[ii], res.ΔCp[ii], comp[1].df/comp[1].d, "", Cag, 0.0, 0.0)
		else
			sub[i] = GasChromatographySimulator.Substance(res.Name[ii], CAS, res.Tchar[ii]+273.15, res.thetachar[ii], res.DeltaCp[ii], comp[1].df/comp[1].d, "", Cag, 0.0, 0.0)
		end
		# this is for phase ratio as set in comp
	end
	if comp[6] == "min"
		a = 60.0
	else
		a = 1.0
	end

	par = Array{GasChromatographySimulator.Parameters}(undef, length(comp[3].measurement))
	pl = Array{DataFrame}(undef, length(comp[3].measurement))
	loss = Array{Float64}(undef, length(comp[3].measurement))
	for i=1:length(comp[3].measurement)
		#ii = findfirst(solute_names[i].==solute_names_compare)
		if "d" in names(res)
			d = mean(res.d) # if d was estimate use the mean value
		elseif @isdefined col_input
			d = col_input.d/1000.0
		else
			d = comp[1].d
		end
		col = GasChromatographySimulator.Column(comp[1].L, d, comp[1].df/comp[1].d*d, comp[1].sp, comp[1].gas)
		par[i] = GasChromatographySimulator.Parameters(col, comp[2][i], sub, opt)
		try
			pl[i] = GasChromatographySimulator.simulate(par[i])[1]
			
		catch
			pl[i] = DataFrame(Name=comp[4], tR=NaN.*ones(length(comp[4])), τR=NaN.*ones(length(comp[4])))
		end
		#CAS = GasChromatographySimulator.CAS_identification(string.(result[j].Name)).CAS
		ΔtR = Array{Float64}(undef, length(comp[4]))
		relΔtR = Array{Float64}(undef, length(comp[4]))
		for k=1:length(comp[4])
			kk = findfirst(pl[i].Name[k].==comp[4])
			tR_compare = Array(comp[3][i, 2:(length(comp[4])+1)]).*a
			ΔtR[k] = pl[i].tR[k] - tR_compare[kk]
			relΔtR[k] = (pl[i].tR[k] - tR_compare[kk])/tR_compare[kk]		
		end
		pl[i][!, :tR] = pl[i].tR./a
		pl[i][!, :τR] = pl[i].τR./a
		pl[i][!, :ΔtR] = ΔtR./a
		pl[i][!, :relΔtR] = relΔtR
		loss[i] = sum(ΔtR.^2)/length(ΔtR)
	end
	return pl, loss, par
end

# ╔═╡ 505f722f-5429-4282-a1c9-d41fc3284c11
function flagged_loss(meas, df, index_flag)
	if meas[6] == "min"
		a = 60.0
	else
		a = 1.0
	end
	tRcalc = Array{Float64}(undef, length(meas[3].measurement), length(index_flag))
	tRmeas = Array{Float64}(undef, length(meas[3].measurement), length(index_flag))
	for j=1:length(index_flag)
		for i=1:length(meas[3].measurement)
			tRcalc[i,j] = RetentionParameterEstimator.tR_calc(df.Tchar[index_flag[j]], df.θchar[index_flag[j]], df.ΔCp[index_flag[j]], meas[1].L, meas[1].d, meas[2][i], meas[1].gas)
			tRmeas[i,j] = meas[3][!, index_flag[j]+1][i]*a
		end
	end
	(tRmeas.-tRcalc).^2, tRcalc, tRmeas
end

# ╔═╡ 72312748-de7d-4b15-8ef8-ea7165d19640
function check_measurement(meas; min_th=0.1, loss_th=1.0)
	Tchar_est, θchar_est, ΔCp_est, Telu_max = RetentionParameterEstimator.estimate_start_parameter(meas[3], meas[1], meas[2]; time_unit=meas[6])
	df = RetentionParameterEstimator.estimate_parameters(meas[3], meas[4], meas[1], meas[2], Tchar_est, θchar_est, ΔCp_est; mode="Kcentric_single", pout=meas[5], time_unit=meas[6])[1]
	index_flag = findall(df.min.>min_th)
	if length(index_flag) == 0
		check = true
		msg = "retention times in normal range"
		df_flag = DataFrame()
	else
		check = false
		msg = "discrapancy of retention times detected"
		loss_flag, tRcalc, tRmeas = flagged_loss(meas, df, index_flag)
		index_loss_flag = findall(loss_flag.>loss_th)
		flag_meas = Array{String}(undef, length(index_loss_flag))
		flag_sub = Array{String}(undef, length(index_loss_flag))
		tRmeas_ = Array{Float64}(undef, length(index_loss_flag))
		tRcalc_ = Array{Float64}(undef, length(index_loss_flag))
		for i=1:length(index_loss_flag)
			flag_meas[i] = meas[3].measurement[index_loss_flag[i][1]]
			flag_sub[i] = meas[4][index_flag[index_loss_flag[i][2]]]
			tRmeas_[i] = tRmeas[index_loss_flag[i][1], index_loss_flag[i][2]]
			tRcalc_[i] = tRcalc[index_loss_flag[i][1], index_loss_flag[i][2]]
		end
		df_flag = DataFrame(measurement=flag_meas, solute=flag_sub, tRmeas=tRmeas_, tRcalc=tRcalc_)
	end
	return check, msg, df_flag, df, Telu_max
end

# ╔═╡ 3b40b0b1-7007-48c7-b47b-dbeaf501b73d
begin
	std_opt = RetentionParameterEstimator.std_opt
	if select_mode == "check measurement for plausibility"
		check, msg, df_flag, res, Telu_max = check_measurement(meas_select; min_th=0.1, loss_th=1.0)
		md"""
		## Results
		
		$(check)

		$(msg)

		$(df_flag)
		
		!!! Hint
		
		Check for the plausibility of measurements, if a solute is listed there, all retention times for all used measurements should be checked, not only the listed flagged ones. An error in one retention time, e.g. transposed digits or typo, could produce deviations from the model in other (correct) measurements
		-> Only list the name of the substance, not the measurements.

		$(res)
		"""
	elseif select_mode == "estimate retention parameters"
		col = GasChromatographySimulator.Column(col_input.L, col_input.d*1e-3, meas_select[1].df, meas_select[1].sp, meas_select[1].gas)
		Tchar_est, θchar_est, ΔCp_est, Telu_max = RetentionParameterEstimator.estimate_start_parameter(meas_select[3], col, meas_select[2]; time_unit=meas_select[6])
		res = RetentionParameterEstimator.estimate_parameters(meas_select[3], meas_select[4], col, meas_select[2], Tchar_est, θchar_est, ΔCp_est; mode="Kcentric_single", pout=meas_select[5], time_unit=meas_select[6], opt=std_opt)[1]
		md"""
		## Results

		L/d ratio: $(col_input[1]/(col_input[2]*1e-3))
		
		$(res)
		"""
	elseif select_mode == "estimate retention parameters and diameter, single"
		tRs = meas_select[3][!,findall((collect(any(ismissing, c) for c in eachcol(meas_select[3]))).==false)]
		solute_names = meas_select[4][findall((collect(any(ismissing, c) for c in eachcol(meas_select[3]))).==false)[2:end].-1]
		Tchar_est, θchar_est, ΔCp_est, Telu_max = RetentionParameterEstimator.estimate_start_parameter(tRs, meas_select[1], meas_select[2]; time_unit=meas_select[6])
		res = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, meas_select[1], meas_select[2], Tchar_est, θchar_est, ΔCp_est; pout=meas_select[5], time_unit=meas_select[6], mode="dKcentric_single", opt=std_opt)[1]
		md"""
		## Results
		
		L/d ratio: $(mean(meas_select[1].L./res.d) ± std(meas_select[1].L./res.d))
		
		$(res)
		"""
	elseif select_mode == "estimate retention parameters and diameter, all"
		tRs = meas_select[3][!,findall((collect(any(ismissing, c) for c in eachcol(meas_select[3]))).==false)]
		solute_names = meas_select[4][findall((collect(any(ismissing, c) for c in eachcol(meas_select[3]))).==false)[2:end].-1]
		Tchar_est, θchar_est, ΔCp_est, Telu_max = RetentionParameterEstimator.estimate_start_parameter(tRs, meas_select[1], meas_select[2]; time_unit=meas_select[6])
		#res_dK_single = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, meas_select[1], meas_select[2], Tchar_est.*(1 + rand((-1.0, 1.0))*rand()/100), θchar_est.*(1 + rand((-1.0, 1.0))*rand()/100), ΔCp_est.*(1 + rand((-1.0, 1.0))*rand()/100); pout=meas_select[5], time_unit=meas_select[6], mode="dKcentric_single")[1]
		res_dK_single = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, meas_select[1], meas_select[2], Tchar_est, θchar_est, ΔCp_est; pout=meas_select[5], time_unit=meas_select[6], mode="dKcentric_single")[1]
	
		#res = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, meas_select[1], meas_select[2], res_dK_single.Tchar, res_dK_single.θchar, res_dK_single.ΔCp; pout=meas_select[5], time_unit=meas_select[6], mode="dKcentric", opt=std_opt)[1]
		
		new_col = GasChromatographySimulator.Column(meas_select[1].L, mean(res_dK_single.d), meas_select[1].df, meas_select[1].sp, meas_select[1].gas)
		res = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, new_col, meas_select[2], res_dK_single.Tchar, res_dK_single.θchar, res_dK_single.ΔCp; pout=meas_select[5], time_unit=meas_select[6], mode="dKcentric", opt=std_opt)[1]
		md"""
		## Results
		
		L/d ratio: $(meas_select[1].L/res.d[1])
		
		$(res)
		"""
	elseif select_mode == "estimate retention parameters and diameter, direct"
		tRs = meas_select[3][!,findall((collect(any(ismissing, c) for c in eachcol(meas_select[3]))).==false)]
		solute_names = meas_select[4][findall((collect(any(ismissing, c) for c in eachcol(meas_select[3]))).==false)[2:end].-1]
		Tchar_est, θchar_est, ΔCp_est, Telu_max = RetentionParameterEstimator.estimate_start_parameter(tRs, meas_select[1], meas_select[2]; time_unit=meas_select[6])
		res = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, meas_select[1], meas_select[2], Tchar_est, θchar_est, ΔCp_est; pout=meas_select[5], time_unit=meas_select[6], mode="dKcentric")[1]
		md"""
		## Results
		
		L/d ratio: $(meas_select[1].L/res.d[1])
		
		$(res)
		"""
	elseif select_mode == "estimate retention parameters and diameter, all, mod"
		tRs = meas_select[3][!,findall((collect(any(ismissing, c) for c in eachcol(meas_select[3]))).==false)]
		solute_names = meas_select[4][findall((collect(any(ismissing, c) for c in eachcol(meas_select[3]))).==false)[2:end].-1]
		Tchar_est, θchar_est, ΔCp_est, Telu_max = RetentionParameterEstimator.estimate_start_parameter(tRs, meas_select[1], meas_select[2]; time_unit=meas_select[6])
		#res_dK_single = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, meas_select[1], meas_select[2], Tchar_est.*(1 + rand((-1.0, 1.0))*rand()/100), θchar_est.*(1 + rand((-1.0, 1.0))*rand()/100), ΔCp_est.*(1 + rand((-1.0, 1.0))*rand()/100); pout=meas_select[5], time_unit=meas_select[6], mode="dKcentric_single")[1]
		res_dK_single = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, meas_select[1], meas_select[2], Tchar_est, θchar_est, ΔCp_est; pout=meas_select[5], time_unit=meas_select[6], mode="dKcentric_single")[1]
	
		#res = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, meas_select[1], meas_select[2], res_dK_single.Tchar, res_dK_single.θchar, res_dK_single.ΔCp; pout=meas_select[5], time_unit=meas_select[6], mode="dKcentric", opt=std_opt)[1]
		
		new_col = GasChromatographySimulator.Column(meas_select[1].L, mean(res_dK_single.d), meas_select[1].df, meas_select[1].sp, meas_select[1].gas)
		res = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, new_col, meas_select[2], res_dK_single.Tchar.+6.0, res_dK_single.θchar.+3.0, res_dK_single.ΔCp.*1.1; pout=meas_select[5], time_unit=meas_select[6], mode="dKcentric", opt=std_opt)[1]
		md"""
		## Results
		
		L/d ratio: $(meas_select[1].L/res.d[1])
		
		$(res)
		"""
	elseif select_mode == "estimate retention parameters and diameter, multistart"
		tRs = meas_select[3][!,findall((collect(any(ismissing, c) for c in eachcol(meas_select[3]))).==false)]
		solute_names = meas_select[4][findall((collect(any(ismissing, c) for c in eachcol(meas_select[3]))).==false)[2:end].-1]
		Tchar_est, θchar_est, ΔCp_est, Telu_max = RetentionParameterEstimator.estimate_start_parameter(tRs, meas_select[1], meas_select[2]; time_unit=meas_select[6])
		
		res_multi = Array{DataFrame}(undef, 10)
		min = Array{Float64}(undef, 10)
		loss_init = Array{Float64}(undef, 10)
		init_Tchar = Array{Array{Float64}}(undef, 10)
		init_θchar = Array{Array{Float64}}(undef, 10)
		init_ΔCp = Array{Array{Float64}}(undef, 10)
		#Threads.@threads for i=1:10
		for i=1:10
			if i==1
				f_Tchar = 1.0
				f_θchar = 1.0
				f_ΔCp = 1.0
			else
				f_Tchar = 1+randn()/20
				f_θchar = 1+randn()/20
				f_ΔCp = 1+randn()/20
			end
			init_Tchar[i] = Tchar_est.*f_Tchar
			init_θchar[i] = θchar_est.*f_θchar
			init_ΔCp[i] = ΔCp_est.*f_ΔCp
			loss_init[i] = RetentionParameterEstimator.loss(Array(tRs[!,2:end]), init_Tchar[i], init_θchar[i], init_ΔCp[i], meas_select[1].L, meas_select[1].d, meas_select[2], meas_select[1].gas)
			if  loss_init[i] < 1
				res_multi[i] = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, meas_select[1], meas_select[2], Tinit_Tchar[i], init_θchar[i], init_ΔCp[i]; pout=meas_select[5], time_unit=meas_select[6], mode="dKcentric", opt=std_opt)[1]
				min[i] = res_multi[i].min[1]
			else
				res_multi[i] = DataFrame()
				min[i] = loss_init[i]
			end
		end
		res = res_multi[findfirst(minimum(min).==min)]
		#=md"""
		## Results
		
		L/d ratio: $(meas_select[1].L/res.d[1])
		
		$(res)
		"""=#
	elseif select_mode == "estimate retention parameters and diameter, 2x single"
		tRs = meas_select[3][!,findall((collect(any(ismissing, c) for c in eachcol(meas_select[3]))).==false)]
		solute_names = meas_select[4][findall((collect(any(ismissing, c) for c in eachcol(meas_select[3]))).==false)[2:end].-1]
		
		Tchar_est, θchar_est, ΔCp_est, Telu_max = RetentionParameterEstimator.estimate_start_parameter(tRs, meas_select[1], meas_select[2]; time_unit=meas_select[6])
		
		res_dKcentric_single = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, meas_select[1], meas_select[2], Tchar_est, θchar_est, ΔCp_est; pout=meas_select[5], time_unit=meas_select[6], mode="dKcentric_single", opt=std_opt)[1]
		
		d = mean(res_dKcentric_single.d) ± std(res_dKcentric_single.d)

		new_col = GasChromatographySimulator.Column(meas_select[1].L, mean(res_dKcentric_single.d), meas_select[1].df, meas_select[1].sp, meas_select[1].gas)
		res = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, new_col, meas_select[2], res_dKcentric_single.Tchar, res_dKcentric_single.θchar, res_dKcentric_single.ΔCp; pout=meas_select[5], time_unit=meas_select[6], mode="Kcentric_single", opt=std_opt)[1]
		#res = RetentionParameterEstimator.estimate_parameters(tRs, solute_names, new_col, meas_select[2], Tchar_est, θchar_est, ΔCp_est; pout=meas_select[5], time_unit=meas_select[6], mode="Kcentric_single", opt=std_opt)[1]

		res[!, :d] = mean(res_dKcentric_single.d).*ones(length(res.Name))
		res[!, :d_std] = std(res_dKcentric_single.d).*ones(length(res.Name))
		md"""
		## Results

		d = $(d)
		
		L/d ratio: $(meas_select[1].L./d)
		
		$(res)
		"""
	end
end

# ╔═╡ 992fae9e-5257-4996-840c-7bf8969d3e02
res.Tchar.-Tchar_est

# ╔═╡ f75dcbf4-b146-47ea-8fa2-1cdf3f772aeb
res.θchar.-θchar_est

# ╔═╡ d09e48ec-1d38-46c8-b8ce-19d9a83bfa51
res.ΔCp.-ΔCp_est

# ╔═╡ 4779a2ed-c827-4d4e-aa5f-614aaca6c173
std_opt # opt

# ╔═╡ 3ac78538-e2ec-466d-bfe4-d55bcc009db0
p = [meas_select[3][!,end].*60.0, meas_select[1].L, meas_select[1].d, meas_select[2], std_opt, meas_select[1].gas, "squared"]

# ╔═╡ f650199a-554d-49ff-8346-f0e5c6819b7c
LF(x) = RetentionParameterEstimator.opt_Kcentric(x, p)

# ╔═╡ 07187629-64b7-4869-9eeb-7d530e43587d
LF([300, 30, 100])

# ╔═╡ 26672fc0-3b99-44e6-8a46-47dbeb3bca38
H(x) = ForwardDiff.hessian(LF, x)

# ╔═╡ deb2179f-06e5-4439-a3a0-1664fc1a5083
xx = [res.Tchar[end], res.θchar[end], res.ΔCp[end]]

# ╔═╡ c3e0b7e9-f8d3-4b63-9e36-fb4995211a6c
LF(xx)

# ╔═╡ 0e6aa86b-1cd5-4281-b0f1-d8e96b2d6e2e
H(xx)

# ╔═╡ a1c0f53a-925c-4800-a920-955e9da97e12
inv(H(xx))

# ╔═╡ 9347878d-66f0-464b-a315-9c6eb46ea43d
sqrt.(abs.(inv(H(xx))))

# ╔═╡ f408ded9-9130-4652-b0aa-219936f4c947
xx[1] ± sqrt.(abs.(inv(H(xx))))[1,1]

# ╔═╡ 6475ab4f-f701-4398-ac25-b65d0844301a
xx[2] ± sqrt.(abs.(inv(H(xx))))[2,2]

# ╔═╡ 20902ac4-f210-4ddd-90f5-a713eb4bd9b5
xx[3] ± sqrt.(abs.(inv(H(xx))))[3,3]

# ╔═╡ 93224376-c99e-4182-b430-3d75a52ee0f1
res.min[end]

# ╔═╡ 47e7621e-58e9-4701-912d-82f646e64e0f
LF(xx)/res.min[1]

# ╔═╡ 0f4c35c4-32f7-4d11-874d-1f23daad7da8
begin
	head = ["file:", file_meas["name"], 
		"selected measurements:", join(selected_measurements, " "),
		"mode:", select_mode
	]
	io = IOBuffer()
	CSV.write(io, DataFrame[], header=head)
	CSV.write(io, res, append=true, writeheader=true)
	#export_str_ = export_str(opt_values, col_values, prog_values, peaklist)
	md"""
	## Export Results
	Filename: $(@bind result_filename TextField((20,1); default="Result.csv"))
	"""
end

# ╔═╡ b4f17579-9994-46e1-a3d0-6030650f0dbe
md"""
$(DownloadButton(io.data, result_filename))
"""

# ╔═╡ 116ccf37-4a18-44ac-ae6e-98932901a8b0
begin
	comp_select = filter_selected_measurements(comp, selected_comparison, selected_solutes)
	pl, loss, par = comparison(res, meas_select, comp_select)
	df_loss = DataFrame(comp=selected_comparison, loss=loss, sqrt_loss=sqrt.(loss))
end

# ╔═╡ 157f24e6-1664-41c8-9079-b9dd0a2c98a9
pl

# ╔═╡ ddadd79d-7b74-4b2c-ad59-a31a5692cf3b
begin
	p_chrom = Array{Plots.Plot}(undef, length(pl))
	for i=1:length(pl)
		p_chrom[i] = plot(xlabel="time in $(meas_select[6])", legend=false)
		plot!(chroms[i][!,1], -1e-6.*chroms[i][!,2])
		i1 = findfirst(minimum(pl[i].tR)*0.95.<chroms[i][!,1])
		i2 = findlast(maximum(pl[i].tR)*1.05.>chroms[i][!,1])
		min_ = minimum(-1e-6.*chroms[i][!,2][i1:i2])*1.05
		max_ = maximum(GasChromatographySimulator.plot_chromatogram(pl[i], (minimum(pl[i].tR)*0.95, maximum(pl[i].tR)*1.05))[3])*1.05
		GasChromatographySimulator.plot_chromatogram!(p_chrom[i], pl[i], (minimum(pl[i].tR)*0.95, maximum(pl[i].tR)*1.05); annotation=false)
		xlims!((minimum(pl[i].tR)*0.95, maximum(pl[i].tR)*1.05))
		ylims!(min_, max_)
		#add marker for measured retention times
		for j=1:length(comp_select[4])
			plot!(p_chrom[i], comp_select[3][i,j+1].*ones(2), [min_, max_], c=:orange)
		end
		plot!(p_chrom[i], title=comp_select[3].measurement[i])
	end
	p_chrom
	plot(p_chrom..., layout=(length(p_chrom),1), size=(800,length(p_chrom)*200), yaxis=nothing, grid=false)
end

# ╔═╡ 3b03f424-57a4-45c6-9838-5b7635e1278a
function plot_lnk!(p, lnk, Tmin, Tmax; lbl="")
	Trange = Tmin:(Tmax-Tmin)/100.0:Tmax
	lnk_calc = Array{Float64}(undef, length(Trange))
	for i=1:length(Trange)
		lnk_calc[i] = lnk(Trange[i])
	end
	plot!(p, Trange, lnk_calc, label=lbl)
	return p
end

# ╔═╡ 4c9fc541-7947-46db-acd7-502cbda831cb
function Tmin(meas) 	
	Tmin_ = Array{Float64}(undef, length(meas[2]))
	for i=1:length(meas[2])
		Tmin_[i] = minimum(meas[2][i].temp_steps)
	end
	return minimum(Tmin_)
end

# ╔═╡ 940baba7-fba9-4b31-838b-eb7a4504246f
function add_min_max_marker!(p, Tmin, Tmax, lnk)
	scatter!(p, (Tmin, lnk(Tmin)), markershape=:rtriangle, markersize=5, c=:black, label = "measured temperature range")
	scatter!(p, (Tmax, lnk(Tmax)), markershape=:ltriangle, markersize=5, c=:black, label = "measured temperature range")
	return p
end

# ╔═╡ 76c80e90-68e6-4cdb-bb72-785d8c7df9e9
function plotbox(df1)
	p_box1 = boxplot(ylabel="relative difference in %", legend=false)
	#for i=1:length(unique(df1.methode))
	#	df1f = filter([:methode] => x -> x .== unique(df1.methode)[i], df1)
		boxplot!(p_box1, ["Tchar"], df1.relΔTchar.*100.0)
		dotplot!(p_box1, ["Tchar"], df1.relΔTchar.*100.0, marker=(:black, stroke(0)))
		boxplot!(p_box1, ["θchar"], df1.relΔθchar.*100.0)
		dotplot!(p_box1, ["θchar"], df1.relΔθchar.*100.0, marker=(:black, stroke(0)))
		boxplot!(p_box1, ["ΔCp"], df1.relΔΔCp.*100.0)
		dotplot!(p_box1, ["ΔCp"], df1.relΔΔCp.*100.0, marker=(:black, stroke(0)))
	#end
	return p_box1
end	

# ╔═╡ 792c4a8b-bd1c-4769-a907-1e562566a168
function difference_estimation_to_alternative_data(res, db)
	ΔTchar = Array{Float64}(undef, length(res.Name))
	Δθchar = Array{Float64}(undef, length(res.Name))
	ΔΔCp = Array{Float64}(undef, length(res.Name))
	relΔTchar = Array{Float64}(undef, length(res.Name))
	relΔθchar = Array{Float64}(undef, length(res.Name))
	relΔΔCp = Array{Float64}(undef, length(res.Name))
	for i=1:length(res.Name)
		ii = findfirst(GasChromatographySimulator.CAS_identification([res.Name[i]]).CAS[1].==db.CAS)
		if !isnothing(ii)
			ΔTchar[i] = res.Tchar[i] - (db.Tchar[ii] + 273.15) 
			Δθchar[i] = res.θchar[i] - db.thetachar[ii]
			ΔΔCp[i] = res.ΔCp[i] - db.DeltaCp[ii]
			relΔTchar[i] = ΔTchar[i]/(db.Tchar[ii] + 273.15)
			relΔθchar[i] = Δθchar[i]/db.thetachar[ii]
			relΔΔCp[i] = ΔΔCp[i]/db.DeltaCp[ii]
		else
			ΔTchar[i] = NaN
			Δθchar[i] = NaN
			ΔΔCp[i] = NaN
			relΔTchar[i] = NaN
			relΔθchar[i] = NaN
			relΔΔCp[i] = NaN
		end
	end
	diff = DataFrame(Name=res.Name, ΔTchar=ΔTchar, Δθchar=Δθchar, ΔΔCp=ΔΔCp, relΔTchar=relΔTchar, relΔθchar=relΔθchar, relΔΔCp=relΔΔCp)
	return diff
end

# ╔═╡ 5123aa1b-b899-4dd6-8909-6be3b82a80d0
begin
	if isnothing(file_db)
		md"""
		Selected file with alternative parameters: _nothing_

		Please select a file with alternative parameters!
		"""
	else
		db = DataFrame(CSV.File(file_db["data"]))
		diff = difference_estimation_to_alternative_data(res, db)
		pbox = plotbox(diff)
		md"""
		Selected file with alternative parameters: $(file_db["name"])

		$(embed_display(db))

		Difference of parameters to estimation:

		$(embed_display(diff))

		Boxplot of relative differences:

		$(embed_display(pbox))
		"""
	end
end

# ╔═╡ c0f0b955-6791-401f-8252-745332c4210f
begin
	gr()
	p_lnk_ = plot(xlabel=L"T \mathrm{\; in \: °C}", ylabel=L"\ln{k}", title=select_solute, minorticks=4, minorgrid=true)
	if @isdefined isolnk
		scatter!(p_lnk_, isolnk.T.-273.15, isolnk[!, select_solute], label="isothermal meas.")
	end
	if @isdefined db
		i = findfirst(GasChromatographySimulator.CAS_identification([select_solute]).CAS[1].==db.CAS)
		if !isnothing(i)
			lnk_db(T) = (db.DeltaCp[i]/8.31446261815324 + (db.Tchar[i]+273.15)/db.thetachar[i])*((db.Tchar[i]+273.15)/(T+273.15)-1) + db.DeltaCp[i]/8.31446261815324*log((T+273.15)/(db.Tchar[i]+273.15))
			plot_lnk!(p_lnk_, lnk_db, Tmin(meas)*0.925, (Telu_max[findfirst(select_solute.==meas[4])]-273.15)*1.025, lbl="alternative database")
		end
	end
	res_f = filter([:Name] => x -> x == select_solute, res)
	lnk(T) = (res_f.ΔCp[1]/8.31446261815324 + res_f.Tchar[1]/res_f.θchar[1])*(res_f.Tchar[1]/(T+273.15)-1) + res_f.ΔCp[1]/8.31446261815324*log((T+273.15)/res_f.Tchar[1])
	
	plot_lnk!(p_lnk_, lnk, Tmin(meas)*0.925, (Telu_max[findfirst(select_solute.==meas[4])]-273.15)*1.025, lbl=select_mode)
	add_min_max_marker!(p_lnk_, Tmin(meas), Telu_max[findfirst(select_solute.==meas[4])]-273.15, lnk)
	p_lnk_
end

# ╔═╡ Cell order:
# ╠═8a1861d4-4acf-11ed-39b0-7b1af92df0ee
# ╠═59dd4f69-4a2a-4fd6-b679-6f807c297948
# ╠═09422105-a747-40ac-9666-591326850d8f
# ╠═5be3a351-f6e0-4cc6-800d-a55f8decf889
# ╠═eb5fc23c-2151-47fa-b56c-5771a4f8b9c5
# ╟─6d4ec54b-01b2-4208-9b9e-fcb70d236c3e
# ╠═3e808e31-081e-48d2-908e-d66726e270e8
# ╠═ebc2a807-4413-4721-930a-6328ae72a1a9
# ╟─eb14e619-82d4-49ac-ab2a-28e56230dbc6
# ╟─d745c22b-1c96-4a96-83da-abb1df91ab87
# ╠═b2c254a2-a5d6-4f18-803a-75d048fc7cdf
# ╠═365ec346-65e4-4bb8-9b8f-824ed0e5e3f2
# ╟─f3ffd4ce-a378-4033-88e9-bc1fb8cc4bbe
# ╟─e98f4b1b-e577-40d0-a7d8-71c53d99ee1b
# ╟─3b40b0b1-7007-48c7-b47b-dbeaf501b73d
# ╠═992fae9e-5257-4996-840c-7bf8969d3e02
# ╠═f75dcbf4-b146-47ea-8fa2-1cdf3f772aeb
# ╠═d09e48ec-1d38-46c8-b8ce-19d9a83bfa51
# ╠═882b9cf0-7320-408d-b7aa-3e20fffb0ff2
# ╠═f650199a-554d-49ff-8346-f0e5c6819b7c
# ╠═31b176dc-b667-4604-b490-046d55e545c9
# ╠═1ed16860-6899-4a17-aa2e-dffd32a4cf82
# ╠═2c418101-9184-49b7-8ad8-ae3be25d7cdb
# ╠═fcc0c0cf-aee7-42bb-9ad9-790325db31aa
# ╠═4779a2ed-c827-4d4e-aa5f-614aaca6c173
# ╠═cdef487b-1a2b-470a-9f13-e08b88e079a2
# ╠═1f34c58b-5f35-422b-a645-4aa9e96a51f6
# ╠═3ac78538-e2ec-466d-bfe4-d55bcc009db0
# ╠═ee641f1d-4840-4250-9eec-0f0166f49ff2
# ╠═deb2179f-06e5-4439-a3a0-1664fc1a5083
# ╠═93224376-c99e-4182-b430-3d75a52ee0f1
# ╠═c3e0b7e9-f8d3-4b63-9e36-fb4995211a6c
# ╠═07187629-64b7-4869-9eeb-7d530e43587d
# ╠═47e7621e-58e9-4701-912d-82f646e64e0f
# ╠═26672fc0-3b99-44e6-8a46-47dbeb3bca38
# ╠═0e6aa86b-1cd5-4281-b0f1-d8e96b2d6e2e
# ╠═a1c0f53a-925c-4800-a920-955e9da97e12
# ╠═9347878d-66f0-464b-a315-9c6eb46ea43d
# ╠═f408ded9-9130-4652-b0aa-219936f4c947
# ╠═6475ab4f-f701-4398-ac25-b65d0844301a
# ╠═20902ac4-f210-4ddd-90f5-a713eb4bd9b5
# ╟─0f4c35c4-32f7-4d11-874d-1f23daad7da8
# ╟─b4f17579-9994-46e1-a3d0-6030650f0dbe
# ╟─0d61fd05-c0c6-4764-9f96-3b867a456ad3
# ╟─762c877d-3f41-49de-a7ea-0eef1456ac11
# ╟─07a7e45a-d73e-4a83-9323-700d3e2a88cc
# ╠═af2aa109-c514-40dc-a5a9-c7b9974dd8f1
# ╟─ae424251-f4f7-48aa-a72b-3be85a193705
# ╠═116ccf37-4a18-44ac-ae6e-98932901a8b0
# ╠═157f24e6-1664-41c8-9079-b9dd0a2c98a9
# ╠═ddadd79d-7b74-4b2c-ad59-a31a5692cf3b
# ╟─f83b26e7-f16d-49ac-ab3b-230533ac9c82
# ╠═6187b9a9-83e4-49e0-be6f-8aea1a09da7c
# ╟─5123aa1b-b899-4dd6-8909-6be3b82a80d0
# ╟─903d12f1-6d6f-4a71-8095-7457cffcafc4
# ╠═347c53d6-9c94-47de-814a-aef43fa5aae2
# ╟─66287dfc-fe77-4fa8-8892-79d2d3de6cb3
# ╠═b4fc2f6e-4500-45da-852e-59fb084e365a
# ╠═c0f0b955-6791-401f-8252-745332c4210f
# ╟─9178967d-26dc-43be-b6e4-f35bbd0b0b04
# ╠═46fab3fe-cf88-4cbf-b1cc-a232bb7520db
# ╠═fb4ccf8c-18af-4b7d-8bff-50028afdcb6a
# ╠═06d5dd5a-6883-4b62-8942-dc19e7b08e4d
# ╠═505f722f-5429-4282-a1c9-d41fc3284c11
# ╠═72312748-de7d-4b15-8ef8-ea7165d19640
# ╠═3b03f424-57a4-45c6-9838-5b7635e1278a
# ╠═4c9fc541-7947-46db-acd7-502cbda831cb
# ╠═940baba7-fba9-4b31-838b-eb7a4504246f
# ╠═76c80e90-68e6-4cdb-bb72-785d8c7df9e9
# ╠═792c4a8b-bd1c-4769-a907-1e562566a168
