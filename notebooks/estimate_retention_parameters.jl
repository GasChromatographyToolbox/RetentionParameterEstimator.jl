### A Pluto.jl notebook ###
# v0.19.41

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

# ╔═╡ 09422105-a747-40ac-9666-591326850d8f
begin 
	# online version
	import Pkg
	version = "0.1.8"
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec(name="RetentionParameterEstimator", version=version)
		])
	using RetentionParameterEstimator
	md"""
	online, Packages, estimate\_retention\_parameters, for RetentionParameterEstimator v$(version)
	"""

#=	# local version (database is still downloaded from github)

	import Pkg
	# activate the shared project environment
	Pkg.activate(Base.current_project())
	using RetentionParameterEstimator
	#Pkg.add([
	#	Pkg.PackageSpec(name="RetentionParameterEstimator", rev="minor_fix_jsb")
	#])
	
	md"""
	local, Packages, estimate\_retention\_parameters.jl, for RetentionParameterEstimator v0.1.6
	"""
=#
end

# ╔═╡ eb5fc23c-2151-47fa-b56c-5771a4f8b9c5
begin
	html"""<style>
	main {
	    max-width: 75%;
	    margin-left: 1%;
	    margin-right: 20% !important;
	}
	"""
end

# ╔═╡ f46b165e-67d9-402f-a225-72d1082007be
TableOfContents()

# ╔═╡ 6d4ec54b-01b2-4208-9b9e-fcb70d236c3e
md"""
# 
$(Resource("https://raw.githubusercontent.com/JanLeppert/RetentionParameterEstimator.jl/main/docs/src/assets/logo_b.svg"))
Estimation of K-centric retention parameters by temperature programmed GC with different temperature programs.
"""

# ╔═╡ ebc2a807-4413-4721-930a-6328ae72a1a9
md"""
## Select measurement data
Measured chromatograms used to **estimate** the parameters by optimization:

Load own data: $(@bind own_data CheckBox(default=false))
"""

# ╔═╡ 51a22a15-24f9-4280-9c91-32e48727003a
if own_data == true
	md"""
	$(@bind file_meas FilePicker([MIME("text/csv")]))
	"""
else
	file_meas = RetentionParameterEstimator.download_data("https://raw.githubusercontent.com/JanLeppert/RetentionParameterEstimator.jl/main/data/meas_df05_Rxi5SilMS.csv");
end

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

# ╔═╡ b2c254a2-a5d6-4f18-803a-75d048fc7cdf
meas_select = RetentionParameterEstimator.filter_selected_measurements(meas, selected_measurements, selected_solutes);

# ╔═╡ f3ffd4ce-a378-4033-88e9-bc1fb8cc4bbe
md"""
## Select mode
* `m1` ... estimate the three retention parameters (``T_{char}``, ``θ_{char}`` and ``ΔC_p``).
* `m1a` ... estimate the three retention parameters (``T_{char}``, ``θ_{char}`` and ``ΔC_p``) and select ``L`` and ``d``.
* `m2` ... estimate the three retention parameters (``T_{char}``, ``θ_{char}`` and ``ΔC_p``) and the column diameter ``d``.
$(@bind select_mode confirm(Select(["m1", "m1a", "m2"])))
"""

# ╔═╡ e98f4b1b-e577-40d0-a7d8-71c53d99ee1b
if select_mode == "m1a"
	#md"""
	#Column parameters:

	#L in m: $(@bind L_input NumberField(0.0:0.01:1000.0; default=meas[1].L))

	#d in mm: $(@bind d_input NumberField(0.0:0.0001:1.0; default=meas[1].d*1000.0))

	@bind col_input confirm(
		PlutoUI.combine() do Child
		md"""
		## Column dimensions
		
		L in m: $(Child("L", NumberField(0.0:0.01:1000.0; default=meas_select[1].L)))

		d in mm: $(Child("d", NumberField(0.0:0.0001:1.0; default=meas_select[1].d*1000.0))) 
		"""
		end
	)
end	

# ╔═╡ 3b40b0b1-7007-48c7-b47b-dbeaf501b73d
begin
	min_th = 0.1
	loss_th = 1.0
	if select_mode == "m1"
		check, msg, df_flag, index_flag, res, Telu_max = RetentionParameterEstimator.check_measurement(meas_select, (L = meas_select[1].L, d = meas_select[1].d*1000); min_th=min_th, loss_th=loss_th, se_col=false)
		md"""
		## Results

		L/d ratio: $(meas_select[1].L/(meas_select[1].d))
		
		$(res)
		"""
	elseif select_mode == "m1a"
		check, msg, df_flag, index_flag, res, Telu_max = RetentionParameterEstimator.check_measurement(meas_select, col_input; min_th=min_th, loss_th=loss_th, se_col=false)
		md"""
		## Results

		L/d ratio: $(col_input[1]/(col_input[2]*1e-3))
		
		$(res)
		"""
	elseif select_mode == "m2"
		res, Telu_max = RetentionParameterEstimator.method_m2(meas_select; se_col=false)
		check = true
		md"""
		## Results

		d = $(1000.0 * res.d[1]) mm
		
		L/d ratio: $(meas_select[1].L./(res.d[1]))
		
		$(res)
		"""
	end
end
# d = $(1000.0 * (res.d[1] ± res.d_std[1])) mm

# ╔═╡ 8cc151a6-a60a-4aba-a813-1a142a073948
begin
	if check == false
		md"""
		## Results
		
		!!! warning
		
		$(msg)

		The found minima for solutes $(meas[4][index_flag]) is above the threshold $(min_th) s².

		$(df_flag)
		
		Check for the plausibility of measurements, if a solute is listed there, all retention times for all used measurements should be checked, not only the listed flagged ones. An error in one retention time, e.g. typo, could produce deviations from the model in other (correct) measurements.

		$(res)
		"""
	end
end

# ╔═╡ 0f4c35c4-32f7-4d11-874d-1f23daad7da8
begin
	head = ["file:", file_meas["name"], 
		"selected measurements:", join(selected_measurements, " "),
		"mode:", select_mode
	]
	res_export = RetentionParameterEstimator.separate_error_columns(res)
	io = IOBuffer()
	CSV.write(io, DataFrame[], header=head)
	CSV.write(io, res_export, append=true, writeheader=true)
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

# ╔═╡ 0d61fd05-c0c6-4764-9f96-3b867a456ad3
md"""
## Select verification data

Measured chromatograms used to **verify** the estimated parameters:
"""

# ╔═╡ 38d1e196-f375-48ac-bc11-80b10472c1cd
if own_data == true
	md"""
	$(@bind file_comp FilePicker([MIME("text/csv")]))
	"""
else
	file_comp = RetentionParameterEstimator.download_data("https://raw.githubusercontent.com/JanLeppert/RetentionParameterEstimator.jl/main/data/comp_df05_Rxi5SilMS.csv");
end

# ╔═╡ 762c877d-3f41-49de-a7ea-0eef1456ac11
begin
	if isnothing(file_comp)
		md"""
		Selected chromatograms for **verification**: _nothing_

		Please select a file of chromatograms for **verification**!
		"""
	else
		if file_meas == file_comp
			comp = RetentionParameterEstimator.load_chromatograms(file_comp);
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

# ╔═╡ ae424251-f4f7-48aa-a72b-3be85a193705
md"""
## Verification

Using the estimated parameters (``T_{char}``, ``θ_{char}``, ``ΔC_p``) resp. (``T_{char}``, ``θ_{char}``, ``ΔC_p``, ``d``) to simulate other temperature programs. Compare the simulated retention times with measured retention times.
"""

# ╔═╡ 116ccf37-4a18-44ac-ae6e-98932901a8b0
begin
	comp_select = RetentionParameterEstimator.filter_selected_measurements(comp, selected_comparison, selected_solutes)
	pl, loss, par = RetentionParameterEstimator.comparison(res, comp_select)
	meanΔtR = Array{Float64}(undef, length(pl))
	meanrelΔtR = Array{Float64}(undef, length(pl))
	for i=1:length(pl)
		meanΔtR[i] = mean(abs.(pl[i].ΔtR))
		meanrelΔtR[i] = mean(abs.(pl[i].relΔtR))
	end
	df_loss = DataFrame(comp=selected_comparison, loss=loss, sqrt_loss=sqrt.(loss), mean_ΔtR=meanΔtR, mean_relΔtR_percent=meanrelΔtR.*100.0)
end

# ╔═╡ 157f24e6-1664-41c8-9079-b9dd0a2c98a9
md"""
Peaklists of the simulated verification programs, including difference to measured retention times.

$(embed_display(pl))
"""

# ╔═╡ f6a55039-8c32-4d21-8119-fc48e3bf0dc2
md"""
**Options for plot**

add labels: $(@bind labels CheckBox(default=true))
"""

# ╔═╡ e16e8b29-0f85-43b6-a405-6ab60b0e3ee0
begin
	if labels==true
		label_db = DataFrame(Name=comp[4])
		md"""
		
		$(embed_display(RetentionParameterEstimator.plot_chromatogram_comparison(pl, meas_select, comp_select; annotation=labels)))
		
		* Simulated chromatograms in **blue**.
		* Measured retention times in **orange**. 
		* Labels: $(embed_display(label_db)) 
		"""
	else
		md"""
		
		$(embed_display(RetentionParameterEstimator.plot_chromatogram_comparison(pl, meas_select, comp_select; annotation=labels)))
		
		* Simulated chromatograms in **blue**.
		* Measured retention times in **orange**. 
		"""
	end
end

# ╔═╡ f83b26e7-f16d-49ac-ab3b-230533ac9c82
md"""
## Alternative parameters

Select file with alternative parameters:
"""

# ╔═╡ 62c014d5-5e57-49d7-98f8-dace2f1aaa32
if own_data == true
	md"""
	$(@bind file_db FilePicker([MIME("text/csv")]))
	"""
else
	file_db = RetentionParameterEstimator.download_data("https://raw.githubusercontent.com/JanLeppert/RetentionParameterEstimator.jl/main/data/database_Rxi5SilMS_beta125.csv");
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
		diff = RetentionParameterEstimator.difference_estimation_to_alternative_data(res, db)
		pbox = RetentionParameterEstimator.plotbox(diff)
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

# ╔═╡ 903d12f1-6d6f-4a71-8095-7457cffcafc4
md"""
## Isothermal ``\ln{k}`` data

Select file with isothermal measured ``\ln{k}`` data:
"""

# ╔═╡ a41148bc-03b0-4bd1-b76a-7406aab63f48
if own_data == true
	md"""
	$(@bind file_isolnk FilePicker([MIME("text/csv")]))
	"""
else
	file_isolnk = RetentionParameterEstimator.download_data("https://raw.githubusercontent.com/JanLeppert/RetentionParameterEstimator.jl/main/data/isothermal_lnk_df05_Rxi5SilMS.csv");
end

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
$(@bind select_solute confirm(Select(meas_select[4]; default=meas_select[4][1])))
"""

# ╔═╡ c0f0b955-6791-401f-8252-745332c4210f
begin
	gr()
	#p_lnk_ = plot(xlabel=L"T \mathrm{\; in \: °C}", ylabel=L"\ln{k}", title=select_solute, minorticks=4, minorgrid=true)
	p_lnk_ = plot(xlabel="T in °C", ylabel="lnk", title=select_solute, minorticks=4, minorgrid=true)
	if @isdefined isolnk
		scatter!(p_lnk_, isolnk.T.-273.15, isolnk[!, select_solute], label="isothermal meas.")
	end
	if @isdefined db
		i = findfirst(RetentionParameterEstimator.GasChromatographySimulator.CAS_identification(select_solute).CAS.==db.CAS)
		if !isnothing(i)
			lnk_db(T) = (db.DeltaCp[i]/8.31446261815324 + (db.Tchar[i]+273.15)/db.thetachar[i])*((db.Tchar[i]+273.15)/(T+273.15)-1) + db.DeltaCp[i]/8.31446261815324*log((T+273.15)/(db.Tchar[i]+273.15))
			RetentionParameterEstimator.plot_lnk!(p_lnk_, lnk_db, RetentionParameterEstimator.Tmin(meas)*0.925, (Telu_max[findfirst(select_solute.==meas[4])]-273.15)*1.025, lbl="alternative database")
		end
	end
	res_f = filter([:Name] => x -> x == select_solute, res)
	lnk(T) = (Measurements.value.(res_f.ΔCp[1])/8.31446261815324 + Measurements.value.(res_f.Tchar[1])/Measurements.value.(res_f.θchar[1]))*(Measurements.value.(res_f.Tchar[1])/(T+273.15)-1) + Measurements.value.(res_f.ΔCp[1])/8.31446261815324*log((T+273.15)/Measurements.value.(res_f.Tchar[1]))
	
	RetentionParameterEstimator.plot_lnk!(p_lnk_, lnk, RetentionParameterEstimator.Tmin(meas)*0.925, (Telu_max[findfirst(select_solute.==meas[4])]-273.15)*1.025, lbl=select_mode)
	RetentionParameterEstimator.add_min_max_marker!(p_lnk_, RetentionParameterEstimator.Tmin(meas), Telu_max[findfirst(select_solute.==meas[4])]-273.15, lnk)
	p_lnk_
end

# ╔═╡ b6c2ad4d-6fc5-4700-80e3-f616c0b9aa91
begin
	gr()
	#p_lnk_all = plot(xlabel=L"T \mathrm{\; in \: °C}", ylabel=L"\ln{k}", title="all", minorticks=4, minorgrid=true, legend=:none)
	p_lnk_all = plot(xlabel="T in °C", ylabel="lnk", title="all", minorticks=4, minorgrid=true, legend=:none)

	for i=1:length(meas_select[4])
		res_ = filter([:Name] => x -> x == meas_select[4][i], res)
		lnk(T) = (Measurements.value.(res_.ΔCp[1])/8.31446261815324 + Measurements.value.(res_.Tchar[1])/Measurements.value.(res_.θchar[1]))*(Measurements.value.(res_.Tchar[1])/(T+273.15)-1) + Measurements.value.(res_.ΔCp[1])/8.31446261815324*log((T+273.15)/Measurements.value.(res_.Tchar[1]))
	
		RetentionParameterEstimator.plot_lnk!(p_lnk_all, lnk, RetentionParameterEstimator.Tmin(meas_select)*0.925, (Telu_max[i]-273.15)*1.025, lbl=meas_select[4][i])
		RetentionParameterEstimator.add_min_max_marker!(p_lnk_all, RetentionParameterEstimator.Tmin(meas_select), Telu_max[i]-273.15, lnk)
	end
	p_lnk_all
end

# ╔═╡ 9178967d-26dc-43be-b6e4-f35bbd0b0b04
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═09422105-a747-40ac-9666-591326850d8f
# ╟─eb5fc23c-2151-47fa-b56c-5771a4f8b9c5
# ╟─f46b165e-67d9-402f-a225-72d1082007be
# ╟─6d4ec54b-01b2-4208-9b9e-fcb70d236c3e
# ╟─ebc2a807-4413-4721-930a-6328ae72a1a9
# ╟─51a22a15-24f9-4280-9c91-32e48727003a
# ╟─eb14e619-82d4-49ac-ab2a-28e56230dbc6
# ╟─d745c22b-1c96-4a96-83da-abb1df91ab87
# ╟─b2c254a2-a5d6-4f18-803a-75d048fc7cdf
# ╟─f3ffd4ce-a378-4033-88e9-bc1fb8cc4bbe
# ╟─e98f4b1b-e577-40d0-a7d8-71c53d99ee1b
# ╟─3b40b0b1-7007-48c7-b47b-dbeaf501b73d
# ╟─8cc151a6-a60a-4aba-a813-1a142a073948
# ╟─0f4c35c4-32f7-4d11-874d-1f23daad7da8
# ╟─b4f17579-9994-46e1-a3d0-6030650f0dbe
# ╟─0d61fd05-c0c6-4764-9f96-3b867a456ad3
# ╟─38d1e196-f375-48ac-bc11-80b10472c1cd
# ╟─762c877d-3f41-49de-a7ea-0eef1456ac11
# ╟─07a7e45a-d73e-4a83-9323-700d3e2a88cc
# ╟─ae424251-f4f7-48aa-a72b-3be85a193705
# ╟─116ccf37-4a18-44ac-ae6e-98932901a8b0
# ╟─157f24e6-1664-41c8-9079-b9dd0a2c98a9
# ╟─f6a55039-8c32-4d21-8119-fc48e3bf0dc2
# ╟─e16e8b29-0f85-43b6-a405-6ab60b0e3ee0
# ╟─f83b26e7-f16d-49ac-ab3b-230533ac9c82
# ╟─62c014d5-5e57-49d7-98f8-dace2f1aaa32
# ╟─5123aa1b-b899-4dd6-8909-6be3b82a80d0
# ╟─903d12f1-6d6f-4a71-8095-7457cffcafc4
# ╟─a41148bc-03b0-4bd1-b76a-7406aab63f48
# ╟─66287dfc-fe77-4fa8-8892-79d2d3de6cb3
# ╟─b4fc2f6e-4500-45da-852e-59fb084e365a
# ╟─c0f0b955-6791-401f-8252-745332c4210f
# ╟─b6c2ad4d-6fc5-4700-80e3-f616c0b9aa91
# ╟─9178967d-26dc-43be-b6e4-f35bbd0b0b04
