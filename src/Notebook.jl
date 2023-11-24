# functions for the notebook

function download_data(url)
	io = IOBuffer();
	download = urldownload(url, save_raw=io);
	return Dict{Any, Any}("data" => io, "name" => split(url, '/')[end])
end

# filter function for selected measurements and selected solutes
function filter_selected_measurements(meas, selected_measurements, selected_solutes)
	index_measurements = Array{Int}(undef, length(selected_measurements))
	for i=1:length(selected_measurements)
		index_measurements[i] = findfirst(selected_measurements[i].==meas[3].measurement)
	end
	index_solutes = Array{Int}(undef, length(selected_solutes))
	for i=1:length(selected_solutes)
		if isnothing(findfirst(selected_solutes[i].==names(meas[3])))
			error("The names of selected analytes are different from the measurement file.")
		else
			index_solutes[i] = findfirst(selected_solutes[i].==names(meas[3]))
		end
	end
	meas_select = (meas[1], meas[2][index_measurements], meas[3][index_measurements, [1; index_solutes]], selected_solutes, meas[5], meas[6])
	return meas_select
end

# comparing predicted retention times using the optimization result `res` with the measured retention times `comp`
function comparison(res, comp)
	opt = GasChromatographySimulator.Options(ng=true, odesys=false)
	#CAS_comp = GasChromatographySimulator.CAS_identification(comp[4]).CAS
	#CAS_meas = GasChromatographySimulator.CAS_identification(meas[4]).CAS
	i_sub = findall(x->x in res.Name, comp[4]) # indices of common elements of CAS_comp in CAS_meas (indeices relative to CAS_meas)
	sub = Array{GasChromatographySimulator.Substance}(undef, length(i_sub))
	for i in i_sub
		id = GasChromatographySimulator.CAS_identification(res.Name[i])
		#if ismissing(id)
		#	id_C15= GasChromatographySimulator.CAS_identification("pentadecane")
		#	Cag = GasChromatographySimulator.diffusivity(id_C15, comp[1].gas) # use C15 value
		#	CAS = id_C15.CAS
		#else
			Cag = GasChromatographySimulator.diffusivity(id, comp[1].gas)
			CAS = id.CAS
		#end
		if "θchar" in names(res) 
			sub[i] = GasChromatographySimulator.Substance(res.Name[i], CAS, Measurements.value(res.Tchar[i]), Measurements.value(res.θchar[i]), Measurements.value(res.ΔCp[i]), comp[1].df/comp[1].d, "", Cag, 0.0, 0.0)
		else
			sub[i] = GasChromatographySimulator.Substance(res.Name[i], CAS, res.Tchar[i]+273.15, res.thetachar[i], res.DeltaCp[i], comp[1].df/comp[1].d, "", Cag, 0.0, 0.0)
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
			d = mean(Measurements.value.(res.d)) # if d was estimate use the mean value
		elseif @isdefined col_input
			d = col_input.d/1000.0
		else
			d = comp[1].d
		end
		col = GasChromatographySimulator.Column(comp[1].L, d, comp[1].df/comp[1].d*d, comp[1].sp, comp[1].gas)
		par[i] = GasChromatographySimulator.Parameters(col, comp[2][i], sub, opt)
		#try
			#pl[i] = GasChromatographySimulator.simulate(par[i])[1]
			sol, peak = GasChromatographySimulator.solve_multithreads(par[i])
			pl[i] = peaklist(sol, peak, par[i])
		#catch
		#	pl[i] = DataFrame(Name=comp[4], tR=NaN.*ones(length(comp[4])), τR=NaN.*ones(length(comp[4])))
		#end
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

# from GasChromatographySimulator
function peaklist(sol, peak, par)
	n = length(par.sub)
    No = Array{Union{Missing, Int64}}(undef, n)
    Name = Array{String}(undef, n)
    CAS = Array{String}(undef, n)
    tR = Array{Float64}(undef, n)
    TR = Array{Float64}(undef, n)
    σR = Array{Float64}(undef, n)
    uR = Array{Float64}(undef, n)
    τR = Array{Float64}(undef, n)
    kR = Array{Float64}(undef, n)
    Res = fill(NaN, n)
    Δs = fill(NaN, n)
    Annotations = Array{String}(undef, n)
    #Threads.@threads for i=1:n
    for i=1:n
        Name[i] = par.sub[i].name
        CAS[i] = par.sub[i].CAS
        if sol[i].t[end]==par.col.L
            tR[i] = sol[i].u[end]
            TR[i] = par.prog.T_itp(par.col.L, tR[i]) - 273.15 
            uR[i] = 1/GasChromatographySimulator.residency(par.col.L, tR[i], par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.df, par.col.gas, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control, k_th=par.opt.k_th)
            τR[i] = sqrt(peak[i].u[end])
            σR[i] = τR[i]*uR[i]
            kR[i] = GasChromatographySimulator.retention_factor(par.col.L, tR[i], par.prog.T_itp, par.col.d, par.col.df, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀; k_th=par.opt.k_th)
        else
            tR[i] = NaN
            TR[i] = NaN
            uR[i] = NaN
            τR[i] = NaN
            σR[i] = NaN
            kR[i] = NaN
        end
        No[i] = try
            parse(Int,split(par.sub[i].ann, ", ")[end])
        catch
            missing
        end
        if ismissing(No[i])
            Annotations[i] = par.sub[i].ann
        else
            Annotations[i] = join(split(par.sub[i].ann, ", ")[1:end-1], ", ")
        end
    end  
    df = sort!(DataFrame(No = No, Name = Name, CAS = CAS, tR = tR, τR = τR, TR=TR, σR = σR, uR = uR, kR = kR, Annotations = Annotations, ), [:tR])
    #Threads.@threads for i=1:n-1
    for i=1:n-1
        Res[i] = (df.tR[i+1] - df.tR[i])/(2*(df.τR[i+1] + df.τR[i]))
        Δs[i] = (df.tR[i+1] - df.tR[i])/(df.τR[i+1] - df.τR[i]) * log(df.τR[i+1]/df.τR[i])
    end
    df[!, :Res] = Res 
    df[!, :Δs] = Δs   
    
    return select(df, [:No, :Name, :CAS, :tR, :τR, :TR, :σR, :uR, :kR, :Res, :Δs, :Annotations])
end

function plot_chromatogram_comparison(pl, meas, comp)
	gr()
	p_chrom = Array{Plots.Plot}(undef, length(pl))
	for i=1:length(pl)
		p_chrom[i] = plot(xlabel="time in $(meas[6])", legend=false)

		max_ = maximum(GasChromatographySimulator.plot_chromatogram(pl[i], (minimum(pl[i].tR)*0.95, maximum(pl[i].tR)*1.05))[3])*1.05

		min_ = - max_/20
		
		GasChromatographySimulator.plot_chromatogram!(p_chrom[i], pl[i], (minimum(pl[i].tR)*0.95, maximum(pl[i].tR)*1.05); annotation=false)
		xlims!((minimum(pl[i].tR)*0.95, maximum(pl[i].tR)*1.05))
		ylims!(min_, max_)
		#add marker for measured retention times
		for j=1:length(comp[4])
			plot!(p_chrom[i], comp[3][i,j+1].*ones(2), [min_, max_], c=:orange)
		end
		plot!(p_chrom[i], title=comp[3].measurement[i])
	end
	
	return plot(p_chrom..., layout=(length(p_chrom),1), size=(800,length(p_chrom)*200), yaxis=nothing, grid=false)
end	

function plot_lnk!(p, lnk, Tmin, Tmax; lbl="")
	Trange = Tmin:(Tmax-Tmin)/100.0:Tmax
	lnk_calc = Array{Float64}(undef, length(Trange))
	for i=1:length(Trange)
		lnk_calc[i] = lnk(Trange[i])
	end
	plot!(p, Trange, lnk_calc, label=lbl)
	return p
end

function Tmin(meas) 	
	Tmin_ = Array{Float64}(undef, length(meas[2]))
	for i=1:length(meas[2])
		Tmin_[i] = minimum(meas[2][i].temp_steps)
	end
	return minimum(Tmin_)
end

function add_min_max_marker!(p, Tmin, Tmax, lnk)
	scatter!(p, (Tmin, lnk(Tmin)), markershape=:rtriangle, markersize=5, c=:black, label = "measured temperature range")
	scatter!(p, (Tmax, lnk(Tmax)), markershape=:ltriangle, markersize=5, c=:black, label = "measured temperature range")
	return p
end

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

function difference_estimation_to_alternative_data(res, db)
	ΔTchar = Array{Float64}(undef, length(res.Name))
	Δθchar = Array{Float64}(undef, length(res.Name))
	ΔΔCp = Array{Float64}(undef, length(res.Name))
	relΔTchar = Array{Float64}(undef, length(res.Name))
	relΔθchar = Array{Float64}(undef, length(res.Name))
	relΔΔCp = Array{Float64}(undef, length(res.Name))
	for i=1:length(res.Name)
		ii = findfirst(GasChromatographySimulator.CAS_identification(res.Name[i]).CAS.==db.CAS)
		if !isnothing(ii)
			ΔTchar[i] = Measurements.value(res.Tchar[i]) - (db.Tchar[ii] + 273.15) 
			Δθchar[i] = Measurements.value(res.θchar[i]) - db.thetachar[ii]
			ΔΔCp[i] = Measurements.value(res.ΔCp[i]) - db.DeltaCp[ii]
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