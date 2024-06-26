# functions for loading the measured retention data and the necessary informations (column definition, programs, ...)

# this function is a modified version from GasChromatographySimulator
# if it works, put this function into GasChromatographySimulator
"""
    Program(TP, FpinP, L; pout="vacuum", time_unit="min")

Construct the structure `Program` with conventional formulation (see function `conventional_program` in `GasChromatographySimulator.jl`) of programs for the case
without a thermal gradient. 

# Arguments
* `TP`: conventional formulation of a temperature program. 
* `FpinP`: conventional formulation of a Flow (in m³/s) resp. inlet pressure (in Pa(a)) program.
* `L`: Length of the capillary measured in m (meter).
* `pout`: Outlet pressure, "vacuum" (default), "atmosphere" or the outlet pressure in Pa(a).
* `time_unit`: unit of time in the programs, "min"` (default) times are measured in minutes, "s" times are measured in seconds.

The argument `L` is used to construct the temperature interpolation `T_itp(x,t)`.

# Examples
```julia
julia> Program((40.0, 1.0, 5.0, 280.0, 2.0, 20.0, 320.0, 2.0),
                (400000.0, 10.0, 5000.0, 500000.0, 20.0),
                10.0)
```
"""
function Program(TP, FpinP, L; pout="vacuum", time_unit="min")
    ts1, Ts = GasChromatographySimulator.conventional_program(TP; time_unit=time_unit)
    ts2, Fps = GasChromatographySimulator.conventional_program(FpinP; time_unit=time_unit)
    # remove additional 0.0 which are not at the first position
    ts1_ = ts1[[1; findall(0.0.!=ts1)]]
	Ts_ = Ts[[1; findall(0.0.!=ts1)]]
	ts2_ = ts2[[1; findall(0.0.!=ts2)]]
	Fps_ = Fps[[1; findall(0.0.!=ts2)]]
    time_steps = GasChromatographySimulator.common_time_steps(ts1_, ts2_)
    temp_steps = GasChromatographySimulator.new_value_steps(Ts_, ts1_, time_steps)
    Fpin_steps = GasChromatographySimulator.new_value_steps(Fps_, ts2_, time_steps)
    if pout == "vacuum"
        pout_steps = zeros(length(time_steps))
    elseif isa(pout, Number)
        pout_steps = pout.*ones(length(time_steps)) 
    else
        pout_steps = pn.*ones(length(time_steps))
    end
    a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) L.*ones(length(time_steps)) zeros(length(time_steps))]
    gf(x) = GasChromatographySimulator.gradient(x, a_gf)
    T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, L)
    Fpin_itp = GasChromatographySimulator.steps_interpolation(time_steps, Fpin_steps)
    pout_itp = GasChromatographySimulator.steps_interpolation(time_steps, pout_steps)
    prog = GasChromatographySimulator.Program(time_steps, temp_steps, Fpin_steps, pout_steps, gf, a_gf, T_itp, Fpin_itp, pout_itp)
    return prog
end

function extract_temperature_and_pressure_programs(TPprog)
    iP = findall(occursin.("p", names(TPprog))) # column index of pressure information
    pamb = TPprog[!, iP[end]]
    #iT = 2:(iP[1]-1) # column index of temperature program information
    TPs = TPprog[!,1:(iP[1]-1)] 
    PPs = DataFrame(measurement = TPs.measurement, p1 = TPprog[!,iP[1]].+pamb, t1 = TPs.t1) # pressure program in Pa(a)
    iTrate = findall(occursin.("RT", names(TPs))) # column index of heating rates 
    heatingrates = Array{Array{Union{Missing, Float64}}}(undef, length(iTrate))
    Tdiff = Array{Array{Union{Missing, Float64}}}(undef, length(iTrate))
    for i=1:length(iTrate)
        heatingrates[i] = TPs[!, iTrate[i]]
        Tdiff[i] = TPs[!, iTrate[i]+1] .- TPs[!, iTrate[i]-2]
    end
    pressurerates = Array{Array{Union{Missing, Float64}}}(undef, length(iTrate))
    for i=1:length(iTrate)
        pressurerates[i] = (TPprog[!,iP[i+1]] .- TPprog[!,iP[i]]) .* heatingrates[i] ./ Tdiff[i]
        if !ismissing(pressurerates[i] != zeros(length(pressurerates[i])))
            PPs[!,"RP$(i)"] = pressurerates[i]
            PPs[!,"p$(i+1)"] = TPprog[!, iP[i+1]].+pamb
            PPs[!,"t$(i+1)"] = TPs[!,"t$(i+1)"]
        else
            PPs[!,"t$(i)"] = TPs[!,"t$(i+1)"] .+ Tdiff[i] ./ heatingrates[i] 
        end
    end 
    PPs[!,"pamb"] = pamb
    return TPs, PPs
end

function extract_measured_program(TPprog, path, L)
    # load the measured programs
    prog = Array{GasChromatographySimulator.Program}(undef, length(TPprog.filename))
    for i=1:length(TPprog.filename)
        meas_prog = DataFrame(CSV.File(joinpath(path, TPprog.filename[i]), silencewarnings=true))
        prog[i] = GasChromatographySimulator.Program(meas_prog.timesteps, meas_prog.tempsteps, meas_prog.pinsteps, meas_prog.poutsteps, L)
    end
    return prog
end

"""
    load_chromatograms(file; filter_missing=true)

Loading of the chromatographic data (column information, GC program information, retention time information, see also "Structure of input data") from a file.

# Arguments
* `file` ... path to the file.
* `filter_missing` ... option to ignore solutes, where some retention times are not given. (`default = true`).

# Output
A tuple of the following quantities:
* `col` ... settings of the column as `GasChromatographySimulator.Column` structure.
* `prog` ... Array of the GC programs as `GasChromatographySimulator.Program` structure.
* `tRs` ... DataFrame of the retention times.
* `solute_names` ... Vector of the solute names.
* `pout` ... outlet pressure (detector pressure), "vacuum" or "atmospheric". 
* `time_unit` ... unit of time scale used in the retention times and GC programs, "min" or "s".
"""
function load_chromatograms(file; filter_missing=true, delim=";") # new version -> check case for `filter_missing=false` and missing values in further functions!!!!
    n = open(f->countlines(f), file)
    #col_df = DataFrame(CSV.File(file, header=1, limit=1, stringtype=String, silencewarnings=true, delim=delim, types=[Float64, Float64, Float64, String, String, String, String]))
    col_df = DataFrame(CSV.File(file, header=1, limit=1, stringtype=String, silencewarnings=true, delim=delim, types=Dict("L" => Float64, "d" => Float64, "df" => Float64)))
	col = GasChromatographySimulator.Column(col_df.L[1], col_df.d[1], col_df.df[1], col_df.sp[1], col_df.gas[1])
    pout = col_df.pout[1]
	time_unit = col_df.time_unit[1]

    n_meas = Int((n - 2 - 2)/2) 
    TPprog = DataFrame(CSV.File(file, header=3, limit=n_meas, stringtype=String, silencewarnings=true, delim=delim))
    #PP = DataFrame(CSV.File(file, header=3+n_meas+1, limit=n_meas, stringtype=String)) # convert pressures from Pa(g) to Pa(a), add p_atm to this data set
    tRs_ = DataFrame(CSV.File(file, header=n-n_meas, stringtype=String, silencewarnings=true, delim=delim))
    solute_names_ = names(tRs_)[2:end] # filter non-solute names out (columnx)
    filter!(x -> !occursin.("Column", x), solute_names_)
    if names(TPprog)[2] == "filename"
        path = dirname(file)
        prog = extract_measured_program(TPprog, path, col.L)
    else 
        TPs, PPs = extract_temperature_and_pressure_programs(TPprog)

        prog = Array{GasChromatographySimulator.Program}(undef, length(TPs.measurement))
        for i=1:length(TPs.measurement)
            if pout == "atmospheric"
                pout_ = PPs[i, end]
            else
                pout_ = "vacuum"
            end
            prog[i] = Program(collect(skipmissing(TPs[i, 2:end])), collect(skipmissing(PPs[i, 2:(end-1)])), col.L; pout=pout_, time_unit=time_unit)
        end
    end

	# filter out substances with retention times missing
	if filter_missing == true
		solute_names = solute_names_[findall((collect(any(ismissing, c) for c in eachcol(tRs_))).==false)[2:end].-1]
		tRs = tRs_[!,findall((collect(any(ismissing, c) for c in eachcol(tRs_))).==false)]
	else
		solute_names = solute_names_ # what to do with missing entries for the optimization????
		tRs = tRs_
	end
	
    return col, prog, tRs[!,1:(length(solute_names)+1)], solute_names, pout, time_unit#, TPs, PPs
end

"""
    load_chromatograms(file::Dict{Any, Any}; filter_missing=true, path=joinpath(dirname(pwd()), "data", "exp_pro"))

Loading of the chromatographic data (column information, GC program information, retention time information, see also "Structure of input data") from a file selected by the FilePicker in a Pluto notebook.

# Arguments
* `file` ... file dictionary from the FilePicker.
* `filter_missing` ... option to ignore solutes, where some retention times are not given. (`default = true`).
* `path` ... if the temperature programs are defined by measured temperatures over time, define the path to these files.

# Output
A tuple of the following quantities:
* `col` ... settings of the column as `GasChromatographySimulator.Column` structure.
* `prog` ... Array of the GC programs as `GasChromatographySimulator.Program` structure.
* `tRs` ... DataFrame of the retention times.
* `solute_names` ... Vector of the solute names.
* `pout` ... outlet pressure (detector pressure), "vacuum" or "atmospheric". 
* `time_unit` ... unit of time scale used in the retention times and GC programs, "min" or "s".
"""
function load_chromatograms(file::Dict{Any, Any}; filter_missing=true, path=joinpath(dirname(pwd()), "data", "exp_pro"), delim=";") # if file is the output of FilePicker()
    n = length(CSV.File(file["data"]; silencewarnings=true, comment=";;"))+1
    col_df = DataFrame(CSV.File(file["data"], header=1, limit=1, stringtype=String, silencewarnings=true, delim=delim))
	col = GasChromatographySimulator.Column(convert(Float64, col_df.L[1]), col_df.d[1], col_df.df[1], col_df.sp[1], col_df.gas[1])
    pout = col_df.pout[1]
	time_unit = col_df.time_unit[1]
 
    n_meas = Int((n - 2 - 2)/2) 
    TPprog = DataFrame(CSV.File(file["data"], header=3, limit=n_meas, stringtype=String, silencewarnings=true, delim=delim))
    #PP = DataFrame(CSV.File(file, header=3+n_meas+1, limit=n_meas, stringtype=String)) # convert pressures from Pa(g) to Pa(a), add p_atm to this data set
    tRs_ = DataFrame(CSV.File(file["data"], header=n-n_meas, stringtype=String, silencewarnings=true, comment=";;", delim=delim))
    solute_names_ = names(tRs_)[2:end] # filter non-solute names out (columnx)
    filter!(x -> !occursin.("Column", x), solute_names_)
    if names(TPprog)[2] == "filename"
        #path = dirname(file)
        prog = extract_measured_program(TPprog, path, col.L)
    else 
        TPs, PPs = extract_temperature_and_pressure_programs(TPprog)

        prog = Array{GasChromatographySimulator.Program}(undef, length(TPs.measurement))
        for i=1:length(TPs.measurement)
            if pout == "vacuum"
                pout_ = "vacuum"
            else # pout="atmospheric"
                pout_ = PPs[i, end]
            end
            prog[i] = Program(collect(skipmissing(TPs[i, 2:end])), collect(skipmissing(PPs[i, 2:(end-1)])), col.L; pout=pout_, time_unit=time_unit)
        end
    end

	# filter out substances with retention times missing
	if filter_missing == true
		solute_names = solute_names_[findall((collect(any(ismissing, c) for c in eachcol(tRs_))).==false)[2:end].-1]
		tRs = tRs_[!,findall((collect(any(ismissing, c) for c in eachcol(tRs_))).==false)]
	else
		solute_names = solute_names_
		tRs = tRs_
	end
	
    return col, prog, tRs[!,1:(length(solute_names)+1)], solute_names, pout, time_unit#, TPs, PPs
end