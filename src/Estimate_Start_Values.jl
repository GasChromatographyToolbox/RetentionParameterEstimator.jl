# functions used to estimate start values of the parameters

"""
    reference_holdup_time(prog, L, d, gas; control="Pressure")

Calculate the reference holdup time for the GC program `prog` for a column with length `L` and diameter `d` and `gas` as mobile phase.
"""
function reference_holdup_time(prog, L, d, gas; control="Pressure")
    Tref = 150.0
    # estimate the time of the temperature program for T=Tref
    t_ = prog.time_steps
    T_ = prog.temp_steps
    interp = interpolate((T_ .- Tref, ), cumsum(t_), Gridded(Linear()))
	tref = interp(0.0)
    # inlet and outlet pressure at time tref
    Fpin_ref = prog.Fpin_itp(tref)
    pout_ref = prog.pout_itp(tref)
    # hold-up time calculated for the time of the program, when T=Tref
    tMref = GasChromatographySimulator.holdup_time(Tref+273.15, Fpin_ref, pout_ref, L, d, gas, control=control)
    return tMref
end

"""
    elution_temperature(tRs, prog)

Calculate the elution temperatures from retention times `tRs` and GC programs `prog`.
"""
function elution_temperature(tRs, prog)
    Telus = Array{Float64}(undef, length(tRs))
    for i=1:length(tRs)
        Telus[i] = prog.T_itp(0.0, tRs[i])
    end
    return Telus
end

#="""
    estimate_start_parameter(tR_meas, TPs, PPs, L, d, gas; pout="vacuum", time_unit="min", control="Pressure")

Estimation of start values for the K-centric parameters `Tchar`, `θchar` and `ΔCp` from measured retention times `tR_meas` for several GC programs 
with temperature programs `TPs` and pressure programs `PPs` for a column with length `L` and diameter `d` and `gas` as mobile phase. The programs `TPs` and `PPs` 
are defined as conventional programs (Vector of "start value", "hold time of start value", "ramp", "end value", "hold time of end value"). **The temperature program
should have only one ramp**.
"""   
function estimate_start_parameter(tR_meas, TPs, PPs, L, d, gas; pout="vacuum", time_unit="min", control="Pressure")
    if time_unit == "min"
        t_conv = 60.0
    else
        t_conv = 1.0
    end
    prog = Array{GasChromatographySimulator.Program}(undef, length(TPs))
    tMref = Array{Float64}(undef, length(TPs))
    RT = Array{Float64}(undef, length(TPs))
    if typeof(size(tR_meas)) == Tuple{Int64, Int64}
        n2 = size(tR_meas)[2]
    else
        n2 = 1
    end 
    Telu_meas = Array{Float64}(undef, size(tR_meas)[1], n2)
    for i=1:length(TPs)
        if pout == "atmospheric"
            pout_ = PPs[i][end]
        else
            pout_ = "vacuum"
        end
        prog[i] = Program(TPs[i], PPs[i][1:(end-1)], L; pout=pout_, time_unit=time_unit)
        tMref[i] = reference_holdup_time(prog[i], L, d, gas; control=control)/t_conv
        RT[i] = TPs[i][3] # single-ramp temperature programs are assumed
        Telu_meas[i,:] = elution_temperature(tR_meas[i,:].*t_conv, prog[i])
    end 
    rT = RT.*tMref./θref
    Telu_max = Array{Float64}(undef, n2)
    Tchar_est = Array{Float64}(undef, n2)
    θchar_est = Array{Float64}(undef, n2)
    ΔCp_est = Array{Float64}(undef, n2)
    for i=1:n2
        interp = interpolate((rT, ), Telu_meas[:,i], Gridded(Linear()))
        #spl = Spline1D(rT, Telu_meas[:,i])
        #Tchar_est[i] = spl(rT_nom)
        Telu_max[i] = maximum(Telu_meas[:,i])
        Tchar_est[i] = interp(rT_nom)
        θchar_est[i] = 22.0*(Tchar_est[i]/Tst)^0.7 # factor of φ?
        ΔCp_est[i] = 100.0
    end
    return Tchar_est, θchar_est, ΔCp_est, Telu_max
end=#

function estimate_start_parameter_single_ramp(tRs::DataFrame, TPs::DataFrame, PPs::DataFrame, L, d, df, gas; pout="vacuum", time_unit="min", control="Pressure")
    if time_unit == "min"
		a = 60.0
	else
		a = 1.0
	end
    tR_meas = Array(tRs[:,2:end]).*a
    if pout == "atmospheric"
		pout = PPs.pamb
	else
		pout = "vacuum"
	end
    nt, ns = size(tR_meas)
    tMref = Array{Float64}(undef, nt)
    RT = Array{Float64}(undef, nt) 
    Telu_meas = Array{Float64}(undef, nt, ns)
    for i=1:nt
        prog = Program(collect(skipmissing(TPs[i, 2:end])), collect(skipmissing(PPs[i, 2:end])), L; pout=pout, time_unit=time_unit)
        tMref[i] = reference_holdup_time(prog, L, d, gas; control=control)/a
        RT[i] = TPs[i, 4] # single-ramp temperature programs are assumed
        Telu_meas[i,:] = elution_temperature(tR_meas[i,:], prog)
    end 
    rT = RT.*tMref./θref
    Telu_max = Array{Float64}(undef, ns)
    Tchar_est = Array{Float64}(undef, ns)
    θchar_est = Array{Float64}(undef, ns)
    ΔCp_est = Array{Float64}(undef, ns)
    for i=1:ns
        interp = interpolate((rT, ), Telu_meas[:,i], Gridded(Linear()))
        Telu_max[i] = maximum(Telu_meas[:,i])
        Tchar_est[i] = interp(rT_nom)
        θchar_est[i] = 22.0*(Tchar_est[i]/Tst)^0.7*(1000*df/d)^0.09
        ΔCp_est[i] = 100.0
    end
    return Tchar_est, θchar_est, ΔCp_est, Telu_max
end

function estimate_start_parameter_mean_elu_temp(tRs::DataFrame, TPs::DataFrame, PPs::DataFrame, L, d, df; pout="vacuum", time_unit="min")
	if time_unit == "min"
		a = 60.0
	else
		a = 1.0
	end
    tR_meas = Array(tRs[:,2:end]).*a
    if pout == "atmospheric"
		pout = PPs.pamb
	else
		pout = "vacuum"
	end
    nt, ns = size(tR_meas)
	Telu_max = Array{Float64}(undef, ns)
	Tchar_elu = Array{Float64}(undef, ns)
	θchar_elu = Array{Float64}(undef, ns)
	for j=1:ns
		Telu = Array{Float64}(undef, nt)
		for i=1:nt
            prog = RetentionParameterEstimator.Program(collect(skipmissing(TPs[i, 2:end])), collect(skipmissing(PPs[i, 2:end])), L; pout=pout, time_unit=time_unit)
			Telu[i] = RetentionParameterEstimator.elution_temperature(tR_meas[i,j], prog)[1]
		end
		Telu_max[j] = maximum(Telu)
		Tchar_elu[j] = mean(Telu)
		θchar_elu[j] = 22.0*(Tchar_elu[j]/273.15)^0.7*(1000*df/d)^0.09
	end
	ΔCp_elu = 100.0.*ones(ns)
	return Tchar_elu, θchar_elu, ΔCp_elu, Telu_max
end

function estimate_start_parameter(tR_meas::DataFrame, TPs::DataFrame, PPs::DataFrame, L, d, df, gas; pout="vacuum", time_unit="min", control="Pressure")
    if size(TPs)[2] == 6
        Tchar_est, θchar_est, ΔCp_est, Telu_max = estimate_start_parameter_single_ramp(tR_meas, TPs::DataFrame, PPs::DataFrame, L, d, df, gas; pout=pout, time_unit=time_unit, control=control)
    else 
        Tchar_est, θchar_est, ΔCp_est, Telu_max = estimate_start_parameter_mean_elu_temp(tR_meas, TPs::DataFrame, PPs::DataFrame, L, d, df; pout=pout, time_unit=time_unit)
    end
    return Tchar_est, θchar_est, ΔCp_est, Telu_max
end

"""
    estimate_start_parameter(Telu_meas, tMref_meas, TPs)

Estimation of start values for the K-centric parameters `Tchar`, `θchar` and `ΔCp` from measured elution temperatures `Telu_meas` and measured reference 
    holdup-times `tMref_meas` at reference temperature of 150°C for several GC programs with temperature programs `TPs`. The programs `TPs`
are defined as conventional programs (Vector of "start value", "hold time of start value", "ramp", "end value", "hold time of end value"). **The temperature program
should have only one ramp. The time units of `tMref_meas` and `TPs` have to be the same.**
"""   
function estimate_start_parameter(Telu_meas, tMref_meas, TPs)
    RT = Array{Float64}(undef, length(TPs))
    for i=1:length(TPs)
        RT[i] = TPs[i][3] # single-ramp temperature programs are assumed
    end
    rT = RT.*tMref_meas./θref
    Tchar_est = Array{Float64}(undef, size(Telu_meas)[2])
    θchar_est = Array{Float64}(undef, size(Telu_meas)[2])
    ΔCp_est = Array{Float64}(undef, size(Telu_meas)[2])
    for i=1:size(Telu_meas)[2]
        interp = interpolate((rT, ), Telu_meas[:,i], Gridded(Linear()))
        #spl = Spline1D(rT, Telu_meas[:,i])
        #Tchar_est[i] = spl(rT_nom)
        Tchar_est[i] = interp(rT_nom)
        θchar_est[i] = 22.0*(Tchar_est[i]/Tst)^0.7 # factor of φ?
        ΔCp_est[i] = 100.0
    end
    return Tchar_est, θchar_est, ΔCp_est
end