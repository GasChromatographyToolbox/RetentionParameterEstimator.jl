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
    spl = Spline1D(cumsum(t_), T_ .- Tref)
    tref = roots(spl)[1]
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

"""
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
        prog[i] = GasChromatographySimulator.Program(TPs[i], PPs[i], L; pout=pout, time_unit=time_unit)
        tMref[i] = reference_holdup_time(prog[i], L, d, gas; control=control)/t_conv
        RT[i] = TPs[i][3] # single-ramp temperature programs are assumed
        Telu_meas[i,:] = elution_temperature(tR_meas[i,:], prog[i])
    end 
    rT = RT.*tMref./θref
    Tchar_est = Array{Float64}(undef, n2)
    θchar_est = Array{Float64}(undef, n2)
    ΔCp_est = Array{Float64}(undef, n2)
    for i=1:n2
        spl = Spline1D(rT, Telu_meas[:,i])
        Tchar_est[i] = spl(rT_nom)
        θchar_est[i] = 22.0*(Tchar_est[i]/Tst)^0.7 # factor of φ?
        ΔCp_est[i] = 100.0
    end
    return Tchar_est, θchar_est, ΔCp_est
end

function estimate_start_parameter(tR_meas, TPs::DataFrame, PPs::DataFrame, L, d, gas; pout="vacuum", time_unit="min", control="Pressure")
    if time_unit == "min"
        t_conv = 60.0
    else
        t_conv = 1.0
    end
    prog = Array{GasChromatographySimulator.Program}(undef, length(TPs.measurement))
    tMref = Array{Float64}(undef, length(TPs.measurement))
    RT = Array{Float64}(undef, length(TPs.measurement)) 
    if typeof(size(tR_meas)) == Tuple{Int64, Int64}
        n2 = size(tR_meas)[2]
    else
        n2 = 1
    end 
    Telu_meas = Array{Float64}(undef, size(tR_meas)[1], n2)
    for i=1:length(TPs.measurement)
        prog[i] = GasChromatographySimulator.Program(collect(skipmissing(TPs[i, 2:end])), collect(skipmissing(PPs[i, 2:end])), L; pout=pout, time_unit=time_unit)
        tMref[i] = reference_holdup_time(prog[i], L, d, gas; control=control)/t_conv
        RT[i] = TPs[i, 4] # single-ramp temperature programs are assumed
        Telu_meas[i,:] = elution_temperature(tR_meas[i,:], prog[i])
    end 
    rT = RT.*tMref./θref
    Tchar_est = Array{Float64}(undef, n2)
    θchar_est = Array{Float64}(undef, n2)
    ΔCp_est = Array{Float64}(undef, n2)
    for i=1:n2
        spl = Spline1D(rT, Telu_meas[:,i])
        Tchar_est[i] = spl(rT_nom)
        θchar_est[i] = 22.0*(Tchar_est[i]/Tst)^0.7 # factor of φ?
        ΔCp_est[i] = 100.0
    end
    return Tchar_est, θchar_est, ΔCp_est
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
        spl = Spline1D(rT, Telu_meas[:,i])
        Tchar_est[i] = spl(rT_nom)
        θchar_est[i] = 22.0*(Tchar_est[i]/Tst)^0.7 # factor of φ?
        ΔCp_est[i] = 100.0
    end
    return Tchar_est, θchar_est, ΔCp_est
end