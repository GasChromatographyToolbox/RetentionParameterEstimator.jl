# functions used to estimate start values of the parameters

"""
    reference_holdup_time(prog, L, d, gas; control="Pressure")

Calculate the reference holdup time for the GC program `prog` for a column with length `L` and diameter `d` and `gas` as mobile phase. The reference holdup time is the holdup time at the reference temperature 150°C.
"""
function reference_holdup_time(col, prog; control="Pressure")
    Tref = 150.0
    # estimate the time of the temperature program for T=Tref
    t_ = prog.time_steps
    T_ = prog.temp_steps
    T_diff = T_ .- Tref
    unique_indices = GasChromatographySimulator.deduplicate_knots!(T_diff; move_knots=true)
    T_diff_unique = T_diff[unique_indices]
    t_cumsum_unique = cumsum(t_)[unique_indices]
    # Use GasChromatographySimulator's linear_interpolation
    interp = GasChromatographySimulator.linear_interpolation((T_diff_unique,), t_cumsum_unique)
    tref = interp(0.0)
    # inlet and outlet pressure at time tref
    Fpin_ref = prog.Fpin_itp(tref)
    pout_ref = prog.pout_itp(tref)
    # hold-up time calculated for the time of the program, when T=Tref
    tMref = GasChromatographySimulator.holdup_time(Tref+273.15, Fpin_ref, pout_ref, col.L, col.d, col.gas, control=control)
    return tMref
end

"""
    elution_temperature(tRs, prog)

Calculate the elution temperatures from retention times `tRs` and GC programs `prog`.
"""
function elution_temperature(tRs, prog)
    Telus = Array{Float64}(undef, length(tRs))
    for i=1:length(tRs)
        if ismissing(tRs[i])
            Telus[i] = NaN
        else
            Telus[i] = prog.T_itp(0.0, tRs[i])
        end
    end
    return Telus
end

"""
    estimate_start_parameter_single_ramp(tRs::DataFrame, col, prog; time_unit="min", control="Pressure")

Estimation of initial parameters for `Tchar`, `θchar` and `ΔCp` based on the elution temperatures calculated from the retention times `tR` and GC programs `prog` for column `col`.
For this function it is assumed, that single ramp heating programs are used. The elution temperatures of all measurements are calculated and than interpolated over the heating rates. For a dimensionless heating rate of 0.69 the elution temperature and the characteristic temperature of a substance are nearly equal.
Based on this estimated `Tchar` estimates for the initial values of `θchar` and `ΔCp` are calculated as
    ``
    \\theta_{char,init} = 22 \\left(\\frac{T_{char,init}}{T_{st}}\\right)^{0.7} \\left(1000\\frac{d_f}{d}\\right)^{0.09} °C
    ``
and
    ``
    \\Delta C_p = (-52 + 0.34 T_{char,init}) \\mathrm{J mol^{-1} K^{-1}}
    ``

# Output
* `Tchar_est` ... estimate for initial guess of the characteristic temperature
* `θchar_est` ... estimate for initial guess of θchar
* `ΔCp_est` ... estimate for initial guess of ΔCp
* `Telu_max` ... the maximum of the calculated elution temperatures of the solutes
"""    
function estimate_start_parameter_single_ramp(tRs::DataFrame, col, prog; time_unit="min", control="Pressure")
    a = time_unit_conversion_factor(time_unit)
    tR_meas = Array(tRs[:,2:end]).*a
    nt, ns = size(tR_meas)
    tMref = Array{Float64}(undef, nt)
    RT = Array{Float64}(undef, nt) 
    Telu_meas = Array{Float64}(undef, nt, ns)
    for i=1:nt
        tMref[i] = reference_holdup_time(col, prog[i]; control=control)/a
        # single-ramp temperature programs with ramp between time_steps 2 and 3 are assumed
        RT[i] = (prog[i].temp_steps[3] - prog[i].temp_steps[2])/prog[i].time_steps[3]*a 
        Telu_meas[i,:] = elution_temperature(tR_meas[i,:], prog[i])
    end 
    rT = RT.*tMref./θref
    Telu_max = Array{Float64}(undef, ns)
    Tchar_est = Array{Float64}(undef, ns)
    θchar_est = Array{Float64}(undef, ns)
    ΔCp_est = Array{Float64}(undef, ns)
    for i=1:ns
		indexExistingRT = Not(findall(ismissing.(tR_meas[:,i])))
        rT_values = rT[indexExistingRT] .- rT_nom
        Telu_values = Telu_meas[indexExistingRT,i]
        # Sort by rT_values for interpolation (GasChromatographySimulator requires sorted x-values)
        sort_indices = sortperm(rT_values)
        rT_sorted = rT_values[sort_indices]
        Telu_sorted = Telu_values[sort_indices]
        # Deduplicate y-values (Telu) to match original Interpolations.jl behavior
        # This modifies Telu_sorted in place, similar to original: knots = Interpolations.deduplicate_knots!(Telu_meas[...])
        unique_indices = GasChromatographySimulator.deduplicate_knots!(Telu_sorted; move_knots=true)
        Telu_max[i] = maximum(Telu_meas[indexExistingRT,i])
        # Use GasChromatographySimulator's linear_interpolation
        interp = GasChromatographySimulator.linear_interpolation((rT_sorted,), Telu_sorted)
        Tchar_est[i] = interp(0.0)
        θchar_est[i] = 22.0*(Tchar_est[i]/Tst)^0.7*(1000*col.df/col.d)^0.09
        ΔCp_est[i] = -52.0 + 0.34*Tchar_est[i]
    end
    return Tchar_est, θchar_est, ΔCp_est, Telu_max
end

"""
    estimate_start_parameter_single_ramp_weighted(tRs::DataFrame, col, prog; time_unit="min", control="Pressure", α=2.0, weighted_fraction=0.7)

Estimation of initial parameters for `Tchar`, `θchar` and `ΔCp` based on the elution temperatures calculated from the retention times `tR` and GC programs `prog` for column `col`.
This function is similar to `estimate_start_parameter_single_ramp`, but uses a weighted approach that gives more weight to measurements with higher heating rates, as these typically provide more accurate Tchar estimates.

For this function it is assumed, that single ramp heating programs are used. The elution temperatures of all measurements are calculated and then interpolated over the heating rates. 

The estimation uses a weighted approach:
1. A weighted average of elution temperatures (default 70% weight), where higher heating rates receive exponentially more weight
2. An interpolation to the nominal heating rate rT_nom = 0.69 (default 30% weight), which provides theoretical correction

The weighting uses exponential scaling: `weight = exp(α * (rT - rT_min))`, where `α` controls the strength of the weighting (default 2.0).

# Arguments
* `tRs` ... DataFrame with retention times
* `col` ... Column object
* `prog` ... Array of GC programs
* `time_unit="min"` ... Time unit of the retention times
* `control="Pressure"` ... Control mode for holdup time calculation
* `α=2.0` ... Weighting factor for exponential weighting (higher values give more weight to high heating rates)
* `weighted_fraction=0.7` ... Fraction of weighted average in the final estimate (remaining fraction uses interpolation)

# Output
* `Tchar_est` ... estimate for initial guess of the characteristic temperature
* `θchar_est` ... estimate for initial guess of θchar
* `ΔCp_est` ... estimate for initial guess of ΔCp
* `Telu_max` ... the maximum of the calculated elution temperatures of the solutes

Based on this estimated `Tchar` estimates for the initial values of `θchar` and `ΔCp` are calculated as
    ``
    \\theta_{char,init} = 22 \\left(\\frac{T_{char,init}}{T_{st}}\\right)^{0.7} \\left(1000\\frac{d_f}{d}\\right)^{0.09} °C
    ``
and
    ``
    \\Delta C_p = (-52 + 0.34 T_{char,init}) \\mathrm{J mol^{-1} K^{-1}}
    ``
"""    
function estimate_start_parameter_single_ramp_weighted(tRs::DataFrame, col, prog; time_unit="min", control="Pressure", α=2.0, weighted_fraction=0.7)
    a = time_unit_conversion_factor(time_unit)
    tR_meas = Array(tRs[:,2:end]).*a
    nt, ns = size(tR_meas)
    tMref = Array{Float64}(undef, nt)
    RT = Array{Float64}(undef, nt) 
    Telu_meas = Array{Float64}(undef, nt, ns)
    for i=1:nt
        tMref[i] = reference_holdup_time(col, prog[i]; control=control)/a
        # single-ramp temperature programs with ramp between time_steps 2 and 3 are assumed
        RT[i] = (prog[i].temp_steps[3] - prog[i].temp_steps[2])/prog[i].time_steps[3]*a 
        Telu_meas[i,:] = elution_temperature(tR_meas[i,:], prog[i])
    end 
    rT = RT.*tMref./θref
    Telu_max = Array{Float64}(undef, ns)
    Tchar_est = Array{Float64}(undef, ns)
    θchar_est = Array{Float64}(undef, ns)
    ΔCp_est = Array{Float64}(undef, ns)
    for i=1:ns
		indexExistingRT = Not(findall(ismissing.(tR_meas[:,i])))
        rT_values = rT[indexExistingRT] .- rT_nom
        Telu_values = Telu_meas[indexExistingRT,i]
        rT_original = rT[indexExistingRT]  # Original rT values (not shifted)
        
        # Sort by rT_values for interpolation (GasChromatographySimulator requires sorted x-values)
        sort_indices = sortperm(rT_values)
        rT_sorted = rT_values[sort_indices]
        Telu_sorted = Telu_values[sort_indices]
        rT_original_sorted = rT_original[sort_indices]
        
        # Deduplicate y-values (Telu) to match original Interpolations.jl behavior
        # This modifies Telu_sorted in place, similar to original: knots = Interpolations.deduplicate_knots!(Telu_meas[...])
        unique_indices = GasChromatographySimulator.deduplicate_knots!(Telu_sorted; move_knots=true)
        # Apply same deduplication to rT arrays
        rT_sorted_unique = rT_sorted[unique_indices]
        Telu_sorted_unique = Telu_sorted[unique_indices]
        rT_original_unique = rT_original_sorted[unique_indices]
        
        Telu_max[i] = maximum(Telu_meas[indexExistingRT,i])
        
        # Weight measurements by heating rate: higher heating rates get more weight
        # This accounts for the observation that higher heating rates give more accurate Tchar estimates
        # Use exponential weighting: weight = exp(α * (rT - rT_min)) where α controls the strength
        # Normalize so weights sum to 1
        if length(rT_original_unique) >= 2
            rT_min = minimum(rT_original_unique)
            # Use exponential weighting: weight increases with rT, with stronger emphasis on higher rates
            # α controls the strength: higher α means more emphasis on high heating rates
            weights = exp.(α .* (rT_original_unique .- rT_min))
            weights = weights ./ sum(weights)
            
            # Use GasChromatographySimulator's linear_interpolation for the base interpolation
            interp = GasChromatographySimulator.linear_interpolation((rT_sorted_unique,), Telu_sorted_unique)
            
            # Calculate weighted Tchar estimate: use interpolation at rT_nom as baseline,
            # but bias towards higher heating rates using weighted average
            Tchar_interp = interp(0.0)  # Interpolation at rT_nom
            
            # Weighted average of Telu values (higher heating rates weighted more)
            Tchar_weighted = sum(weights .* Telu_sorted_unique)
            
            # Combine: use weighted average but bias towards interpolation result
            # The interpolation provides theoretical correction, weighted average provides empirical trend
            Tchar_est[i] = weighted_fraction * Tchar_weighted + (1.0 - weighted_fraction) * Tchar_interp
        else
            # Fallback: if only one measurement, use it directly
            Tchar_est[i] = Telu_sorted_unique[1]
        end
        
        θchar_est[i] = 22.0*(Tchar_est[i]/Tst)^0.7*(1000*col.df/col.d)^0.09
        ΔCp_est[i] = -52.0 + 0.34*Tchar_est[i]
    end
    return Tchar_est, θchar_est, ΔCp_est, Telu_max
end


"""
    estimate_start_parameter_mean_elu_temp(tRs::DataFrame, col, prog; time_unit="min", control="Pressure")

Estimation of initial parameters for `Tchar`, `θchar` and `ΔCp` based on the elution temperatures calculated from the retention times `tR` and GC programs `prog` for column `col`.
This function is used, if the temperature program is not a single ramp heating program. The elution temperatures of all measurements are calculated and the mean value of the elution temperatures is used as the initial characteristic temperature of a substance.
Based on this estimated `Tchar` estimates for the initial values of `θchar` and `ΔCp` are calculated as
    ``
    \\theta_{char,init} = 22 \\left(\\frac{T_{char,init}}{T_{st}}\\right)^{0.7} \\left(1000\\frac{d_f}{d}\\right)^{0.09} °C
    ``
and
    ``
    \\Delta C_p = (-52 + 0.34 T_{char,init}) \\mathrm{J mol^{-1} K^{-1}}
    ``

# Output
* `Tchar_est` ... estimate for initial guess of the characteristic temperature
* `θchar_est` ... estimate for initial guess of θchar
* `ΔCp_est` ... estimate for initial guess of ΔCp
* `Telu_max` ... the maximum of the calculated elution temperatures of the solutes
""" 
function estimate_start_parameter_mean_elu_temp(tRs::DataFrame, col, prog; time_unit="min")
	a = time_unit_conversion_factor(time_unit)
    tR_meas = Array(tRs[:,2:end]).*a
    
    nt, ns = size(tR_meas)
	Telu_meas = Array{Float64}(undef, nt, ns)
    for i=1:nt
        Telu_meas[i,:] = elution_temperature(tR_meas[i,:], prog[i])
    end
	Telu_max = Array{Float64}(undef, ns)
	Tchar_est = Array{Float64}(undef, ns)
	θchar_est = Array{Float64}(undef, ns)
    ΔCp_est = Array{Float64}(undef, ns)
	for i=1:ns
		indexExistingRT = Not(findall(ismissing.(tR_meas[:,i])))
		Telu_max[i] = maximum(Telu_meas[indexExistingRT,i])
		Tchar_est[i] = mean(Telu_meas[indexExistingRT,i])
		θchar_est[i] = 22.0*(Tchar_est[i]/273.15)^0.7*(1000*col.df/col.d)^0.09
        ΔCp_est[i] = -52.0 + 0.34*Tchar_est[i]
	end
	return Tchar_est, θchar_est, ΔCp_est, Telu_max
end

"""
    estimate_start_parameter(tRs::DataFrame, col, prog; time_unit="min", control="Pressure")

Estimation of initial parameters for `Tchar`, `θchar` and `ΔCp` based on the elution temperatures calculated from the retention times `tR` and GC programs `prog` for column `col`.
The initial value of `Tchar` is estimated from the elution temperatures of the measurements.
Based on this estimated `Tchar` estimates for the initial values of `θchar` and `ΔCp` are calculated as
    ``
    \\theta_{char,init} = 22 \\left(\\frac{T_{char,init}}{T_{st}}\\right)^{0.7} \\left(1000\\frac{d_f}{d}\\right)^{0.09} °C
    ``
and
    ``
    \\Delta C_p = (-52 + 0.34 T_{char,init}) \\mathrm{J mol^{-1} K^{-1}}
    ``

# Output
* `Tchar_est` ... estimate for initial guess of the characteristic temperature
* `θchar_est` ... estimate for initial guess of θchar
* `ΔCp_est` ... estimate for initial guess of ΔCp
* `Telu_max` ... the maximum of the calculated elution temperatures of the solutes
""" 
function estimate_start_parameter(tR_meas::DataFrame, col, prog; time_unit="min", control="Pressure")
    Tchar_est, θchar_est, ΔCp_est, Telu_max = try
        estimate_start_parameter_single_ramp(tR_meas, col, prog; time_unit=time_unit, control=control)
    catch 
        estimate_start_parameter_mean_elu_temp(tR_meas, col, prog; time_unit=time_unit)
    end
    return Tchar_est, θchar_est, ΔCp_est, Telu_max
end