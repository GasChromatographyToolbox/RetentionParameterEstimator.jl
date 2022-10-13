# functions used for the optimization of the loss-function 

# this function is a modified version from GasChromatographySimulator
"""
    Program(TP, FpinP, L; pout="vacuum", time_unit="min")

Construct the structure `Program` with conventional formulation (see [`conventional_program`](@ref)) of programs for the case
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
    time_steps = GasChromatographySimulator.common_time_steps(ts1, ts2)
    temp_steps = GasChromatographySimulator.new_value_steps(Ts, ts1, time_steps)
    Fpin_steps = GasChromatographySimulator.new_value_steps(Fps, ts2, time_steps)
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

"""
    opt_Kcentric(x_Kcentric, p)

Function used for optimization of the loss-function in regards to the three K-centric parameters.

# Arguments
* `x_Kentric` ... 3n-vector of the three K-centric parameters of n solutes. Elements 1:n are Tchar, n+1:2n are θchar and 2n+1:3n are ΔCp values.
* `p` ... vector containing the fixed parameters:
    * `tR = p[1]` ... mxn-array of the measured retention times in seconds.
    * `L = p[2]` ... number of the length of the column in m.
    * `d = p[3]` ... number of the diameters of the column in m.
    * `prog = p[4]` ... m-array of structure GasChromatographySimulator.Programs containing the definition of the GC-programs.
    * `opt = p[5]` ... struture GasChromatographySimulator.Options containing the settings for options for the simulation.
    * `gas = p[6]` ... string of name of the mobile phase gas. 

# Output
* `sum((tR.-tRcalc).^2)` ... sum of the squared residuals over m GC-programs and n solutes.
"""
function opt_Kcentric(x_Kcentric, p)
	tR = p[1]
	L = p[2]
	d = p[3]
	prog = p[4]
	opt = p[5]
    gas = p[6]
    if length(size(tR)) == 1
		ns = 1
	else
		ns = size(tR)[2]
	end
	Tchar = x_Kcentric[1:ns] # Array length = number solutes
	θchar = x_Kcentric[ns+1:2*ns] # Array length = number solutes
	ΔCp = x_Kcentric[2*ns+1:3*ns] # Array length = number solutes
    return loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas)[1]
end

function opt_dL(x_dL, p)
	tR = p[1]
	prog = p[2]
	opt = p[3]
    gas = p[4]
	Tchar = p[5] # Array length = number solutes
	θchar = p[6] # Array length = number solutes
	ΔCp = p[7] # Array length = number solutes
    #if length(size(tR)) == 1
	#	ns = 1
	#else
	#	ns = size(tR)[2]
	#end
	L = x_dL[1]
	d = x_dL[2]
    return loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas)[1]
end

function opt_dLKcentric(x_Kcentric, p)
	tR = p[1]
	prog = p[2]
	opt = p[3]
    gas = p[4]
    if length(size(tR)) == 1
		ns = 1
	else
		ns = size(tR)[2]
	end
	L = x_dLKcentric[1]
	d = x_dLKcentric[2]
	Tchar = x_dLKcentric[3:ns+2] # Array length = number solutes
	θchar = x_dLKcentric[ns+2+1:2*(ns+2)] # Array length = number solutes
	ΔCp = x_dLKcentric[2*(ns+2)+1:3*(ns+2)] # Array length = number solutes
    return loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas)[1]
end

# optimize every solute separatly, tR is a 2D-array with RT of different programs in the first dimension and different solutes in the second dimension  
function optimize_Kcentric_single(tR, L, d, gas, prog, opt, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=10000)
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
                BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
	
    if typeof(size(tR)) == Tuple{Int64, Int64}
        n2 = size(tR)[2]
    else
        n2 = 1
    end 
    opt_sol = Array{Any}(undef, n2)
	for i=1:n2
		p = [tR[:,i], L, d, prog, opt, gas]
		x0 = [Tchar_e[i], θchar_e[i], ΔCp_e[i]]
		lb = [lb_Tchar[i], lb_θchar[i], lb_ΔCp[i]]
		ub = [ub_Tchar[i], ub_θchar[i], ub_ΔCp[i]]
		if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton())
			prob = OptimizationProblem(optf, x0, p)
		else
			prob = OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
		end
        if method in optimisers
            opt_sol[i] = solve(prob, method, maxiters=maxiters)
        #elseif method in bbos
        #    opt_sol[i] = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
        else
		    opt_sol[i] = solve(prob, method) #-> :u (Array of the optimized parameters), :minimum (minima of the optimization function) , :retcode (Boolean, successful?)
        end
    end
	return opt_sol
end

# optimization of only one solute
function optimize_Kcentric_single(tR, L, d, gas, prog, opt, Tchar_e::Float64, θchar_e::Float64, ΔCp_e::Float64, lb_Tchar::Float64, lb_θchar::Float64, lb_ΔCp::Float64, ub_Tchar::Float64, ub_θchar::Float64, ub_ΔCp::Float64, method; maxiters=10000)
	opt_sol[1] = optimize_Kcentric_single(tR, L, d, gas, prog, opt, [Tchar_e], [θchar_e], [ΔCp_e], [lb_Tchar], [lb_θchar], [lb_ΔCp], [ub_Tchar], [ub_θchar], [ub_ΔCp], method; maxiters=maxiters)
end

# optimize all solutes together
function optimize_Kcentric_all(tR, L, d, gas, prog, opt, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=10000)
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
                BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    p = [tR, L, d, prog, opt, gas]
	x0 = [Tchar_e; θchar_e; ΔCp_e]
	lb = [lb_Tchar; lb_θchar; lb_ΔCp]
	ub = [ub_Tchar; ub_θchar; ub_ΔCp]
	optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
	if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton())
		prob = OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
	else
		prob = OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
	end
	if method in optimisers
        opt_sol = solve(prob, method, maxiters=maxiters)
    elseif method in bbos
        opt_sol = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
    else
        opt_sol = solve(prob, method) #-> :u (Array of the optimized parameters), :minimum (minima of the optimization function) , :retcode (Boolean, successful?)
    end
	return opt_sol
end



# optimization for retention parameters, every solute separatly
function optimize(tR_meas, solute_names, column, options, TPs, PPs, method, method_short_name; maxiters=10000, relbound=0.5, mode="single")
	L = column.L[1]
    d = column.d[1]
    gas = string(column.gas[1])
    pout = string(column.pout[1])
    time_unit = string(column.time_unit[1])

    if length(size(tR_meas)) == 1
        ns = 1
    else
        ns = size(tR_meas)[2]
    end

	prog = Array{GasChromatographySimulator.Program}(undef, length(TPs.measurement))
    for i=1:length(TPs.measurement)
        prog[i] = GasChromatographySimulator.Program(collect(skipmissing(TPs[i, 2:end])), collect(skipmissing(PPs[i, 2:end])), column.L[1]; pout=string(column.pout[1]), time_unit=string(column.time_unit[1]))
    end
    # estimation of start parameter only works with single ramp programs
	Tchar_est, θchar_est, ΔCp_est = estimate_start_parameter(tR_meas, TPs, PPs, L, d, gas; pout=pout, time_unit=time_unit, control=options.control)

    method_short_names = Array{String}(undef, ns)
    Tchar = Array{Float64}(undef, ns)
	θchar = Array{Float64}(undef, ns)
	ΔCp = Array{Float64}(undef, ns)
    min = Array{Float64}(undef, ns)
    retcode = Array{Any}(undef, ns)
    if mode == "single"
	    sol = optimize_Kcentric_single(tR_meas, L, d, gas, prog, options, Tchar_est, θchar_est, ΔCp_est, relbound.*Tchar_est, relbound.*θchar_est, relbound.*ΔCp_est, (1.0+relbound).*Tchar_est, (1.0+relbound).*θchar_est, (1.0+relbound).*ΔCp_est, method; maxiters=maxiters)
        for j=1:ns
            method_short_names[j] = method_short_name
            Tchar[j] = sol[j][1]
            θchar[j] = sol[j][2]
            ΔCp[j] = sol[j][3]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
    else
        sol = optimize_Kcentric_all(tR_meas, L, d, gas, prog, options, Tchar_est, θchar_est, ΔCp_est, relbound.*Tchar_est, relbound.*θchar_est, relbound.*ΔCp_est, (1.0+relbound).*Tchar_est, (1.0+relbound).*θchar_est, (1.0+relbound).*ΔCp_est, method; maxiters=maxiters)
        Tchar = sol[1:ns] # Array length = number solutes
        θchar = sol[ns+1:2*ns] # Array length = number solutes
        ΔCp = sol[2*ns+1:3*ns] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
            method_short_names[j] = method_short_name
        end
    end
	df = DataFrame(Name=solute_names, Methods=method_short_names, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
    #CSV.write(file, df)
	return df, sol
end

# optimization for retention parameters, every solute separatly and write results in file
function optimize(file, tR_meas, solute_names, column, options, TPs, PPs, method, method_short_name; maxiters=10000, relbound=0.5, mode="single")
	df, sol = optimize(tR_meas, solute_names, column, options, TPs, PPs, method, method_short_name; maxiters=maxiters, relbound=relbound, mode=mode)
    CSV.write(file, df)
	return df, sol
end

# fixed initial value
function optimize(tR_meas, solute_names, column, options, TPs, PPs, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method, method_short_name; maxiters=10000, mode="single")
	L = column.L[1]
    d = column.d[1]
    gas = string(column.gas[1])
    pout = string(column.pout[1])
    time_unit = string(column.time_unit[1])

    if length(size(tR_meas)) == 1
        ns = 1
    else
        ns = size(tR_meas)[2]
    end

	prog = Array{GasChromatographySimulator.Program}(undef, length(TPs.measurement))
    for i=1:length(TPs.measurement)
        prog[i] = GasChromatographySimulator.Program(collect(skipmissing(TPs[i, 2:end])), collect(skipmissing(PPs[i, 2:end])), column.L[1]; pout=string(column.pout[1]), time_unit=string(column.time_unit[1]))
    end
    method_short_names = Array{String}(undef, ns)
    Tchar = Array{Float64}(undef, ns)
	θchar = Array{Float64}(undef, ns)
	ΔCp = Array{Float64}(undef, ns)
    min = Array{Float64}(undef, ns)
    retcode = Array{Any}(undef, ns)
    if mode == "single"
	    sol = optimize_Kcentric_single(tR_meas, L, d, gas, prog, options, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        for j=1:ns
            method_short_names[j] = method_short_name
            Tchar[j] = sol[j][1]
            θchar[j] = sol[j][2]
            ΔCp[j] = sol[j][3]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
    else
        sol = optimize_Kcentric_all(tR_meas, L, d, gas, prog, options, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        Tchar = sol[1:ns] # Array length = number solutes
        θchar = sol[ns+1:2*ns] # Array length = number solutes
        ΔCp = sol[2*ns+1:3*ns] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
            method_short_names[j] = method_short_name
        end
    end
	df = DataFrame(Name=solute_names, Methods=method_short_names, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
	return df, sol
end

function optimize(file, tR_meas, solute_names, column, options, TPs, PPs, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method, method_short_name; maxiters=10000, mode="single")
	df, sol = optimize(tR_meas, solute_names, column, options, TPs, PPs, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method, method_short_name; maxiters=maxiters, mode=mode)
    CSV.write(file, df)
	return df, sol
end

# benchmark for one solute, using optimization with estimated initial values
# benchmark time measurement -> function, extract :times, :memory, allocs / input: parameters for BenchmarkTools, meassured RT, system parameters as dictionary
#function benchmark_optimize(file, tR_meas, solute_names, column, options, TPs, PPs, method, method_short_name)
#	#times = Array{Float64}(undef, length(solute_names))
#	#memory = Array{Float64}(undef, length(solute_names))
#	#allocs = Array{Float64}(undef, length(solute_names))
#	b = @benchmark optimize($tR_meas, $solute_names, $column, $options, $TPs, $PPs, $method, $method_short_name) seconds = 600
#	times = sum(b.times)/length(b.times)
#	memory = b.memory
#	allocs = b.allocs
#	df = DataFrame(Methods=method_short_name, Times=times, Memory=memory, Allocs=allocs)
#   CSV.write(file, df; append=true)
#	return df
#end

## functions to estimate L and d while the retention parameters are known
function optimize_dL(tR, gas, prog, opt, Tchar::Array{Number}, θchar::Array{Number}, ΔCp::Array{Number}, L_e, d_e, lb_L, lb_d, ub_L, ub_d, method; maxiters=10000)
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
                BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    p = [tR, prog, opt, gas, Tchar, θchar, ΔCp]
	x0 = [L_e, d_e]
	lb = [lb_L, lb_d]
	ub = [ub_L, ub_d]
	optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
	if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton())
		prob = OptimizationProblem(optf, x0, p)
	else
		prob = OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
	end
	if method in optimisers
        opt_sol = solve(prob, method, maxiters=maxiters)
    elseif method in bbos
        opt_sol = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
    else
        opt_sol = solve(prob, method) #-> :u (Array of the optimized parameters), :minimum (minima of the optimization function) , :retcode (Boolean, successful?)
    end
	return opt_sol
end