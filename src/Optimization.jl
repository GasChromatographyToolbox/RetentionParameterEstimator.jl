# functions used for the optimization of the loss-function 

# this function is a modified version from GasChromatographySimulator
# if it works, put this function into GasChromatographySimulator
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

#------------------
# Optimization only for the K-centric retention parameters
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

function opt_Kcentric_(x_Kcentric, p)
	tR = p[1]
	L = p[2]
	d = p[3]
    df = p[4]
	prog = p[5]
	opt = p[6]
    gas = p[7]
    if length(size(tR)) == 1
		ns = 1
	else
		ns = size(tR)[2]
	end
	Tchar = x_Kcentric[1:ns] # Array length = number solutes
	θchar = x_Kcentric[ns+1:2*ns] # Array length = number solutes
	ΔCp = x_Kcentric[2*ns+1:3*ns] # Array length = number solutes
    return loss_(tR, Tchar, θchar, ΔCp, L, d, df, prog, opt, gas)[1]
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
		if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
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

function optimize_Kcentric_single_(tR, L, d, df, gas, prog, opt, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=10000)
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
                BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    optf = OptimizationFunction(opt_Kcentric_, Optimization.AutoForwardDiff())
	
    if typeof(size(tR)) == Tuple{Int64, Int64}
        n2 = size(tR)[2]
    else
        n2 = 1
    end 
    opt_sol = Array{Any}(undef, n2)
	for i=1:n2
		p = [tR[:,i], L, d, df, prog, opt, gas]
		x0 = [Tchar_e[i], θchar_e[i], ΔCp_e[i]]
		lb = [lb_Tchar[i], lb_θchar[i], lb_ΔCp[i]]
		ub = [ub_Tchar[i], ub_θchar[i], ub_ΔCp[i]]
		if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
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
	if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
		prob = OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
	else
		prob = OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
	end
	if method in optimisers
        opt_sol = solve(prob, method, maxiters=maxiters)
    #elseif method in bbos
    #    opt_sol = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
    else
        opt_sol = solve(prob, method)
    end
	return opt_sol
end

# with given initial value
# DELETE
function optimize(tR_meas, solute_names, column, options, TPs, PPs, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=10000, mode="single")
    if column[:time_unit] == "min"
        a = 60.0
    else
        a = 1.0
    end

    if length(size(tR_meas)) == 1
        ns = 1
    else
        ns = size(tR_meas)[2]
    end

	prog = Array{GasChromatographySimulator.Program}(undef, length(TPs.measurement))
    for i=1:length(TPs.measurement)
        if column[:pout] == "atmospheric"
            pout = PPs[i, end]
        else
            pout = "vacuum"
        end
        prog[i] = Program(collect(skipmissing(TPs[i, 2:end])), collect(skipmissing(PPs[i, 2:(end-1)])), column[:L]; pout=pout, time_unit=column[:time_unit])
    end
    #method_short_names = Array{String}(undef, ns)
    Tchar = Array{Float64}(undef, ns)
	θchar = Array{Float64}(undef, ns)
	ΔCp = Array{Float64}(undef, ns)
    min = Array{Float64}(undef, ns)
    retcode = Array{Any}(undef, ns)
    if mode == "single"
	    sol = optimize_Kcentric_single(tR_meas.*a, column[:L], column[:d], column[:gas], prog, options, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        for j=1:ns
            #method_short_names[j] = method_short_name
            Tchar[j] = sol[j][1]
            θchar[j] = sol[j][2]
            ΔCp[j] = sol[j][3]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
    else
        sol = optimize_Kcentric_all(tR_meas.*a, column[:L], column[:d], column[:gas], prog, options, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        Tchar = sol[1:ns] # Array length = number solutes
        θchar = sol[ns+1:2*ns] # Array length = number solutes
        ΔCp = sol[2*ns+1:3*ns] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
            #method_short_names[j] = method_short_name
        end
    end
	df = DataFrame(Name=solute_names, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
	return df, sol
end

# optimization for retention parameters
# DELETE
function optimize(tR_meas, solute_names, column, options, TPs, PPs, method; maxiters=10000, relbound=0.5, mode="single")
    if column[:time_unit] == "min"
        a = 60.0
    else
        a = 1.0
    end
    Tchar_e, θchar_e, ΔCp_e = estimate_start_parameter(tR_meas, TPs, PPs, column[:L], column[:d], column[:gas]; pout=column[:pout], time_unit=column[:time_unit], control=options.control)
    df, sol = optimize(tR_meas, solute_names, column, options, TPs, PPs,
                        Tchar_e, θchar_e, ΔCp_e,
                        Tchar_e.*relbound, θchar_e.*relbound, ΔCp_e.*relbound,
                        Tchar_e.*(1+relbound), θchar_e.*(1+relbound), ΔCp_e.*(1+relbound),
                        method; maxiters=maxiters, mode=mode)
    return df, sol
end

#------------------
# Optimization for L, d and K-centric retention parameters together
function opt_LdKcentric(x_LdKcentric, p)
	tR = p[1]
	prog = p[2]
	opt = p[3]
    gas = p[4]
    if length(size(tR)) == 1
		ns = 1
	else
		ns = size(tR)[2]
	end
	L = x_LdKcentric[1]
	d = x_LdKcentric[2]
	Tchar = x_LdKcentric[3:ns+2] # Array length = number solutes
	θchar = x_LdKcentric[ns+2+1:2*ns+2] # Array length = number solutes
	ΔCp = x_LdKcentric[2*ns+2+1:3*ns+2] # Array length = number solutes
    return loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas)[1]
end

function opt_LKcentric(x_LKcentric, p)
	tR = p[1]
	d = p[2]
	prog = p[3]
	opt = p[4]
    gas = p[5]
    if length(size(tR)) == 1
		ns = 1
	else
		ns = size(tR)[2]
	end
	L = x_LKcentric[1]
	Tchar = x_LKcentric[2:ns+1] # Array length = number solutes
	θchar = x_LKcentric[ns+1+1:2*ns+1] # Array length = number solutes
	ΔCp = x_LKcentric[2*ns+1+1:3*ns+1] # Array length = number solutes
    return RetentionParameterEstimator.loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas)[1]
end

function opt_dKcentric(x_dKcentric, p)
	tR = p[1]
	L = p[2]
	prog = p[3]
	opt = p[4]
    gas = p[5]
    if length(size(tR)) == 1
		ns = 1
	else
		ns = size(tR)[2]
	end
	d = x_dKcentric[1]
	Tchar = x_dKcentric[2:ns+1] # Array length = number solutes
	θchar = x_dKcentric[ns+1+1:2*ns+1] # Array length = number solutes
	ΔCp = x_dKcentric[2*ns+1+1:3*ns+1] # Array length = number solutes
    return RetentionParameterEstimator.loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas)[1]
end

function optimize_LdKcentric(tR, gas, prog, opt, L_e, d_e, Tchar_e, θchar_e, ΔCp_e, lb_L, lb_d, lb_Tchar, lb_θchar, lb_ΔCp, ub_L, ub_d, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=10000)
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
                BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    p = [tR, prog, opt, gas]
	x0 = [L_e; d_e; Tchar_e; θchar_e; ΔCp_e]
	lb = [lb_L; lb_d; lb_Tchar; lb_θchar; lb_ΔCp]
	ub = [ub_L; ub_d; ub_Tchar; ub_θchar; ub_ΔCp]
	optf = OptimizationFunction(opt_LdKcentric, Optimization.AutoForwardDiff())
	if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
		prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
	else
		prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
	end
	if method in optimisers
        opt_sol = solve(prob, method, maxiters=maxiters)
    #elseif method in bbos
    #    opt_sol = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
    else
        opt_sol = solve(prob, method)
    end
	return opt_sol
end

function optimize_LKcentric(tR, d, gas, prog, opt, L_e, Tchar_e, θchar_e, ΔCp_e, lb_L, lb_Tchar, lb_θchar, lb_ΔCp, ub_L, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=10000)
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
                BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    p = [tR, d, prog, opt, gas]
	x0 = [L_e; Tchar_e; θchar_e; ΔCp_e]
	lb = [lb_L; lb_Tchar; lb_θchar; lb_ΔCp]
	ub = [ub_L; ub_Tchar; ub_θchar; ub_ΔCp]
	optf = OptimizationFunction(opt_LKcentric, Optimization.AutoForwardDiff())
	if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
		prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
	else
		prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
	end
	if method in optimisers
        opt_sol = solve(prob, method, maxiters=maxiters)
    #elseif method in bbos
    #    opt_sol = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
    else
        opt_sol = solve(prob, method)
    end
	return opt_sol
end

function optimize_dKcentric(tR, L, gas, prog, opt, d_e, Tchar_e, θchar_e, ΔCp_e, lb_d, lb_Tchar, lb_θchar, lb_ΔCp, ub_d, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=10000)
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
                BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    p = [tR, L, prog, opt, gas]
	x0 = [d_e; Tchar_e; θchar_e; ΔCp_e]
	lb = [lb_d; lb_Tchar; lb_θchar; lb_ΔCp]
	ub = [ub_d; ub_Tchar; ub_θchar; ub_ΔCp]
	optf = OptimizationFunction(opt_dKcentric, Optimization.AutoForwardDiff())
	if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
		prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
	else
		prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
	end
	if method in optimisers
        opt_sol = solve(prob, method, maxiters=maxiters)
    #elseif method in bbos
    #    opt_sol = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
    else
        opt_sol = solve(prob, method)
    end
	return opt_sol
end

# rename this function later
function optimize_all(tR_meas, solute_names, column, options, TPs, PPs, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=10000, mode="LdKcentric")
    # mode = "LdKcentric", "LKcentric", "dKcentric", "Kcentric", "Kcentric_single"
    if column[:time_unit] == "min"
        a = 60.0
    else
        a = 1.0
    end

    if length(size(tR_meas)) == 1
        ns = 1
    else
        ns = size(tR_meas)[2]
    end

	prog = Array{GasChromatographySimulator.Program}(undef, length(TPs.measurement))
    for i=1:length(TPs.measurement)
        if column[:pout] == "atmospheric"
            pout = PPs[i, end]
        else
            pout = "vacuum"
        end
        prog[i] = Program(collect(skipmissing(TPs[i, 2:end])), collect(skipmissing(PPs[i, 2:(end-1)])), column[:L]; pout=pout, time_unit=column[:time_unit])
    end
    
    Tchar = Array{Float64}(undef, ns)
	θchar = Array{Float64}(undef, ns)
	ΔCp = Array{Float64}(undef, ns)
    min = Array{Float64}(undef, ns)
    retcode = Array{Any}(undef, ns)
    if mode == "Kcentric_single"
	    sol = optimize_Kcentric_single(tR_meas.*a, column[:L], column[:d], column[:gas], prog, options, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        for j=1:ns
            Tchar[j] = sol[j][1]
            θchar[j] = sol[j][2]
            ΔCp[j] = sol[j][3]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
    elseif mode == "Kcentric"
        sol = optimize_Kcentric_all(tR_meas.*a, column[:L], column[:d], column[:gas], prog, options, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        Tchar = sol[1:ns] # Array length = number solutes
        θchar = sol[ns+1:2*ns] # Array length = number solutes
        ΔCp = sol[2*ns+1:3*ns] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
    elseif mode == "LdKcentric"
        L_e = column[:L]
        d_e = column[:d]
        lb_L = L_e/100
        ub_L = L_e*100
        lb_d = d_e/100
        ub_d = d_e*100
        sol = optimize_LdKcentric(tR_meas.*a, column[:gas], prog, options, L_e, d_e, Tchar_e, θchar_e, ΔCp_e, lb_L, lb_d, lb_Tchar, lb_θchar, lb_ΔCp, ub_L, ub_d, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        L = sol[1].*ones(ns)
        d = sol[2].*ones(ns)
        Tchar = sol[3:ns+2] # Array length = number solutes
        θchar = sol[ns+2+1:2*ns+2] # Array length = number solutes
        ΔCp = sol[2*ns+2+1:3*ns+2] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, L=L, d=d, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
    elseif mode == "LKcentric"
        L_e = column[:L]
        lb_L = L_e/100
        ub_L = L_e*100
        sol = optimize_LKcentric(tR_meas.*a, column[:d], column[:gas], prog, options, L_e, Tchar_e, θchar_e, ΔCp_e, lb_L, lb_Tchar, lb_θchar, lb_ΔCp, ub_L, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        L = sol[1].*ones(ns)
        Tchar = sol[2:ns+1] # Array length = number solutes
        θchar = sol[ns+1+1:2*ns+1] # Array length = number solutes
        ΔCp = sol[2*ns+1+1:3*ns+1] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, L=L, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
    elseif mode == "dKcentric"
        d_e = column[:d]
        lb_d = d_e/100
        ub_d = d_e*100
        sol = optimize_dKcentric(tR_meas.*a, column[:L], column[:gas], prog, options, d_e, Tchar_e, θchar_e, ΔCp_e, lb_d, lb_Tchar, lb_θchar, lb_ΔCp, ub_d, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        d = sol[1].*ones(ns)
        Tchar = sol[2:ns+1] # Array length = number solutes
        θchar = sol[ns+1+1:2*ns+1] # Array length = number solutes
        ΔCp = sol[2*ns+1+1:3*ns+1] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, d=d, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
    end
	
	return df, sol
end

function optimize_all_(tR_meas, solute_names, column, options, TPs, PPs, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=10000, mode="LdKcentric")
    # mode = "LdKcentric", "LKcentric", "dKcentric", "Kcentric", "Kcentric_single"
    if column[:time_unit] == "min"
        a = 60.0
    else
        a = 1.0
    end

    if length(size(tR_meas)) == 1
        ns = 1
    else
        ns = size(tR_meas)[2]
    end

	prog = Array{GasChromatographySimulator.Program}(undef, length(TPs.measurement))
    for i=1:length(TPs.measurement)
        if column[:pout] == "atmospheric"
            pout = PPs[i, end]
        else
            pout = "vacuum"
        end
        prog[i] = Program(collect(skipmissing(TPs[i, 2:end])), collect(skipmissing(PPs[i, 2:(end-1)])), column[:L]; pout=pout, time_unit=column[:time_unit])
    end
    
    Tchar = Array{Float64}(undef, ns)
	θchar = Array{Float64}(undef, ns)
	ΔCp = Array{Float64}(undef, ns)
    min = Array{Float64}(undef, ns)
    retcode = Array{Any}(undef, ns)
    if mode == "Kcentric_single"
	    sol = optimize_Kcentric_single_(tR_meas.*a, column[:L], column[:d], column[:df], column[:gas], prog, options, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        for j=1:ns
            Tchar[j] = sol[j][1]
            θchar[j] = sol[j][2]
            ΔCp[j] = sol[j][3]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
    elseif mode == "Kcentric"
        sol = optimize_Kcentric_all_(tR_meas.*a, column[:L], column[:d], column[:gas], prog, options, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        Tchar = sol[1:ns] # Array length = number solutes
        θchar = sol[ns+1:2*ns] # Array length = number solutes
        ΔCp = sol[2*ns+1:3*ns] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
    elseif mode == "LdKcentric"
        L_e = column[:L]
        d_e = column[:d]
        lb_L = L_e/100
        ub_L = L_e*100
        lb_d = d_e/100
        ub_d = d_e*100
        sol = optimize_LdKcentric_(tR_meas.*a, column[:gas], prog, options, L_e, d_e, Tchar_e, θchar_e, ΔCp_e, lb_L, lb_d, lb_Tchar, lb_θchar, lb_ΔCp, ub_L, ub_d, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        L = sol[1].*ones(ns)
        d = sol[2].*ones(ns)
        Tchar = sol[3:ns+2] # Array length = number solutes
        θchar = sol[ns+2+1:2*ns+2] # Array length = number solutes
        ΔCp = sol[2*ns+2+1:3*ns+2] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, L=L, d=d, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
    elseif mode == "LKcentric"
        L_e = column[:L]
        lb_L = L_e/100
        ub_L = L_e*100
        sol = optimize_LKcentric_(tR_meas.*a, column[:d], column[:gas], prog, options, L_e, Tchar_e, θchar_e, ΔCp_e, lb_L, lb_Tchar, lb_θchar, lb_ΔCp, ub_L, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        L = sol[1].*ones(ns)
        Tchar = sol[2:ns+1] # Array length = number solutes
        θchar = sol[ns+1+1:2*ns+1] # Array length = number solutes
        ΔCp = sol[2*ns+1+1:3*ns+1] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, L=L, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
    elseif mode == "dKcentric"
        d_e = column[:d]
        lb_d = d_e/100
        ub_d = d_e*100
        sol = optimize_dKcentric_(tR_meas.*a, column[:L], column[:gas], prog, options, d_e, Tchar_e, θchar_e, ΔCp_e, lb_d, lb_Tchar, lb_θchar, lb_ΔCp, ub_d, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=maxiters)
        d = sol[1].*ones(ns)
        Tchar = sol[2:ns+1] # Array length = number solutes
        θchar = sol[ns+1+1:2*ns+1] # Array length = number solutes
        ΔCp = sol[2*ns+1+1:3*ns+1] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, d=d, Tchar=Tchar, θchar=θchar, ΔCp=ΔCp, min=min, retcode=retcode)
    end
	
	return df, sol
end