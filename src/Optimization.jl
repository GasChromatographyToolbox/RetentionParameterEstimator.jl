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
function opt_Kcentric(x, p)
	tR = p[1]
	L = p[2]
	d = p[3]
	prog = p[4]
	opt = p[5]
    gas = p[6]
    metric = p[7]
    if length(size(tR)) == 1
		ns = 1
	else
		ns = size(tR)[2]
	end
	Tchar = x[1:ns] # Array length = number solutes
	θchar = x[ns+1:2*ns] # Array length = number solutes
	ΔCp = x[2*ns+1:3*ns] # Array length = number solutes
    return loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas; metric=metric)[1]
end

function opt_dKcentric(x, p)
    tR = p[1]
    L = p[2]
    prog = p[3]
    opt = p[4]
    gas = p[5]
    metric = p[6]
    if length(size(tR)) == 1
        ns = 1
    else
        ns = size(tR)[2]
    end
    d = x[1]
    Tchar = x[2:ns+1] # Array length = number solutes
    θchar = x[ns+1+1:2*ns+1] # Array length = number solutes
    ΔCp = x[2*ns+1+1:3*ns+1] # Array length = number solutes
    return loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas; metric=metric)
end

function opt_d(x, p)
    tR = p[1]
    L = p[2]
    Tchar = p[3]
    θchar = p[4]
    ΔCp = p[5]
    prog = p[6]
    opt = p[7]
    gas = p[8]
    metric = p[9]
    if length(size(tR)) == 1
        ns = 1
    else
        ns = size(tR)[2]
    end
    d = x[1]
    return loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas; metric=metric)
end

function opt_ABC(x, p)
	tR = p[1]
	L = p[2]
	d = p[3]
    β = p[4]
	prog = p[5]
	opt = p[6]
    gas = p[7]
    metric = p[8]
    if length(size(tR)) == 1
		ns = 1
	else
		ns = size(tR)[2]
	end
	A = x[1:ns] # Array length = number solutes
	B = x[ns+1:2*ns] # Array length = number solutes
	C = x[2*ns+1:3*ns] # Array length = number solutes
    return loss(tR, A, B, C, L, d, β, prog, opt, gas; metric=metric)[1]
end

function opt_dABC(x, p)
    tR = p[1]
    L = p[2]
    β = p[3]
    prog = p[4]
    opt = p[5]
    gas = p[6]
    metric = p[7]
    if length(size(tR)) == 1
        ns = 1
    else
        ns = size(tR)[2]
    end
    d = x[1]
    A = x[2:ns+1] # Array length = number solutes
    B = x[ns+1+1:2*ns+1] # Array length = number solutes
    C = x[2*ns+1+1:3*ns+1] # Array length = number solutes
    return loss(tR, A, B, C, L, d, β, prog, opt, gas; metric=metric)[1]
end

function opt_βABC(x, p)
    tR = p[1]
    L = p[2]
    d = p[3]
    prog = p[4]
    opt = p[5]
    gas = p[6]
    metric = p[7]
    if length(size(tR)) == 1
        ns = 1
    else
        ns = size(tR)[2]
    end
    β = x[1]
    A = x[2:ns+1] # Array length = number solutes
    B = x[ns+1+1:2*ns+1] # Array length = number solutes
    C = x[2*ns+1+1:3*ns+1] # Array length = number solutes
    return loss(tR, A, B, C, L, d, β, prog, opt, gas; metric=metric)[1]
end

function opt_β(x, p)
    tR = p[1]
    L = p[2]
    d = p[3]
    A = p[4]
    B = p[5]
    C = p[6]
    prog = p[7]
    opt = p[8]
    gas = p[9]
    metric = p[10]
    if length(size(tR)) == 1
        ns = 1
    else
        ns = size(tR)[2]
    end
    β = x[1]
    return loss(tR, A, B, C, L, d, β, prog, opt, gas; metric=metric)[1]
end

function opt_dβABC(x, p)
    tR = p[1]
    L = p[2]
    prog = p[3]
    opt = p[4]
    gas = p[5]
    metric = p[6]
    if length(size(tR)) == 1
        ns = 1
    else
        ns = size(tR)[2]
    end
    d = x[1]
    β = x[2]
    A = x[3:ns+2] # Array length = number solutes
    B = x[ns+2+1:2*ns+2] # Array length = number solutes
    C = x[2*ns+2+1:3*ns+2] # Array length = number solutes
    return loss(tR, A, B, C, L, d, β, prog, opt, gas; metric=metric)[1]
end

function opt_dβ(x, p)
    tR = p[1]
    L = p[2]
    A = p[3]
    B = p[4]
    C = p[5]
    prog = p[6]
    opt = p[7]
    gas = p[8]
    metric = p[9]
    if length(size(tR)) == 1
        ns = 1
    else
        ns = size(tR)[2]
    end
    d = x[1]
    β = x[2]
    return loss(tR, A, B, C, L, d, β, prog, opt, gas; metric=metric)[1]
end

# optimize every solute separatly, tR is a 2D-array with RT of different programs in the first dimension and different solutes in the second dimension  
#function optimize_Kcentric_single(tR, L, d, gas, prog, opt, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=10000)
#	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
#                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
#                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
#    bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
#                BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
#    
#    optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
#	
#    if typeof(size(tR)) == Tuple{Int64, Int64}
#        n2 = size(tR)[2]
#    else
#        n2 = 1
#    end 
#    opt_sol = Array{Any}(undef, n2)
#	for i=1:n2
#		p = [tR[:,i], L, d, prog, opt, gas]
#		x0 = [Tchar_e[i], θchar_e[i], ΔCp_e[i]]
#		lb = [lb_Tchar[i], lb_θchar[i], lb_ΔCp[i]]
#		ub = [ub_Tchar[i], ub_θchar[i], ub_ΔCp[i]]
#		if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
#			prob = OptimizationProblem(optf, x0, p)
#		else
#			prob = OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
#		end
#        if method in optimisers
#            opt_sol[i] = solve(prob, method, maxiters=maxiters)
#        #elseif method in bbos
#        #    opt_sol[i] = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
#        else
#		    opt_sol[i] = solve(prob, method) #-> :u (Array of the optimized parameters), :minimum (minima of the optimization function) , :retcode (Boolean, successful?)
#        end
#    end
#	return opt_sol
#end

#function optimize_Kcentric_single(tR, L, d, gas, prog, opt, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=10000, metric="quadratic")
#	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
#                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
#                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
#    #bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
#    #            BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
#    
#    optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
#	
#    if typeof(size(tR)) == Tuple{Int64, Int64}
#        n2 = size(tR)[2]
#    else
#        n2 = 1
#    end 
#    opt_sol = Array{Any}(undef, n2)
#	for i=1:n2
#		p = [tR[:,i], L, d, prog, opt, gas, metric]
#		x0 = [Tchar_e[i], θchar_e[i], ΔCp_e[i]]
#		lb = [lb_Tchar[i], lb_θchar[i], lb_ΔCp[i]]
#		ub = [ub_Tchar[i], ub_θchar[i], ub_ΔCp[i]]
#		#if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
#		if method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton())
#            prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
#        else
#            prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
#        end
#        #if method in optimisers
#        #    opt_sol[i] = solve(prob, method, maxiters=maxiters)
#        #elseif method in bbos
#        #    opt_sol[i] = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
#        #else#
#		opt_sol[i] = solve(prob, method, maxiters=maxiters) #-> :u (Array of the optimized parameters), :minimum (minima of the optimization function) , :retcode (Boolean, successful?)
#        #end
#    end
#	return opt_sol
#end

#=function optimize_Kcentric_single(tR, L, d, gas, prog, opt, Tchar_e, θchar_e, ΔCp_e, method; maxiters=10000, metric="quadratic")
    
    optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
	
    if typeof(size(tR)) == Tuple{Int64, Int64}
        ns = size(tR)[2]
    else
        ns = 1
    end 
    opt_sol = Array{SciMLBase.OptimizationSolution}(undef, ns)
	for i=1:ns
		p = [tR[:,i], L, d, prog, opt, gas, metric]
		x0 = [Tchar_e[i], θchar_e[i], ΔCp_e[i]]
		
		prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
		opt_sol[i] = solve(prob, method, maxiters=maxiters)
    end
	return opt_sol
end=#

#=
function optimize_ABC_single(tR, L, λ, β, gas, prog, opt, Tchar_e, θchar_e, ΔCp_e, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp, method; maxiters=10000, metric="quadratic")
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    #bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
    #            BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    optf = OptimizationFunction(opt_ABC, Optimization.AutoForwardDiff())
	
    if typeof(size(tR)) == Tuple{Int64, Int64}
        n2 = size(tR)[2]
    else
        n2 = 1
    end 
    opt_sol = Array{Any}(undef, n2)
	for i=1:n2
		p = [tR[:,i], L, λ, β, prog, opt, gas, metric]
		x0 = [Tchar_e[i], θchar_e[i], ΔCp_e[i]]
		lb = [lb_Tchar[i], lb_θchar[i], lb_ΔCp[i]]
		ub = [ub_Tchar[i], ub_θchar[i], ub_ΔCp[i]]
		#if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
		if method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton())
            prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
        else
            prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
        end
        #if method in optimisers
        #    opt_sol[i] = solve(prob, method, maxiters=maxiters)
        #elseif method in bbos
        #    opt_sol[i] = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
        #else
		opt_sol[i] = solve(prob, method, maxiters=maxiters) #-> :u (Array of the optimized parameters), :minimum (minima of the optimization function) , :retcode (Boolean, successful?)
        #end
    end
	return opt_sol
end
=#

# optimization of only one solute
#function optimize_Kcentric_single(tR, L, d, gas, prog, opt, Tchar_e::Float64, θchar_e::Float64, ΔCp_e::Float64, lb_Tchar::Float64, lb_θchar::Float64, lb_ΔCp::Float64, ub_Tchar::Float64, ub_θchar::Float64, ub_ΔCp::Float64, method; maxiters=10000)
#	opt_sol[1] = optimize_Kcentric_single(tR, L, d, gas, prog, opt, [Tchar_e], [θchar_e], [ΔCp_e], [lb_Tchar], [lb_θchar], [lb_ΔCp], [ub_Tchar], [ub_θchar], [ub_ΔCp], method; maxiters=maxiters)
#end

# optimize all solutes together
#=
function optimize_Kcentric_all(tR, L, d, gas, prog, opt, A_e, B_e, C_e, lb_A, lb_B, lb_C, ub_A, ub_B, ub_C, method; maxiters=10000, metric="quadratic")
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    #bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
    #            BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    p = [tR, L, d, prog, opt, gas, metric]
	x0 = [A_e; B_e; C_e]
	lb = [lb_A; lb_B; lb_C]
	ub = [ub_A; ub_B; ub_C]
	optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
	#if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
	if method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton())
		prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
	else
		prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
	end
	#if method in optimisers
    #    opt_sol = solve(prob, method, maxiters=maxiters)
    #elseif method in bbos
    #    opt_sol = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
    #else
    opt_sol = solve(prob, method, maxiters=maxiters)
    #end
	return opt_sol
end
=#

"""
    optimize_Kcentric(tR, L, gas, prog, opt, d_e, Tchar_e, θchar_e, ΔCp_e, method; maxiters=10000, metric="quadratic") 

Optimization regarding the estimization of the retention parameters `Tchar`, `θchar` and `ΔCp`. The initial guess is a vector (`Tchar_e`, `θchar_e` and `ΔCp_e`) for optimization
algorithms, which do not need lower/upper bounds. It should be a matrix (`Tchar_e`, `θchar_e` and `ΔCp_e`), with the first column (`[1,:]`)
beeing the initial guess, the second column (`[2,:]`) the lower bound and the third column (`[3,:]`) the upper bound.
"""
function optimize_Kcentric(tR, L, d, gas, prog, opt, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}, method; maxiters=10000, metric="quadratic", g_tol=1e-4) where T<:Number
    p = [tR, L, d, prog, opt, gas, metric]
	x0 = [Tchar_e; θchar_e; ΔCp_e]
	optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters, g_tol=g_tol)
    opt_sol = solve(prob, method, maxiters=maxiters)
	return opt_sol
end

function optimize_Kcentric(tR, L, d, gas, prog, opt, Tchar_e::Matrix{T}, θchar_e::Matrix{T}, ΔCp_e::Matrix{T}, method; maxiters=10000, metric="quadratic", g_tol=1e-4) where T<:Number
    p = [tR, L, d, prog, opt, gas, metric]
	x0 = [Tchar_e[1,:]; θchar_e[1,:]; ΔCp_e[1,:]]
    lb = [Tchar_e[2,:]; θchar_e[2,:]; ΔCp_e[2,:]]
    ub = [Tchar_e[3,:]; θchar_e[3,:]; ΔCp_e[3,:]]
	optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub, f_calls_limit=maxiters, g_tol=g_tol)
    opt_sol = solve(prob, method, maxiters=maxiters)
	return opt_sol
end

#=
function optimize_ABC_all(tR, L, λ, β, gas, prog, opt, A_e, B_e, C_e, lb_A, lb_B, lb_C, ub_A, ub_B, ub_C, method; maxiters=10000, metric="quadratic")
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    #bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
    #            BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    p = [tR, L, λ, β, prog, opt, gas, metric]
	x0 = [A_e; B_e; C_e]
	lb = [lb_A; lb_B; lb_C]
	ub = [ub_A; ub_B; ub_C]
	optf = OptimizationFunction(opt_ABC, Optimization.AutoForwardDiff())
	#if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
	if method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton())
		prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
	else
		prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
	end
	#if method in optimisers
    #    opt_sol = solve(prob, method, maxiters=maxiters)
    #elseif method in bbos
    #    opt_sol = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
    #else
    opt_sol = solve(prob, method, maxiters=maxiters)
    #end
	return opt_sol
end
=#


#------------------

"""
    optimize_dKcentric(tR, L, gas, prog, opt, d_e, Tchar_e, θchar_e, ΔCp_e, method; maxiters=10000, metric="quadratic") 

Optimization regarding the estimization of the column diameter `d` and the retention parameters `Tchar`, `θchar` and `ΔCp`. The initial guess is a number (`d_e`) or a vector (`Tchar_e`, `θchar_e` and `ΔCp_e`) for optimization
algorithms, which do not need lower/upper bounds. It should be a vector of length 3 (`d_e`) or a matrix (`Tchar_e`, `θchar_e` and `ΔCp_e`), with the first element/column 
beeing the initial guess, the second element/column the lower bound and the third element/column the upper bound.
"""
function optimize_dKcentric(tR, L, gas, prog, opt, d_e::Number, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}, method; maxiters=10000, metric="quadratic") where T<:Number
    p = [tR, L, prog, opt, gas, metric]
	x0 = [d_e; Tchar_e; θchar_e; ΔCp_e]
	optf = OptimizationFunction(opt_dKcentric, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters)
	return opt_sol
end

function optimize_dKcentric(tR, L, gas, prog, opt, d_e::Vector{T}, Tchar_e::Matrix{T}, θchar_e::Matrix{T}, ΔCp_e::Matrix{T}, method; maxiters=10000, metric="quadratic") where T<:Number    
    p = [tR, L, prog, opt, gas, metric]
	x0 = [d_e[1]; Tchar_e[1,:]; θchar_e[1,:]; ΔCp_e[1,:]]
    lb = [d_e[2]; Tchar_e[2,:]; θchar_e[2,:]; ΔCp_e[2,:]]
    ub = [d_e[3]; Tchar_e[3,:]; θchar_e[3,:]; ΔCp_e[3,:]]
	optf = OptimizationFunction(opt_dKcentric, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters)
	return opt_sol
end

"""
    optimize_d(tR, L, Tchar, θchar, ΔCp, gas, prog, opt, d_e, method; maxiters=10000, metric="squared")

Optimization regarding the estimization of the column diameter `d`. The initial guess `d_e`, is a number for optimization algorithms, which do not need lower/upper bounds. It should be a vector of length 3, with the first element 
beeing the initial guess, the second element the lower bound and the third element the upper bound for `d`.
"""
function optimize_d(tR, L, Tchar, θchar, ΔCp, gas, prog, opt, d_e::Number, method; maxiters=10000, metric="squared")    
    p = [tR, L, Tchar, θchar, ΔCp, prog, opt, gas, metric]
	x0 = [d_e]
	optf = OptimizationFunction(opt_d, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters)
	return opt_sol
end

function optimize_d(tR, L, Tchar, θchar, ΔCp, gas, prog, opt, d_e::Vector{T}, method; maxiters=10000, metric="squared") where T<:Number
    p = [tR, L, Tchar, θchar, ΔCp, prog, opt, gas, metric]
	x0 = [d_e[1]]
    lb = [d_e[2]]
    ub = [d_e[3]]
	optf = OptimizationFunction(opt_d, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters)
	return opt_sol
end

#=
function optimize_λABC(tR, L, β, gas, prog, opt, λ_e, A_e, B_e, C_e, lb_λ, lb_A, lb_B, lb_C, ub_λ, ub_A, ub_B, ub_C, method; maxiters=10000, metric="quadratic")
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    #bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
    #            BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    p = [tR, L, β, prog, opt, gas, metric]
	x0 = [λ_e; A_e; B_e; C_e]
	lb = [lb_λ; lb_A; lb_B; lb_C]
	ub = [ub_λ; ub_A; ub_B; ub_C]
	optf = OptimizationFunction(opt_λABC, Optimization.AutoForwardDiff())
	if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
		prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
	else
		prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
	end
	#if method in optimisers
    #    opt_sol = solve(prob, method, maxiters=maxiters)
    #elseif method in bbos
    #    opt_sol = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
    #else
    opt_sol = solve(prob, method, maxiters=maxiters)
    #end
	return opt_sol
end
=#

#=
function optimize_βABC(tR, L, λ, gas, prog, opt, β_e, A_e, B_e, C_e, lb_β, lb_A, lb_B, lb_C, ub_β, ub_A, ub_B, ub_C, method; maxiters=10000, metric="quadratic")
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    #bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
    #            BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    p = [tR, L, λ, prog, opt, gas, metric]
	x0 = [β_e; A_e; B_e; C_e]
	lb = [lb_β; lb_A; lb_B; lb_C]
	ub = [ub_β; ub_A; ub_B; ub_C]
	optf = OptimizationFunction(opt_βABC, Optimization.AutoForwardDiff())
	if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
		prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
	else
		prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
	end
	#if method in optimisers
    #    opt_sol = solve(prob, method, maxiters=maxiters)
    #elseif method in bbos
    #    opt_sol = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
    #else
    opt_sol = solve(prob, method, maxiters=maxiters)
    #end
	return opt_sol
end
=#

#=
function optimize_λβABC(tR, L, gas, prog, opt, λ_e, β_e, A_e, B_e, C_e, lb_λ, lb_β, lb_A, lb_B, lb_C, ub_λ, ub_β, ub_A, ub_B, ub_C, method; maxiters=10000, metric="quadratic")
	optimisers = [ Optimisers.Descent(), Optimisers.Momentum(), Optimisers.Nesterov(), Optimisers.RMSProp(), Optimisers.Adam(),
                    Optimisers.RAdam(), Optimisers.OAdam(), Optimisers.AdaMax(), Optimisers.ADAGrad(), Optimisers.ADADelta(),
                    Optimisers.AMSGrad(), Optimisers.NAdam(), Optimisers.AdamW()]
    #bbos = [BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_separable_nes(), BBO_xnes(), BBO_dxnes(), BBO_adaptive_de_rand_1_bin(), BBO_de_rand_1_bin(),
    #            BBO_de_rand_1_bin_radiuslimited(), BBO_de_rand_2_bin(), BBO_de_rand_2_bin_radiuslimited()]
    
    p = [tR, L, prog, opt, gas, metric]
	x0 = [λ_e; β_e; A_e; B_e; C_e]
	lb = [lb_λ; lb_β; lb_A; lb_B; lb_C]
	ub = [ub_λ; ub_β; ub_A; ub_B; ub_C]
	optf = OptimizationFunction(opt_λβABC, Optimization.AutoForwardDiff())
	if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton()) || method in optimisers
		prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
	else
		prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
	end
	#if method in optimisers
    #    opt_sol = solve(prob, method, maxiters=maxiters)
    #elseif method in bbos
    #    opt_sol = solve(prob, method, maxiters=maxiters, TraceMode=:silent)
    #else
    opt_sol = solve(prob, method, maxiters=maxiters)
    #end
	return opt_sol
end
=#


# rename this function later
function estimate_parameters(tR_meas, solute_names, column, options, TPs, PPs, rp1_e, rp2_e, rp3_e, method; maxiters=10000, mode="dKcentric", metric="quadratic", φ₀=1e-3)
    # mode = "Kcentric", "Kcentric_single", "dKcentric", "dKcentric_single", "d"
    # later add here modes for ABC-Parameters and ABC+df
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
    
    d_e = column[:d]

    rp1 = Array{Float64}(undef, ns)
	rp2 = Array{Float64}(undef, ns)
	rp3 = Array{Float64}(undef, ns)
    min = Array{Float64}(undef, ns)
    retcode = Array{Any}(undef, ns)
    if mode == "Kcentric_single"
        sol = Array{SciMLBase.OptimizationSolution}(undef, ns)
        for j=1:ns
            sol[j] = optimize_Kcentric(tR_meas[:,j].*a, column[:L], column[:d], column[:gas], prog, options, rp1_e[j,:], rp2_e[j,:], rp3_e[j,:], method; maxiters=maxiters, metric=metric)
            rp1[j] = sol[j][1]
            rp2[j] = sol[j][2]
            rp3[j] = sol[j][3]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "Kcentric"
        sol = optimize_Kcentric(tR_meas.*a, column[:L], column[:d], column[:gas], prog, options, rp1_e, rp2_e, rp3_e, method; maxiters=maxiters, metric=metric)
        rp1 = sol[1:ns] # Array length = number solutes
        rp2 = sol[ns+1:2*ns] # Array length = number solutes
        rp3 = sol[2*ns+1:3*ns] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "dKcentric"
        sol = optimize_dKcentric(tR_meas.*a, column[:L], column[:gas], prog, options, d_e, rp1_e, rp2_e, rp3_e, method; maxiters=maxiters, metric=metric)
        d = sol[1].*ones(ns)
        rp1 = sol[2:ns+1] # Array length = number solutes
        rp2 = sol[ns+1+1:2*ns+1] # Array length = number solutes
        rp3 = sol[2*ns+1+1:3*ns+1] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, d=d, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "dKcentric_single"
        sol = Array{SciMLBase.OptimizationSolution}(undef, ns)
        d = Array{Float64}(undef, ns)
        for j=1:ns
            sol[j] = optimize_dKcentric(tR_meas[:,j].*a, column[:L], column[:gas], prog, options, d_e, rp1_e[j,:], rp2_e[j,:], rp3_e[j,:], method; maxiters=maxiters, metric=metric)
            d[j] = sol[j][1]
            rp1[j] = sol[j][2]
            rp2[j] = sol[j][3]
            rp3[j] = sol[j][4]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, d=d, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "d"
        sol = optimize_d(tR_meas.*a, column[:L], rp1_e, rp2_e, rp3_e, column[:gas], prog, options, d_e, method; maxiters=maxiters, metric=metric)
        φ = sol[1].*ones(ns)
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, d=d, min=min, retcode=retcode)
    end
	
	return df, sol
end

#=
function estimate_parameters(tR_meas, solute_names, column, options, TPs, PPs, rp1_e, rp2_e, rp3_e, lb_rp1, lb_rp2, lb_rp3, ub_rp1, ub_rp2, ub_rp3, method; maxiters=10000, mode="LdKcentric", metric="quadratic", φ₀=1e-3)
    # mode = "Kcentric", "Kcentric_single", "λKcentric", "λKcentric_single", "φKcentric", "φKcentric_single", "λφKcentric", "λφKcentric_single", "λφ", "ddfKcentric", "ddfKcentric_single", "ABC", "ABC_single", "λABC", "dfABC", "λdfABC"
    # later add here modes for ABC-Parameters and ABC+df
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
    
    λ_e = column[:L]/column[:d]
    lb_λ = λ_e/10
    ub_λ = λ_e*10
    φ_e = column[:df]/column[:d]
    lb_φ = φ_e/10
    ub_φ = φ_e*10
    d_e = column[:d]
    lb_d = d_e/10
    ub_d = d_e*10
    df_e = column[:df]
    lb_df = df_e/10
    ub_df = df_e*10
    β_e = column[:d]/(4*column[:df])
    lb_β = β_e/10
    ub_β = β_e*10

    rp1 = Array{Float64}(undef, ns)
	rp2 = Array{Float64}(undef, ns)
	rp3 = Array{Float64}(undef, ns)
    min = Array{Float64}(undef, ns)
    retcode = Array{Any}(undef, ns)
    if mode == "Kcentric_single"
	    sol = optimize_Kcentric_single(tR_meas.*a, column[:L], column[:d], column[:gas], prog, options, rp1_e, rp2_e, rp3_e, lb_rp1, lb_rp2, lb_rp3, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        for j=1:ns
            rp1[j] = sol[j][1]
            rp2[j] = sol[j][2]
            rp3[j] = sol[j][3]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "Kcentric"
        sol = optimize_Kcentric_all(tR_meas.*a, column[:L], column[:d], column[:gas], prog, options, rp1_e, rp2_e, rp3_e, lb_rp1, lb_rp2, lb_rp3, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        rp1 = sol[1:ns] # Array length = number solutes
        rp2 = sol[ns+1:2*ns] # Array length = number solutes
        rp3 = sol[2*ns+1:3*ns] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "λKcentric"
        #λ_e = column[:L]/column[:d]
        #lb_λ = λ_e/10
        #ub_λ = λ_e*10
        sol = optimize_λKcentric(tR_meas.*a, column[:L], column[:gas], prog, options, λ_e, rp1_e, rp2_e, rp3_e, lb_λ, lb_rp1, lb_rp2, lb_rp3, ub_λ, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        λ = sol[1].*ones(ns)
        rp1 = sol[2:ns+1] # Array length = number solutes
        rp2 = sol[ns+1+1:2*ns+1] # Array length = number solutes
        rp3 = sol[2*ns+1+1:3*ns+1] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, λ=λ, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "λKcentric_single"
        #λ_e = column[:L]/column[:d]
        #lb_λ = λ_e/10
        #ub_λ = λ_e*10
        sol = optimize_λKcentric_single(tR_meas.*a, column[:L], column[:gas], prog, options, λ_e, rp1_e, rp2_e, rp3_e, lb_λ, lb_rp1, lb_rp2, lb_rp3, ub_λ, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        λ = Array{Float64}(undef, ns)
        for j=1:ns
            λ[j] = sol[j][1]
            rp1[j] = sol[j][2]
            rp2[j] = sol[j][3]
            rp3[j] = sol[j][4]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, λ=λ, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "φKcentric"
        #φ_e = column[:df]/column[:d]
        #lb_φ = φ_e/10
        #ub_φ = φ_e*10
        sol = optimize_φKcentric(tR_meas.*a, column[:L], column[:d], φ₀, column[:gas], prog, options, φ_e, rp1_e, rp2_e, rp3_e, lb_φ, lb_rp1, lb_rp2, lb_rp3, ub_φ, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        φ = sol[1].*ones(ns)
        rp1 = sol[2:ns+1] # Array length = number solutes
        rp2 = sol[ns+1+1:2*ns+1] # Array length = number solutes
        rp3 = sol[2*ns+1+1:3*ns+1] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, φ=φ, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "φKcentric_single"
        #φ_e = column[:df]/column[:d]
        #lb_φ = φ_e/10
        #ub_φ = φ_e*10
        sol = optimize_φKcentric_single(tR_meas.*a, column[:L], column[:d], φ₀, column[:gas], prog, options, φ_e, rp1_e, rp2_e, rp3_e, lb_φ, lb_rp1, lb_rp2, lb_rp3, ub_φ, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        φ = Array{Float64}(undef, ns)
        for j=1:ns
            φ[j] = sol[j][1]
            rp1[j] = sol[j][2]
            rp2[j] = sol[j][3]
            rp3[j] = sol[j][4]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, φ=φ, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "λφKcentric"
        #λ_e = column[:L]/column[:d]
        #lb_λ = λ_e/10
        #ub_λ = λ_e*10
        #φ_e = column[:df]/column[:d]
        #lb_φ = φ_e/10
        #ub_φ = φ_e*10
        sol = optimize_λφKcentric(tR_meas.*a, column[:L], φ₀, column[:gas], prog, options, λ_e, φ_e, rp1_e, rp2_e, rp3_e, lb_λ, lb_φ, lb_rp1, lb_rp2, lb_rp3, ub_λ, ub_φ, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        λ = sol[1].*ones(ns)
        φ = sol[2].*ones(ns)
        rp1 = sol[3:ns+2] # Array length = number solutes
        rp2 = sol[ns+3:2*ns+2] # Array length = number solutes
        rp3 = sol[2*ns+3:3*ns+2] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, λ=λ, φ=φ, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "λφKcentric_single"
        #λ_e = column[:L]/column[:d]
        #lb_λ = λ_e/10
        #ub_λ = λ_e*10
        #φ_e = column[:df]/column[:d]
        #lb_φ = φ_e/10
        #ub_φ = φ_e*10
        sol = optimize_λφKcentric_single(tR_meas.*a, column[:L], φ₀, column[:gas], prog, options, λ_e, φ_e, rp1_e, rp2_e, rp3_e, lb_λ, lb_φ, lb_rp1, lb_rp2, lb_rp3, ub_λ, ub_φ, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        λ = Array{Float64}(undef, ns)
        φ = Array{Float64}(undef, ns)
        for j=1:ns
            λ[j] = sol[j][1]
            φ[j] = sol[j][2]
            rp1[j] = sol[j][3]
            rp2[j] = sol[j][4]
            rp3[j] = sol[j][5]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, λ=λ, φ=φ, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "λφ"
        #λ_e = column[:L]/column[:d]
        #lb_λ = λ_e/10
        #ub_λ = λ_e*10
        #φ_e = column[:df]/column[:d]
        #lb_φ = φ_e/10
        #ub_φ = φ_e*10
        sol = optimize_λφ(tR_meas.*a, column[:L], φ₀, rp1_e, rp2_e, rp3_e, column[:gas], prog, options, λ_e, φ_e, lb_λ, lb_φ, ub_λ, ub_φ, method; maxiters=maxiters, metric=metric)
        λ = sol[1].*ones(ns)
        φ = sol[2].*ones(ns)
        #rp1 = sol[3:ns+2] # Array length = number solutes
        #rp2 = sol[ns+3:2*ns+2] # Array length = number solutes
        #rp3 = sol[2*ns+3:3*ns+2] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, λ=λ, φ=φ, Tchar=rp1_e, θchar=rp2_e, ΔCp=rp3_e, min=min, retcode=retcode)
    elseif mode == "ddfKcentric"
        #d_e = column[:d]
        #lb_d = d_e/10
        #ub_d = d_e*10
        #df_e = column[:df]
        #lb_df = df_e/10
        #ub_df = df_e*10
        sol = optimize_λφKcentric(tR_meas.*a, column[:L], φ₀, column[:gas], prog, options, d_e, df_e, rp1_e, rp2_e, rp3_e, lb_d, lb_df, lb_rp1, lb_rp2, lb_rp3, ub_d, ub_df, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        d = sol[1].*ones(ns)
        df = sol[2].*ones(ns)
        rp1 = sol[3:ns+2] # Array length = number solutes
        rp2 = sol[ns+3:2*ns+2] # Array length = number solutes
        rp3 = sol[2*ns+3:3*ns+2] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, d=d, df=df, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "ddfKcentric_single"
        #d_e = column[:d]
        #lb_d = d_e/10
        #ub_d = d_e*10
        #df_e = column[:df]
        #lb_df = df_e/10
        #ub_df = df_e*10
        sol = optimize_ddfKcentric_single(tR_meas.*a, column[:L], φ₀, column[:gas], prog, options, d_e, df_e, rp1_e, rp2_e, rp3_e, lb_d, lb_df, lb_rp1, lb_rp2, lb_rp3, ub_d, ub_df, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        d = Array{Float64}(undef, ns)
        df = Array{Float64}(undef, ns)
        for j=1:ns
            d[j] = sol[j][1]
            df[j] = sol[j][2]
            rp1[j] = sol[j][3]
            rp2[j] = sol[j][4]
            rp3[j] = sol[j][5]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, d=d, df=df, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min, retcode=retcode)
    elseif mode == "ABC_single"
        sol = optimize_ABC_single(tR_meas.*a, column[:L], λ_e, β_e, column[:gas], prog, options, rp1_e, rp2_e, rp3_e, lb_rp1, lb_rp2, lb_rp3, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        for j=1:ns
            rp1[j] = sol[j][1]
            rp2[j] = sol[j][2]
            rp3[j] = sol[j][3]
            min[j] = sol[j].minimum
            retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, A=rp1, B=rp2, C=rp3, min=min, retcode=retcode)
    elseif mode == "ABC"
        sol = optimize_ABC(tR_meas.*a, column[:L], λ_e, β_e, column[:gas], prog, options, rp1_e, rp2_e, rp3_e, lb_rp1, lb_rp2, lb_rp3, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        rp1 = sol[1:ns] # Array length = number solutes
        rp2 = sol[ns+1:2*ns] # Array length = number solutes
        rp3 = sol[2*ns+1:3*ns]
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, A=rp1, B=rp2, C=rp3, min=min, retcode=retcode)
    elseif mode == "λABC"
        #λ_e = column[:L]/column[:d]
        #lb_λ = λ_e/100
        #ub_λ = λ_e*100
        sol = optimize_λABC(tR_meas.*a, column[:L], β_e, column[:gas], prog, options, λ_e, rp1_e, rp2_e, rp3_e, lb_λ, lb_rp1, lb_rp2, lb_rp3, ub_λ, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        λ = sol[1].*ones(ns)
        rp1 = sol[2:ns+1] # Array length = number solutes
        rp2 = sol[ns+2:2*ns+1] # Array length = number solutes
        rp3 = sol[2*ns+2:3*ns+1]
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, λ=λ, A=rp1, B=rp2, C=rp3, min=min, retcode=retcode)
    elseif mode == "βABC"
        #df_e = column[:df]
        #lb_df = df_e/10
        #ub_df = df_e*100
        sol = optimize_βABC(tR_meas.*a, column[:L], λ_e, column[:gas], prog, options, β_e, rp1_e, rp2_e, rp3_e, lb_β, lb_rp1, lb_rp2, lb_rp3, ub_β, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        β = sol[1].*ones(ns)
        rp1 = sol[2:ns+1] # Array length = number solutes
        rp2 = sol[ns+2:2*ns+1] # Array length = number solutes
        rp3 = sol[2*ns+2:3*ns+1]
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, β=β, A=rp1, B=rp2, C=rp3, min=min, retcode=retcode)
    elseif mode == "λβABC"
        #λ_e = column[:L]/column[:d]
        #lb_λ = λ_e/100
        #ub_λ = λ_e*100
        #df_e = column[:df]
        #lb_df = df_e/10
        #ub_df = df_e*100
        sol = optimize_λβABC(tR_meas.*a, column[:L], column[:gas], prog, options, λ_e, β_e, rp1_e, rp2_e, rp3_e, lb_λ, lb_β, lb_rp1, lb_rp2, lb_rp3, ub_λ, ub_β, ub_rp1, ub_rp2, ub_rp3, method; maxiters=maxiters, metric=metric)
        λ = sol[1].*ones(ns)
        β = sol[2].*ones(ns)
        rp1 = sol[3:ns+2] # Array length = number solutes
        rp2 = sol[ns+3:2*ns+2] # Array length = number solutes
        rp3 = sol[2*ns+3:3*ns+2]
        for j=1:ns
            min[j] = sol.minimum
            retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, λ=λ, β=β, A=rp1, B=rp2, C=rp3, min=min, retcode=retcode)
    end
	
	return df, sol
end
=#