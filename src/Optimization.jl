# functions used for the optimization of the loss-function 




#------------------
# Optimization only for the K-centric retention parameters
"""
    opt_Kcentric(x, p)

Function used for optimization of the loss-function in regards to the three K-centric parameters.

# Arguments
* `x` ... 3n-vector of the three K-centric parameters of n solutes. Elements 1:n are Tchar, n+1:2n are θchar and 2n+1:3n are ΔCp values.
* `p` ... vector containing the fixed parameters:
    * `tR = p[1]` ... vector of the measured retention times in seconds.
    * `substance_list = p[2]` ... vector of the names of the solutes related to `tR` and `prog`.
    * `L = p[3]` ... number of the length of the column in m.
    * `d = p[4]` ... number of the diameters of the column in m.
    * `prog = p[5]` ... vector of structure GasChromatographySimulator.Programs containing the definition of the GC-programs.
    * `opt = p[6]` ... struture GasChromatographySimulator.Options containing the settings for options for the simulation.
    * `gas = p[7]` ... string of name of the mobile phase gas. 
    * `metric = p[8]` ... string of the metric used for the loss function (`squared` or `abs`). 

# Output
* `sum((tR.-tRcalc).^2)` ... sum of the squared residuals over m GC-programs and n solutes.
"""
function opt_Kcentric(x, p)
    # for 2-parameter version `x` is shorter, the part for `ΔCp = x[2*ns+1+1:3*ns+1]` will be missing
    # instead set `ΔCp = zeros(ns)`
	tR = p[1]
	substance_list = p[2]
	L = p[3]
	d = p[4]
	prog = p[5]
	opt = p[6]
    gas = p[7]
    metric = p[8]

	ns = length(unique(substance_list))
		
	Tchar = x[1:ns] # Array length = number solutes
	θchar = x[ns+1:2*ns] # Array length = number solutes
	ΔCp = x[2*ns+1:3*ns] # Array length = number solutes
    return loss(tR, Tchar, θchar, ΔCp, substance_list, L, d, prog, gas; opt=opt, metric=metric)[1]
end

#=function opt_Kcentric_(x, p)
    # for 2-parameter version `x` is shorter, the part for `ΔCp = x[2*ns+1+1:3*ns+1]` will be missing
    # instead set `ΔCp = zeros(ns)`
	tR = p[1]
    φ₀ = p[2]
	L = p[3]
	d = p[4]
    df = p[5]
	prog = p[6]
	opt = p[7]
    gas = p[8]
    metric = p[9]
    if length(size(tR)) == 1
		ns = 1
	else
		ns = size(tR)[2]
	end
	Tchar = x[1:ns] # Array length = number solutes
	θchar = x[ns+1:2*ns] # Array length = number solutes
	ΔCp = x[2*ns+1:3*ns] # Array length = number solutes
    return loss(tR, Tchar, θchar, ΔCp, φ₀, L, d, df, prog, gas; opt=opt, metric=metric)[1]
end=#

#="""
    optimize_Kcentric(tR, col, prog, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, metric="squared") 

Optimization regarding the estimization of the retention parameters `Tchar`, `θchar` and `ΔCp`. The initial guess is a vector (`Tchar_e`, `θchar_e` and `ΔCp_e`) for optimization
algorithms, which do not need lower/upper bounds. If a method is used, which needs lower/upper bounds the initial parameters should be a matrix, with the first column (`[1,:]`)
beeing the initial guess, the second column (`[2,:]`) the lower bound and the third column (`[3,:]`) the upper bound.
"""
function optimize_Kcentric(tR, col, prog, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number
    p = (tR, col.L, col.d, prog, opt, col.gas, metric)
	x0 = [Tchar_e; θchar_e; ΔCp_e]
	optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end=#

"""
    optimize_Kcentric(tR::Vector{T}, substance_list, col, prog, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=RetentionParameterEstimator.NewtonTrustRegion(), opt=RetentionParameterEstimator.std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number

Optimize the K-centric retention parameters for a given set of substances, measured retention times and used programs.

# Arguments
- `tR::Vector{T}`: Vector of retention times.
- `substance_list`: List of substances to be optimized.
- `col`: Column characteristics including length `L`, diameter `d`, and gas type `gas`.
- `prog`: Program conditions.
- `Tchar_e::Vector{T}`: Initial estimates for characteristic temperatures.
- `θchar_e::Vector{T}`: Initial estimates for characteristic phase ratios.
- `ΔCp_e::Vector{T}`: Initial estimates for heat capacity changes.
- `method`: Optimization method to be used (default: `RetentionParameterEstimator.NewtonTrustRegion()`).
- `opt`: Optimization options (default: `RetentionParameterEstimator.std_opt`).
- `maxiters`: Maximum number of iterations (default: 10000).
- `maxtime`: Maximum time allowed for optimization in seconds (default: 600.0).
- `metric`: Metric to be used for optimization (default: "squared").

# Returns
- `opt_sol`: The solution of the optimization problem.
"""
function optimize_Kcentric(tR::Vector{T}, substance_list, col, prog, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=RetentionParameterEstimator.NewtonTrustRegion(), opt=RetentionParameterEstimator.std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number
    p = (tR, substance_list, col.L, col.d, prog, opt, col.gas, metric)
	x0 = [Tchar_e; θchar_e; ΔCp_e]
	optf = OptimizationFunction(opt_Kcentric, RetentionParameterEstimator.Optimization.AutoForwardDiff())
	prob = RetentionParameterEstimator.Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end

#=function optimize_Kcentric(tR, col, prog, Tchar_e::Matrix{T}, θchar_e::Matrix{T}, ΔCp_e::Matrix{T}; method=BBO_adaptive_de_rand_1_bin_radiuslimited(), opt=std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number # here default method should be one which needs bounds
    p = (tR, col.L, col.d, prog, opt, col.gas, metric)
	x0 = [Tchar_e[1,:]; θchar_e[1,:]; ΔCp_e[1,:]]
    lb = [Tchar_e[2,:]; θchar_e[2,:]; ΔCp_e[2,:]]
    ub = [Tchar_e[3,:]; θchar_e[3,:]; ΔCp_e[3,:]]
	optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end

function optimize_Kcentric_(tR, col, prog, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number
    p = (tR, col.df/col.d, col.L, col.d, col.df, prog, opt, col.gas, metric)
	x0 = [Tchar_e; θchar_e; ΔCp_e]
	optf = OptimizationFunction(opt_Kcentric_, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end=#

#="""
    opt_dKcentric(x, p)

Function used for optimization of the loss-function in regards to the three K-centric parameters and the column diameter.

# Arguments
* `x` ... (3n+1)-vector of the column diameter and the three K-centric parameters of n solutes. The first element represents the column diameterd `d` and elements 2:n+1 are `Tchar`, n+2:2n+1 are `θchar` and 2n+2:3n+1 are `ΔCp` values.
* `p` ... vector containing the fixed parameters:
    * `tR = p[1]` ... mxn-array of the measured retention times in seconds.
    * `L = p[2]` ... number of the length of the column in m.
    * `prog = p[3]` ... m-array of structure GasChromatographySimulator.Programs containing the definition of the GC-programs.
    * `opt = p[4]` ... struture GasChromatographySimulator.Options containing the settings for options for the simulation.
    * `gas = p[5]` ... string of name of the mobile phase gas. 
    * `metric = p[6]` ... string of the metric used for the loss function (`squared` or `abs`). 

# Output
* `sum((tR.-tRcalc).^2)` ... sum of the squared residuals over m GC-programs and n solutes.
"""
function opt_dKcentric(x, p)
    # for 2-parameter version `x` is shorter, the part for `ΔCp = x[2*ns+1+1:3*ns+1]` will be missing
    # instead set `ΔCp = zeros(ns)`
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

    return loss(tR, Tchar, θchar, ΔCp, L, d, prog, gas; opt=opt, metric=metric)
end=#

"""
    opt_dKcentric(x, p)

Function used for optimization of the loss-function in regards to the colum diameter and the three K-centric parameters.

# Arguments
* `x` ... 1+3n-vector of column diameter and the three K-centric parameters of n solutes. Element 1 is the diameter, 2:n+1 are Tchar, n+2:2n+1 are θchar and 2n+2:3n+1 are ΔCp values.
* `p` ... vector containing the fixed parameters:
    * `tR = p[1]` ... vector of the measured retention times in seconds.
    * `substance_list = p[2]` ... vector of the names of the solutes related to `tR` and `prog`.
    * `L = p[3]` ... number of the length of the column in m.
    * `prog = p[4]` ... vector of structure GasChromatographySimulator.Programs containing the definition of the GC-programs.
    * `opt = p[5]` ... struture GasChromatographySimulator.Options containing the settings for options for the simulation.
    * `gas = p[6]` ... string of name of the mobile phase gas. 
    * `metric = p[7]` ... string of the metric used for the loss function (`squared` or `abs`). 

# Output
* `sum((tR.-tRcalc).^2)` ... sum of the squared residuals over m GC-programs and n solutes.
"""
function opt_dKcentric(x, p)
    # for 2-parameter version `x` is shorter, the part for `ΔCp = x[2*ns+1+1:3*ns+1]` will be missing
    # instead set `ΔCp = zeros(ns)`
    tR = p[1]
	substance_list = p[2]
    L = p[3]
    prog = p[4]
    opt = p[5]
    gas = p[6]
    metric = p[7]
	
    ns = length(unique(substance_list))
	
    d = x[1]
    Tchar = x[2:ns+1] # Array length = number solutes
    θchar = x[ns+1+1:2*ns+1] # Array length = number solutes
    ΔCp = x[2*ns+1+1:3*ns+1] # Array length = number solutes

    return loss(tR, Tchar, θchar, ΔCp, substance_list, L, d, prog, gas; opt=opt, metric=metric)
end

#=function opt_dKcentric_(x, p)
    # for 2-parameter version `x` is shorter, the part for `ΔCp = x[2*ns+1+1:3*ns+1]` will be missing
    # instead set `ΔCp = zeros(ns)`
    tR = p[1]
    φ₀ = p[2]
    L = p[3]
    df = p[4]
    prog = p[5]
    opt = p[6]
    gas = p[7]
    metric = p[8]
    if length(size(tR)) == 1
        ns = 1
    else
        ns = size(tR)[2]
    end
    d = x[1]
    Tchar = x[2:ns+1] # Array length = number solutes
    θchar = x[ns+1+1:2*ns+1] # Array length = number solutes
    ΔCp = x[2*ns+1+1:3*ns+1] # Array length = number solutes

    return loss(tR, Tchar, θchar, ΔCp, φ₀, L, d, df, prog, gas; opt=opt, metric=metric)
end=#

#="""
    optimize_dKcentric(tR, col, prog, d_e, Tchar_e, θchar_e, ΔCp_e; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, metric="squared") 

Optimization regarding the estimization of the column diameter `d` and the retention parameters `Tchar`, `θchar` and `ΔCp`. The initial guess is a number (`d_e`) or a vector (`Tchar_e`, `θchar_e` and `ΔCp_e`) for optimization
algorithms, which do not need lower/upper bounds. If a method is used, which needs lower/upper bounds the initial parameters should be a vector of length 3 (`d_e`) or a matrix (`Tchar_e`, `θchar_e` and `ΔCp_e`), with the first element/column 
beeing the initial guess, the second element/column the lower bound and the third element/column the upper bound.
"""
function optimize_dKcentric(tR, col, prog, d_e::Number, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number
    p = (tR, col.L, prog, opt, col.gas, metric)
	x0 = [d_e; Tchar_e; θchar_e; ΔCp_e]
	optf = OptimizationFunction(opt_dKcentric, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end=#

"""
    optimize_dKcentric(tR::Vector{T}, substance_list, col, prog, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=RetentionParameterEstimator.NewtonTrustRegion(), opt=RetentionParameterEstimator.std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number

Optimize the column diameter and K-centric retention parameters for a given set of substances, measured retention times and used programs.

# Arguments
- `tR::Vector{T}`: Vector of retention times.
- `substance_list`: List of substances to be optimized.
- `col`: Column characteristics including length `L` and gas type `gas`.
- `prog`: Program conditions.
- `Tchar_e::Vector{T}`: Initial estimates for characteristic temperatures.
- `θchar_e::Vector{T}`: Initial estimates for characteristic phase ratios.
- `ΔCp_e::Vector{T}`: Initial estimates for heat capacity changes.
- `method`: Optimization method to be used (default: `RetentionParameterEstimator.NewtonTrustRegion()`).
- `opt`: Optimization options (default: `RetentionParameterEstimator.std_opt`).
- `maxiters`: Maximum number of iterations (default: 10000).
- `maxtime`: Maximum time allowed for optimization in seconds (default: 600.0).
- `metric`: Metric to be used for optimization (default: "squared").

# Returns
- `opt_sol`: The solution of the optimization problem.
"""
function optimize_dKcentric(tR::Vector{T}, substance_list, col, prog, d_e::Number, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number
    p = (tR, substance_list, col.L, prog, opt, col.gas, metric)
	x0 = [d_e; Tchar_e; θchar_e; ΔCp_e]
	optf = OptimizationFunction(opt_dKcentric, RetentionParameterEstimator.Optimization.AutoForwardDiff())
	prob = RetentionParameterEstimator.Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end

"""
    opt_d_only(x, p)

Function used for optimization of the loss-function with respect to only the column diameter `d`, while keeping the K-centric parameters fixed.

# Arguments
* `x` ... 1-element vector containing the column diameter `d`.
* `p` ... vector containing the fixed parameters:
    * `tR = p[1]` ... vector of the measured retention times in seconds.
    * `substance_list = p[2]` ... vector of the names of the solutes related to `tR` and `prog`.
    * `Tchar = p[3]` ... vector of characteristic temperatures (fixed).
    * `θchar = p[4]` ... vector of characteristic constants (fixed).
    * `ΔCp = p[5]` ... vector of heat capacity changes (fixed).
    * `L = p[6]` ... number of the length of the column in m.
    * `prog = p[7]` ... vector of structure GasChromatographySimulator.Programs containing the definition of the GC-programs.
    * `opt = p[8]` ... structure GasChromatographySimulator.Options containing the settings for options for the simulation.
    * `gas = p[9]` ... string of name of the mobile phase gas.
    * `metric = p[10]` ... string of the metric used for the loss function (`squared` or `abs`).

# Output
* Loss value from the loss function.
"""
function opt_d_only(x, p)
    d = x[1]
    tR = p[1]
    substance_list = p[2]
    Tchar = p[3]  # fixed
    θchar = p[4]  # fixed
    ΔCp = p[5]    # fixed
    L = p[6]
    prog = p[7]
    opt = p[8]
    gas = p[9]
    metric = p[10]
    
    return loss(tR, Tchar, θchar, ΔCp, substance_list, L, d, prog, gas; opt=opt, metric=metric)
end

"""
    optimize_d_only(tR::Vector{T}, substance_list, Tchar, θchar, ΔCp, L, prog, gas, d_e; method=RetentionParameterEstimator.NewtonTrustRegion(), opt=RetentionParameterEstimator.std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number

Optimize only the column diameter `d` while keeping the K-centric retention parameters fixed.

# Arguments
- `tR::Vector{T}`: Vector of retention times.
- `substance_list`: List of substances.
- `Tchar`: Vector of characteristic temperatures (fixed).
- `θchar`: Vector of characteristic phase ratios (fixed).
- `ΔCp`: Vector of heat capacity changes (fixed).
- `L`: Length of the column in m.
- `prog`: Program conditions.
- `gas`: Gas type.
- `d_e`: Initial estimate for column diameter.
- `method`: Optimization method to be used (default: `RetentionParameterEstimator.NewtonTrustRegion()`).
- `opt`: Optimization options (default: `RetentionParameterEstimator.std_opt`).
- `maxiters`: Maximum number of iterations (default: 10000).
- `maxtime`: Maximum time allowed for optimization in seconds (default: 600.0).
- `metric`: Metric to be used for optimization (default: "squared").

# Returns
- `d_opt`: The optimized column diameter value.
"""
function optimize_d_only(tR::Vector{T}, substance_list, Tchar, θchar, ΔCp, L, prog, gas, d_e; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number
    p = (tR, substance_list, Tchar, θchar, ΔCp, L, prog, opt, gas, metric)
    x0 = [d_e]
    optf = OptimizationFunction(opt_d_only, RetentionParameterEstimator.Optimization.AutoForwardDiff())
    prob = RetentionParameterEstimator.Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
    return opt_sol[1]  # return just the d value
end


#=function optimize_dKcentric(tR, col, prog, d_e::Vector{T}, Tchar_e::Matrix{T}, θchar_e::Matrix{T}, ΔCp_e::Matrix{T}; method=BBO_adaptive_de_rand_1_bin_radiuslimited(), opt=opt_std, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number  # here default method should be one which needs bounds  
    p = (tR, col.L, prog, opt, col.gas, metric)
	x0 = [d_e[1]; Tchar_e[1,:]; θchar_e[1,:]; ΔCp_e[1,:]]
    lb = [d_e[2]; Tchar_e[2,:]; θchar_e[2,:]; ΔCp_e[2,:]]
    ub = [d_e[3]; Tchar_e[3,:]; θchar_e[3,:]; ΔCp_e[3,:]]
	optf = OptimizationFunction(opt_dKcentric, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end

function optimize_dKcentric_(tR, col, prog, d_e::Number, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number
    p = (tR, col.df/col.d, col.L, col.df, prog, opt, col.gas, metric)
	x0 = [d_e; Tchar_e; θchar_e; ΔCp_e]
	optf = OptimizationFunction(opt_dKcentric_, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end=#

#=function estimate_parameters(tRs, solute_names, col, prog, rp1_e, rp2_e, rp3_e; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, mode="dKcentric", metric="squared", pout="vacuum", time_unit="min")
    # mode = "Kcentric", "Kcentric_single", "dKcentric", "dKcentric_single"
    
    # add the case for 2-parameter model, where rp3 === 0.0 always
    # -> alternative versions of the different `optimize_` functions (without the third retention parameter)
    # -> similar, alternative versions for `opt_` functions needed
    a = time_unit_conversion_factor(time_unit)
	tR_meas = Array(tRs[:,2:end]).*a
    if length(size(tR_meas)) == 1
        ns = 1
    else
        ns = size(tR_meas)[2]
    end
    
    d_e = col.d

    rp1 = Array{Float64}(undef, ns)
	rp2 = Array{Float64}(undef, ns)
	rp3 = Array{Float64}(undef, ns)
    min = Array{Float64}(undef, ns)
    #retcode = Array{Any}(undef, ns)
    if mode == "Kcentric_single"
        sol = Array{SciMLBase.OptimizationSolution}(undef, ns)
        for j=1:ns
            sol[j] = optimize_Kcentric(tR_meas[:,j], col, prog, rp1_e[j,:], rp2_e[j,:], rp3_e[j,:]; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
            rp1[j] = sol[j][1]
            rp2[j] = sol[j][2]
            rp3[j] = sol[j][3]
            min[j] = sol[j].minimum
            #retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
    elseif mode == "Kcentric"
        sol = optimize_Kcentric(tR_meas, col, prog, rp1_e, rp2_e, rp3_e; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
        rp1 = sol[1:ns] # Array length = number solutes
        rp2 = sol[ns+1:2*ns] # Array length = number solutes
        rp3 = sol[2*ns+1:3*ns] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            #retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
    elseif mode == "dKcentric"
        sol = optimize_dKcentric(tR_meas, col, prog, d_e, rp1_e, rp2_e, rp3_e; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
        d = sol[1].*ones(ns)
        rp1 = sol[2:ns+1] # Array length = number solutes
        rp2 = sol[ns+1+1:2*ns+1] # Array length = number solutes
        rp3 = sol[2*ns+1+1:3*ns+1] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
            #retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, d=d, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
    elseif mode == "dKcentric_single"
        sol = Array{SciMLBase.OptimizationSolution}(undef, ns)
        d = Array{Float64}(undef, ns)
        for j=1:ns
            sol[j] = optimize_dKcentric(tR_meas[:,j], col, prog, d_e, rp1_e[j,:], rp2_e[j,:], rp3_e[j,:]; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
            d[j] = sol[j][1]
            rp1[j] = sol[j][2]
            rp2[j] = sol[j][3]
            rp3[j] = sol[j][4]
            min[j] = sol[j].minimum
            #retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, d=d, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
    end
	
	return df, sol
end=#

function estimate_parameters(tRs, solute_names, col, prog, rp1_e, rp2_e, rp3_e; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, mode="dKcentric", metric="squared", pout="vacuum", time_unit="min", parallel=false)
    # mode = "Kcentric", "Kcentric_single", "dKcentric", "dKcentric_single", "d_only"
    
    # add the case for 2-parameter model, where rp3 === 0.0 always
    # -> alternative versions of the different `optimize_` functions (without the third retention parameter)
    # -> similar, alternative versions for `opt_` functions needed
	a = time_unit_conversion_factor(time_unit)
	tR_meas = Array(tRs[:,2:end]).*a
	
    if length(size(tR_meas)) == 1
        ns = 1
    else
        ns = size(tR_meas)[2]
    end
    
    d_e = col.d

    rp1 = Array{Float64}(undef, ns)
	rp2 = Array{Float64}(undef, ns)
	rp3 = Array{Float64}(undef, ns)
    min = Array{Float64}(undef, ns)
    #retcode = Array{Any}(undef, ns)
    if mode == "Kcentric_single" #-> new version for missing values ok -> check with no missing		
        sol = Array{SciMLBase.OptimizationSolution}(undef, ns)
        if parallel
            Base.Threads.@threads for j=1:ns
                # filter-out missing values:
                tRs_, prog_, subst_list_ = prepare_single_substance_data(tR_meas[:,j], prog, solute_names[j])
                    
                sol[j] = optimize_Kcentric(tRs_, subst_list_, col, prog_, rp1_e[j,:], rp2_e[j,:], rp3_e[j,:]; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
                
                rp1[j] = sol[j][1]
                rp2[j] = sol[j][2]
                rp3[j] = sol[j][3]
                min[j] = sol[j].objective
            end
        else
        for j=1:ns
			# filter-out missing values:
                tRs_, prog_, subst_list_ = prepare_single_substance_data(tR_meas[:,j], prog, solute_names[j])
				
            sol[j] = optimize_Kcentric(tRs_, subst_list_, col, prog_, rp1_e[j,:], rp2_e[j,:], rp3_e[j,:]; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
			
            rp1[j] = sol[j][1]
            rp2[j] = sol[j][2]
            rp3[j] = sol[j][3]
            min[j] = sol[j].objective
        end
        end
        
        df = DataFrame(Name=solute_names, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
        
        # IMPORTANT: When parallelization is used, ensure DataFrame is sorted by Name to match solute_names order
        # This prevents any potential ordering issues from parallel execution
        if parallel
            # Reorder DataFrame to match solute_names order explicitly
            # Create a dictionary for O(1) lookup instead of O(n) findfirst
            name_to_index = Dict(name => idx for (idx, name) in enumerate(df.Name))
            sorted_indices = [name_to_index[name] for name in solute_names]
            df = df[sorted_indices, :]
        end
    elseif mode == "Kcentric" #-> new version for missing values ok -> check with no missing
		# vectorizing and filter-out missing values:
		tRs_, prog_, subst_list_ = prepare_optimization_data(tRs, solute_names, prog, time_unit)
		
        sol = optimize_Kcentric(tRs_, subst_list_, col, prog_, rp1_e, rp2_e, rp3_e; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
        rp1 = sol[1:ns] # Array length = number solutes
        rp2 = sol[ns+1:2*ns] # Array length = number solutes
        rp3 = sol[2*ns+1:3*ns] # Array length = number solutes
        for j=1:ns
            min[j] = sol.objective
            #retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
    elseif mode == "dKcentric"
		# vectorizing and filter-out missing values:
		tRs_, prog_, subst_list_ = prepare_optimization_data(tRs, solute_names, prog, time_unit)
		
        sol = optimize_dKcentric(tRs_, subst_list_, col, prog_, d_e, rp1_e, rp2_e, rp3_e; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
        d = sol[1].*ones(ns)
        rp1 = sol[2:ns+1] # Array length = number solutes
        rp2 = sol[ns+1+1:2*ns+1] # Array length = number solutes
        rp3 = sol[2*ns+1+1:3*ns+1] # Array length = number solutes
        for j=1:ns
            min[j] = sol.objective
            #retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, d=d, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
    elseif mode == "dKcentric_single"#-> new version for missing values ok -> check with no missing
        sol = Array{SciMLBase.OptimizationSolution}(undef, ns)
        d = Array{Float64}(undef, ns)
        if parallel
            Base.Threads.@threads for j=1:ns
                # filter-out missing values:
                tRs_, prog_, subst_list_ = prepare_single_substance_data(tR_meas[:,j], prog, solute_names[j])
                
                sol[j] = optimize_dKcentric(tRs_, subst_list_, col, prog_, d_e, rp1_e[j,:], rp2_e[j,:], rp3_e[j,:]; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
                d[j] = sol[j][1]
                rp1[j] = sol[j][2]
                rp2[j] = sol[j][3]
                rp3[j] = sol[j][4]
                min[j] = sol[j].objective
            end
        else
        for j=1:ns
			# filter-out missing values:
                tRs_, prog_, subst_list_ = prepare_single_substance_data(tR_meas[:,j], prog, solute_names[j])
			
            sol[j] = optimize_dKcentric(tRs_, subst_list_, col, prog_, d_e, rp1_e[j,:], rp2_e[j,:], rp3_e[j,:]; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
            d[j] = sol[j][1]
            rp1[j] = sol[j][2]
            rp2[j] = sol[j][3]
            rp3[j] = sol[j][4]
            min[j] = sol[j].objective
        end
        end
        
        df = DataFrame(Name=solute_names, d=d, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
        
        # IMPORTANT: When parallelization is used, ensure DataFrame is sorted by Name to match solute_names order
        # This prevents any potential ordering issues from parallel execution
        if parallel
            # Reorder DataFrame to match solute_names order explicitly
            # Create a dictionary for O(1) lookup instead of O(n) findfirst
            name_to_index = Dict(name => idx for (idx, name) in enumerate(df.Name))
            sorted_indices = [name_to_index[name] for name in solute_names]
            df = df[sorted_indices, :]
        end
    elseif mode == "d_only"
        # Optimize only d while keeping substance parameters fixed
        # rp1_e, rp2_e, rp3_e should be vectors (fixed values) or matrices (use first row as fixed values)
        # Extract fixed values (handle both vector and matrix inputs)
        if isa(rp1_e, Matrix)
            Tchar_fixed_input = rp1_e[1,:]  # Use first row if matrix
        else
            Tchar_fixed_input = copy(rp1_e)  # Already a vector - make a copy
        end
        if isa(rp2_e, Matrix)
            θchar_fixed_input = rp2_e[1,:]
        else
            θchar_fixed_input = copy(rp2_e)
        end
        if isa(rp3_e, Matrix)
            ΔCp_fixed_input = rp3_e[1,:]
        else
            ΔCp_fixed_input = copy(rp3_e)
        end
        
        # Ensure vectors have correct length
        if length(Tchar_fixed_input) != ns
            error("Tchar_fixed length ($(length(Tchar_fixed_input))) does not match number of substances ($ns)")
        end
        
        # Prepare vectorized data (filter out missing values)
        tRs_, prog_, subst_list_ = prepare_optimization_data(tRs, solute_names, prog, time_unit)
        
        # IMPORTANT: The loss function uses unique(substance_list) to map parameters
        # So we need to reorder Tchar_fixed, θchar_fixed, ΔCp_fixed to match unique(subst_list_) order
        unique_subst = unique(subst_list_)
        Tchar_fixed = Array{Float64}(undef, length(unique_subst))
        θchar_fixed = Array{Float64}(undef, length(unique_subst))
        ΔCp_fixed = Array{Float64}(undef, length(unique_subst))
        
        for (idx, subst_name) in enumerate(unique_subst)
            orig_idx = findfirst(solute_names .== subst_name)
            if isnothing(orig_idx)
                error("Substance $subst_name not found in solute_names")
            end
            Tchar_fixed[idx] = Tchar_fixed_input[orig_idx]
            θchar_fixed[idx] = θchar_fixed_input[orig_idx]
            ΔCp_fixed[idx] = ΔCp_fixed_input[orig_idx]
        end
        
        # Optimize d (1D optimization, shared across all substances)
        d_opt = optimize_d_only(tRs_, subst_list_, Tchar_fixed, θchar_fixed, ΔCp_fixed, 
                                col.L, prog_, col.gas, d_e; 
                                method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
        
        # Calculate loss value for the optimized d
        loss_val = loss(tRs_, Tchar_fixed, θchar_fixed, ΔCp_fixed, subst_list_, col.L, d_opt, prog_, col.gas; opt=opt, metric=metric)[1]
        
        # Return DataFrame with same d for all substances, fixed substance parameters
        # Use original solute_names order (not unique_subst order)
        d = d_opt.*ones(ns)
        rp1 = Tchar_fixed_input  # Use original order for DataFrame
        rp2 = θchar_fixed_input
        rp3 = ΔCp_fixed_input
        min = loss_val.*ones(ns)  # Same loss for all since it's joint optimization
        
        df = DataFrame(Name=solute_names, d=d, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)
        sol = nothing  # d_only doesn't return individual solutions
    end
	
	return df, sol
end

"""
    estimate_parameters(chrom; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, mode="dKcentric", metric="squared")

Calculate the estimates for the K-centric parameters and (optional) the column diameter.

# Arguments
* `chrom` ... Tuple of the loaded chromatogram, see [`load_chromatograms`](@ref)

# Options
* `method=NewtonTrustRegion()` ... used optimization method
* `opt=std_opt` ... general options, `std_opt = GasChromatographySimulator.Options(abstol=1e-8, reltol=1e-5, ng=true, odesys=false)`
* `maxiters=10000` ... maximum number of iterations for every single optimization
* `maxtime=600.0` ... maximum time for every single optimization
* `mode="dKcentric"` ... mode of the estimation. 
    Possible options: 
    * "Kcentric_single" ... optimization for the three K-centric retention parameters separatly for every solute
    * "Kcentric" ... optimization for the three K-centric retention parameters together for all solutes
    * "dKcentric_single" ... optimization for the column diameter and the three K-centric retention parameters separatly for every solute
    * "dKcentric" ... optimization for the column diameter and the three K-centric retention parameters together for all solutes
    * "d_only" ... optimization for only the column diameter while keeping the K-centric retention parameters fixed (requires fixed parameter values as vectors)
    * "d_only" ... optimization for only the column diameter while keeping the K-centric retention parameters fixed (requires fixed parameter values as vectors)
* `metric="squared"` ... used metric for the loss function ("squared" or "abs")

# Output
* `df` ... DataFrame with the columns `Name` (solute names), `d` (estimated column diameter, optional), `Tchar` (estimated Tchar), `θchar` (estimated θchar), `ΔCp` (estimated ΔCp) and `min` (value of the loss function at the found optima)
* `sol` ... Array of `SciMLBase.OptimizationSolution` with the results of the optimization with some additional informations.
"""    
function estimate_parameters(chrom; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, mode="dKcentric", metric="squared", parallel=false)
	Tchar_est, θchar_est, ΔCp_est, Telu_max = estimate_start_parameter(chrom[3], chrom[1], chrom[2]; time_unit=chrom[6])
	return estimate_parameters(chrom[3], chrom[4], chrom[1], chrom[2], Tchar_est, θchar_est, ΔCp_est; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, mode=mode, metric=metric, pout=chrom[5], time_unit=chrom[6], parallel=parallel)
end

function estimate_parameters_(tRs, solute_names, col, prog, rp1_e, rp2_e, rp3_e; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, mode="dKcentric", metric="squared", pout="vacuum", time_unit="min")
    # mode = "Kcentric", "Kcentric_single", "dKcentric", "dKcentric_single"
    a = time_unit_conversion_factor(time_unit)
	tR_meas = Array(tRs[:,2:end]).*a
    if length(size(tR_meas)) == 1
        ns = 1
    else
        ns = size(tR_meas)[2]
    end
    
    d_e = col.d

    rp1 = Array{Float64}(undef, ns)
	rp2 = Array{Float64}(undef, ns)
	rp3 = Array{Float64}(undef, ns)
    min = Array{Float64}(undef, ns)
    #retcode = Array{Any}(undef, ns)
    if mode == "Kcentric_single"
        sol = Array{SciMLBase.OptimizationSolution}(undef, ns)
        for j=1:ns
            sol[j] = optimize_Kcentric_(tR_meas[:,j], col, prog, rp1_e[j,:], rp2_e[j,:], rp3_e[j,:]; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
            rp1[j] = sol[j][1]
            rp2[j] = sol[j][2]
            rp3[j] = sol[j][3]
            min[j] = sol[j].objective
            #retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
    elseif mode == "Kcentric"
        sol = optimize_Kcentric_(tR_meas, col, prog, rp1_e, rp2_e, rp3_e; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
        rp1 = sol[1:ns] # Array length = number solutes
        rp2 = sol[ns+1:2*ns] # Array length = number solutes
        rp3 = sol[2*ns+1:3*ns] # Array length = number solutes
        for j=1:ns
            min[j] = sol.objective
            #retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
    elseif mode == "dKcentric"
        sol = optimize_dKcentric_(tR_meas, col, prog, d_e, rp1_e, rp2_e, rp3_e; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
        d = sol[1].*ones(ns)
        rp1 = sol[2:ns+1] # Array length = number solutes
        rp2 = sol[ns+1+1:2*ns+1] # Array length = number solutes
        rp3 = sol[2*ns+1+1:3*ns+1] # Array length = number solutes
        for j=1:ns
            min[j] = sol.objective
            #retcode[j] = sol.retcode
        end
        df = DataFrame(Name=solute_names, d=d, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
    elseif mode == "dKcentric_single"
        sol = Array{SciMLBase.OptimizationSolution}(undef, ns)
        d = Array{Float64}(undef, ns)
        for j=1:ns
            sol[j] = optimize_dKcentric_(tR_meas[:,j], col, prog, d_e, rp1_e[j,:], rp2_e[j,:], rp3_e[j,:]; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
            d[j] = sol[j][1]
            rp1[j] = sol[j][2]
            rp2[j] = sol[j][3]
            rp3[j] = sol[j][4]
            min[j] = sol[j].objective
            #retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, d=d, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
    end
	
	return df, sol
end

# full methods

"""
    check_measurement(meas, col_input; min_th=0.1, loss_th=1.0, se_col=true)

Similar to `method_m1` ([`method_m1`](@ref)) estimate the three retention parameters ``T_{char}``, ``θ_{char}`` and ``ΔC_p`` including standard errors, see [`RetentionParameterEstimator.stderror`](@ref).
In addition, if the found optimized minima is above a threshold `min_th`, it is flagged and the squared differences of single measured retention times and calculated retention times above
another threshold `loss_th` are recorded.   

# Arguments
* `meas` ... Tuple with the loaded measurement data, see [`load_chromatograms`](@ref).
* `col_input` ... Named tuple with `col_input.L` the column length in m and `col_input.d` the column diameter in mm. If this parameter is not gicen, than these parameters are taken from `meas`. 
* `se_col=true` ... If `true` the standard errors (from the Hessian matrix, see [`RetentionParameterEstimator.stderror`](@ref)) of the estimated parameters are added as separate columns to the result dataframe. If `false` the standard errors are added to the values as `Masurement` type.  

# Output 
* `check` ... Boolean. `true` if all values are below the thresholds, `false` if not.
* `msg` ... String. Description of `check`
* `df_flag` ... Dataframe containing name of the flagged measurements, solutes and the corresponding mesured and calculated retention times.
* `index_flag` ... Indices of flagged results
* `res` ... Dataframe with the optimized parameters and the found minima.
* `Telu_max` ... The maximum of elution temperatures every solute experiences in the measured programs.
"""  
function check_measurement(meas, col_input; min_th=0.1, loss_th=1.0, se_col=true, method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, parallel=false)
	col = GasChromatographySimulator.Column(col_input.L, col_input.d*1e-3, meas[1].df, meas[1].sp, meas[1].gas)
	Tchar_est, θchar_est, ΔCp_est, Telu_max = RetentionParameterEstimator.estimate_start_parameter(meas[3], col, meas[2]; time_unit=meas[6])
	df = estimate_parameters(meas[3], meas[4], col, meas[2], Tchar_est, θchar_est, ΔCp_est; mode="Kcentric_single", pout=meas[5], time_unit=meas[6], method=method, opt=std_opt, maxiters=maxiters, maxtime=maxtime, parallel=parallel)[1]
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
    # calculate the standard errors of the 3 parameters using the hessian matrix
    stderrors = stderror(meas, df, col_input; opt=opt, parallel=parallel)[1]
	# output dataframe
    res = if se_col == true
	    DataFrame(Name=df.Name, min=df.min, Tchar=df.Tchar, Tchar_std=stderrors.sd_Tchar, θchar=df.θchar, θchar_std=stderrors.sd_θchar, ΔCp=df.ΔCp, ΔCp_std=stderrors.sd_ΔCp)
    else
	    DataFrame(Name=df.Name, min=df.min, Tchar=df.Tchar.±stderrors.sd_Tchar, θchar=df.θchar.±stderrors.sd_θchar, ΔCp=df.ΔCp.±stderrors.sd_ΔCp)
    end
	return check, msg, df_flag, index_flag, res, Telu_max
end

function flagged_loss(meas, df, index_flag)
	a = time_unit_conversion_factor(meas[6])
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

"""
    method_m1(meas, col_input; se_col=true, parallel=false)

Estimation of the three retention parameters ``T_{char}``, ``θ_{char}`` and ``ΔC_p`` including standard errors, see [`RetentionParameterEstimator.stderror`](@ref).

# Arguments
* `meas` ... Tuple with the loaded measurement data, see [`load_chromatograms`](@ref).
* `col_input` ... Named tuple with `col_input.L` the column length in m and `col_input.d` the column diameter in mm. 
* `se_col=true` ... If `true` the standard errors (from the Hessian matrix, see [`RetentionParameterEstimator.stderror`](@ref)) of the estimated parameters are added as separate columns to the result dataframe. If `false` the standard errors are added to the values as `Masurement` type.  
* `parallel=false` ... If `true`, use parallelization for per-substance optimizations and standard error calculations (requires Julia to be started with multiple threads, e.g., `julia -t 4`).

# Output 
* `res` ... Dataframe with the optimized parameters and the found minima.
* `Telu_max` ... The maximum of elution temperatures every solute experiences in the measured programs.
"""   
function method_m1(meas, col_input; se_col=true, method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, parallel=false)
	# definition of the column
	col = GasChromatographySimulator.Column(col_input.L, col_input.d*1e-3, meas[1].df, meas[1].sp, meas[1].gas)
	# calculate start parameters	
	Tchar_est, θchar_est, ΔCp_est, Telu_max = estimate_start_parameter(meas[3], col, meas[2]; time_unit=meas[6])
	# optimize every solute separatly for the 3 remaining parameters `Tchar`, `θchar`, `ΔCp`
	res_ = estimate_parameters(meas[3], meas[4], col, meas[2], Tchar_est, θchar_est, ΔCp_est; mode="Kcentric_single", pout=meas[5], time_unit=meas[6], method=method, opt=std_opt, maxiters=maxiters, maxtime=maxtime, parallel=parallel)[1]
	# calculate the standard errors of the 3 parameters using the hessian matrix
	stderrors = stderror(meas, res_, col_input; opt=opt, parallel=parallel)[1]
	# Match standard errors by name only if parallel (non-parallel: already in correct order)
	if parallel
		sd_Tchar, sd_θchar, sd_ΔCp = match_stderrors_by_name(res_, stderrors)
	else
		sd_Tchar = stderrors.sd_Tchar
		sd_θchar = stderrors.sd_θchar
		sd_ΔCp = stderrors.sd_ΔCp
	end
	# output dataframe
    res = if se_col == true
	    DataFrame(Name=res_.Name, min=res_.min, Tchar=res_.Tchar, Tchar_std=sd_Tchar, θchar=res_.θchar, θchar_std=sd_θchar, ΔCp=res_.ΔCp, ΔCp_std=sd_ΔCp)
    else
	    DataFrame(Name=res_.Name, min=res_.min, Tchar=res_.Tchar.±sd_Tchar, θchar=res_.θchar.±sd_θchar, ΔCp=res_.ΔCp.±sd_ΔCp)
	end
    return res, Telu_max
end

"""
    method_m2(meas; se_col=true, parallel=false)

Estimation of the column diameter ``d`` and three retention parameters ``T_{char}``, ``θ_{char}`` and ``Δ C_p`` including standard errors, see [`RetentionParameterEstimator.stderror`](@ref).
In a first run all four parameters are estimated for every substance separatly, resulting in different optimized column diameters. The mean value of the column diameter is used for 
a second optimization using this mean diameter and optimize the remainig thre retention parameters ``T_{char}``, ``θ_{char}`` and ``Δ C_p``.

# Arguments
* `meas` ... Tuple with the loaded measurement data, see [`load_chromatograms`](@ref).
* `se_col=true` ... If `true` the standard errors (from the Hessian matrix, see [`RetentionParameterEstimator.stderror`](@ref)) of the estimated parameters are added as separate columns to the result dataframe. If `false` the standard errors are added to the values as `Masurement` type.  
* `parallel=false` ... If `true`, use parallelization for per-substance optimizations and standard error calculations (requires Julia to be started with multiple threads, e.g., `julia -t 4`).

# Output 
* `res` ... Dataframe with the optimized parameters and the found minima.
* `Telu_max` ... The maximum of elution temperatures every solute experiences in the measured programs.
"""   
function method_m2(meas; se_col=true, method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, parallel=false)
	# calculate start parameters	
	Tchar_est, θchar_est, ΔCp_est, Telu_max = estimate_start_parameter(meas[3], meas[1], meas[2]; time_unit=meas[6])
	# optimize every solute separatly for the 4 parameters `Tchar`, `θchar`, `ΔCp` and `d`	
	res_dKcentric_single = estimate_parameters(meas[3], meas[4], meas[1], meas[2], Tchar_est, θchar_est, ΔCp_est; pout=meas[5], time_unit=meas[6], mode="dKcentric_single", method=method, opt=std_opt, maxiters=maxiters, maxtime=maxtime, parallel=parallel)[1]
	# define a new column with the mean value of the estimated `d` over all solutes
	new_col = GasChromatographySimulator.Column(meas[1].L, mean(res_dKcentric_single.d), meas[1].df, meas[1].sp, meas[1].gas)
	
	# IMPORTANT: Extract parameters by name to ensure correct ordering (only needed for parallel)
	Tchar_from_res = Array{Float64}(undef, length(meas[4]))
	θchar_from_res = Array{Float64}(undef, length(meas[4]))
	ΔCp_from_res = Array{Float64}(undef, length(meas[4]))
	if parallel
		# Create a dictionary for O(1) lookup instead of O(n) findfirst in loop
		name_to_idx = Dict(name => idx for (idx, name) in enumerate(res_dKcentric_single.Name))
		for (idx, subst_name) in enumerate(meas[4])
			res_idx = get(name_to_idx, subst_name, nothing)
			if isnothing(res_idx)
				error("Substance $subst_name not found in res_dKcentric_single")
			end
			Tchar_from_res[idx] = res_dKcentric_single.Tchar[res_idx]
			θchar_from_res[idx] = res_dKcentric_single.θchar[res_idx]
			ΔCp_from_res[idx] = res_dKcentric_single.ΔCp[res_idx]
		end
	else
		# Non-parallel: DataFrame is already in correct order, use direct indexing
		for idx=1:length(meas[4])
			Tchar_from_res[idx] = res_dKcentric_single.Tchar[idx]
			θchar_from_res[idx] = res_dKcentric_single.θchar[idx]
			ΔCp_from_res[idx] = res_dKcentric_single.ΔCp[idx]
		end
	end
	
	# optimize every solute separatly for the 3 remaining parameters `Tchar`, `θchar`, `ΔCp`
	res_ = estimate_parameters(meas[3], meas[4], new_col, meas[2], Tchar_from_res, θchar_from_res, ΔCp_from_res; pout=meas[5], time_unit=meas[6], mode="Kcentric_single", method=method, opt=std_opt, maxiters=maxiters, maxtime=maxtime, parallel=parallel)[1]

	res_[!, :d] = mean(res_dKcentric_single.d).*ones(length(res_.Name))
	
	res_[!, :d_std] = std(res_dKcentric_single.d).*ones(length(res_.Name))
	# calculate the standard errors of the 3 parameters using the hessian matrix
    res_col = (L=new_col.L, d=res_.d[1])
	stderrors = stderror(meas, res_, res_col; opt=opt, parallel=parallel)[1]
	# Match standard errors by name only if parallel (non-parallel: already in correct order)
	if parallel
		sd_Tchar, sd_θchar, sd_ΔCp = match_stderrors_by_name(res_, stderrors)
	else
		sd_Tchar = stderrors.sd_Tchar
		sd_θchar = stderrors.sd_θchar
		sd_ΔCp = stderrors.sd_ΔCp
	end
	# output dataframe
    res = if se_col == true
	    DataFrame(Name=res_.Name, min=res_.min, Tchar=res_.Tchar, Tchar_std=sd_Tchar, θchar=res_.θchar, θchar_std=sd_θchar, ΔCp=res_.ΔCp, ΔCp_std=sd_ΔCp, d=res_.d, d_std=res_.d_std)
    else
	    DataFrame(Name=res_.Name, min=res_.min, Tchar=res_.Tchar.±sd_Tchar, θchar=res_.θchar.±sd_θchar, ΔCp=res_.ΔCp.±sd_ΔCp, d=res_.d.±res_.d_std)
    end
	return res, Telu_max
end

"""
    method_m3(meas; se_col=true, parallel=false)

Estimation of the column diameter ``d`` and three retention parameters ``T_{char}``, ``θ_{char}`` and ``Δ C_p`` including standard errors, see [`RetentionParameterEstimator.stderror`](@ref).
Brute-force method, where all parameters (`3n+1` for `n` substances) are estimate in one optimization. 

ATTENTION: This method can take long time to finish. The more substances, the longer it takes.

# Arguments
* `meas` ... Tuple with the loaded measurement data, see [`load_chromatograms`](@ref).
* `se_col=true` ... If `true` the standard errors (from the Hessian matrix, see [`RetentionParameterEstimator.stderror`](@ref)) of the estimated parameters are added as separate columns to the result dataframe. If `false` the standard errors are added to the values as `Masurement` type.  
* `parallel=false` ... If `true`, use parallelization for standard error calculations (requires Julia to be started with multiple threads, e.g., `julia -t 4`). Note: The main optimization cannot be parallelized as it's a single joint optimization.

# Output 
* `res` ... Dataframe with the optimized parameters and the found minima.
* `Telu_max` ... The maximum of elution temperatures every solute experiences in the measured programs.
"""  
function method_m3(meas; se_col=true, method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, parallel=false)
	# definition of the column, the given diameter is used as start value
	col = meas[1]
	# calculate start parameters	
	Tchar_est, θchar_est, ΔCp_est, Telu_max = estimate_start_parameter(meas[3], col, meas[2]; time_unit=meas[6])
	# optimize every solute separatly for the 3 remaining parameters `Tchar`, `θchar`, `ΔCp`
	res_ = estimate_parameters(meas[3], meas[4], col, meas[2], Tchar_est, θchar_est, ΔCp_est; mode="dKcentric", pout=meas[5], time_unit=meas[6], method=method, opt=std_opt, maxiters=maxiters, maxtime=maxtime, parallel=parallel)[1]
	# calculate the standard errors of the 3 parameters using the hessian matrix
	stderrors, hessian = stderror_m3(meas, res_; opt=opt, parallel=parallel)
	# Match standard errors by name only if parallel (non-parallel: already in correct order)
	if parallel
		sd_d, sd_Tchar, sd_θchar, sd_ΔCp = match_stderrors_m3_by_name(res_, stderrors)
	else
		sd_d = stderrors.sd_d
		sd_Tchar = stderrors.sd_Tchar
		sd_θchar = stderrors.sd_θchar
		sd_ΔCp = stderrors.sd_ΔCp
	end
	# output dataframe
    res = if se_col == true
	    DataFrame(Name=res_.Name, min=res_.min, Tchar=res_.Tchar, Tchar_std=sd_Tchar, θchar=res_.θchar, θchar_std=sd_θchar, ΔCp=res_.ΔCp, ΔCp_std=sd_ΔCp, d=res_.d, d_std=sd_d)
    else
	    DataFrame(Name=res_.Name, min=res_.min, Tchar=res_.Tchar.±sd_Tchar, θchar=res_.θchar.±sd_θchar, ΔCp=res_.ΔCp.±sd_ΔCp, d=res_.d.±sd_d)
	end
    return res, Telu_max
end

"""
    method_m4(meas; se_col=true, method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, max_alternating_iters=10, tol=1e-6, parallel=false)

Estimation of the column diameter ``d`` and three retention parameters ``T_{char}``, ``θ_{char}`` and ``Δ C_p`` including standard errors, see [`RetentionParameterEstimator.stderror`](@ref).
Uses alternating/block coordinate descent optimization:
1. Initialize `d` from a quick estimate (mean of `dKcentric_single` results with reduced iterations)
2. **Block 1**: Optimize `d` (1D optimization) while fixing substance parameters
3. **Block 2**: Optimize substance parameters (parallelizable) while fixing `d`
4. Iterate steps 2-3 until convergence

This approach properly enforces that `d` is the same for all substances and is more efficient than joint optimization for many substances (>10).

# Arguments
* `meas` ... Tuple with the loaded measurement data, see [`load_chromatograms`](@ref).
* `se_col=true` ... If `true` the standard errors (from the Hessian matrix, see [`RetentionParameterEstimator.stderror`](@ref)) of the estimated parameters are added as separate columns to the result dataframe. If `false` the standard errors are added to the values as `Masurement` type.
* `method=NewtonTrustRegion()` ... Optimization method to use.
* `opt=std_opt` ... Options for the ODE solver.
* `maxiters=10000` ... Maximum number of iterations for each optimization.
* `maxtime=600.0` ... Maximum time for each optimization in seconds.
* `max_alternating_iters=10` ... Maximum number of alternating iterations.
* `tol=1e-6` ... Convergence tolerance for `d` and substance parameters (iteration stops when relative changes in `d`, `Tchar`, `θchar`, and `ΔCp` are all less than `tol`).
* `parallel=false` ... If `true`, use parallelization for per-substance optimizations (Block 2) and standard error calculations (requires Julia to be started with multiple threads, e.g., `julia -t 4`). This can provide significant speedup for many substances.

# Output 
* `res` ... Dataframe with the optimized parameters and the found minima.
* `Telu_max` ... The maximum of elution temperatures every solute experiences in the measured programs.
"""   
function method_m4(meas; se_col=true, method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, max_alternating_iters=10, tol=1e-6, parallel=false)
	# calculate start parameters	
	Tchar_est, θchar_est, ΔCp_est, Telu_max = estimate_start_parameter(meas[3], meas[1], meas[2]; time_unit=meas[6])
	
	# Initialize d: use mean from a quick dKcentric_single pass (or use col.d as fallback)
	# For initialization, we can use a subset or all substances with reduced iterations
	# Here we use all substances but this could be optimized to use a subset
	res_dKcentric_single_init = estimate_parameters(meas[3], meas[4], meas[1], meas[2], Tchar_est, θchar_est, ΔCp_est; 
	                                                pout=meas[5], time_unit=meas[6], mode="dKcentric_single", 
	                                                method=method, opt=std_opt, maxiters=min(maxiters, 1000), 
	                                                maxtime=min(maxtime, 60.0), parallel=parallel)[1]
	d_current = mean(res_dKcentric_single_init.d)
	
	# Use initial estimates from the initialization pass
	# IMPORTANT: Extract parameters by name to ensure correct ordering (only needed for parallel)
	Tchar_current = Array{Float64}(undef, length(meas[4]))
	θchar_current = Array{Float64}(undef, length(meas[4]))
	ΔCp_current = Array{Float64}(undef, length(meas[4]))
	if parallel
		# Create a dictionary for O(1) lookup instead of O(n) findfirst in loop
		name_to_idx = Dict(name => idx for (idx, name) in enumerate(res_dKcentric_single_init.Name))
		for (idx, subst_name) in enumerate(meas[4])
			res_idx = get(name_to_idx, subst_name, nothing)
			if isnothing(res_idx)
				error("Substance $subst_name not found in res_dKcentric_single_init")
			end
			Tchar_current[idx] = res_dKcentric_single_init.Tchar[res_idx]
			θchar_current[idx] = res_dKcentric_single_init.θchar[res_idx]
			ΔCp_current[idx] = res_dKcentric_single_init.ΔCp[res_idx]
		end
	else
		# Non-parallel: DataFrame is already in correct order, use direct indexing
		for idx=1:length(meas[4])
			Tchar_current[idx] = res_dKcentric_single_init.Tchar[idx]
			θchar_current[idx] = res_dKcentric_single_init.θchar[idx]
			ΔCp_current[idx] = res_dKcentric_single_init.ΔCp[idx]
		end
	end
	
	# Alternating optimization
	res_ = nothing
	new_col = nothing
	# Store previous substance parameters for convergence checking
	Tchar_prev = copy(Tchar_current)
	θchar_prev = copy(θchar_current)
	ΔCp_prev = copy(ΔCp_current)
	
	for iter = 1:max_alternating_iters
		d_prev = d_current
		
		# Block 1: Optimize d while fixing substance parameters
		# Use estimate_parameters with mode="d_only" for uniform API
		res_d = estimate_parameters(meas[3], meas[4], meas[1], meas[2], Tchar_current, θchar_current, ΔCp_current; 
		                           pout=meas[5], time_unit=meas[6], mode="d_only", method=method, 
		                           opt=opt, maxiters=maxiters, maxtime=maxtime)[1]
		d_current = res_d.d[1]  # All d values are the same, take first one
		
		# Block 2: Optimize substance parameters while fixing d (can be parallelized)
		new_col = GasChromatographySimulator.Column(meas[1].L, d_current, meas[1].df, meas[1].sp, meas[1].gas)
		res_ = estimate_parameters(meas[3], meas[4], new_col, meas[2], Tchar_current, θchar_current, ΔCp_current; 
		                           pout=meas[5], time_unit=meas[6], mode="Kcentric_single", method=method, 
		                           opt=opt, maxiters=maxiters, maxtime=maxtime, parallel=parallel)[1]
		
		# IMPORTANT: Extract parameters by name to ensure correct ordering (only needed for parallel)
		if parallel
			# Create a dictionary for O(1) lookup instead of O(n) findfirst in loop
			name_to_idx = Dict(name => idx for (idx, name) in enumerate(res_.Name))
			for (idx, subst_name) in enumerate(meas[4])
				res_idx = get(name_to_idx, subst_name, nothing)
				if isnothing(res_idx)
					error("Substance $subst_name not found in res_")
				end
				Tchar_current[idx] = res_.Tchar[res_idx]
				θchar_current[idx] = res_.θchar[res_idx]
				ΔCp_current[idx] = res_.ΔCp[res_idx]
			end
		else
			# Non-parallel: DataFrame is already in correct order, use direct indexing
			for idx=1:length(meas[4])
				Tchar_current[idx] = res_.Tchar[idx]
				θchar_current[idx] = res_.θchar[idx]
				ΔCp_current[idx] = res_.ΔCp[idx]
			end
		end
		
		# Check convergence: both d and substance parameters must have converged
		d_converged = abs(d_current - d_prev) < tol
		
		# Check substance parameter convergence using relative tolerance
		# Use maximum relative change across all substances
		max_tchar_rel_change = maximum(abs.((Tchar_current .- Tchar_prev) ./ (Tchar_prev .+ 1e-10)))  # Add small epsilon to avoid division by zero
		max_thetachar_rel_change = maximum(abs.((θchar_current .- θchar_prev) ./ (θchar_prev .+ 1e-10)))
		max_deltacp_rel_change = maximum(abs.((ΔCp_current .- ΔCp_prev) ./ (abs.(ΔCp_prev) .+ 1e-10)))
		
		# Use same tolerance for relative changes (or could use separate parameter)
		substance_params_converged = (max_tchar_rel_change < tol) && (max_thetachar_rel_change < tol) && (max_deltacp_rel_change < tol)
		
		if d_converged && substance_params_converged
			break
		end
		
		# Update previous values for next iteration
		Tchar_prev = copy(Tchar_current)
		θchar_prev = copy(θchar_current)
		ΔCp_prev = copy(ΔCp_current)
	end
	
	# Ensure new_col is set (should always be set from loop, but just in case)
	if new_col === nothing
		new_col = GasChromatographySimulator.Column(meas[1].L, d_current, meas[1].df, meas[1].sp, meas[1].gas)
	end
	
	# Add d column to results
	res_[!, :d] = d_current.*ones(length(res_.Name))
	
	# Calculate standard error for d using the spread from initialization (or could use Hessian)
	# For now, use std from initialization as approximation
	res_[!, :d_std] = std(res_dKcentric_single_init.d).*ones(length(res_.Name))
	
	# Calculate the standard errors of the 3 parameters using the hessian matrix
    res_col = (L=new_col.L, d=d_current)
	stderrors = stderror(meas, res_, res_col; opt=opt, parallel=parallel)[1]
	# Match standard errors by name only if parallel (non-parallel: already in correct order)
	if parallel
		sd_Tchar, sd_θchar, sd_ΔCp = match_stderrors_by_name(res_, stderrors)
	else
		sd_Tchar = stderrors.sd_Tchar
		sd_θchar = stderrors.sd_θchar
		sd_ΔCp = stderrors.sd_ΔCp
	end
	
	# Output dataframe
    res = if se_col == true
	    DataFrame(Name=res_.Name, min=res_.min, Tchar=res_.Tchar, Tchar_std=sd_Tchar, 
	              θchar=res_.θchar, θchar_std=sd_θchar, ΔCp=res_.ΔCp, ΔCp_std=sd_ΔCp, 
	              d=res_.d, d_std=res_.d_std)
    else
	    DataFrame(Name=res_.Name, min=res_.min, Tchar=res_.Tchar.±sd_Tchar, 
	              θchar=res_.θchar.±sd_θchar, ΔCp=res_.ΔCp.±sd_ΔCp, 
	              d=res_.d.±res_.d_std)
    end
	return res, Telu_max
end


"""
    stderror(meas, res, col_input; opt=std_opt, metric="squared", parallel=false)

Calculation of the standard error of the found optimized parameters using the hessian matrix at the optima.

# Attention
The used loss-function is hard coded in the function `opt_Kcentric` and has to be changed if another loss-function is used.

# Arguments
* `meas` ... Tuple with the loaded measurement data, see [`load_chromatograms`](@ref).
* `res` ... Dataframe with the result of the optimization, see [`estimate_parameters`](@ref).
* `col_input` ... Named tuple with `col_input.L` the column length in m and `col_input.d` the column diameter in mm.

Optional parameters:
* `opt=std_opt` ... Options for the ODE solver.
* `metric="squared"` ... Metric used for the loss function (`"squared"` or `"abs"`).
* `parallel=false` ... If `true`, use parallelization for standard error calculations (requires Julia to be started with multiple threads, e.g., `julia -t 4`).

# Output
* `stderrors` ... Dataframe with the standard errors of the optimized parameters.
* `Hessian` ... The hessian matrix at the found optima. 
"""

function stderror(meas, res, col_input; opt=std_opt, metric="squared", parallel=false)
	sdTchar = Array{Float64}(undef, size(res)[1])
	sdθchar = Array{Float64}(undef, size(res)[1])
	sdΔCp = Array{Float64}(undef, size(res)[1])
	Hessian = Array{Any}(undef, size(res)[1])

    a = time_unit_conversion_factor(meas[6])
    
    # Create a dictionary for O(1) lookup instead of O(n) findfirst in loop
    meas4_to_idx = Dict(name => idx for (idx, name) in enumerate(meas[4]))
    
    if parallel
        Base.Threads.@threads for i=1:size(res)[1]
            # IMPORTANT: Match by name, not by index
            subst_name = res.Name[i]
            meas_idx = get(meas4_to_idx, subst_name, nothing)
            if isnothing(meas_idx)
                error("Substance $subst_name not found in meas[4]")
            end
            
            # filter-out missing values:
            tR_meas = meas[3][!,meas_idx+1].*a
            tRs_, prog_, subst_list_ = prepare_single_substance_data(tR_meas, meas[2], subst_name)
            
            p = (tRs_, subst_list_, col_input.L, col_input.d, prog_, opt, meas[1].gas, metric) 
        # the loss-function used in the optimization !!! HAS TO BE CHANGED !!!
		LF(x) = opt_Kcentric(x, p)
		# the hessian matrix of the loss-function, calculated with ForwardDiff.jl
		H(x) = ForwardDiff.hessian(LF, x)
		# the hessian matrix at the found optima
		Hessian[i] = H([res.Tchar[i], res.θchar[i], res.ΔCp[i]])
		# the calculated standard errors of the parameters
		sdTchar[i] = sqrt.(abs.(inv(Hessian[i])))[1,1]
		sdθchar[i] = sqrt.(abs.(inv(Hessian[i])))[2,2]
		sdΔCp[i] = sqrt.(abs.(inv(Hessian[i])))[3,3]
	end
    else
	for i=1:size(res)[1]
        # IMPORTANT: Match by name, not by index
        subst_name = res.Name[i]
        meas_idx = get(meas4_to_idx, subst_name, nothing)
        if isnothing(meas_idx)
            error("Substance $subst_name not found in meas[4]")
        end
		
        # filter-out missing values:
        tR_meas = meas[3][!,meas_idx+1].*a
            tRs_, prog_, subst_list_ = prepare_single_substance_data(tR_meas, meas[2], subst_name)
        
		p = (tRs_, subst_list_, col_input.L, col_input.d, prog_, opt, meas[1].gas, metric) 
        # the loss-function used in the optimization !!! HAS TO BE CHANGED !!!
		LF(x) = opt_Kcentric(x, p)
		# the hessian matrix of the loss-function, calculated with ForwardDiff.jl
		H(x) = ForwardDiff.hessian(LF, x)
		# the hessian matrix at the found optima
		Hessian[i] = H([res.Tchar[i], res.θchar[i], res.ΔCp[i]])
		# the calculated standard errors of the parameters
		sdTchar[i] = sqrt.(abs.(inv(Hessian[i])))[1,1]
		sdθchar[i] = sqrt.(abs.(inv(Hessian[i])))[2,2]
		sdΔCp[i] = sqrt.(abs.(inv(Hessian[i])))[3,3]
        end
	end
	stderrors = DataFrame(Name=res.Name, sd_Tchar=sdTchar, sd_θchar=sdθchar, sd_ΔCp=sdΔCp)
	return stderrors, Hessian
end

"""
    match_stderrors_by_name(res, stderrors)

Match standard errors to results by substance name, ensuring correct association regardless of DataFrame order.
Returns arrays of standard errors in the same order as `res.Name`.
"""
function match_stderrors_by_name(res, stderrors)
	sd_Tchar = Array{Float64}(undef, length(res.Name))
	sd_θchar = Array{Float64}(undef, length(res.Name))
	sd_ΔCp = Array{Float64}(undef, length(res.Name))
	
	# Create a dictionary for O(1) lookup instead of O(n) findfirst in loop
	name_to_idx = Dict(name => idx for (idx, name) in enumerate(stderrors.Name))
	
	for (idx, subst_name) in enumerate(res.Name)
		stderr_idx = get(name_to_idx, subst_name, nothing)
		if isnothing(stderr_idx)
			error("Substance $subst_name not found in stderrors")
		end
		sd_Tchar[idx] = stderrors.sd_Tchar[stderr_idx]
		sd_θchar[idx] = stderrors.sd_θchar[stderr_idx]
		sd_ΔCp[idx] = stderrors.sd_ΔCp[stderr_idx]
	end
	
	return sd_Tchar, sd_θchar, sd_ΔCp
end

"""
    stderror_m3(meas, res; opt=std_opt, metric="squared", parallel=false)

Calculation of the standard error of the found optimized parameters (including column diameter `d`) using the hessian matrix at the optima.
Used for `method_m3` which optimizes all parameters jointly.

# Attention
The used loss-function is hard coded in the function `opt_dKcentric` and has to be changed if another loss-function is used.

# Arguments
* `meas` ... Tuple with the loaded measurement data, see [`load_chromatograms`](@ref).
* `res` ... Dataframe with the result of the optimization, see [`estimate_parameters`](@ref).

Optional parameters:
* `opt=std_opt` ... Options for the ODE solver.
* `metric="squared"` ... Metric used for the loss function (`"squared"` or `"abs"`).
* `parallel=false` ... If `true`, use parallelization for standard error calculations (requires Julia to be started with multiple threads, e.g., `julia -t 4`).

# Output
* `stderrors` ... Dataframe with the standard errors of the optimized parameters (including `d`).
* `Hessian` ... The hessian matrix at the found optima. 
"""
function stderror_m3(meas, res; opt=std_opt, metric="squared", parallel=false)
	sd_d = Array{Float64}(undef, size(res)[1])
    sdTchar = Array{Float64}(undef, size(res)[1])
	sdθchar = Array{Float64}(undef, size(res)[1])
	sdΔCp = Array{Float64}(undef, size(res)[1])
	Hessian = Array{Any}(undef, size(res)[1])

    a = time_unit_conversion_factor(meas[6])
    
    if parallel
        # Create a dictionary for O(1) lookup instead of O(n) findfirst in loop
        meas4_to_idx = Dict(name => idx for (idx, name) in enumerate(meas[4]))
        Base.Threads.@threads for i=1:size(res)[1]
            # IMPORTANT: Match by name, not by index
            subst_name = res.Name[i]
            meas_idx = get(meas4_to_idx, subst_name, nothing)
            if isnothing(meas_idx)
                error("Substance $subst_name not found in meas[4]")
            end
            
            # filter-out missing values:
            tR_meas = meas[3][!,meas_idx+1].*a
            tRs_, prog_, subst_list_ = prepare_single_substance_data(tR_meas, meas[2], subst_name)
            
            p = (tRs_, subst_list_, meas[1].L, prog_, opt, meas[1].gas, metric) 
            # the loss-function used in the optimization !!! HAS TO BE CHANGED !!!
            LF(x) = opt_dKcentric(x, p)
            # the hessian matrix of the loss-function, calculated with ForwardDiff.jl
            H(x) = ForwardDiff.hessian(LF, x)
            # the hessian matrix at the found optima
            Hessian[i] = H([res.d[i], res.Tchar[i], res.θchar[i], res.ΔCp[i]])
            # the calculated standard errors of the parameters
            sd_d[i] = sqrt.(abs.(inv(Hessian[i])))[1,1]
            sdTchar[i] = sqrt.(abs.(inv(Hessian[i])))[2,2]
            sdθchar[i] = sqrt.(abs.(inv(Hessian[i])))[3,3]
            sdΔCp[i] = sqrt.(abs.(inv(Hessian[i])))[4,4]
        end
    else
	for i=1:size(res)[1]
        # Non-parallel: res is in same order as meas[4], use direct indexing
        # filter-out missing values:
        tR_meas = meas[3][!,i+1].*a
            tRs_, prog_, subst_list_ = prepare_single_substance_data(tR_meas, meas[2], res.Name[i])
        
		p = (tRs_, subst_list_, meas[1].L, prog_, opt, meas[1].gas, metric) 
        # the loss-function used in the optimization !!! HAS TO BE CHANGED !!!
		LF(x) = opt_dKcentric(x, p)
		# the hessian matrix of the loss-function, calculated with ForwardDiff.jl
		H(x) = ForwardDiff.hessian(LF, x)
		# the hessian matrix at the found optima
		Hessian[i] = H([res.d[i], res.Tchar[i], res.θchar[i], res.ΔCp[i]])
		# the calculated standard errors of the parameters
        sd_d[i] = sqrt.(abs.(inv(Hessian[i])))[1,1]
		sdTchar[i] = sqrt.(abs.(inv(Hessian[i])))[2,2]
		sdθchar[i] = sqrt.(abs.(inv(Hessian[i])))[3,3]
		sdΔCp[i] = sqrt.(abs.(inv(Hessian[i])))[4,4]
        end
	end
	stderrors = DataFrame(Name=res.Name, sd_d=sd_d, sd_Tchar=sdTchar, sd_θchar=sdθchar, sd_ΔCp=sdΔCp)
	return stderrors, Hessian
end

"""
    match_stderrors_m3_by_name(res, stderrors)

Match standard errors from stderror_m3 to results by substance name, ensuring correct association regardless of DataFrame order.
Returns arrays of standard errors in the same order as `res.Name`.
"""
function match_stderrors_m3_by_name(res, stderrors)
	sd_d = Array{Float64}(undef, length(res.Name))
	sd_Tchar = Array{Float64}(undef, length(res.Name))
	sd_θchar = Array{Float64}(undef, length(res.Name))
	sd_ΔCp = Array{Float64}(undef, length(res.Name))
	
	# Create a dictionary for O(1) lookup instead of O(n) findfirst in loop
	name_to_idx = Dict(name => idx for (idx, name) in enumerate(stderrors.Name))
	
	for (idx, subst_name) in enumerate(res.Name)
		stderr_idx = get(name_to_idx, subst_name, nothing)
		if isnothing(stderr_idx)
			error("Substance $subst_name not found in stderrors")
		end
		sd_d[idx] = stderrors.sd_d[stderr_idx]
		sd_Tchar[idx] = stderrors.sd_Tchar[stderr_idx]
		sd_θchar[idx] = stderrors.sd_θchar[stderr_idx]
		sd_ΔCp[idx] = stderrors.sd_ΔCp[stderr_idx]
	end
	
	return sd_d, sd_Tchar, sd_θchar, sd_ΔCp
end