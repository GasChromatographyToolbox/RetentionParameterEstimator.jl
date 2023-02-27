# functions used for the optimization of the loss-function 




#------------------
# Optimization only for the K-centric retention parameters
"""
    opt_Kcentric(x, p)

Function used for optimization of the loss-function in regards to the three K-centric parameters.

# Arguments
* `x` ... 3n-vector of the three K-centric parameters of n solutes. Elements 1:n are Tchar, n+1:2n are θchar and 2n+1:3n are ΔCp values.
* `p` ... vector containing the fixed parameters:
    * `tR = p[1]` ... mxn-array of the measured retention times in seconds.
    * `L = p[2]` ... number of the length of the column in m.
    * `d = p[3]` ... number of the diameters of the column in m.
    * `prog = p[4]` ... m-array of structure GasChromatographySimulator.Programs containing the definition of the GC-programs.
    * `opt = p[5]` ... struture GasChromatographySimulator.Options containing the settings for options for the simulation.
    * `gas = p[6]` ... string of name of the mobile phase gas. 
    * `metric = p[7]` ... string of the metric used for the loss function (`squared` or `abs`). 

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
    return loss(tR, Tchar, θchar, ΔCp, L, d, prog, gas; opt=opt, metric=metric)[1]
end

function opt_Kcentric_(x, p)
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
end

"""
    optimize_Kcentric(tR, col, prog, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, metric="squared") 

Optimization regarding the estimization of the retention parameters `Tchar`, `θchar` and `ΔCp`. The initial guess is a vector (`Tchar_e`, `θchar_e` and `ΔCp_e`) for optimization
algorithms, which do not need lower/upper bounds. If a method is used, which needs lower/upper bounds the initial parameters should be a matrix, with the first column (`[1,:]`)
beeing the initial guess, the second column (`[2,:]`) the lower bound and the third column (`[3,:]`) the upper bound.
"""
function optimize_Kcentric(tR, col, prog, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number
    p = [tR, col.L, col.d, prog, opt, col.gas, metric]
	x0 = [Tchar_e; θchar_e; ΔCp_e]
	optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end

function optimize_Kcentric(tR, col, prog, Tchar_e::Matrix{T}, θchar_e::Matrix{T}, ΔCp_e::Matrix{T}; method=BBO_adaptive_de_rand_1_bin_radiuslimited(), opt=std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number # here default method should be one which needs bounds
    p = [tR, col.L, col.d, prog, opt, col.gas, metric]
	x0 = [Tchar_e[1,:]; θchar_e[1,:]; ΔCp_e[1,:]]
    lb = [Tchar_e[2,:]; θchar_e[2,:]; ΔCp_e[2,:]]
    ub = [Tchar_e[3,:]; θchar_e[3,:]; ΔCp_e[3,:]]
	optf = OptimizationFunction(opt_Kcentric, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end

function optimize_Kcentric_(tR, col, prog, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number
    p = [tR, col.df/col.d, col.L, col.d, col.df, prog, opt, col.gas, metric]
	x0 = [Tchar_e; θchar_e; ΔCp_e]
	optf = OptimizationFunction(opt_Kcentric_, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end

"""
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
end

function opt_dKcentric_(x, p)
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
end

"""
    optimize_dKcentric(tR, col, prog, d_e, Tchar_e, θchar_e, ΔCp_e; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, metric="squared") 

Optimization regarding the estimization of the column diameter `d` and the retention parameters `Tchar`, `θchar` and `ΔCp`. The initial guess is a number (`d_e`) or a vector (`Tchar_e`, `θchar_e` and `ΔCp_e`) for optimization
algorithms, which do not need lower/upper bounds. If a method is used, which needs lower/upper bounds the initial parameters should be a vector of length 3 (`d_e`) or a matrix (`Tchar_e`, `θchar_e` and `ΔCp_e`), with the first element/column 
beeing the initial guess, the second element/column the lower bound and the third element/column the upper bound.
"""
function optimize_dKcentric(tR, col, prog, d_e::Number, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number
    p = [tR, col.L, prog, opt, col.gas, metric]
	x0 = [d_e; Tchar_e; θchar_e; ΔCp_e]
	optf = OptimizationFunction(opt_dKcentric, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end


function optimize_dKcentric(tR, col, prog, d_e::Vector{T}, Tchar_e::Matrix{T}, θchar_e::Matrix{T}, ΔCp_e::Matrix{T}; method=BBO_adaptive_de_rand_1_bin_radiuslimited(), opt=opt_std, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number  # here default method should be one which needs bounds  
    p = [tR, col.L, prog, opt, col.gas, metric]
	x0 = [d_e[1]; Tchar_e[1,:]; θchar_e[1,:]; ΔCp_e[1,:]]
    lb = [d_e[2]; Tchar_e[2,:]; θchar_e[2,:]; ΔCp_e[2,:]]
    ub = [d_e[3]; Tchar_e[3,:]; θchar_e[3,:]; ΔCp_e[3,:]]
	optf = OptimizationFunction(opt_dKcentric, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, lb=lb, ub=ub, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end

function optimize_dKcentric_(tR, col, prog, d_e::Number, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, metric="squared") where T<:Number
    p = [tR, col.df/col.d, col.L, col.df, prog, opt, col.gas, metric]
	x0 = [d_e; Tchar_e; θchar_e; ΔCp_e]
	optf = OptimizationFunction(opt_dKcentric_, Optimization.AutoForwardDiff())
	prob = Optimization.OptimizationProblem(optf, x0, p, f_calls_limit=maxiters)
    opt_sol = solve(prob, method, maxiters=maxiters, maxtime=maxtime)
	return opt_sol
end

function estimate_parameters(tRs, solute_names, col, prog, rp1_e, rp2_e, rp3_e; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, mode="dKcentric", metric="squared", pout="vacuum", time_unit="min")
    # mode = "Kcentric", "Kcentric_single", "dKcentric", "dKcentric_single"
    if time_unit == "min"
        a = 60.0
    else
        a = 1.0
    end
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
end

"""
    estimate_parameters()

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
* `metric="squared"` ... used metric for the loss function ("squared" or "abs")

# Output
* `df` ... DataFrame with the columns `Name` (solute names), `d` (estimated column diameter, optional), `Tchar` (estimated Tchar), `θchar` (estimated θchar), `ΔCp` (estimated ΔCp) and `min` (value of the loss function at the found optima)
* `sol` ... Array of `SciMLBase.OptimizationSolution` with the results of the optimization with some additional informations.
"""    
function estimate_parameters(chrom; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, mode="dKcentric", metric="squared")
	Tchar_est, θchar_est, ΔCp_est, Telu_max = estimate_start_parameter(chrom[3], chrom[1], chrom[2]; time_unit=chrom[6])
	return estimate_parameters(chrom[3], chrom[4], chrom[1], chrom[2], Tchar_est, θchar_est, ΔCp_est; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, mode=mode, metric=metric, pout=chrom[5], time_unit=chrom[6])
end

function estimate_parameters_(tRs, solute_names, col, prog, rp1_e, rp2_e, rp3_e; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, maxtime=600.0, mode="dKcentric", metric="squared", pout="vacuum", time_unit="min")
    # mode = "Kcentric", "Kcentric_single", "dKcentric", "dKcentric_single"
    if time_unit == "min"
        a = 60.0
    else
        a = 1.0
    end
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
            min[j] = sol[j].minimum
            #retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
    elseif mode == "Kcentric"
        sol = optimize_Kcentric_(tR_meas, col, prog, rp1_e, rp2_e, rp3_e; method=method, opt=opt, maxiters=maxiters, maxtime=maxtime, metric=metric)
        rp1 = sol[1:ns] # Array length = number solutes
        rp2 = sol[ns+1:2*ns] # Array length = number solutes
        rp3 = sol[2*ns+1:3*ns] # Array length = number solutes
        for j=1:ns
            min[j] = sol.minimum
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
            min[j] = sol.minimum
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
            min[j] = sol[j].minimum
            #retcode[j] = sol[j].retcode
        end
        df = DataFrame(Name=solute_names, d=d, Tchar=rp1, θchar=rp2, ΔCp=rp3, min=min)#, retcode=retcode)
    end
	
	return df, sol
end