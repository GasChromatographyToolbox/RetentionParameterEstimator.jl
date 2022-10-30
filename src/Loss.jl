# functions defining the loss-function

"""
    tR_calc(Tchar, θchar, ΔCp, φ₀, L, d, df, prog, opt, gas)

Calculates the retention time tR for a solute with the K-centric parameters `Tchar` `θchar` and `ΔCp` and the corresponding dimensionless film 
thickness `φ₀` for a column with length `L`, internal diameter `d` and film thickness `df`, the (conventional) program `prog`, options `opt` and mobile phase `gas`.
For this calculation only the ODE for the migration of a solute in a GC column is solved, using the function `GasChromatographySimulator.solving_migration`.
"""
function tR_calc(Tchar, θchar, ΔCp, φ₀, L, d, df, prog, opt, gas)
	solution = GasChromatographySimulator.solving_migration(Tchar, θchar, ΔCp, φ₀, L, d, df, prog, opt, gas)
	tR = solution.u[end]
	return tR
end
# make a version, where all functions are defined here and instead d -> λ and instead df -> φ

function tR_calc_(Tchar, θchar, ΔCp, λ, φ, L, φ₀, prog, opt, gas)
	solution = GasChromatographySimulator.solving_migration(Tchar, θchar, ΔCp, φ₀, L, L/λ, L/λ*φ, prog, opt, gas)
	tR = solution.u[end]
	return tR
end

function tR_calc(A, B, C, L, λ, β, prog, opt, gas)
	k(x,t,A,B,C,β) = exp(A + B/prog.T_itp(x,t) + C*log(prog.T_itp(x,t)) - log(β))
	#rM(x,t,L,d) = GasChromatographySimulator.mobile_phase_residency(x,t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, L, d, gas; ng=opt.ng, vis=opt.vis, control=opt.control)
	rM(x,t,L,λ) = 64 * sqrt(prog.pin_itp(t)^2 - x/L*(prog.pin_itp(t)^2-prog.pout_itp(t)^2)) / (prog.pin_itp(t)^2-prog.pout_itp(t)^2) * λ^2/L * GasChromatoggraphySimulator.viscosity(x, t, T_itp, gas, vis=vis) * prog.T_itp(x, t)
	r(t,p,x) = (1+k(x,t,p[1],p[2],p[3], p[6]))*rM(x,t,p[4],p[5])
	t₀ = 0.0
	xspan = (0.0, L)
	p = [A, B, C, L, λ, β]
	prob = ODEProblem(r, t₀, xspan, p)
	solution = solve(prob, alg=opt.alg, abstol=opt.abstol, reltol=opt.reltol)
	tR = solution.u[end]
	return tR
end

"""
    loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas; metric="quadratic")

Loss function as sum of squares of the residuals between the measured and calculated retention times.
    
# Arguments
* `tR` ... mxn-array of the measured retention times in seconds.
* `Tchar` ... n-array of characteristic temperatures in K.
* `θchar` ... n-array of characteristic constants in °C.
* `ΔCp` ... n-array of the change of adiabatic heat capacity in J mol^-1 K^-1.
* `L` ... number of the length of the column in m.
* `d` ... number of the diameters of the column in m.
* `prog` ... m-array of structure GasChromatographySimulator.Programs containing the definition of the GC-programs.
* `opt` ... struture GasChromatographySimulator.Options containing the settings for options for the simulation.
* `gas` ... string of name of the mobile phase gas. 

# Output
The output is a tuple of the following quantites:
* `sum((tR.-tRcalc).^2)` ... sum of the squared residuals over m GC-programs and n solutes.
* `tRcalc` ... mxn-array of the calculated retention times
"""
function loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas; metric="quadratic")
	if length(size(tR)) == 1
		ns = 1 # number of solutes
		np = size(tR)[1] # number of programs
	elseif length(size(tR)) == 0
		ns = 1
		np = 1
	else
		ns = size(tR)[2]
		np = size(tR)[1]
	end
	tRcalc = Array{Any}(undef, np, ns)
	for j=1:ns
		for i=1:np
			tRcalc[i,j] = tR_calc(Tchar[j], θchar[j], ΔCp[j], 1e-3, L, d, d*1e-3, prog[i], opt, gas)
		end
	end
	if metric == "abs"
		l = sum(abs.(tR.-tRcalc))/(ns*np)
	elseif metric == "quadratic"
		l = sum((tR.-tRcalc).^2)/(ns*np)
	else
		l = sum((tR.-tRcalc).^2)/(ns*np)
	end
	return l, tRcalc
end

function loss(tR, Tchar, θchar, ΔCp, φ₀, L, d, df, prog, opt, gas; metric="quadratic")
	if length(size(tR)) == 1
		ns = 1 # number of solutes
		np = size(tR)[1] # number of programs
	elseif length(size(tR)) == 0
		ns = 1
		np = 1
	else
		ns = size(tR)[2]
		np = size(tR)[1]
	end
	tRcalc = Array{Any}(undef, np, ns)
	for j=1:ns
		for i=1:np
			tRcalc[i,j] = tR_calc(Tchar[j], θchar[j], ΔCp[j], φ₀, L, d, df, prog[i], opt, gas)
		end
	end
	if metric == "abs"
		l = sum(abs.(tR.-tRcalc))/(ns*np)
	elseif metric == "quadratic"
		l = sum((tR.-tRcalc).^2)/(ns*np)
	else
		l = sum((tR.-tRcalc).^2)/(ns*np)
	end
	return l, tRcalc
end

function loss(tR, A, B, C, L, d, df, prog, opt, gas; metric="quadratic")
	# loss function as sum over all programs and solutes
	if length(size(tR)) == 1
		ns = 1 # number of solutes
		np = size(tR)[1] # number of programs
	elseif length(size(tR)) == 0
		ns = 1
		np = 1
	else
		ns = size(tR)[2]
		np = size(tR)[1]
	end
	tRcalc = Array{Any}(undef, np, ns)
	for j=1:ns
		for i=1:np
			tRcalc[i,j] = tR_calc(A[j], B[j], C[j], L, d, df, prog[i], opt, gas)
		end
	end
	if metric == "abs"
		l = sum(abs.(tR.-tRcalc))/(ns*np)
	elseif metric == "quadratic"
		l = sum((tR.-tRcalc).^2)/(ns*np)
	else
		l = sum((tR.-tRcalc).^2)/(ns*np)
	end
	return l, tRcalc
end	

