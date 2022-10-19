# functions defining the loss-function

"""
    tR_calc(Tchar, θchar, ΔCp, L, d, prog, opt, gas)

Calculates the retention time tR for a solute with the K-centric parameters `Tchar` `θchar` and `ΔCp` for a column with length `L` and diameter `d`, the (conventional) program `prog`, options `opt` and mobile phase `gas`.
    
"""
function tR_calc(Tchar, θchar, ΔCp, L, d, prog, opt, gas)
	# df has no influence on the result (is hidden in Tchar, θchar, ΔCp)
	R = 8.3145
	k(x,t,Tchar,θchar,ΔCp) = exp((ΔCp/R + Tchar/θchar)*(Tchar/prog.T_itp(x,t)-1) + ΔCp/R*log(prog.T_itp(x,t)/Tchar))
	rM(x,t,L,d) = GasChromatographySimulator.mobile_phase_residency(x,t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, L, d, gas; ng=opt.ng, vis=opt.vis, control=opt.control)
	r(t,p,x) = (1+k(x,t,p[1],p[2],p[3]))*rM(x,t,p[4],p[5])
	t₀ = 0.0
	xspan = (0.0, L)
	p = [Tchar, θchar, ΔCp, L, d]
	prob = ODEProblem(r, t₀, xspan, p)
	solution = solve(prob, alg=opt.alg, abstol=opt.abstol, reltol=opt.reltol)
	tR = solution.u[end]
	return tR
end

"""
    loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas)

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
function loss(tR, Tchar, θchar, ΔCp, L, d, prog, opt, gas)
	# loss function as sum over all programs and solutes, fixed L and d
	#gas = "He"
	#tR = p[1]
	#L = p[2]
	#d = p[3]
	#prog = p[4]
	#opt = p[5]

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
	#Tchar = p_Kcentric[1:ns] # Array length = number solutes
	#θchar = p_Kcentric[ns+1:2*ns] # Array length = number solutes
	#ΔCp = p_Kcentric[2*ns+1:3*ns] # Array length = number solutes
	tRcalc = Array{Any}(undef, np, ns)
	for j=1:ns
		for i=1:np
			tRcalc[i,j] = tR_calc(Tchar[j], θchar[j], ΔCp[j], L, d, prog[i], opt, gas)
		end
	end
	return sum((tR.-tRcalc).^2), tRcalc
end	