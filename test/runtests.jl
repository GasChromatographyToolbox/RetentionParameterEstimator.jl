using Test, RetentionParameterEstimator, GasChromatographySimulator

L = 30.0
d = 0.25e-3
df = 0.25e-6
sp = "SLB5ms"
gas = "He"
TP = [40.0, 3.0, 15.0, 340.0, 5.0]
PP = [150000.0, 3.0, 5000.0, 250000.0, 5.0]
solutes = ["Decane", "2-Octanol", "Pentadecane"]
db_path = "./data"
db_file = "Database_SLB5ms.csv"
pout = "vacuum"
time_unit = "min"

# Simulation_Test.jl
    par = RetentionParameterEstimator.conventional_GC(L, d, df, sp, gas, TP, PP, solutes, db_path, db_file; pout=pout, time_unit=time_unit)
    @test par.prog.time_steps[2] == 3.0*60.0

    tR_sim, par_sim = RetentionParameterEstimator.sim_test_chrom(L, d, df, sp, gas, [TP, TP], [PP, PP], solutes, db_path, db_file; pout=pout, time_unit=time_unit)
    @test tR_sim[1,1] == tR_sim[2,1]
    # more meaningful test?

    # compare GasChromatographySimulator.simulate() and RetentionParameterEstimator.tR_calc() 

    pl, sol = GasChromatographySimulator.simulate(par)

    Tchar = par.sub[3].Tchar
    θchar = par.sub[3].θchar
    ΔCp = par.sub[3].ΔCp
    opt = par.opt

    tR = RetentionParameterEstimator.tR_calc(Tchar, θchar, ΔCp, df/d, L, d, df, par.prog, gas)

    pl.tR[3] == tR # false, because of odesys=true for GasChromatographySimulator
        # and in RetentionParameterEstimator only the migration ODE is used 
    @test isapprox(pl.tR[3], tR, atol=1e-4)

    #opt_ = GasChromatographySimulator.Options(ng=true)
    #par_ = GasChromatographySimulator.Parameters(par.col, par.prog, par.sub, opt_)
    #pl_, sol_ = GasChromatographySimulator.simulate(par_)

    #opt__ = GasChromatographySimulator.Options(ng=true, odesys=false)
    #par__ = GasChromatographySimulator.Parameters(par.col, par.prog, par.sub, opt__)
    #pl__, sol__ = GasChromatographySimulator.simulate(par__)

    #opt___ = GasChromatographySimulator.Options(ng=false, odesys=false)
    #par___ = GasChromatographySimulator.Parameters(par.col, par.prog, par.sub, opt___)
    #pl___, sol___ = GasChromatographySimulator.simulate(par___)

    #[pl.tR[3], pl_.tR[3], pl__.tR[3], pl___.tR[3]].-tR
    # lowest difference for opt__ and opt___ to tR ≈ 1e-5

# Load.jl

# Loss.jl

# Estimate_Start_Values.jl

# Optimization.jl