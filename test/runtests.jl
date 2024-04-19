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

#    pl, sol = GasChromatographySimulator.simulate(par)

#    Tchar = par.sub[3].Tchar
#    θchar = par.sub[3].θchar
#    ΔCp = par.sub[3].ΔCp
#    opt = par.opt

    ##tR = RetentionParameterEstimator.tR_calc(Tchar, θchar, ΔCp, df/d, L, d, df, par.prog, gas)

    ##pl.tR[3] == tR # false, because of odesys=true for GasChromatographySimulator
        # and in RetentionParameterEstimator only the migration ODE is used 
    ##@test isapprox(pl.tR[3], tR, atol=1e-4)

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
file = "./data/meas_test.csv"
meas = RetentionParameterEstimator.load_chromatograms(file; filter_missing=true)
@test meas[4][1] == "2-Octanone"

meas_select = RetentionParameterEstimator.filter_selected_measurements(meas, ["meas3", "meas4", "meas5"], ["2-Octanone"]);

file = "./data/meas_test_2.csv" # from Email from mathijs.ruiter@go-jsb.com -> Issue #29
meas_ = RetentionParameterEstimator.load_chromatograms(file; filter_missing=true)
@test meas_[2][1].Fpin_steps[3] == meas_[2][1].Fpin_steps[3]
# Loss.jl

# Estimate_Start_Values.jl

# Optimization.jl
#df, sol = RetentionParameterEstimator.estimate_parameters(meas; mode="Kcentric_single")
#@test isapprox(df.Tchar[1], 400.15; atol = 0.01)
#@test isapprox(sol[2].minimum, 0.009; atol = 0.0001)

col_input = (L = meas_select[1].L, d = meas_select[1].d*1000)
check, msg, df_flag, index_flag, res, Telu_max = RetentionParameterEstimator.check_measurement(meas_select, col_input; min_th=0.1, loss_th=1.0, se_col=false)
@test check == true
@test isapprox(res.Tchar[1], 400.0; atol = 1.0)
#@test isapprox(res.min[2], 0.009; atol = 0.001)

#col_input_ = (L = meas[1].L, d = meas[1].d*1000*1.1)
#check_, msg_, df_flag_, index_flag_, res_, Telu_max_ = RetentionParameterEstimator.check_measurement(meas, col_input_; min_th=0.1, loss_th=1.0)
#@test msg_ == "discrapancy of retention times detected"
#@test isapprox(res_.Tchar[1], 415.5; atol = 0.1)
#@test isapprox(res_.min[2], 0.21; atol = 0.01)

res_m1, Telu_max_m1 = RetentionParameterEstimator.method_m1(meas_select, col_input, se_col=false)
@test res_m1.Tchar == res.Tchar
@test Telu_max_m1 == Telu_max

res_m2, Telu_max_m2 = RetentionParameterEstimator.method_m2(meas_select, se_col=true)
@test isapprox(res_m2.d[1], 0.00024, atol=0.00001)  