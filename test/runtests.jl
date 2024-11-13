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
@testset "Simulation Test" begin
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
    end

# Load.jl
@testset "Loading Data" begin
    file = "./data/meas_test.csv"
    meas = RetentionParameterEstimator.load_chromatograms(file)
    @test meas[4][1] == "2-Octanone"

    meas_select = RetentionParameterEstimator.filter_selected_measurements(meas, ["meas3", "meas4", "meas5"], ["2-Octanone"]);
    @test length(meas_select[4]) == 1

    file = "./data/meas_test_2.csv" # from Email from mathijs.ruiter@go-jsb.com -> Issue #29
    meas_2 = RetentionParameterEstimator.load_chromatograms(file)
    @test meas_2[2][1].Fpin_steps[3] == meas_2[2][1].Fpin_steps[3]

    file = "./data/meas_test_missing.csv" # from some tR are missing
    meas_missing = RetentionParameterEstimator.load_chromatograms(file)
    @test any(ismissing, meas_missing[3][:,2]) # tR of first substance should contain missing values 
    @test typeof(meas_missing[2]) == Array{RetentionParameterEstimator.GasChromatographySimulator.Program, 1} # correct type for the programs
    @test all(x -> typeof(x) == String, meas_missing[4])  # All solute names should be strings
    # TODO: these verifications should be done while loading the data and give an error/warning if it is not the case
    @test length(unique(meas_missing[4])) == length(meas_missing[4]) # No duplicate solutes 
    @test all(r -> any(!ismissing, r), eachrow(meas_missing[3]))  # Each solute has some measurements
    @test all(c -> any(!ismissing, c), eachcol(meas_missing[3]))  # Each measurement has some solutes
end

# Loss.jl

# Estimate_Start_Values.jl

# Optimization.jl
@testset "Optimization" begin
    file = "./data/meas_test.csv"
    meas = RetentionParameterEstimator.load_chromatograms(file)
    meas_select = RetentionParameterEstimator.filter_selected_measurements(meas, ["meas3", "meas4", "meas5"], ["2-Octanone"]);

    col_input = (L = meas_select[1].L, d = meas_select[1].d*1000)
    check, msg, df_flag, index_flag, res, Telu_max = RetentionParameterEstimator.check_measurement(meas_select, col_input; min_th=0.1, loss_th=1.0, se_col=false)
    @test check == true
    @test isapprox(res.Tchar[1], 400.0; atol = 1.0)

    res_m1, Telu_max_m1 = RetentionParameterEstimator.method_m1(meas_select, col_input, se_col=false)
    @test res_m1.Tchar == res.Tchar
    @test Telu_max_m1 == Telu_max

    res_m2, Telu_max_m2 = RetentionParameterEstimator.method_m2(meas, se_col=true)
    @test isapprox(res_m2.d[1], 0.00025, atol=1e-5)  

    res_m3, Telu_missing_m3 = RetentionParameterEstimator.method_m3(meas)
    @test isapprox(res_m3.d[1], res_m2.d[1], atol=1e-5)
end

@testset "Optimization with Missing Values" begin
    file = "./data/meas_test_missing.csv" # from some tR are missing
    meas_missing = RetentionParameterEstimator.load_chromatograms(file)
    col_input = (L = meas_missing[1].L, d = meas_missing[1].d*1000)
    
    # Test parameter estimation with partially missing data
    res_missing_m1 = RetentionParameterEstimator.method_m1(meas_missing, col_input, se_col=false)[1]
    @test !ismissing(res_missing_m1.Tchar[1])
    @test !isnothing(res_missing_m1.θchar[1])

    res_missing_m2 = RetentionParameterEstimator.method_m2(meas_missing, se_col=false)[1]
    @test isapprox(res_missing_m2.Tchar[2], res_missing_m1.Tchar[2], atol=1.0)

    res_missing_m3 = RetentionParameterEstimator.method_m3(meas_missing)[1]
    @test isapprox(res_missing_m3.d[1], res_missing_m2.d[1], atol=1e-5)
    
    # Test with all measurements missing for one solute
    # TODO: this case is not covered yet
    #meas_all_missing = deepcopy(meas_missing)
    #meas_all_missing[3][:,2] .= missing
    #res_all_missing = RetentionParameterEstimator.method_m1(meas_all_missing, col_input, se_col=false)
    #@test length(res_all_missing.Tchar) == length(meas_all_missing[4]) - 1  # Should skip fully missing solute
end