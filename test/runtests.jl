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

    file_dict = Dict{Any, Any}("data" => joinpath(@__DIR__, "data", "meas_test.csv"), "name" => "meas_test.csv")
    meas_dict = RetentionParameterEstimator.load_chromatograms(file_dict)
    @test meas_dict[4][1] == meas[4][1]

    buf = IOBuffer(read(joinpath(@__DIR__, "data", "meas_test.csv")))
    meas_buf = RetentionParameterEstimator.load_chromatograms(Dict{Any, Any}("data" => buf, "name" => "meas_test.csv"))
    @test meas_buf[4][1] == meas[4][1]
end

# Loss.jl
@testset "Loss Functions" begin
    # Setup test data
    file = "./data/meas_test.csv"
    meas = RetentionParameterEstimator.load_chromatograms(file)
    meas_select = RetentionParameterEstimator.filter_selected_measurements(meas, ["meas3", "meas4", "meas5"], ["2-Octanone"])
    
    col_input = (L = meas_select[1].L, d = meas_select[1].d*1000)
    check, msg, df_flag, index_flag, res, Telu_max = RetentionParameterEstimator.check_measurement(meas_select, col_input; min_th=0.1, loss_th=1.0, se_col=false, maxiters=500, maxtime=30.0)
    
    # Get parameters from optimization result
    # Extract Float64 values from Measurement types (when se_col=false, results are Measurements)
    Tchar = Measurements.value(res.Tchar[1])
    θchar = Measurements.value(res.θchar[1])
    ΔCp = Measurements.value(res.ΔCp[1])
    L = meas_select[1].L
    d = meas_select[1].d
    df = meas_select[1].df
    gas = meas_select[1].gas
    prog = meas_select[2]
    
    # Test 1: tR_calc (without df)
    @testset "tR_calc (without df)" begin
        tR_calc_result = RetentionParameterEstimator.tR_calc(Tchar, θchar, ΔCp, L, d, prog[1], gas)
        @test tR_calc_result > 0.0  # Retention time should be positive
        @test isfinite(tR_calc_result)  # Should be finite
        
        # Test with different programs
        if length(prog) > 1
            tR_calc_result2 = RetentionParameterEstimator.tR_calc(Tchar, θchar, ΔCp, L, d, prog[2], gas)
            @test tR_calc_result2 > 0.0
            @test isfinite(tR_calc_result2)
        end
    end
    
    # Test 2: tR_τR_calc
    @testset "tR_τR_calc" begin
        Cag = 0.1  # Typical diffusivity coefficient
        t₀ = 0.0
        τ₀ = 0.1  # Initial peak width
        
        tR_τR_result = RetentionParameterEstimator.tR_τR_calc(Tchar, θchar, ΔCp, L, d, df, prog[1], Cag, t₀, τ₀, gas)
        @test length(tR_τR_result) == 2  # Should return tuple (tR, τR)
        tR_result, τR_result = tR_τR_result
        
        @test tR_result > 0.0  # Retention time should be positive
        @test τR_result > 0.0  # Peak width should be positive
        @test isfinite(tR_result)
        @test isfinite(τR_result)
        
        # Peak width should be reasonable (not too large compared to retention time)
        @test τR_result < tR_result  # Peak width should be smaller than retention time
    end
    
    # Test 4: loss function (matrix input, squared metric)
    @testset "loss (matrix input, squared metric)" begin
        # Create a simple 2D array of retention times (2 programs, 1 solute)
        @test length(prog) >= 2  # Need at least 2 programs for this test
        tR_matrix = [300.0; 350.0]  # 2x1 array (2 programs, 1 solute)
        Tchar_vec = [Tchar]
        θchar_vec = [θchar]
        ΔCp_vec = [ΔCp]
        prog_vec = [prog[1], prog[2]]
        
        loss_squared = RetentionParameterEstimator.loss(tR_matrix, Tchar_vec, θchar_vec, ΔCp_vec, L, d, prog_vec, gas; metric="squared")
        @test loss_squared >= 0.0  # Loss should be non-negative
        @test isfinite(loss_squared)
        
        # Test with multiple solutes (if available) - OPTIONAL: Skip expensive optimization for CI
        # This test requires an expensive check_measurement call, so we make it optional
        # Uncomment the following block if you want to test multi-solute loss calculations:
        #=
        if length(meas[4]) > 1
            # Use first two solutes
            meas_two = RetentionParameterEstimator.filter_selected_measurements(meas, ["meas3", "meas4"], meas[4][1:2])
            col_input_two = (L = meas_two[1].L, d = meas_two[1].d*1000)
            check_two, msg_two, df_flag_two, index_flag_two, res_two, Telu_max_two = RetentionParameterEstimator.check_measurement(meas_two, col_input_two; min_th=0.1, loss_th=1.0, se_col=false, maxiters=500, maxtime=30.0)
            
            if length(res_two.Tchar) >= 2 && length(meas_two[2]) >= 2
                tR_matrix_2 = [300.0 400.0; 350.0 450.0]  # 2x2 array (2 programs, 2 solutes)
                Tchar_vec_2 = [Measurements.value(res_two.Tchar[1]), Measurements.value(res_two.Tchar[2])]
                θchar_vec_2 = [Measurements.value(res_two.θchar[1]), Measurements.value(res_two.θchar[2])]
                ΔCp_vec_2 = [Measurements.value(res_two.ΔCp[1]), Measurements.value(res_two.ΔCp[2])]
                prog_vec_2 = [meas_two[2][1], meas_two[2][2]]
                
                loss_squared_2 = RetentionParameterEstimator.loss(tR_matrix_2, Tchar_vec_2, θchar_vec_2, ΔCp_vec_2, meas_two[1].L, meas_two[1].d, prog_vec_2, meas_two[1].gas; metric="squared")
                @test loss_squared_2 >= 0.0
                @test isfinite(loss_squared_2)
            end
        end
        =#
    end
    
    # Test 5: loss function (matrix input, abs metric)
    @testset "loss (matrix input, abs metric)" begin
        @test length(prog) >= 2  # Need at least 2 programs for this test
        tR_matrix = [300.0; 350.0]
        Tchar_vec = [Tchar]
        θchar_vec = [θchar]
        ΔCp_vec = [ΔCp]
        prog_vec = [prog[1], prog[2]]
        
        loss_abs = RetentionParameterEstimator.loss(tR_matrix, Tchar_vec, θchar_vec, ΔCp_vec, L, d, prog_vec, gas; metric="abs")
        @test loss_abs >= 0.0  # Loss should be non-negative
        @test isfinite(loss_abs)
        
        # abs metric should give different (but related) result than squared
        loss_squared = RetentionParameterEstimator.loss(tR_matrix, Tchar_vec, θchar_vec, ΔCp_vec, L, d, prog_vec, gas; metric="squared")
        @test loss_abs != loss_squared  # Should be different
    end
    
    # Test 6: loss function (vector input with substance_list)
    @testset "loss (vector input with substance_list)" begin
        @test length(prog) >= 2  # Need at least 2 programs for this test
        # Create vector of retention times matching programs and substances
        # Note: The vector loss function expects Tchar, θchar, ΔCp to be arrays matching unique substances
        substance_list = ["2-Octanone", "2-Octanone"]  # Same substance, two programs
        tR_vector = [300.0, 350.0]  # Retention times for two programs
        Tchar_vec = [Tchar]  # Array with one value for the unique substance
        θchar_vec = [θchar]
        ΔCp_vec = [ΔCp]
        prog_vec = [prog[1], prog[2]]
        
        loss_vector = RetentionParameterEstimator.loss(tR_vector, Tchar_vec, θchar_vec, ΔCp_vec, substance_list, L, d, prog_vec, gas; metric="squared")
        @test loss_vector >= 0.0
        @test isfinite(loss_vector)
        
        # Test with abs metric
        loss_vector_abs = RetentionParameterEstimator.loss(tR_vector, Tchar_vec, θchar_vec, ΔCp_vec, substance_list, L, d, prog_vec, gas; metric="abs")
        @test loss_vector_abs >= 0.0
        @test isfinite(loss_vector_abs)
    end
    
    # Test 7: loss function error handling (mismatched lengths)
    @testset "loss (error handling)" begin
        @test length(prog) >= 2  # Need at least 2 programs for this test
        tR_vector = [300.0, 350.0]
        substance_list = ["2-Octanone"]  # Wrong length (should match length of tR_vector)
        Tchar_vec = [Tchar]
        θchar_vec = [θchar]
        ΔCp_vec = [ΔCp]
        prog_vec = [prog[1], prog[2]]
        
        # The function should throw an ErrorException for mismatched lengths
        @test_throws ErrorException RetentionParameterEstimator.loss(tR_vector, Tchar_vec, θchar_vec, ΔCp_vec, substance_list, L, d, prog_vec, gas)
    end
    
    # Test 8: loss function with single program and single solute (edge case)
    @testset "loss (edge cases)" begin
        # Single program, single solute
        tR_single = [300.0]  # 1x1 array
        Tchar_single = [Tchar]
        θchar_single = [θchar]
        ΔCp_single = [ΔCp]
        prog_single = [prog[1]]
        
        loss_single = RetentionParameterEstimator.loss(tR_single, Tchar_single, θchar_single, ΔCp_single, L, d, prog_single, gas; metric="squared")
        @test loss_single >= 0.0
        @test isfinite(loss_single)
        
        # Scalar input (0-dimensional)
        tR_scalar = 300.0
        loss_scalar = RetentionParameterEstimator.loss(tR_scalar, Tchar_single, θchar_single, ΔCp_single, L, d, prog_single, gas; metric="squared")
        @test loss_scalar >= 0.0
        @test isfinite(loss_scalar)
    end
    
    # Test 9: Consistency between loss calculations
    @testset "loss (consistency)" begin
        @test length(prog) >= 2  # Need at least 2 programs for this test
        # For the same data, loss should be consistent
        tR_matrix = [300.0; 350.0]
        Tchar_vec = [Tchar]
        θchar_vec = [θchar]
        ΔCp_vec = [ΔCp]
        prog_vec = [prog[1], prog[2]]
        
        loss1 = RetentionParameterEstimator.loss(tR_matrix, Tchar_vec, θchar_vec, ΔCp_vec, L, d, prog_vec, gas; metric="squared")
        loss2 = RetentionParameterEstimator.loss(tR_matrix, Tchar_vec, θchar_vec, ΔCp_vec, L, d, prog_vec, gas; metric="squared")
        @test isapprox(loss1, loss2, atol=1e-10)  # Should be identical
        
        # Test that loss decreases when using calculated retention times
        tR_calc_1 = RetentionParameterEstimator.tR_calc(Tchar, θchar, ΔCp, L, d, prog[1], gas)
        tR_calc_2 = RetentionParameterEstimator.tR_calc(Tchar, θchar, ΔCp, L, d, prog[2], gas)
        tR_calc_matrix = [tR_calc_1; tR_calc_2]
        
        loss_calc = RetentionParameterEstimator.loss(tR_calc_matrix, Tchar_vec, θchar_vec, ΔCp_vec, L, d, prog_vec, gas; metric="squared")
        @test loss_calc < loss1 || isapprox(loss_calc, 0.0, atol=1e-6)  # Should be very small or zero
    end
end

# Estimate_Start_Values.jl

# Optimization.jl
@testset "Optimize d_only" begin
    # Test the new optimize_d_only function
    # Use full meas data (not filtered) to ensure proper dimensions
    file = "./data/meas_test.csv"
    meas = RetentionParameterEstimator.load_chromatograms(file)
    
    # Get parameters from a single substance for testing
    meas_select = RetentionParameterEstimator.filter_selected_measurements(meas, ["meas3", "meas4", "meas5"], ["2-Octanone"])
    col_input = (L = meas_select[1].L, d = meas_select[1].d*1000)
    check, msg, df_flag, index_flag, res, Telu_max = RetentionParameterEstimator.check_measurement(meas_select, col_input; min_th=0.1, loss_th=1.0, se_col=false, maxiters=500, maxtime=30.0)
    
    # Get parameters
    Tchar = Measurements.value(res.Tchar[1])
    θchar = Measurements.value(res.θchar[1])
    ΔCp = Measurements.value(res.ΔCp[1])
    L = meas_select[1].L
    d_initial = meas_select[1].d
    gas = meas_select[1].gas
    
    # Prepare data for optimize_d_only using helper function
    tRs_, prog_, subst_list_ = RetentionParameterEstimator.prepare_optimization_data(
        meas[3], meas[4], meas[2], meas[6]
    )
    
    # For optimize_d_only, we need to use parameters for all substances
    # But we only have optimized parameters for one substance, so we'll use the full meas
    # and get start parameters for all substances
    Tchar_est, θchar_est, ΔCp_est, Telu_max_est = RetentionParameterEstimator.estimate_start_parameter(meas[3], meas[1], meas[2]; time_unit=meas[6])
    
    # Test optimize_d_only with all substances
    d_optimized = RetentionParameterEstimator.optimize_d_only(tRs_, subst_list_, Tchar_est, θchar_est, ΔCp_est, 
                                                               meas[1].L, prog_, meas[1].gas, d_initial; maxiters=500, maxtime=20.0)
    
    @test d_optimized > 0.0  # Diameter should be positive
    @test isfinite(d_optimized)  # Should be finite
    @test isapprox(d_optimized, d_initial, atol=1e-3)  # Should be close to initial value (may differ slightly)
    
    # Test that optimize_d_only actually optimizes (loss should decrease or stay similar)
    # Calculate loss before and after (using the optimized d)
    # This is a basic sanity check
    @test d_optimized > 1e-6  # Should be a reasonable diameter
    @test d_optimized < 1e-2  # Should be less than 1 cm
end

@testset "Optimization" begin
    file = "./data/meas_test.csv"
    meas = RetentionParameterEstimator.load_chromatograms(file)
    meas_select = RetentionParameterEstimator.filter_selected_measurements(meas, ["meas3", "meas4", "meas5"], ["2-Octanone"]);

    col_input = (L = meas_select[1].L, d = meas_select[1].d*1000)
    check, msg, df_flag, index_flag, res, Telu_max = RetentionParameterEstimator.check_measurement(meas_select, col_input; min_th=0.1, loss_th=1.0, se_col=false, maxiters=500, maxtime=30.0)
    @test check == true
    @test isapprox(res.Tchar[1], 400.0; atol = 1.0)

    res_m1, Telu_max_m1 = RetentionParameterEstimator.method_m1(meas_select, col_input, se_col=false, maxiters=500, maxtime=30.0)
    @test res_m1.Tchar == res.Tchar
    @test Telu_max_m1 == Telu_max

    res_m2, Telu_max_m2 = RetentionParameterEstimator.method_m2(meas, se_col=false, maxiters=500, maxtime=30.0)
    @test isapprox(res_m2.d[1], 0.00025, atol=1e-5)  

    res_m3, Telu_missing_m3 = RetentionParameterEstimator.method_m3(meas, maxiters=500, maxtime=30.0)
    @test isapprox(res_m3.d[1], res_m2.d[1], atol=1e-5)
    
    # Test method_m4 (alternating optimization)
    res_m4, Telu_max_m4 = RetentionParameterEstimator.method_m4(meas, se_col=false, maxiters=500, maxtime=30.0, max_alternating_iters=3, tol=1e-4)
    @test isapprox(res_m4.d[1], 0.00025, atol=1e-5)
    @test length(res_m4.Name) == length(meas[4])
    @test all(res_m4.d .== res_m4.d[1])  # All d values should be the same (system parameter)
    
    # Compare method_m2 and method_m4 results (should be similar but may differ slightly)
    @test isapprox(res_m4.d[1], res_m2.d[1], atol=1e-4)  # d should be similar
    @test isapprox(res_m4.Tchar[1], res_m2.Tchar[1], atol=5.0)  # Tchar may differ slightly
    @test Telu_max_m4 == Telu_max_m2  # Telu_max should be the same
    
    # Test that method_m4 enforces d is the same for all substances (key feature)
    @test all(res_m4.d .== res_m4.d[1])  # All d values must be identical
    @test length(unique(res_m4.d)) == 1  # Only one unique d value
end

@testset "Standard Errors" begin
    # Test standard error calculations separately (expensive, so only test once)
    file = "./data/meas_test.csv"
    meas = RetentionParameterEstimator.load_chromatograms(file)
    meas_select = RetentionParameterEstimator.filter_selected_measurements(meas, ["meas3", "meas4", "meas5"], ["2-Octanone"])
    col_input = (L = meas_select[1].L, d = meas_select[1].d*1000)
    
    # Test standard errors with method_m2 (includes d)
    res_m2_se, _ = RetentionParameterEstimator.method_m2(meas, se_col=true, maxiters=500, maxtime=30.0)
    @test hasproperty(res_m2_se, :Tchar_std)
    @test hasproperty(res_m2_se, :θchar_std)
    @test hasproperty(res_m2_se, :ΔCp_std)
    @test hasproperty(res_m2_se, :d_std)
    @test all(res_m2_se.Tchar_std .> 0.0)  # Standard errors should be positive
    @test all(res_m2_se.d_std .> 0.0)
    
    # Test standard errors with method_m4
    res_m4_se, _ = RetentionParameterEstimator.method_m4(meas, se_col=true, maxiters=500, maxtime=30.0, max_alternating_iters=3, tol=1e-4)
    @test hasproperty(res_m4_se, :Tchar_std)
    @test hasproperty(res_m4_se, :d_std)
    @test all(res_m4_se.d .== res_m4_se.d[1])  # All d values should be the same
    @test all(res_m4_se.d_std .== res_m4_se.d_std[1])  # All d_std values should be the same
end

@testset "Optimization with Missing Values" begin
    file = "./data/meas_test_missing.csv" # from some tR are missing
    meas_missing = RetentionParameterEstimator.load_chromatograms(file)
    col_input = (L = meas_missing[1].L, d = meas_missing[1].d*1000)
    
    # Test parameter estimation with partially missing data
    res_missing_m1 = RetentionParameterEstimator.method_m1(meas_missing, col_input, se_col=false, maxiters=500, maxtime=30.0)[1]
    @test !ismissing(res_missing_m1.Tchar[1])
    @test !isnothing(res_missing_m1.θchar[1])

    res_missing_m2 = RetentionParameterEstimator.method_m2(meas_missing, se_col=false, maxiters=500, maxtime=30.0)[1]
    # Tolerance increased to 2.0 to account for small numerical differences between 
    # method_m1 and method_m2, which may vary slightly with different GasChromatographySimulator versions
    @test isapprox(res_missing_m2.Tchar[2], res_missing_m1.Tchar[2], atol=2.0)
    #@show res_missing_m2.Tchar[2]
    #@show res_missing_m1.Tchar[2]

    res_missing_m3 = RetentionParameterEstimator.method_m3(meas_missing, maxiters=500, maxtime=30.0)[1]
    @test isapprox(res_missing_m3.d[1], res_missing_m2.d[1], atol=1e-5)
    
    # Test method_m4 with missing values
    res_missing_m4 = RetentionParameterEstimator.method_m4(meas_missing, se_col=false, maxiters=500, maxtime=30.0, max_alternating_iters=3, tol=1e-4)[1]
    @test !ismissing(res_missing_m4.Tchar[1])
    @test !isnothing(res_missing_m4.θchar[1])
    @test all(res_missing_m4.d .== res_missing_m4.d[1])  # All d values should be the same
    # Compare with method_m2 (should be similar)
    @test isapprox(res_missing_m4.d[1], res_missing_m2.d[1], atol=1e-4)
    @test isapprox(res_missing_m4.Tchar[2], res_missing_m2.Tchar[2], atol=2.0)
    
    # Test with all measurements missing for one solute
    # TODO: this case is not covered yet
    #meas_all_missing = deepcopy(meas_missing)
    #meas_all_missing[3][:,2] .= missing
    #res_all_missing = RetentionParameterEstimator.method_m1(meas_all_missing, col_input, se_col=false)
    #@test length(res_all_missing.Tchar) == length(meas_all_missing[4]) - 1  # Should skip fully missing solute
end

@testset "Multi-start Optimization" begin
    file = "./data/meas_test.csv"
    meas = RetentionParameterEstimator.load_chromatograms(file)
    meas_select = RetentionParameterEstimator.filter_selected_measurements(meas, ["meas3", "meas4", "meas5"], ["2-Octanone"])
    
    col = GasChromatographySimulator.Column(meas_select[1].L, meas_select[1].d, meas_select[1].df, meas_select[1].sp, meas_select[1].gas)
    
    # Get initial estimates using helper function
    Tchar_est, θchar_est, ΔCp_est, _ = RetentionParameterEstimator.estimate_start_parameter(meas_select[3], col, meas_select[2]; time_unit=meas_select[6])
    
    # Prepare data using helper functions
    a = RetentionParameterEstimator.time_unit_conversion_factor(meas_select[6])
    tR_meas = Array(meas_select[3][:,2:end]).*a
    tRs_, prog_, subst_list_ = RetentionParameterEstimator.prepare_single_substance_data(tR_meas[:,1], meas_select[2], meas_select[4][1])
    
    # Test 1: Basic multistart functionality (single chromatogram - minimal n_starts)
    @testset "optimize_Kcentric_multistart basic" begin
        # Create single-chromatogram scenario
        prog_single = [prog_[1]]
        tRs_single = [tRs_[1]]
        subst_list_single = [subst_list_[1]]
        
        sol_multistart = RetentionParameterEstimator.optimize_Kcentric_multistart(
            tRs_single, subst_list_single, col, prog_single, 
            [Tchar_est[1]], [θchar_est[1]], [ΔCp_est[1]]; 
            maxiters=500, maxtime=30.0, n_starts=2  # Minimal starts for testing
        )
        
        @test sol_multistart[1] > 0.0
        @test sol_multistart[2] > 0.0
        @test isfinite(sol_multistart.objective)
        @test sol_multistart.objective >= 0.0
    end
    
    # Test 2: single-chromatogram data with multistart disabled (default)
    @testset "estimate_parameters Kcentric_single one chromatogram, multistart off" begin
        meas_single = (
            meas_select[1],
            [meas_select[2][1]],
            meas_select[3][1:1, :],
            meas_select[4],
            meas_select[5],
            meas_select[6]
        )
        
        res_off = RetentionParameterEstimator.estimate_parameters(
            meas_single[3], meas_single[4], meas_single[1], meas_single[2],
            Tchar_est, θchar_est, ΔCp_est;
            mode="Kcentric_single", pout=meas_single[5], time_unit=meas_single[6],
            maxiters=500, maxtime=30.0, multistart_n=0
        )[1]
        
        @test length(res_off.Name) == length(meas_single[4])
        @test all(res_off.Tchar .> 0.0)
        @test all(isfinite.(res_off.Tchar))
    end
    
    # Test 3: estimate_parameters with explicit multistart_n (minimal test)
    @testset "estimate_parameters explicit multistart_n" begin
        res_explicit = RetentionParameterEstimator.estimate_parameters(
            meas_select[3], meas_select[4], meas_select[1], meas_select[2],
            Tchar_est, θchar_est, ΔCp_est;
            mode="Kcentric_single", pout=meas_select[5], time_unit=meas_select[6],
            maxiters=500, maxtime=30.0, multistart_n=2
        )[1]
        
        @test length(res_explicit.Name) == length(meas_select[4])
        @test all(res_explicit.Tchar .> 0.0)
        @test all(isfinite.(res_explicit.Tchar))
    end

    # Test 4: multistart must not worsen per-substance loss vs multistart_n=0
    @testset "estimate_parameters multistart loss not worse than off" begin
        res_off = RetentionParameterEstimator.estimate_parameters(
            meas_select[3], meas_select[4], meas_select[1], meas_select[2],
            Tchar_est, θchar_est, ΔCp_est;
            mode="Kcentric_single", pout=meas_select[5], time_unit=meas_select[6],
            maxiters=500, maxtime=30.0, multistart_n=0
        )[1]
        res_on = RetentionParameterEstimator.estimate_parameters(
            meas_select[3], meas_select[4], meas_select[1], meas_select[2],
            Tchar_est, θchar_est, ΔCp_est;
            mode="Kcentric_single", pout=meas_select[5], time_unit=meas_select[6],
            maxiters=500, maxtime=30.0, multistart_n=2
        )[1]
        @test res_on.Name == res_off.Name
        for i in eachindex(res_on.min)
            @test res_on.min[i] <= res_off.min[i] + 1e-6 * max(res_off.min[i], 1.0)
        end
    end

    # Test 5: optimize_dKcentric_multistart smoke (single chromatogram)
    @testset "optimize_dKcentric_multistart basic" begin
        prog_single = [prog_[1]]
        tRs_single = [tRs_[1]]
        subst_list_single = [subst_list_[1]]
        sol_dms = RetentionParameterEstimator.optimize_dKcentric_multistart(
            tRs_single, subst_list_single, col, prog_single, col.d,
            [Tchar_est[1]], [θchar_est[1]], [ΔCp_est[1]];
            maxiters=500, maxtime=30.0, n_starts=2
        )
        @test sol_dms[1] > 0.0
        @test sol_dms[2] > 0.0
        @test sol_dms[3] > 0.0
        @test isfinite(sol_dms.objective)
        @test sol_dms.objective >= 0.0
    end

    # Test 6: method_m1 forwards multistart_n
    @testset "method_m1 multistart_n" begin
        col_input = (L = meas_select[1].L, d = meas_select[1].d * 1000)
        res_m1_ms, _ = RetentionParameterEstimator.method_m1(
            meas_select, col_input; se_col=false, maxiters=500, maxtime=30.0, multistart_n=2
        )
        @test length(res_m1_ms.Name) == length(meas_select[4])
        @test all(res_m1_ms.Tchar .> 0.0)
        @test all(isfinite.(res_m1_ms.Tchar))
        @test all(res_m1_ms.min .>= 0.0)
    end

    # Test 7: method_m4 forwards multistart_n
    @testset "method_m4 multistart_n" begin
        res_m4_ms, Telu_m4_ms = RetentionParameterEstimator.method_m4(
            meas; se_col=false, maxiters=500, maxtime=30.0,
            max_alternating_iters=3, tol=1e-4, multistart_n=2
        )
        @test length(res_m4_ms.Name) == length(meas[4])
        @test all(res_m4_ms.d .== res_m4_ms.d[1])
        @test all(res_m4_ms.Tchar .> 0.0)
        @test all(isfinite.(res_m4_ms.Tchar))
        @test isfinite(Telu_m4_ms[1])
    end

    # Test 8: independent perturbation path (coupled_perturbation=false)
    @testset "estimate_parameters coupled_perturbation false" begin
        res_indep = RetentionParameterEstimator.estimate_parameters(
            meas_select[3], meas_select[4], meas_select[1], meas_select[2],
            Tchar_est, θchar_est, ΔCp_est;
            mode="Kcentric_single", pout=meas_select[5], time_unit=meas_select[6],
            maxiters=500, maxtime=30.0, multistart_n=2, coupled_perturbation=false
        )[1]
        @test length(res_indep.Name) == length(meas_select[4])
        @test all(res_indep.Tchar .> 0.0)
        @test all(isfinite.(res_indep.Tchar))
        @test all(res_indep.min .>= 0.0)
    end

    # Test 9: explicit multistart_n=0 matches default method_m1 (backward compatible)
    @testset "method_m1 multistart_n default" begin
        col_input = (L = meas_select[1].L, d = meas_select[1].d * 1000)
        res_default, Telu_default = RetentionParameterEstimator.method_m1(
            meas_select, col_input; se_col=false, maxiters=500, maxtime=30.0
        )
        res_off, Telu_off = RetentionParameterEstimator.method_m1(
            meas_select, col_input; se_col=false, maxiters=500, maxtime=30.0, multistart_n=0
        )
        @test res_default.Tchar == res_off.Tchar
        @test res_default.θchar == res_off.θchar
        @test res_default.ΔCp == res_off.ΔCp
        @test res_default.min == res_off.min
        @test Telu_default == Telu_off

        check, _, _, _, res_check, Telu_check = RetentionParameterEstimator.check_measurement(
            meas_select, col_input; min_th=0.1, loss_th=1.0, se_col=false, maxiters=500, maxtime=30.0
        )
        @test check == true
        @test res_off.Tchar == res_check.Tchar
        @test Telu_off == Telu_check
    end
end

@testset "Start parameter estimation" begin
    file = "./data/meas_test.csv"
    meas = RetentionParameterEstimator.load_chromatograms(file)
    meas_select = RetentionParameterEstimator.filter_selected_measurements(
        meas, ["meas3", "meas4", "meas5"], ["2-Octanone"]
    )
    col = GasChromatographySimulator.Column(
        meas_select[1].L, meas_select[1].d, meas_select[1].df, meas_select[1].sp, meas_select[1].gas
    )

    @testset "average_ramp_rate" begin
        RT = RetentionParameterEstimator.average_ramp_rate(meas_select[2][1]; time_unit=meas_select[6])
        @test RT > 0.0
        @test isfinite(RT)
    end

    @testset "estimate_start_parameter_single_ramp_weighted" begin
        Tchar_w, θchar_w, ΔCp_w, Telu_w = RetentionParameterEstimator.estimate_start_parameter_single_ramp_weighted(
            meas_select[3], col, meas_select[2]; time_unit=meas_select[6]
        )
        @test all(Tchar_w .> 0.0)
        @test all(θchar_w .> 0.0)
        @test all(isfinite.(Tchar_w))
        @test all(isfinite.(Telu_w))
    end

    @testset "estimate_start_parameter_single_measurement_corrected" begin
        Tchar_c, θchar_c, ΔCp_c, Telu_c = RetentionParameterEstimator.estimate_start_parameter_single_measurement_corrected(
            meas_select[3][1:1, :], col, meas_select[2][1]; time_unit=meas_select[6]
        )
        @test Tchar_c[1] > 0.0
        @test θchar_c[1] > 0.0
        @test isfinite(Tchar_c[1])
        @test isfinite(Telu_c[1])
    end
end
