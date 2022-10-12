# Functions used to simulate chromatograms to use as test values for the optimization

"""
conventional_GC(L, d, df, sp, gas, TP, PP, solutes, db_path, db_file)

Description.
"""
function conventional_GC(L, d, df, sp, gas, TP, PP, solutes, db_path, db_file; pout="vacuum", time_unit="min")
opt = GasChromatographySimulator.Options()
col = GasChromatographySimulator.Column(L, d, df, sp, gas)
prog = GasChromatographySimulator.Program(TP, PP, L; pout=pout, time_unit=time_unit)
sub = GasChromatographySimulator.load_solute_database(db_path, db_file, sp, gas, solutes, zeros(length(solutes)), zeros(length(solutes)))
par = GasChromatographySimulator.Parameters(col, prog, sub, opt)
return par
end

"""
    sim_test_chrom(L, d, df, sp, gas, TPs, PPs, solutes, db_path, db_file)

Description.
"""
function sim_test_chrom(L, d, df, sp, gas, TPs, PPs, solutes, db_path, db_file; pout="vacuum", time_unit="min")
    par_meas = Array{GasChromatographySimulator.Parameters}(undef, length(TPs))
    for i=1:length(TPs)
        par_meas[i] = conventional_GC(L, d, df, sp, gas, TPs[i], PPs[i], solutes, db_path, db_file; pout=pout, time_unit=time_unit)
    end

    #pl_meas = Array{DataFrame}(undef, length(TPs))
    tR_meas = Array{Float64}(undef, length(TPs), length(solutes))
    tR_meas_randn = Array{Float64}(undef, length(TPs), length(solutes))
    for i=1:length(TPs)
        pl_meas = GasChromatographySimulator.simulate(par_meas[i])[1]
        for j=1:length(solutes)
            jj = findfirst(pl_meas.Name.==solutes[j])
            tR_meas[i,j] = pl_meas.tR[jj]
            tR_meas_randn[i,j] = pl_meas.tR[jj]*(1+randn()*0.005) # add here a random ± time
        end
    end
    return tR_meas, tR_meas_randn, par_meas
end

function convert_vector_of_vector_to_2d_array(TPs)
    nTP = Array{Int}(undef, length(TPs))
    for i=1:length(TPs)
        nTP[i] = length(TPs[i])
    end
    TPm = Array{Union{Float64,Missing}}(undef, length(TPs), maximum(nTP))
    for i=1:length(TPs)
        for j=1:maximum(nTP)
            if j>length(TPs[i])
                TPm[i,j] = missing
            else
                TPm[i,j] = TPs[i][j]
            end
        end
    end
    return TPm
end

function program_header!(df_P, c)
    header = Array{Symbol}(undef, size(df_P)[2])
    header[1] = :measurement
    i1 = 1
    i2 = 1
    i3 = 1
    for i=2:size(df_P)[2]
        if i in collect(2:3:size(df_P)[2])
            header[i] = Symbol(string(c, "$(i1)"))
            i1 = i1 + 1
        elseif i in collect(3:3:size(df_P)[2])
            header[i] = Symbol("t$(i2)")
            i2 = i2 + 1
        else
            header[i] = Symbol(string("R", c, "$(i3)"))
            i3 = i3 + 1
        end
    end
    return rename!(df_P, header)
end

function save_simulation(file, column, measurement_name, TPs, PPs, solutes, tR_meas)
    # save the results:
    #
    # system information in header
    # L, d, df, sp, gas, pout, time_unit
    CSV.write(file, DataFrame(column))
    # temperature Program
    # measurment_name, TP
    df_TP = DataFrame([measurement_name convert_vector_of_vector_to_2d_array(TPs)], :auto)
    program_header!(df_TP, "T")
    CSV.write(file, df_TP, append=true, writeheader=true)
    # pressure Program
    # measurement_name, PP 
    df_PP = DataFrame([measurement_name convert_vector_of_vector_to_2d_array(PPs)], :auto)
    program_header!(df_PP, "p")
    CSV.write(file, df_PP, append=true, writeheader=true)
    # retention times
    # measurement_name, tR 
    df_tR = DataFrame(measurement=measurement_name)
    for i=1:length(solutes)
        df_tR[!, solutes[i]] = tR_meas[:,i]
    end
    CSV.write(file, df_tR, append=true, writeheader=true)
end