# functions for loading the measured retention data and the necessary informations (column definition, programs, ...)

function load_chromatograms_(file)
    n = open(f->countlines(f), file)
    col_df = DataFrame(CSV.File(file, header=1, limit=1, stringtype=String))
    column = Dict(   :L => col_df.L[1],
                        :d => col_df.d[1],
                        :df => col_df.df[1],
                        :gas => col_df.gas[1],
                        :pout => col_df.pout[1],
                        :sp => col_df.sp[1],
                        :time_unit => col_df.time_unit[1]
                    )
    n_meas = Int((n - 2 - 3)/3) 
    TP = DataFrame(CSV.File(file, header=3, limit=n_meas, stringtype=String))
    PP = DataFrame(CSV.File(file, header=3+n_meas+1, limit=n_meas, stringtype=String)) # convert pressures from Pa(g) to Pa(a), add p_atm to this data set
    tRs = DataFrame(CSV.File(file, header=n-n_meas, stringtype=String))
    solute_names = names(tRs)[2:end] # filter non-solute names out (columnx)
    filter!(x -> !occursin.("Column", x), solute_names)
    return column, TP, PP, tRs, solute_names
end

function extract_temperature_and_pressure_programs(TPprog)
    iP = findall(occursin.("p", names(TPprog))) # column index of pressure information
    #iT = 2:(iP[1]-1) # column index of temperature program information
    TPs = TPprog[!,1:(iP[1]-1)] 
    PPs = DataFrame(measurement = TPs.measurement, p1 = TPprog[!,iP[1]].+TPprog[!,iP[end]], t1 = TPs.t1)
    iTrate = findall(occursin.("RT", names(TPs))) # column index of heating rates 
    heatingrates = Array{Array{Float64}}(undef, length(iTrate))
    Tdiff = Array{Array{Float64}}(undef, length(iTrate))
    for i=1:length(iTrate)
        heatingrates[i] = TPs[!, iTrate[i]]
        Tdiff[i] = TPs[!, iTrate[i]+1] .- TPs[!, iTrate[i]-2]
    end
    pressurerates = Array{Array{Float64}}(undef, length(iTrate))
    for i=1:length(iTrate)
        pressurerates[i] = (TPprog[!,iP[i+1]] .- TPprog[!,iP[i]]) .* heatingrates[i] ./ Tdiff[i]
        PPs[!,"RP$(i)"] = pressurerates[i]
        PPs[!,"p$(i+1)"] = TPprog[!, iP[i+1]].+TPprog[!, iP[end]]
        PPs[!,"t$(i+1)"] = TPs[!,"t$(i+1)"]
    end 
    pamb = TPprog[!, iP[end]]
    return TPs, PPs, pamb
end

function load_chromatograms__(file)
    n = open(f->countlines(f), file)
    col_df = DataFrame(CSV.File(file, header=1, limit=1, stringtype=String))
    column = Dict(   :L => col_df.L[1],
                        :d => col_df.d[1],
                        :df => col_df.df[1],
                        :gas => col_df.gas[1],
                        :pout => col_df.pout[1],
                        :sp => col_df.sp[1],
                        :time_unit => col_df.time_unit[1]
                    )
    n_meas = Int((n - 2 - 2)/2) 
    TPprog = DataFrame(CSV.File(file, header=3, limit=n_meas, stringtype=String))
    #PP = DataFrame(CSV.File(file, header=3+n_meas+1, limit=n_meas, stringtype=String)) # convert pressures from Pa(g) to Pa(a), add p_atm to this data set
    tRs = DataFrame(CSV.File(file, header=n-n_meas, stringtype=String))
    solute_names = names(tRs)[2:end] # filter non-solute names out (columnx)
    filter!(x -> !occursin.("Column", x), solute_names)
    TPs, PPs, pamb = extract_temperature_and_pressure_programs(TPprog)
    return column, TPs, PPs, tRs, solute_names, pamb
end

function load_chromatograms(file; version="with ambient pressure")
    if version == "with ambient pressure" # standard version, pressure program listed together with temperature program
        column, TPs, PPs, tRs, solute_names, pamb = load_chromatograms__(file)
        return column, TPs, PPs, tRs, solute_names, pamb
    elseif version == "with pressure program" # old version, where pressure program was listed separatly
        column, TP, PP, tRs, solute_names = load_chromatograms_(file)
        return column, TP, PP, tRs, solute_names
    else
        column, TPs, PPs, tRs, solute_names, pamb = load_chromatograms__(file)
        return column, TPs, PPs, tRs, solute_names, pamb
    end
end