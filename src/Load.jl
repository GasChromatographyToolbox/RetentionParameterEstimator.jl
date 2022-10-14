# functions for loading the measured retention data and the necessary informations (column definition, programs, ...)

function load_chromatograms(file)
    n = open(f->countlines(f), file)
    col_df = DataFrame(CSV.File(file, header=1, limit=1))
    column = Dict(   :L => col_df.L[1],
                        :d => col_df.d[1],
                        :df => col_df.df[1],
                        :gas => col_df.gas[1],
                        :pout => col_df.pout[1],
                        :sp => col_df.sp[1],
                        :time_unit => col_df.time_unit[1]
                    )
    n_meas = Int((n - 2 - 3)/3) 
    TP = DataFrame(CSV.File(file, header=3, limit=n_meas))
    PP = DataFrame(CSV.File(file, header=3+n_meas+1, limit=n_meas))
    tRs = DataFrame(CSV.File(file, header=n-n_meas))
    solute_names = names(tRs)[2:end] # filter non-solute names out (columnx)
    filter!(x -> !occursin.("Column", x), solute_names)
    return column, TP, PP, tRs, solute_names
end