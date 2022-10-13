# functions for loading the measured retention data and the necessary informations (column definition, programs, ...)

function load_chromatograms(file)
    n = open(f->countlines(f), file)
    column = DataFrame(CSV.File(file, header=1, limit=1))
    n_meas = Int((n - 2 - 3)/3) 
    TP = DataFrame(CSV.File(file, header=3, limit=n_meas))
    PP = DataFrame(CSV.File(file, header=3+n_meas+1, limit=n_meas))
    tRs = DataFrame(CSV.File(file, header=n-n_meas))
    solute_names = names(tRs)[2:end]
    return column, TP, PP, tRs, solute_names
end