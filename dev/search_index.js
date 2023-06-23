var documenterSearchIndex = {"docs":
[{"location":"#RetentionParameterEstimator.jl","page":"Home","title":"RetentionParameterEstimator.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for RetentionParameterEstimator.jl","category":"page"},{"location":"","page":"Home","title":"Home","text":"RetentionParameterEstimator.jl","category":"page"},{"location":"docstrings/","page":"Docstrings","title":"Docstrings","text":"RetentionParameterEstimator","category":"page"},{"location":"docstrings/#Module-Index","page":"Docstrings","title":"Module Index","text":"","category":"section"},{"location":"docstrings/","page":"Docstrings","title":"Docstrings","text":"Modules = [RetentionParameterEstimator]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"docstrings/#Detailed-API","page":"Docstrings","title":"Detailed API","text":"","category":"section"},{"location":"docstrings/","page":"Docstrings","title":"Docstrings","text":"Modules = [RetentionParameterEstimator]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"docstrings/#RetentionParameterEstimator.Program-Tuple{Any, Any, Any}","page":"Docstrings","title":"RetentionParameterEstimator.Program","text":"Program(TP, FpinP, L; pout=\"vacuum\", time_unit=\"min\")\n\nConstruct the structure Program with conventional formulation (see conventional_program) of programs for the case without a thermal gradient. \n\nArguments\n\nTP: conventional formulation of a temperature program. \nFpinP: conventional formulation of a Flow (in m³/s) resp. inlet pressure (in Pa(a)) program.\nL: Length of the capillary measured in m (meter).\npout: Outlet pressure, \"vacuum\" (default), \"atmosphere\" or the outlet pressure in Pa(a).\ntime_unit: unit of time in the programs, \"min\"` (default) times are measured in minutes, \"s\" times are measured in seconds.\n\nThe argument L is used to construct the temperature interpolation T_itp(x,t).\n\nExamples\n\njulia> Program((40.0, 1.0, 5.0, 280.0, 2.0, 20.0, 320.0, 2.0),\n                (400000.0, 10.0, 5000.0, 500000.0, 20.0),\n                10.0)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.check_measurement-Tuple{Any, Any}","page":"Docstrings","title":"RetentionParameterEstimator.check_measurement","text":"check_measurement(meas, col_input; min_th=0.1, loss_th=1.0)\n\nSimilar to method_m1 (method_m1) estimate the three retention parameters T_char, θ_char and ΔC_p including standard errors, see stderror. In addition, if the found optimized minima is above a threshold min_th, it is flagged and the squared differences of single measured retention times and calculated retention times above another threshold loss_th are recorded.   \n\nArguments\n\nmeas ... Tuple with the loaded measurement data, see load_chromatograms.\ncol_input ... Named tuple with col_input.L the column length in m and col_input.d the column diameter in mm. If this parameter is not gicen, than these parameters are taken from meas. \n\nOutput\n\ncheck ... Boolean. true if all values are below the thresholds, false if not.\nmsg ... String. Description of check\ndf_flag ... Dataframe containing name of the flagged measurements, solutes and the corresponding mesured and calculated retention times.\nindex_flag ... Indices of flagged results\nres ... Dataframe with the optimized parameters and the found minima.\nTelu_max ... The maximum of elution temperatures every solute experiences in the measured programs.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.elution_temperature-Tuple{Any, Any}","page":"Docstrings","title":"RetentionParameterEstimator.elution_temperature","text":"elution_temperature(tRs, prog)\n\nCalculate the elution temperatures from retention times tRs and GC programs prog.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.estimate_parameters-Tuple{Any}","page":"Docstrings","title":"RetentionParameterEstimator.estimate_parameters","text":"estimate_parameters()\n\nCalculate the estimates for the K-centric parameters and (optional) the column diameter.\n\nArguments\n\nchrom ... Tuple of the loaded chromatogram, see load_chromatograms\n\nOptions\n\nmethod=NewtonTrustRegion() ... used optimization method\nopt=std_opt ... general options, std_opt = GasChromatographySimulator.Options(abstol=1e-8, reltol=1e-5, ng=true, odesys=false)\nmaxiters=10000 ... maximum number of iterations for every single optimization\nmaxtime=600.0 ... maximum time for every single optimization\nmode=\"dKcentric\" ... mode of the estimation.    Possible options: \n\"Kcentric_single\" ... optimization for the three K-centric retention parameters separatly for every solute\n\"Kcentric\" ... optimization for the three K-centric retention parameters together for all solutes\n\"dKcentric_single\" ... optimization for the column diameter and the three K-centric retention parameters separatly for every solute\n\"dKcentric\" ... optimization for the column diameter and the three K-centric retention parameters together for all solutes\nmetric=\"squared\" ... used metric for the loss function (\"squared\" or \"abs\")\n\nOutput\n\ndf ... DataFrame with the columns Name (solute names), d (estimated column diameter, optional), Tchar (estimated Tchar), θchar (estimated θchar), ΔCp (estimated ΔCp) and min (value of the loss function at the found optima)\nsol ... Array of SciMLBase.OptimizationSolution with the results of the optimization with some additional informations.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.estimate_start_parameter-Tuple{DataFrame, Any, Any}","page":"Docstrings","title":"RetentionParameterEstimator.estimate_start_parameter","text":"estimate_start_parameter(tRs::DataFrame, col, prog; time_unit=\"min\", control=\"Pressure\")\n\nEstimation of initial parameters for Tchar, θchar and ΔCp based on the elution temperatures calculated from the retention times tR and GC programs prog for column col. The initial value of Tchar is estimated from the elution temperatures of the measurements. Based on this estimated Tchar estimates for the initial values of θchar and ΔCp are calculated as     theta_charinit = 22 left(fracT_charinitT_stright)^07 left(1000fracd_fdright)^009 C and     Delta C_p = (-52 + 034 T_charinit) mathrmJ mol^-1 K^-1\n\nOutput\n\nTchar_est ... estimate for initial guess of the characteristic temperature\nθchar_est ... estimate for initial guess of θchar\nΔCp_est ... estimate for initial guess of ΔCp\nTelu_max ... the maximum of the calculated elution temperatures of the solutes\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.estimate_start_parameter_mean_elu_temp-Tuple{DataFrame, Any, Any}","page":"Docstrings","title":"RetentionParameterEstimator.estimate_start_parameter_mean_elu_temp","text":"estimate_start_parameter_mean_elu_temp(tRs::DataFrame, col, prog; time_unit=\"min\", control=\"Pressure\")\n\nEstimation of initial parameters for Tchar, θchar and ΔCp based on the elution temperatures calculated from the retention times tR and GC programs prog for column col. This function is used, if the temperature program is not a single ramp heating program. The elution temperatures of all measurements are calculated and the mean value of the elution temperatures is used as the initial characteristic temperature of a substance. Based on this estimated Tchar estimates for the initial values of θchar and ΔCp are calculated as     theta_charinit = 22 left(fracT_charinitT_stright)^07 left(1000fracd_fdright)^009 C and     Delta C_p = (-52 + 034 T_charinit) mathrmJ mol^-1 K^-1\n\nOutput\n\nTchar_est ... estimate for initial guess of the characteristic temperature\nθchar_est ... estimate for initial guess of θchar\nΔCp_est ... estimate for initial guess of ΔCp\nTelu_max ... the maximum of the calculated elution temperatures of the solutes\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.estimate_start_parameter_single_ramp-Tuple{DataFrame, Any, Any}","page":"Docstrings","title":"RetentionParameterEstimator.estimate_start_parameter_single_ramp","text":"estimate_start_parameter_single_ramp(tRs::DataFrame, col, prog; time_unit=\"min\", control=\"Pressure\")\n\nEstimation of initial parameters for Tchar, θchar and ΔCp based on the elution temperatures calculated from the retention times tR and GC programs prog for column col. For this function it is assumed, that single ramp heating programs are used. The elution temperatures of all measurements are calculated and than interpolated over the heating rates. For a dimensionless heating rate of 0.6 the elution temperature and the characteristic temperature of a substance are nearly equal. Based on this estimated Tchar estimates for the initial values of θchar and ΔCp are calculated as     theta_charinit = 22 left(fracT_charinitT_stright)^07 left(1000fracd_fdright)^009 C and     Delta C_p = (-180 + 063 T_charinit) mathrmJ mol^-1 K^-1\n\nOutput\n\nTchar_est ... estimate for initial guess of the characteristic temperature\nθchar_est ... estimate for initial guess of θchar\nΔCp_est ... estimate for initial guess of ΔCp\nTelu_max ... the maximum of the calculated elution temperatures of the solutes\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.load_chromatograms-Tuple{Any}","page":"Docstrings","title":"RetentionParameterEstimator.load_chromatograms","text":"load_chromatograms(file; filter_missing=true)\n\nLoading of the chromatographic data (column information, GC program information, retention time information, see also \"Structure of input data\") from a file.\n\nArguments\n\nfile ... path to the file.\nfilter_missing ... option to ignore solutes, where some retention times are not given. (default = true).\n\nOutput\n\nA tuple of the following quantities:\n\ncol ... settings of the column as GasChromatographySimulator.Column structure.\nprog ... Array of the GC programs as GasChromatographySimulator.Program structure.\ntRs ... DataFrame of the retention times.\nsolute_names ... Vector of the solute names.\npout ... outlet pressure (detector pressure), \"vacuum\" or \"atmospheric\". \ntime_unit ... unit of time scale used in the retention times and GC programs, \"min\" or \"s\".\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.load_chromatograms-Tuple{Dict{Any, Any}}","page":"Docstrings","title":"RetentionParameterEstimator.load_chromatograms","text":"load_chromatograms(file::Dict{Any, Any}; filter_missing=true, path=joinpath(dirname(pwd()), \"data\", \"exp_pro\"))\n\nLoading of the chromatographic data (column information, GC program information, retention time information, see also \"Structure of input data\") from a file selected by the FilePicker in a Pluto notebook.\n\nArguments\n\nfile ... file dictionary from the FilePicker.\nfilter_missing ... option to ignore solutes, where some retention times are not given. (default = true).\npath ... if the temperature programs are defined by measured temperatures over time, define the path to these files.\n\nOutput\n\nA tuple of the following quantities:\n\ncol ... settings of the column as GasChromatographySimulator.Column structure.\nprog ... Array of the GC programs as GasChromatographySimulator.Program structure.\ntRs ... DataFrame of the retention times.\nsolute_names ... Vector of the solute names.\npout ... outlet pressure (detector pressure), \"vacuum\" or \"atmospheric\". \ntime_unit ... unit of time scale used in the retention times and GC programs, \"min\" or \"s\".\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.loss-NTuple{8, Any}","page":"Docstrings","title":"RetentionParameterEstimator.loss","text":"loss(tR, Tchar, θchar, ΔCp, L, d, prog, gas; opt=std_opt, metric=\"squared\")\n\nLoss function as sum of squares of the residuals between the measured and calculated retention times.\n\nArguments\n\ntR ... mxn-array of the measured retention times in seconds.\nTchar ... n-array of characteristic temperatures in K.\nθchar ... n-array of characteristic constants in °C.\nΔCp ... n-array of the change of adiabatic heat capacity in J mol^-1 K^-1.\nL ... number of the length of the column in m.\nd ... number of the diameters of the column in m.\nprog ... m-array of structure GasChromatographySimulator.Programs containing the definition of the GC-programs.\ngas ... string of name of the mobile phase gas. \n\nOutput\n\nThe output is a tuple of the following quantites:\n\nsum((tR.-tRcalc).^2) ... sum of the squared residuals over m GC-programs and n solutes.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.method_m1-Tuple{Any, Any}","page":"Docstrings","title":"RetentionParameterEstimator.method_m1","text":"method_m1(meas, col_input)\n\nEstimation of the three retention parameters T_char, θ_char and ΔC_p including standard errors, see stderror.\n\nArguments\n\nmeas ... Tuple with the loaded measurement data, see load_chromatograms.\ncol_input ... Named tuple with col_input.L the column length in m and col_input.d the column diameter in mm. If this parameter is not gicen, than these parameters are taken from meas. \n\nOutput\n\nres ... Dataframe with the optimized parameters and the found minima.\nTelu_max ... The maximum of elution temperatures every solute experiences in the measured programs.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.method_m2-Tuple{Any}","page":"Docstrings","title":"RetentionParameterEstimator.method_m2","text":"method_m2(meas)\n\nEstimation of the column diameter d and three retention parameters T_char, θ_char and Δ C_p including standard errors, see stderror. In a first run all four parameters are estimated for every substance separatly, resulting in different optimized column diameters. The mean value of the column diameter is used for  a second optimization using this mean diameter and optimize the remainig thre retention parameters T_char, θ_char and Δ C_p.\n\nArguments\n\nmeas ... Tuple with the loaded measurement data, see load_chromatograms.\ncol_input ... Named tuple with col_input.L the column length in m and col_input.d the column diameter in mm. If this parameter is not gicen, than these parameters are taken from meas. \n\nOutput\n\nres ... Dataframe with the optimized parameters and the found minima.\nTelu_max ... The maximum of elution temperatures every solute experiences in the measured programs.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.opt_Kcentric-Tuple{Any, Any}","page":"Docstrings","title":"RetentionParameterEstimator.opt_Kcentric","text":"opt_Kcentric(x, p)\n\nFunction used for optimization of the loss-function in regards to the three K-centric parameters.\n\nArguments\n\nx ... 3n-vector of the three K-centric parameters of n solutes. Elements 1:n are Tchar, n+1:2n are θchar and 2n+1:3n are ΔCp values.\np ... vector containing the fixed parameters:\ntR = p[1] ... mxn-array of the measured retention times in seconds.\nL = p[2] ... number of the length of the column in m.\nd = p[3] ... number of the diameters of the column in m.\nprog = p[4] ... m-array of structure GasChromatographySimulator.Programs containing the definition of the GC-programs.\nopt = p[5] ... struture GasChromatographySimulator.Options containing the settings for options for the simulation.\ngas = p[6] ... string of name of the mobile phase gas. \nmetric = p[7] ... string of the metric used for the loss function (squared or abs). \n\nOutput\n\nsum((tR.-tRcalc).^2) ... sum of the squared residuals over m GC-programs and n solutes.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.opt_dKcentric-Tuple{Any, Any}","page":"Docstrings","title":"RetentionParameterEstimator.opt_dKcentric","text":"opt_dKcentric(x, p)\n\nFunction used for optimization of the loss-function in regards to the three K-centric parameters and the column diameter.\n\nArguments\n\nx ... (3n+1)-vector of the column diameter and the three K-centric parameters of n solutes. The first element represents the column diameterd d and elements 2:n+1 are Tchar, n+2:2n+1 are θchar and 2n+2:3n+1 are ΔCp values.\np ... vector containing the fixed parameters:\ntR = p[1] ... mxn-array of the measured retention times in seconds.\nL = p[2] ... number of the length of the column in m.\nprog = p[3] ... m-array of structure GasChromatographySimulator.Programs containing the definition of the GC-programs.\nopt = p[4] ... struture GasChromatographySimulator.Options containing the settings for options for the simulation.\ngas = p[5] ... string of name of the mobile phase gas. \nmetric = p[6] ... string of the metric used for the loss function (squared or abs). \n\nOutput\n\nsum((tR.-tRcalc).^2) ... sum of the squared residuals over m GC-programs and n solutes.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.optimize_Kcentric-Union{Tuple{T}, Tuple{Any, Any, Any, Vector{T}, Vector{T}, Vector{T}}} where T<:Number","page":"Docstrings","title":"RetentionParameterEstimator.optimize_Kcentric","text":"optimize_Kcentric(tR, col, prog, Tchar_e::Vector{T}, θchar_e::Vector{T}, ΔCp_e::Vector{T}; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, metric=\"squared\")\n\nOptimization regarding the estimization of the retention parameters Tchar, θchar and ΔCp. The initial guess is a vector (Tchar_e, θchar_e and ΔCp_e) for optimization algorithms, which do not need lower/upper bounds. If a method is used, which needs lower/upper bounds the initial parameters should be a matrix, with the first column ([1,:]) beeing the initial guess, the second column ([2,:]) the lower bound and the third column ([3,:]) the upper bound.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.optimize_dKcentric-Union{Tuple{T}, Tuple{Any, Any, Any, Number, Vector{T}, Vector{T}, Vector{T}}} where T<:Number","page":"Docstrings","title":"RetentionParameterEstimator.optimize_dKcentric","text":"optimize_dKcentric(tR, col, prog, d_e, Tchar_e, θchar_e, ΔCp_e; method=NewtonTrustRegion(), opt=std_opt, maxiters=10000, metric=\"squared\")\n\nOptimization regarding the estimization of the column diameter d and the retention parameters Tchar, θchar and ΔCp. The initial guess is a number (d_e) or a vector (Tchar_e, θchar_e and ΔCp_e) for optimization algorithms, which do not need lower/upper bounds. If a method is used, which needs lower/upper bounds the initial parameters should be a vector of length 3 (d_e) or a matrix (Tchar_e, θchar_e and ΔCp_e), with the first element/column  beeing the initial guess, the second element/column the lower bound and the third element/column the upper bound.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.reference_holdup_time-Tuple{Any, Any}","page":"Docstrings","title":"RetentionParameterEstimator.reference_holdup_time","text":"reference_holdup_time(prog, L, d, gas; control=\"Pressure\")\n\nCalculate the reference holdup time for the GC program prog for a column with length L and diameter d and gas as mobile phase. The reference holdup time is the holdup time at the reference temperature 150°C.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.stderror-Tuple{Any, Any}","page":"Docstrings","title":"RetentionParameterEstimator.stderror","text":"stderror(meas, res)\n\nCalculation of the standard error of the found optimized parameters using the hessian matrix at the optima.\n\nArguments\n\nmeas ... Tuple with the loaded measurement data, see load_chromatograms.\nres ... Dataframe with the result of the optimization, see estimate_parameters.\n\nOptional parameters\n\ncol_input ... Named tuple with col_input.L the column length in m and col_input.d the column diameter in mm. If this parameter is not gicen, than these parameters are taken from meas. \n\nOutput\n\nstderrors ... Dataframe with the standard errors of the optimized parameters.\nHessian ... The hessian matrix at the found optims. \n\n\n\n\n\n","category":"method"},{"location":"docstrings/#RetentionParameterEstimator.tR_calc-NTuple{7, Any}","page":"Docstrings","title":"RetentionParameterEstimator.tR_calc","text":"tR_calc(Tchar, θchar, ΔCp, L, d, prog, gas; opt=std_opt)\n\nCalculates the retention time tR for a solute with the K-centric parameters Tchar θchar and ΔCp for a column with length L, internal diameter d, the (conventional) program prog, options opt and mobile phase gas. For this calculation only the ODE for the migration of a solute in a GC column is solved, using the function GasChromatographySimulator.solving_migration.\n\n\n\n\n\n","category":"method"}]
}
