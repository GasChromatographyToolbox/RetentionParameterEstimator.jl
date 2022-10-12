# RetentionParameterEstimator
 
Estimation of thermodynamic parameters for the interaction of analytes with a stationary phase in Gas Chromatography (GC).

The thermodynamic parameters are estimated from a set of temperature programmed GC runs. The GC simulation ['GasChromatographySimulator.jl'](https://github.com/JanLeppert/GasChromatographySimulator.jl) is used to compute the retention times with several sets of estimated thermodynamic parameters and compare these computed retention times with the measured retention times. An optimization process is used to minimize the difference between computed and measured retention times. The thermodynamic parameters resulting in this minimized difference are the final result.  
