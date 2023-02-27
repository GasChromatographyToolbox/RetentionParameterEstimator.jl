# RetentionParameterEstimator.jl

[![DOI](https://zenodo.org/badge/550339258.svg)](https://zenodo.org/badge/latestdoi/550339258)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JanLeppert.github.io/RetentionParameterEstimator.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JanLeppert.github.io/RetentionParameterEstimator.jl/dev)

Estimation of retention parameters for the interaction of analytes with a stationary phase in Gas Chromatography (GC).

The retention parameters are estimated from a set of temperature programmed GC runs. The GC simulation ['GasChromatographySimulator.jl'](https://github.com/JanLeppert/GasChromatographySimulator.jl) is used to compute the retention times with several sets of estimated retention parameters and compare these computed retention times with the measured retention times. An optimization process is used to minimize the difference between computed and measured retention times. The retention parameters resulting in this minimized difference are the final result. In addition it is also possible to estimate the column diameter _d_.

## Installation

To install the package type:

```julia
julia> ] add RetentionParameterEstimator
```

To use the package type:

```julia
julia> using RetentionParameterEstimator
```

## Documentation

Please read the [documentation page](...) for more information.

## Notebooks

In the folder [notebooks](https://github.com/JanLeppert/RetentionParameterEstimator/tree/main/notebooks) notebooks, using [Pluto.jl](https://github.com/fonsp/Pluto.jl), for the estimation of retention parameters from temperature programmed GC measurements are available. 

To use these notebooks [Julia, v1.6 or above,](https://julialang.org/downloads/#current_stable_release) must be installed and **Pluto** must be added:

```julia
julia> ]
(v1.7) pkg> add Pluto
```

To run Pluto, use the following commands:

```julia
julia> using Pluto
julia> Pluto.run()
```

Pluto will open your browser. In the field `Open from file` the URL of a notebook or the path to a locally downloaded notebook can be insert and the notebook will open and load the necessary packages. 

### Overview of notebooks



## Contribution

Please open an issue if you:
- want to report a bug 
- have problems using the package (please first look at the documentation)
- have ideas for new features or ways to improve the usage of this package 

You can contribute (e.g. fix bugs, add new features, add to the documentation) to this package by Pull Request: 
- first discuss your contributions in a new issue
- ensure that all tests pass locally before starting the pull request
- new features should be included in `runtests.jl`
- add description to the pull request, link to corresponding issues by `#` and issue number
- the pull request will be reviewed

## Citation

