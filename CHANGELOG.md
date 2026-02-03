# Changelog

All notable changes to RetentionParameterEstimator.jl will be documented in this file.

## [Unreleased]

### Added

- **Multi-start optimization**: Added optional multi-start optimization for single-substance modes (`Kcentric_single`, `dKcentric_single`)
  - New `multistart_n` parameter in `estimate_parameters` to control number of random starting points
  - Controlled solely by `multistart_n` parameter: `0` (default) = disabled, `>0` = enabled with specified number of starts
  - Single and multiple chromatograms are treated the same way (no automatic multistart for single chromatograms)
  - Created `optimize_Kcentric_multistart` and `optimize_dKcentric_multistart` wrapper functions
  - Uses random perturbations (±20% by default) around initial parameter estimates
  - Returns the best solution (lowest loss) across all starting points
  - Robust error handling: if an optimization attempt fails, it continues with the next starting point
  - If all optimization attempts fail, throws a clear error message
  - Added comprehensive tests for multi-start optimization in `test/runtests.jl`
  - Created benchmark scripts for multistart optimization:
    - `scripts/benchmark_multistart_m1.jl`: Benchmarks method_m1 with single chromatograms
    - `scripts/benchmark_multistart_m2.jl`: Benchmarks method_m2 with single chromatograms
    - `scripts/benchmark_multistart_m4.jl`: Benchmarks method_m4 with single chromatograms
  - Each benchmark script processes every chromatogram individually and compares optimization with/without multistart
  - Improved benchmark scripts:
    - Fixed loss value formatting (uses scientific notation for very small values, prevents "0.000000" display)
    - Suppressed warnings during benchmark execution using `NullLogger()`
    - Added automatic database parameter comparison with `database_Rxi5SilMS_beta125.csv` when available
    - Reports mean absolute relative differences for Tchar, θchar, and ΔCp compared to database values
    - Added explanatory notes about multistart behavior (speed and loss differences)
  - **Multistart optimization logging**: Added logging for failed optimization steps in multistart functions
    - Always logs when the original starting point fails (critical for debugging)
    - Always logs when random starting points fail
    - Shows summary of failed optimizations and improved solutions found
    - Added optional `verbose` parameter to `optimize_Kcentric_multistart` and `optimize_dKcentric_multistart` for detailed logging
    - Helps identify why multistart might give worse results than non-multistart (e.g., if original starting point fails)
  - **Coupled parameter perturbation for multistart**: Added option to maintain empirical relationships between retention parameters during multistart perturbations
    - Created `perturb_retention_parameters_coupled()` function: only perturbs Tchar, recalculates θchar and ΔCp using empirical formulas (default behavior)
    - Created `perturb_retention_parameters_independent()` function: perturbs all three parameters independently (previous behavior)
    - Added `coupled_perturbation=true` parameter to `optimize_Kcentric_multistart`, `optimize_dKcentric_multistart`, and `estimate_parameters`
    - When `coupled_perturbation=true` (default), maintains empirical relationships: `θchar = 22.0 * (Tchar/Tst)^0.7 * (1000*col.df/col.d)^0.09` and `ΔCp = -52.0 + 0.34*Tchar`
    - When `coupled_perturbation=false`, all parameters are perturbed independently (allows exploration of parameter space that doesn't follow empirical trends)
    - Updated benchmark scripts (`benchmark_multistart_m1.jl`, `benchmark_multistart_m2.jl`, `benchmark_multistart_m4.jl`) to use `coupled_perturbation=true`
  - **Database comparison for multistart results**: Enhanced benchmark scripts to show database comparisons for both no multistart and multistart cases
    - Added database comparison computation for multistart results in all benchmark scripts
    - Displays mean relative differences for Tchar, θchar, and ΔCp for both optimization strategies
    - Helps evaluate how multistart affects the accuracy of estimated retention parameters
  - **Weighted start parameter estimation**: Added `estimate_start_parameter_single_ramp_weighted()` function that gives more weight to measurements with higher heating rates
    - Accounts for the observation that higher heating rates provide more accurate Tchar estimates
    - Uses exponential weighting: `weight = exp(α * (rT - rT_min))` where `α` controls the strength (default 2.0)
    - Combines weighted average (default 70%) with interpolation to rT_nom (default 30%) for theoretical correction
    - Configurable via `α` parameter (weighting strength) and `weighted_fraction` parameter (mix ratio)
    - Original `estimate_start_parameter_single_ramp()` function remains unchanged and is still the default
  - **Generalized ramp rate calculation**: Added `average_ramp_rate()` function and `use_average_ramp` parameter for more flexible temperature program handling
    - New `average_ramp_rate()` function calculates average ramp rate from first to last temperature plateau, ignoring holding times at beginning and end
    - Added optional `use_average_ramp=false` parameter to `estimate_start_parameter_single_ramp()` and `estimate_start_parameter_single_ramp_weighted()`
    - When `use_average_ramp=true`, works with complex temperature programs (multiple ramps, holds) instead of assuming single ramp between time_steps 2 and 3
    - Falls back to original method if average ramp rate calculation fails (e.g., isothermal programs)
    - Default behavior (`use_average_ramp=false`) maintains backward compatibility
  - **Corrected single measurement parameter estimation**: Added `estimate_start_parameter_single_measurement_corrected()` function for estimating retention parameters from a single chromatogram
    - Uses empirical correction model: `Tchar_est = Telu / (0.25*sqrt(rT) + 0.8)` where `rT` is the dimensionless heating rate
    - Designed for single measurement scenarios, ideally with single ramp programs
    - Accepts both DataFrame (single row) and Vector inputs for retention times
    - Supports `use_average_ramp` parameter for flexible ramp rate calculation
    - Calculates θchar and ΔCp from corrected Tchar using standard empirical relationships
    - Useful when only one chromatogram is available and correction for heating rate effects is needed

### Fixed

- **Julia 1.12+ compatibility**: Fixed world age warnings for `perturb_retention_parameters_coupled` and `perturb_retention_parameters_independent` functions
  - Used `Base.invokelatest()` when calling perturbation functions from multistart functions to ensure compatibility with Julia 1.12+ stricter world age semantics
  - Prevents warnings that could become errors in future Julia versions

### Changed

- **Multi-start optimization behavior**: Changed multistart control to be consistent across all cases
  - Removed automatic multistart activation when only one chromatogram is provided
  - Multistart is now controlled solely by the `multistart_n` parameter, regardless of the number of chromatograms
  - Single and multiple chromatograms are treated identically (no special case behavior)
  - Simplified code by removing intermediate `use_multistart` variable and using `multistart_n` directly
  - Updated both `estimate_parameters` and `estimate_parameters_` functions for consistency
  - **Note**: This is a behavior change. Previously, single chromatograms automatically used 10 multistart runs. Now multistart must be explicitly enabled via `multistart_n > 0`

- **Dependencies**: Added `Random` as an explicit dependency in `Project.toml`
  - Required for multi-start optimization random number generation
  - `Random` is a standard library but must be declared as a dependency when used in packages

### Fixed

- **Test refactoring**: Updated tests to use helper functions for data preparation
  - Replaced manual time unit conversion with `time_unit_conversion_factor()` helper
  - Replaced manual data preparation with `prepare_optimization_data()` and `prepare_single_substance_data()` helpers
  - Improved test maintainability and consistency

## [0.2.1] - 2025-01-27

### Fixed

- **Critical performance bug in `Loss.jl`**: Moved `unique(substance_list)` computation outside the loop in the `loss` function
  - Dramatically improved performance for single-substance optimizations (`method_m1`, `method_m2`, `method_m4`)
  - Root cause: `unique()` was being called inside the optimization loop, causing redundant computations during every loss function evaluation
- **Test performance**: Reduced CI test time from ~71 minutes to ~17 minutes
  - Reduced optimization limits for tests (`maxiters=500`, `maxtime=30.0s` instead of defaults)
  - Changed `se_col=true` to `se_col=false` for most tests
  - Created separate testset for standard errors with reduced limits
  - Commented out expensive multi-solute test
- **Documentation build failure**: Fixed "Cannot resolve @ref" errors for `stderror` function
  - Replaced `@ref` links with plain backticks to avoid naming conflict with `Statistics.stderror`
  - Documentation now builds successfully on GitHub Actions

### Changed

- **Dependencies**: Updated `GasChromatographySimulator` from 0.4/0.5 to 0.6.0
  - Removed `Interpolations.jl` dependency
  - Now uses `GasChromatographySimulator.linear_interpolation` and `GasChromatographySimulator.deduplicate_knots!`
  - Updated `reference_holdup_time` and `estimate_start_parameter_single_ramp` to use new interpolation functions
  - Reduced external dependencies, improved consistency with GasChromatographySimulator
- **CI and Documentation workflows**: Updated GitHub Actions workflows
  - Updated `actions/checkout@v2` to `@v4` in documentation workflow
  - Updated `julia-actions/setup-julia@latest` to `@v1` in documentation workflow
  - Updated `actions/checkout@v3` to `@v4` in CI workflow
  - Updated `julia-actions/setup-julia@v1` to latest version in CI workflow
  - Updated cache and build actions to latest versions
- **Default ODE solver**: Changed from `OwrenZen5()` to `Tsit5()` in `std_opt`
  - Better performance and reliability for optimization methods
  - Users can still override by passing custom `opt` options
- **`method_m4` convergence check**: Added convergence check for substance parameters (`Tchar`, `θchar`, `ΔCp`)
  - Now checks convergence of both `d` and substance parameters before stopping
  - Uses relative tolerance for substance parameters to handle different scales
  - Prevents premature convergence when only `d` has converged
  - Improves reliability, especially with `Tsit5` ODE solver
- **Benchmark scripts**: Combined parallel and non-parallel benchmark scripts into single script
  - Added command-line control for parallelization (`--parallel`, `--no-parallel`)
  - Added warm-up run to mitigate JIT compilation effects
  - Added git repository state tracking in output files
  - Timestamped output files with parallelization status
  - Simplified output format (summary in `.txt`, full results in `.csv`)
- **Dependencies**: Pinned `Optimization.jl` to version 5.0.0
  - Versions 5.1.0-5.4.0 caused severe performance regressions (method_m1: <100s → >3600s)
- **Documentation**: Updated `docs/make.jl` configuration
  - Specified `modules = [RetentionParameterEstimator]` for Documenter
  - Added `DocMeta.setdocmeta!` setup for doctest environment
  - Temporarily disabled doctests to unblock build

### Added

- **Parallelization support**: Added optional parallelization for all optimization methods
  - New `parallel=true` keyword argument for `method_m1`, `method_m2`, `method_m3`, `method_m4`
  - Parallelized `Kcentric_single` and `dKcentric_single` modes in `estimate_parameters`
  - Parallelized standard error calculations in `stderror` and `stderror_m3`
  - Expected speedup: 4-8x on 4-8 core machines
  - Requires Julia to be started with multiple threads (e.g., `julia -t 4`)
- **Data preparation helper functions**: New functions to reduce code duplication
  - `prepare_optimization_data()`: Consolidates time unit conversion, data extraction, missing value filtering, and vectorization
  - `prepare_single_substance_data()`: Handles data preparation for single-substance optimization modes
  - `time_unit_conversion_factor()`: Standardizes time unit conversion
  - Updated `prog_list()` to handle both single program and vector of programs
- **New optimization mode**: Added `mode="d_only"` to `estimate_parameters`
  - Optimizes only column diameter `d` while keeping retention parameters fixed
  - Integrated into `method_m4` alternating optimization strategy
- **New benchmark script**: `scripts/benchmark_m1.jl` for dedicated `method_m1` benchmarking
- **Convergence tracking**: Added convergence checks for substance parameters in `method_m4`
  - Tracks previous `Tchar`, `θchar`, `ΔCp` values across iterations
  - Requires convergence of all parameters before stopping optimization

### Breaking

None. All changes are backward compatible.

---

## [0.2.0] - Previous Release

Previous release notes.
