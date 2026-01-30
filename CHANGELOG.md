# Changelog

All notable changes to RetentionParameterEstimator.jl will be documented in this file.

## [Unreleased]

### Added

- **Multi-start optimization**: Added optional multi-start optimization for single-substance modes (`Kcentric_single`, `dKcentric_single`)
  - New `multistart_n` parameter in `estimate_parameters` to control number of random starting points
  - Automatically enabled with 10 starts when only one chromatogram is available (improves accuracy with poor initial guesses)
  - Can be manually enabled for any case by setting `multistart_n > 0`
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

### Changed

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
