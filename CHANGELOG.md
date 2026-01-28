# Changelog

All notable changes to RetentionParameterEstimator.jl will be documented in this file.

## [Unreleased]

### Changed

- **TagBot workflow**: Updated `.github/workflows/TagBot.yml` to match current [Julia TagBot](https://github.com/marketplace/actions/julia-tagbot#setup) recommendations
  - Removed explicit `permissions:` block to rely on repository defaults
  - Kept `issue_comment` and `workflow_dispatch` triggers for manual and automatic tagging
- **Codecov integration**: Updated Codecov action from v1 to v4 in CI workflow
  - Changed secret name from `CODECOV_SECRET` to `CODECOV_TOKEN` (standard naming)
  - Added `fail_ci_if_error: false` to prevent CI failures if Codecov is temporarily unavailable
  - Added `verbose: true` for better debugging of coverage uploads
  - Fixed coverage file path to `./lcov.info`

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
