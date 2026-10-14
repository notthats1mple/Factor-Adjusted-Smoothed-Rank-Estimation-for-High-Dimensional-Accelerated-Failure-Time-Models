# Reproducibility

## Prepare the R environment

The preserved [renv.lock](../renv.lock) records R **4.5.2** and package versions. The main external solver dependency is `lpSolve` **5.6.23**; configuration and tests use `jsonlite` **2.0.0** and `testthat` **3.3.2**. The lock also records their dependency tree and `renv` **1.2.4**.

From a fresh checkout, run:

```sh
Rscript --vanilla scripts/restore_dependencies.R
Rscript --vanilla scripts/check_environment.R
```

The restoration script prepares the locked package environment using [`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html). It fixes `renv` at **1.2.4**, uses `retry = FALSE`, and checks the installed package versions against the lock. Preserve installation errors and resolve the reported dependency or toolchain issue while retaining the lock file. A prepared environment can go directly to the environment check.

For comparable numerical records, use a single BLAS thread when launching R:

```sh
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
Rscript --vanilla scripts/check_environment.R
```

Record the R executable, `.libPaths()`, package versions, platform, and BLAS/LAPACK reported by the check alongside any computation. Floating-point details can differ across R builds, BLAS implementations, and processors.

## Verify this checkout

```sh
python3 scripts/check_repository.py
Rscript --vanilla scripts/run_checks.R
Rscript --vanilla examples/fixed_array_demo.R
```

The Python repository check uses the standard library for file and SHA256 verification. R checks cover the fixed arithmetic and bounded integration cases included in `tests/`. The example runs a small fit on fixed arrays and demonstrates the frozen prediction map. These commands exercise the published implementation with bundled fixtures.

The fixed checks address preprocessing identities, stable loss arithmetic, partial penalization, full KKT calculations, active-face refinement, and reconstruction of final adapter fields. Treat a failed check as a reproducible result: retain its output and environment record before changing code or controls.

## Load and use the code

This is a source-based R project. Load its dependency-ordered entry point:

```r
source("/path/to/checkout/scripts/load_project.R")
config <- aft_public_config()
```

The loader resolves paths relative to the checkout. The code accepts user-supplied training arrays through the interfaces in [Algorithms](ALGORITHMS.md). Save input hashes, the configuration, initial coefficients, returned fit, status, and evaluation counts for each application. For the common adapter, save both raw solver records and the final `standard` record.

## Source identity and scope

[provenance/algorithm-source.json](../provenance/algorithm-source.json) records the origin and SHA256 of preserved algorithm files from the authors' local archive. Portable loading, public examples, and verification entry points accompany these archived implementations.

The public evidence consists of code, locked dependencies, fixed fixtures, and executable checks. Full study cohorts, selected fits, saved trajectories, and manuscript materials remain in the local archive. The published checks validate the specified computational cases; reproducing the complete study requires its archived training inputs and execution records.
