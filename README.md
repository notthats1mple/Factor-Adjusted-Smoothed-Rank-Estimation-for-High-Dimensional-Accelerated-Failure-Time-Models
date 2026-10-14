# Factor-Adjusted Smoothed Rank Estimation for AFT Models

Research software for sparse estimation in high-dimensional, right-censored accelerated failure time (AFT) models. The R implementation combines training-only PCA, a Gaussian-smoothed Gehan objective with partial L1 penalization, full-coordinate KKT diagnostics, and frozen prediction maps.

**Authors:** Weineng Zhu and Yang Ju · **Code version:** 0.17.0 · [中文说明](README_zh.md)

## Start here

Use R **4.5.2**, the version recorded in [renv.lock](renv.lock). From a fresh checkout:

```sh
git clone https://github.com/notthats1mple/Factor-Adjusted-Smoothed-Rank-Estimation-for-High-Dimensional-Accelerated-Failure-Time-Models.git
cd Factor-Adjusted-Smoothed-Rank-Estimation-for-High-Dimensional-Accelerated-Failure-Time-Models
Rscript --vanilla scripts/restore_dependencies.R
Rscript --vanilla scripts/check_environment.R
Rscript --vanilla scripts/run_checks.R
Rscript --vanilla examples/fixed_array_demo.R
```

Run dependency restoration once when preparing an environment. An existing environment with the locked versions can start with `check_environment.R`. The demonstration uses fixed arrays and a bounded fit; it requires no downloaded data. See [Reproducibility](docs/REPRODUCIBILITY.md) for the environment and verification workflow.

## Implementation

| Component | Entry point | Purpose |
|---|---|---|
| Training preprocessing | `aft_preprocess()` | Center covariates, estimate factors, and retain the training transformation |
| Objective and derivatives | `aft_loss()`, `aft_hvp()` | Evaluate the event-anchored pairwise loss, full gradient, and Hessian-vector products |
| Baseline solver | `aft_solve()` | Proximal gradient with partial penalization |
| Accelerated solver | `aft_solve_fista_experimental()` | Non-restarted backtracking FISTA with a nondecreasing inverse step |
| Local refinement | `aft_polish_active_face()` | Newton updates on a fixed active sign face |
| Common solve adapter | `aft_solve_cv_v012()` | Verify the endpoint and conditionally invoke refinement |
| Prediction | `aft_predict_map()` | Apply the frozen training map to new covariates |
| Training validation | `p3a_prepare_split()`, `p3a_validation_blocked()` | Fold-local preprocessing and within-fold validation arithmetic |

[Algorithms and interfaces](docs/ALGORITHMS.md) specifies the objective, coefficient ordering, stopping criteria, and result fields. The [review guide](docs/REVIEW_GUIDE.md) connects these contracts to source and tests.

Load the implementation into an R session from any working directory:

```r
source("/path/to/checkout/scripts/load_project.R")
config <- aft_public_config()
```

For the common adapter, use `result$standard$coefficients` and its reconstructed fields for downstream work. After a successful refinement, these fields contain the refined coefficients, gradient, sparse `beta`, deployment slope `b_pred`, and fitted scores. Raw solver records remain available in `result$fista` and `result$polish`.

## Repository contents

```text
R/                         Core preprocessing, arithmetic, PG solver, and path helpers
optim/                     FISTA and fixed-active-face Newton implementations
diagnostics/               Independent endpoint arithmetic
tuning/phase3a/             Common adapter and training-validation components
config/                    Portable, frozen solver controls
examples/                  Fixed-array demonstration
tests/                     Arithmetic, integration, and regression checks
scripts/                   Loading, environment, and verification entry points
provenance/                Archived-source identities and SHA256 records
docs/                      Algorithm, reproducibility, and review documentation
```

This code release supports implementation review and experiments with user-supplied training data. The manuscript, study cohorts, fitted objects, and research trajectories are retained in the authors' local archive. The bundled checks reproduce fixed computational cases; recreating the full study requires those archived inputs and execution records.

## Citation and contributions

Use [CITATION.cff](CITATION.cff) to cite this software version. Report reproducible issues and proposed changes using [CONTRIBUTING.md](CONTRIBUTING.md).

The authors have not assigned a reuse license to this code. All rights are reserved; contact the authors through the repository for reuse permission. Dependency packages retain their respective licenses.
