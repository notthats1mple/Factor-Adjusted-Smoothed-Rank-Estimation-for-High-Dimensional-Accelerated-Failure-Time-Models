# Code review guide

Start with the fixed demonstration, then follow an input through preprocessing, optimization, endpoint verification, and prediction:

```sh
Rscript --vanilla scripts/check_environment.R
Rscript --vanilla scripts/run_checks.R
Rscript --vanilla examples/fixed_array_demo.R
```

## Review map

| Question | Primary source | What to inspect |
|---|---|---|
| Which coordinates are penalized? | [R/preprocess.R](../R/preprocess.R), [R/diagnostics.R](../R/diagnostics.R) | `(theta, gamma, beta)` ordering, nuisance dimension `q`, partial proximal map |
| Is the pairwise objective implemented consistently? | [R/loss.R](../R/loss.R), [diagnostics/phase2b/endpoint.R](../diagnostics/phase2b/endpoint.R) | Event anchoring, off-diagonal normalization, tail arithmetic, independent endpoint calculation |
| Does FISTA follow its stated contract? | [optim/phase2c/fista.R](../optim/phase2c/fista.R) | Extrapolated gradient, nondecreasing `L`, whole-segment acceptance, returned-point KKT |
| Does Newton remain on the specified face? | [optim/active_face_polish/polish.R](../optim/active_face_polish/polish.R), [optim/active_face_polish/arithmetic.R](../optim/active_face_polish/arithmetic.R) | Hessian cross-check, signs, exact zeros, inactive inequalities, recorded trials |
| Are final prediction fields consistent after refinement? | [tuning/phase3a/adapter.R](../tuning/phase3a/adapter.R) | Rebuilt gradient, `beta`, `b_pred`, and fitted scores at final coefficients |
| Is validation isolated from preprocessing? | [tuning/phase3a/data_validation.R](../tuning/phase3a/data_validation.R), [pilot/schema.R](../pilot/schema.R) | Training-complement fit, frozen map, column schema, within-fold pairs, eligibility rules |
| Are source and environment versions identifiable? | [provenance/algorithm-source.json](../provenance/algorithm-source.json), [renv.lock](../renv.lock) | Source hashes, locked package versions, local environment output |

[Algorithms](ALGORITHMS.md) provides formulas and field definitions. [Reproducibility](REPRODUCIBILITY.md) explains how to run the checks.

## Numerical outcomes

Read the returned status together with full KKT and the tracked numerical budget. Review successful endpoints, explicit failures, iteration limits, and gate rejections as distinct outcomes. Counts of accepted updates and gradient evaluations describe different aspects of computational cost; preserve both when comparing runs.

The fixed arrays support implementation-level review. Claims about the full research study should be assessed against its archived training inputs, fit records, and execution protocol.

## Report a finding

Include the commit, R and package versions, command, minimal input, expected behavior, actual behavior, and relevant status or trace fields. For a formula discrepancy, identify the coordinate or pair and show the two computed quantities. For a documentation correction, identify the function and the proposed wording. Keep private data and credentials out of issue attachments.
