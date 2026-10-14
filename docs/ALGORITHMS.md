# Algorithms and interfaces

## Inputs and coefficient coordinates

Load `scripts/load_project.R` before using the functions below. The implementation uses finite numeric arrays with a common row order:

- `Y`: observed log-time, a vector of length `n`.
- `delta`: event indicator, coded `1` for an observed event and `0` for right censoring.
- `X`: an `n × p` matrix of high-dimensional covariates.
- `Z`: optional `n × d` unpenalized covariates.
- `K`: a prespecified nonnegative number of factors.
- `h`: a positive smoothing bandwidth; `lambda`: a positive penalty level for the penalized solvers.

`aft_preprocess(X, Z, K)` estimates the transformation using these training rows. It centers columns without rescaling them and stores the means, factor scores `Fhat`, loadings `Bhat`, and selected eigenvalues `Lambdahat`. The implicit design is

$$W = [Z_c,\widehat F,X_c-\widehat F\widehat B^\top].$$

Coefficients are ordered as `(theta, gamma, beta)`. The first `q = d + K` coordinates are unpenalized; only the final `p` coordinates receive the L1 penalty. `K = 0` gives the unadjusted design. For a supplied numeric design, `aft_design_matrix(W, q)` preserves `W` exactly and uses the same partition convention.

The pairwise objective identifies relative scores. Predictions are centered linear scores in log-time coordinates; an absolute survival-time prediction would require a separately specified location/calibration procedure.

## Objective, gradient, and KKT

For `c = (theta, gamma, beta)`, write `e_i = Y_i - w_i^T c`. The code evaluates

$$Q(c)=\frac{1}{n(n-1)}\sum_{i\ne j}\delta_i\ell_h(e_j-e_i)+\lambda\lVert\beta\rVert_1,$$

$$\ell_h(r)=r\Phi(r/h)+h\phi(r/h),$$

where `Phi` and `phi` are the standard normal distribution and density. Pairs are ordered, diagonal terms are excluded, and the event indicator belongs to the anchor `i`. `aft_loss()` returns the smooth loss, full gradient, pair counts, and numerical metadata. The penalty is added by the solvers.

For the smooth gradient `g`, `aft_kkt(c, g, q, lambda)` computes each coordinate residual:

$$\rho_j=\begin{cases}|g_j|,&j\le q,\\|g_j+\lambda\operatorname{sign}(c_j)|,&j>q,\ c_j\ne0,\\\max(|g_j|-\lambda,0),&j>q,\ c_j=0.\end{cases}$$

The complete residual is `max(rho)`. Active coordinates are defined by exact nonzero coefficients. The accelerated solver and adapter use the complete residual plus the tracked gradient-underflow budget, with the released adapter threshold fixed at `1e-7`. The budget consumes the tolerance. Numerical bounds are floating-point evaluations; result records retain `global_optimum_certified = FALSE` and the corresponding interval-certification metadata.

## Solver entry points

```r
aft_solve(Y, delta, design, h, lambda, initial = NULL, control = list())

aft_solve_fista_experimental(
  Y, delta, design, h, lambda,
  initial = NULL, control = list(), checkpoint = NULL, observer = NULL
)

aft_polish_active_face(
  Y, delta, design, initial, h, lambda,
  control = list(), checkpoint = NULL
)

aft_solve_cv_v012(
  Y, delta, design, initial, h, lambda,
  config = aft_public_config()
)
```

The archived baseline `aft_solve()` remains a proximal-gradient implementation. Its own defaults are defined in `R/solver.R`.

### Backtracking FISTA

`optim/phase2c/fista.R` evaluates the gradient at the extrapolated point `y`. The partial proximal step leaves the first `q` coordinates unthresholded and soft-thresholds the `beta` block. Backtracking starts at the preceding accepted inverse step `L`, doubling it when necessary. The implementation fixes `L_initial = 1`, `L_growth = 2`, `omega = 0.5`, and uses the standard non-restarted momentum recurrence.

For `d = z - y`, acceptance requires the whole-segment smooth-loss remainder bound to be at most `omega * L * ||d||^2 / 2`. Full objective, gradient, and KKT arithmetic are then evaluated at the accepted point `z`. The accelerated objective may vary between accepted points; the acceptance rule is based on the segment from `y` to `z`.

`max_iter` counts accepted updates and is bounded by 30,000. The released configuration allows up to 80 backtracking increases for an update. Each trace row corresponds to an accepted update; evaluation counters separately record gradient evaluations, trials, and remainder evaluations. Periodic checkpoints occur every 1,000 accepted updates. Failure records retain the last accepted point and diagnostic state.

### Fixed-active-face Newton refinement

`optim/active_face_polish/polish.R` fixes all initially zero penalized coordinates at zero and retains the signs of initially nonzero penalized coordinates. All nuisance coordinates remain free. An already stationary starting point returns after zero updates.

For an eligible face, the active Hessian is evaluated by independent pairwise and Laplacian constructions. Positive-definiteness and linear-solve accuracy checks precede the Newton step. The line search uses a fraction-to-sign-boundary step, a whole-segment Armijo sufficient condition including the gradient budget, and complete candidate arithmetic. Every accepted point preserves the original signs and zeros and satisfies all inactive-coordinate inequalities. Objective differences, trial decisions, directions, and numerical diagnostics are retained.

The released adapter permits at most ten Newton updates in one refinement call. The face, objective, and penalty scale remain fixed; failed curvature or line-search conditions produce an explicit status.

### Common solve adapter and final fields

`tuning/phase3a/adapter.R` calls FISTA once, independently verifies its returned endpoint, and invokes refinement only when FISTA returns `iteration_limit` and the fixed gate passes. The gate requires five distinct saved terminal checkpoints with identical exact `beta` support and signs, a gradient budget at most `kkt_tol / 10`, and a strictly positive inactive-coordinate margin. A stationary FISTA endpoint completes without refinement.

After each solver return, the adapter reconstructs the objective, full gradient, KKT, coefficient blocks, and predictions at the returned coefficient vector. It compares these fields with independent pairwise arithmetic. The standard result is:

```r
result$status
result$success
result$standard$coefficients
result$standard$gradient
result$standard$kkt
result$standard$beta
result$standard$b_pred
result$standard$fitted_scores
result$fista
result$polish
result$verification
```

Check `result$success` and `result$status` before using a fit as a converged candidate. `stationary_to_tolerance` denotes a verified endpoint; an iteration limit, gate rejection, or numerical failure remains a recorded unsuccessful outcome. The adapter includes call counts, original statuses, checkpoints, and timing records, with automatic retries disabled.

## Frozen prediction

For a preprocessing object from `aft_preprocess()`, the deployment slope is

$$b_{\rm pred}=\beta+\widehat B\widehat\Lambda^{-1}(\gamma-\widehat B^\top\beta).$$

`aft_predict_map(design, coefficients, Xnew, Znew)` centers new rows using training means and applies `b_pred` together with the clinical coefficient `theta`. The sparse estimand `beta` and the generally dense deployment slope `b_pred` serve different roles. The training PCA object remains fixed during prediction.

## Training-validation components

`tuning/phase3a/data_validation.R` contains reusable components for a specified training-only selection procedure. Its `phase2a_schema()` contract requires `X` columns named `gene_0001`, `gene_0002`, … and `Z` columns named `clin_01`, `clin_02`, …, in that exact order. Validation checks names and order without rearranging columns. The lower-level `aft_preprocess()` checks numeric inputs and dimensions and accepts other column names.

1. `p3a_prepare_split()` fits centering and PCA on the training complement and records the column schema.
2. `p3a_predict_frozen()` validates new columns and applies that saved transformation.
3. `p3a_validation_blocked()` evaluates the unsmoothed, baseline-subtracted pairwise score within validation rows, using compensated summation. It normalizes by `n_valid * (n_valid - 1)`.
4. `p3a_aggregate_folds()` pools fold numerators and denominators after every required fold succeeds.
5. `p3a_select_lambda()` selects the lowest eligible score, resolving numerical ties toward the larger penalty. It records a boundary selection.

These components and the solver adapter are included for review and reuse with a specified data split. The original six-scenario execution scheduler and its archived inputs remain part of the local study archive.
