# Section 5 eq:comp-lambda-zero: fit the nuisance block before the null threshold.
# The threshold is usable only at the recorded absolute full-gradient precision.
aft_lambda_zero <- function(Y, delta, design, h, control = list()) {
  A <- aft_nuisance_matrix(design)
  nuisance <- aft_design_matrix(A, design$q)
  null <- aft_solve(Y, delta, nuisance, h, lambda = 1, control = control)
  precision <- list(absolute_tolerance = null$tolerances$kkt_tol,
                    lost_gradient_fraction = .1,
                    lost_gradient_log_budget = log(null$tolerances$kkt_tol) - log(10),
                    source = 'nuisance control kkt_tol, also applied to the full design',
                    scope = 'tracked gradient underflow plus computed nuisance residual; not an interval-arithmetic certificate')
  if (!null$success) return(list(status = null$status, success = FALSE, threshold_usable = FALSE,
    null_fit = null, lambda_zero = NA_real_, precision = precision, numerics = NULL,
    approximate_threshold = TRUE, global_optimum_certified = FALSE))
  x <- c(null$coefficients, numeric(design$p))
  ev <- tryCatch(aft_loss(Y, delta, design, x, h, null$tolerances$block_size), error = identity)
  if (inherits(ev, 'error')) return(list(status = 'numerical_failure', success = FALSE,
    threshold_usable = FALSE, null_fit = null, lambda_zero = NA_real_,
    precision = precision, numerics = NULL, approximate_threshold = TRUE,
    global_optimum_certified = FALSE, message = conditionMessage(ev)))
  ib <- if (design$q < design$P) seq.int(design$q + 1L, design$P) else integer()
  computed <- aft_max_abs(ev$gradient[ib])
  nuisance_residual <- aft_max_abs(ev$gradient[seq_len(design$q)])
  lost <- ev$numerics$lost_gradient_log_upper_bound
  metadata_valid <- is.numeric(lost) && length(lost) == 1L &&
    !is.na(lost) && (is.finite(lost) || identical(lost, -Inf))
  lost_bound <- if (metadata_valid) exp(lost) else NA_real_
  positive_bound_underflow <- metadata_valid && is.finite(lost) && lost_bound == 0
  # Only the reported error bound is rounded upward when it cannot be represented.
  # This never floors lambda_zero or a requested lambda. The log bound remains primary.
  if (positive_bound_underflow) lost_bound <- .Machine$double.xmin * .Machine$double.eps
  combined <- nuisance_residual + lost_bound
  reliable <- metadata_valid && lost <= precision$lost_gradient_log_budget &&
    is.finite(combined) && combined <= precision$absolute_tolerance
  # A computable zero can be usable at an absolute error budget without proving
  # the mathematical threshold is exactly zero. Keep the raw value on refusal.
  list(status = if (reliable) 'stationary_to_tolerance' else 'numerical_uncertainty',
       success = reliable, threshold_usable = reliable, null_fit = null,
       coefficients = x, lambda_zero = if (reliable) computed else NA_real_,
       computed_lambda_zero = computed,
       threshold_underflow_diagnostic_range = if (metadata_valid)
         c(lower = max(0, computed - lost_bound), upper = computed + lost_bound) else c(lower = NA_real_, upper = NA_real_),
       diagnostic_range_is_interval_certificate = FALSE,
       positive_error_bound_underflow = positive_bound_underflow,
       nuisance_residual = nuisance_residual,
       nuisance_residual_with_underflow_budget = combined,
       lost_gradient_log_upper_bound = lost, lost_gradient_upper_bound = lost_bound,
       numerics = ev$numerics, precision = precision,
       gradient = ev$gradient, approximate_threshold = TRUE,
       global_optimum_certified = FALSE,
       message = if (reliable) NULL else 'Full-design gradient uncertainty or nuisance residual exceeds the declared precision budget')
}

aft_path <- function(Y, delta, design, h, lambdas, control = list(), null_control = list()) {
  if (!is.numeric(lambdas) || !length(lambdas) || any(!is.finite(lambdas)) || any(lambdas <= 0))
    stop('Supply explicit strictly positive lambda levels')
  null <- aft_lambda_zero(Y, delta, design, h, null_control)
  if (!null$success) return(list(status = null$status, success = FALSE, null = null, fits = list(),
                                global_optimum_certified = FALSE))
  ordered <- sort(unique(lambdas), decreasing = TRUE)
  x <- null$coefficients; fits <- vector('list', length(ordered))
  for (j in seq_along(ordered)) {
    fit <- aft_solve(Y, delta, design, h, ordered[j], initial = x, control = control)
    fits[[j]] <- fit
    if (fit$success) x <- fit$coefficients
  }
  list(status = if (all(vapply(fits, function(x) x$success, logical(1)))) 'all_stationary_to_tolerance' else 'candidate_failures',
       success = all(vapply(fits, function(x) x$success, logical(1))), null = null,
       lambdas = ordered, requested_lambdas = lambdas, fits = fits,
       global_optimum_certified = FALSE)
}
