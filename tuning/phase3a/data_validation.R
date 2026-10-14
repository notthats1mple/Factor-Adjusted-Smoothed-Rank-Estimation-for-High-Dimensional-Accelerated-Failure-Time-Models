# Phase 3a v0.12: observation-only fold preparation and validation arithmetic.
# The entry point extracts X/Z/Y/delta before calling these functions. This file
# performs no file access or optimization. Only p3a_make_folds uses RNG.

.p3a_data_stop <- function(code, message, ...) {
  stop(structure(c(list(message = message, call = NULL, code = code), list(...)),
    class = c("p3a_data_error", "phase3a_data_error", "error", "condition")))
}

.p3a_vector <- function(x, n, name, finite = TRUE) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || length(x) != n)
    .p3a_data_stop("invalid_observations", paste(name, "must be a numeric vector of length", n))
  if (finite && any(!is.finite(x)))
    .p3a_data_stop("nonfinite_observations", paste(name, "must be finite"))
  as.numeric(x)
}

.p3a_row_ids <- function(row_id, n) {
  row_id <- .p3a_vector(row_id, n, "row_id")
  if (any(row_id < 1 | row_id != floor(row_id)) || anyDuplicated(row_id))
    .p3a_data_stop("invalid_row_id", "row_id must contain unique positive original integer row identifiers")
  row_id
}

.p3a_delta <- function(delta, n) {
  delta <- .p3a_vector(delta, n, "delta")
  if (any(!delta %in% c(0, 1)))
    .p3a_data_stop("invalid_event_indicator", "delta must contain only zero and one")
  delta
}

p3a_observations <- function(X, Z, Y, delta, row_id = seq_len(nrow(X))) {
  if (!is.matrix(X) || !is.numeric(X) || is.complex(X) || any(!is.finite(X)))
    .p3a_data_stop("invalid_observations", "X must be a finite real numeric matrix")
  n <- nrow(X)
  if (n < 2L) .p3a_data_stop("invalid_observations", "At least two observation rows are required")
  if (is.null(Z)) Z <- matrix(numeric(), n, 0L)
  if (!is.matrix(Z) || !is.numeric(Z) || is.complex(Z) || nrow(Z) != n || any(!is.finite(Z)))
    .p3a_data_stop("invalid_observations", "Z must be a finite real numeric matrix with the X row count")
  list(X = X, Z = Z, Y = .p3a_vector(Y, n, "Y"),
    delta = .p3a_delta(delta, n), row_id = .p3a_row_ids(row_id, n))
}

.p3a_rng_snapshot <- function() {
  present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  list(kind = RNGkind(), seed_exists = present,
    seed = if (present) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL)
}

p3a_make_folds <- function(delta, row_id, seed, folds = 3L,
                           rng_kind = c("Mersenne-Twister", "Inversion", "Rejection")) {
  delta <- .p3a_delta(delta, length(delta))
  row_id <- .p3a_row_ids(row_id, length(delta))
  if (length(folds) != 1L || !is.numeric(folds) || !is.finite(folds) || folds < 2 || folds != floor(folds))
    .p3a_data_stop("invalid_fold_control", "folds must be an integer at least two")
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) || seed != floor(seed) || abs(seed) > .Machine$integer.max)
    .p3a_data_stop("invalid_fold_control", "seed must be a representable integer")
  required_rng <- c("Mersenne-Twister", "Inversion", "Rejection")
  if (!identical(as.character(rng_kind), required_rng))
    .p3a_data_stop("invalid_fold_control", "The v0.12 RNG kind is fixed")
  folds <- as.integer(folds)
  before <- .p3a_rng_snapshot()
  do.call(RNGkind, as.list(required_rng))
  set.seed(as.integer(seed))
  event_order <- sort(row_id[delta == 1])
  censored_order <- sort(row_id[delta == 0])
  if (length(event_order) > 1L) event_order <- event_order[sample.int(length(event_order))]
  if (length(censored_order) > 1L) censored_order <- censored_order[sample.int(length(censored_order))]
  permutation <- c(event_order, censored_order)
  allocated <- rep(seq_len(folds), length.out = length(permutation))
  labels <- as.integer(allocated[match(row_id, permutation)])
  after <- .p3a_rng_snapshot()
  checks <- do.call(rbind, lapply(seq_len(folds), function(k) {
    v <- labels == k
    data.frame(fold = k, validation_n = sum(v), training_n = sum(!v),
      validation_events = sum(delta[v]), training_events = sum(delta[!v]),
      valid = sum(v) >= 2L && sum(!v) >= 2L && sum(delta[v]) >= 1 && sum(delta[!v]) >= 1)
  }))
  ok <- all(checks$valid)
  list(status = if (ok) "valid" else "invalid_fold_plan", success = ok, valid = ok,
    reason = if (ok) NULL else paste("Partition size/event requirements failed for folds", paste(checks$fold[!checks$valid], collapse = ",")),
    labels = labels, assignment = data.frame(row_id = row_id, delta = delta, fold = labels),
    row_id = row_id, event_order = event_order, censored_order = censored_order,
    permutation = permutation, fold_checks = checks, seed = as.integer(seed), folds = folds,
    rng_before = before, rng_after = after,
    algorithm = "Ascending original row IDs; independently sample events then censored; concatenate; continuous cyclic labels",
    automatic_retry = FALSE)
}

p3a_prepare_split <- function(X, Z, Y, delta, row_id, train_rows, valid_rows, K) {
  obs <- p3a_observations(X, Z, Y, delta, row_id)
  train_rows <- .p3a_row_ids(train_rows, length(train_rows))
  valid_rows <- .p3a_row_ids(valid_rows, length(valid_rows))
  if (length(intersect(train_rows, valid_rows)) ||
      !setequal(c(train_rows, valid_rows), obs$row_id))
    .p3a_data_stop("invalid_split", "Training and validation rows must be disjoint and together cover the observation rows")
  if (length(train_rows) < 2L || sum(obs$delta[match(train_rows, obs$row_id)]) < 1 ||
      (length(valid_rows) && (length(valid_rows) < 2L || sum(obs$delta[match(valid_rows, obs$row_id)]) < 1)))
    .p3a_data_stop("invalid_fold_plan", "Each nonempty partition requires at least two patients and one event")
  take <- function(ids) {
    ix <- match(sort(ids), obs$row_id)
    list(X = obs$X[ix, , drop = FALSE], Z = obs$Z[ix, , drop = FALSE],
      Y = obs$Y[ix], delta = obs$delta[ix], row_id = obs$row_id[ix])
  }
  train <- take(train_rows)
  valid <- take(valid_rows)
  schema <- phase2a_schema(train$X, train$Z)
  phase2a_validate_schema(schema, valid$X, valid$Z)
  design <- aft_preprocess(train$X, train$Z, K = K)
  phase2a_validate_schema(schema, design$Xc, design$Zc)
  list(status = "prepared", success = TRUE, design = design, schema = schema,
    train = train, valid = valid, train_row_id = train$row_id, valid_row_id = valid$row_id,
    preprocessing_scope = if (length(valid_rows)) "training_complement_only" else "full_training_refit",
    preprocessing_calls = 1L, scaled_columns = FALSE)
}

p3a_predict_frozen <- function(prepared, coefficients, Xvalid, Zvalid) {
  if (!is.list(prepared) || !inherits(prepared$design, "aft_design") ||
      !identical(prepared$design$kind, "pca"))
    .p3a_data_stop("invalid_prediction_record", "A prepared fold with frozen PCA design is required")
  design <- prepared$design
  schema <- prepared$schema
  if (design$n != schema$training_rows)
    .p3a_data_stop("training_schema_mismatch", "Frozen design and training schema row counts differ")
  phase2a_validate_schema(schema, design$Xc, design$Zc)
  phase2a_validate_schema(schema, Xvalid, Zvalid)
  coefficients <- .p3a_vector(coefficients, design$P, "coefficients")
  before <- serialize(list(design = design, schema = schema), NULL, version = 3L)
  scores <- aft_predict_map(design, coefficients, Xvalid, Zvalid)
  if (!identical(before, serialize(list(design = design, schema = schema), NULL, version = 3L)))
    .p3a_data_stop("frozen_mapping_changed", "Prediction changed the frozen preprocessing or schema")
  scores
}

p3a_D <- function(a, b) {
  a <- .p3a_vector(a, length(a), "a", finite = FALSE)
  b <- .p3a_vector(b, length(a), "b", finite = FALSE)
  if (any(!is.finite(a)) || any(!is.finite(b)))
    .p3a_data_stop("validation_numerical_failure", "Pair differences a and b must be finite")
  out <- numeric(length(a))
  positive <- b >= 0
  ix <- positive & a >= 0
  out[ix] <- b[ix]
  ix <- positive & a < 0 & a > -b
  out[ix] <- a[ix] + b[ix]
  ix <- !positive & a > 0 & a >= -b
  out[ix] <- b[ix]
  ix <- !positive & a > 0 & a < -b
  out[ix] <- -a[ix]
  if (any(!is.finite(out)) || any(abs(out) > abs(b)))
    .p3a_data_stop("validation_numerical_failure", "The finite pair difference or its magnitude identity failed")
  out
}

.p3a_neumaier <- function(s, correction, value) {
  t <- s + value
  if (!is.finite(t)) .p3a_data_stop("validation_numerical_failure", "Validation running sum became nonfinite")
  increment <- if (abs(s) >= abs(value)) (s - t) + value else (value - t) + s
  correction <- correction + increment
  if (!is.finite(correction)) .p3a_data_stop("validation_numerical_failure", "Validation compensation became nonfinite")
  list(s = t, correction = correction)
}

.p3a_validation <- function(Y, delta, scores, row_id, block_size) {
  n <- length(Y)
  Y <- .p3a_vector(Y, n, "Y", finite = FALSE)
  scores <- .p3a_vector(scores, n, "scores", finite = FALSE)
  delta <- .p3a_delta(delta, n)
  row_id <- .p3a_row_ids(row_id, n)
  base <- list(success = FALSE, numerator = NA_real_, denominator = n * (n - 1),
    score = NA_real_, fold_score = NA_real_, patient_count = n,
    pairs = n * (n - 1), event_pairs = sum(delta) * (n - 1),
    row_id = sort(row_id), pair_order = "ascending original row_id; i outer, j inner; ordered non-diagonal pairs",
    summation = "Neumaier", block_size = as.integer(block_size))
  if (n < 2L || !any(delta == 1))
    return(c(list(status = "invalid_fold_plan", reason = "Validation requires at least two patients and one event"), base))
  if (any(!is.finite(Y)) || any(!is.finite(scores)))
    return(c(list(status = "validation_numerical_failure", reason = "Nonfinite validation outcome or prediction"), base))
  ord <- order(row_id)
  Y <- Y[ord]; scores <- scores[ord]; delta <- delta[ord]; row_id <- row_id[ord]
  s <- 0; correction <- 0; completed_pairs <- 0L; at_i <- at_j <- NA_integer_
  computation <- tryCatch({
    for (i in seq_len(n)) {
      at_i <- i
      for (first in seq.int(1L, n, by = block_size)) {
        for (j in seq.int(first, min(n, first + block_size - 1L))) {
          if (i == j) next
          at_j <- j
          a <- Y[j] - Y[i]
          b <- scores[i] - scores[j]
          value <- delta[i] * p3a_D(a, b)
          acc <- .p3a_neumaier(s, correction, value)
          s <- acc$s; correction <- acc$correction
          completed_pairs <- completed_pairs + 1L
        }
      }
    }
    numerator <- s + correction
    score <- numerator / base$denominator
    if (!is.finite(numerator) || !is.finite(score))
      .p3a_data_stop("validation_numerical_failure", "Final compensated validation value became nonfinite")
    list(numerator = numerator, score = score)
  }, p3a_data_error = function(e) e)
  if (inherits(computation, "error")) {
    if (!identical(computation$code, "validation_numerical_failure")) stop(computation)
    return(c(list(status = "validation_numerical_failure", reason = conditionMessage(computation),
      completed_pairs = completed_pairs, failure_i = row_id[at_i], failure_j = row_id[at_j],
      sum_part = s, correction = correction), base))
  }
  base$success <- TRUE
  base$numerator <- computation$numerator
  base$score <- base$fold_score <- computation$score
  c(list(status = "valid", completed_pairs = completed_pairs, sum_part = s, correction = correction), base)
}

p3a_validation_scalar <- function(Y, delta, scores, row_id = seq_along(Y)) {
  .p3a_validation(Y, delta, scores, row_id, block_size = 1L)
}

p3a_validation_blocked <- function(Y, delta, scores, row_id = seq_along(Y), block_size = 64L) {
  if (length(block_size) != 1L || !is.numeric(block_size) || !is.finite(block_size) ||
      block_size < 1 || block_size != floor(block_size) || block_size > .Machine$integer.max)
    .p3a_data_stop("invalid_validation_control", "block_size must be a positive integer")
  .p3a_validation(Y, delta, scores, row_id, as.integer(block_size))
}

p3a_aggregate_folds <- function(fold_results, required_folds = 3L) {
  if (!is.list(fold_results) || length(required_folds) != 1L || !is.numeric(required_folds) ||
      !is.finite(required_folds) || required_folds < 1 || required_folds != floor(required_folds))
    .p3a_data_stop("invalid_aggregation_input", "Provide a list of fold records and a positive required fold count")
  reasons <- character()
  if (length(fold_results) != required_folds) reasons <- c(reasons, "required_fold_count_not_met")
  scalar_finite <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  for (k in seq_along(fold_results)) {
    r <- fold_results[[k]]
    validation_status <- if (is.list(r)) {
      if (!is.null(r$validation_status)) r$validation_status else r$status
    } else NULL
    ok <- is.list(r) && isTRUE(r$success) && identical(validation_status, "valid") &&
      scalar_finite(r$numerator) && scalar_finite(r$denominator) && r$denominator > 0
    if (!ok) {
      status <- if (length(validation_status) == 1L) as.character(validation_status) else "missing_validation"
      detail <- if (is.list(r) && length(r$reason) == 1L) paste0(":", r$reason) else ""
      reasons <- c(reasons, paste0("fold_", k, ":", status, detail))
    }
  }
  result <- list(status = "ineligible", success = FALSE, eligible = FALSE,
    CV = NA_real_, numerator = NA_real_, denominator = NA_real_,
    reasons = reasons, folds = fold_results, required_folds = as.integer(required_folds), cross_fold_pairs = FALSE)
  if (length(reasons)) return(result)
  totals <- tryCatch({
    s <- 0; correction <- 0; denominator <- 0
    for (r in fold_results) {
      acc <- .p3a_neumaier(s, correction, r$numerator)
      s <- acc$s; correction <- acc$correction
      denominator <- denominator + r$denominator
    }
    numerator <- s + correction
    CV <- numerator / denominator
    if (any(!is.finite(c(numerator, denominator, CV))))
      .p3a_data_stop("validation_numerical_failure", "Pooled validation arithmetic became nonfinite")
    list(numerator = numerator, denominator = denominator, CV = CV)
  }, p3a_data_error = function(e) e)
  if (inherits(totals, "error")) {
    result$reasons <- conditionMessage(totals)
    return(result)
  }
  result$status <- "eligible"; result$success <- result$eligible <- TRUE
  result$CV <- totals$CV; result$numerator <- totals$numerator; result$denominator <- totals$denominator
  result
}

p3a_select_lambda <- function(aggregate_rows, lambdas, atol = 1e-10, rtol = 1e-10) {
  lambdas <- .p3a_vector(lambdas, length(aggregate_rows), "lambdas")
  if (!length(lambdas) || any(lambdas <= 0) || anyDuplicated(lambdas))
    .p3a_data_stop("invalid_selection_input", "The fixed lambda grid must be nonempty, positive, and unique")
  if (!is.list(aggregate_rows) || length(atol) != 1L || length(rtol) != 1L ||
      !is.numeric(atol) || !is.numeric(rtol) || any(!is.finite(c(atol, rtol))) || atol < 0 || rtol < 0)
    .p3a_data_stop("invalid_selection_input", "Invalid aggregate records or tie tolerances")
  CV <- rep(NA_real_, length(lambdas))
  eligible <- rep(FALSE, length(lambdas))
  for (k in seq_along(lambdas)) {
    r <- aggregate_rows[[k]]
    if (!is.list(r)) next
    if (!is.null(r$lambda) && (length(r$lambda) != 1L || !is.numeric(r$lambda) ||
        !is.finite(r$lambda) || r$lambda != lambdas[k]))
      .p3a_data_stop("lambda_record_mismatch", "Aggregate lambda does not match its fixed grid position")
    declared <- if (!is.null(r$eligible)) isTRUE(r$eligible) else identical(r$status, "eligible")
    eligible[k] <- declared && is.numeric(r$CV) && length(r$CV) == 1L && is.finite(r$CV)
    if (eligible[k]) CV[k] <- r$CV
  }
  result <- list(status = "tuning_failed", success = FALSE, selected_lambda = NA_real_,
    selected_index = NA_integer_, selected_CV = NA_real_, best_CV = NA_real_,
    tie_tolerance = NA_real_, tied_lambdas = numeric(), boundary_selected = NA,
    boundary_side = NA_character_, eligible = eligible, lambdas = lambdas, CV = CV,
    automatic_grid_expansion = FALSE)
  if (!any(eligible)) return(result)
  best <- min(CV[eligible])
  tolerance <- atol + rtol * max(1, abs(best))
  tied <- which(eligible & CV <= best + tolerance)
  selected <- tied[which.max(lambdas[tied])]
  result$status <- "selected"; result$success <- TRUE
  result$selected_lambda <- lambdas[selected]; result$selected_index <- selected
  result$selected_CV <- CV[selected]; result$best_CV <- best; result$tie_tolerance <- tolerance
  result$tied_lambdas <- lambdas[tied]
  result$boundary_selected <- lambdas[selected] %in% range(lambdas)
  result$boundary_side <- if (length(lambdas) == 1L) "both" else
    if (lambdas[selected] == max(lambdas)) "upper" else if (lambdas[selected] == min(lambdas)) "lower" else "interior"
  result
}

p3a_warm_start_new <- function(initial, scope_id) {
  initial <- .p3a_vector(initial, length(initial), "initial")
  if (!length(initial) || !is.character(scope_id) || length(scope_id) != 1L || is.na(scope_id) || !nzchar(scope_id))
    .p3a_data_stop("invalid_warm_start", "A finite initial vector and one explicit method/fold/h scope are required")
  list(scope_id = scope_id, initial = initial, coefficients = initial, current_coefficients = initial,
    last_success_coefficients = NULL, last_success_lambda = NA_real_, source = "nuisance_initialization",
    source_id = paste0(scope_id, "/initializer"), attempts = 0L, successful_candidates = 0L)
}

p3a_warm_start_update <- function(state, result, scope_id) {
  if (!is.list(state) || !identical(state$scope_id, scope_id))
    .p3a_data_stop("warm_start_scope_mismatch", "Warm starts remain within their original method/fold/h scope")
  if (!is.list(result)) .p3a_data_stop("invalid_warm_start", "A saved candidate result is required")
  state$attempts <- state$attempts + 1L
  state$last_attempt_status <- result$status
  if (isTRUE(result$success) && identical(result$status, "stationary_to_tolerance")) {
    coefficients <- .p3a_vector(result$coefficients, length(state$initial), "successful coefficients")
    state$coefficients <- state$current_coefficients <- state$last_success_coefficients <- coefficients
    state$last_success_lambda <- if (is.numeric(result$lambda) && length(result$lambda) == 1L) result$lambda else NA_real_
    state$source <- "last_successful_candidate"
    state$source_id <- if (is.character(result$id) && length(result$id) == 1L) result$id else "last_successful_candidate"
    state$successful_candidates <- state$successful_candidates + 1L
  }
  state
}

p3a_refit_decision <- function(selection, refit_result = NULL) {
  selected <- is.list(selection) && identical(selection$status, "selected") &&
    is.numeric(selection$selected_lambda) && length(selection$selected_lambda) == 1L &&
    is.finite(selection$selected_lambda) && selection$selected_lambda > 0
  if (!selected) return(list(status = "tuning_failed", attempt_refit = FALSE,
    selected_lambda = NA_real_, automatic_retry = FALSE, fallback = FALSE))
  result <- list(status = "refit_planned", attempt_refit = is.null(refit_result),
    selected_lambda = selection$selected_lambda, automatic_retry = FALSE, fallback = FALSE,
    initialization = "new_full_training_nuisance_initialization")
  if (!is.null(refit_result)) {
    if (!is.list(refit_result)) .p3a_data_stop("invalid_refit_result", "A saved refit record is required")
    result$success <- isTRUE(refit_result$success) && identical(refit_result$status, "stationary_to_tolerance")
    result$status <- if (result$success) "final_refit_success" else "final_refit_failed"
    result$fit_status <- refit_result$status
  }
  result
}
