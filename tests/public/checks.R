# Deterministic checks: fixed arrays, direct arithmetic, and isolated solver mocks.
# Loaded by scripts/run_checks.R; check() records each assertion.
public_checks <- function(check) {
  close <- function(x, y, tol = 1e-10)
    length(x) == length(y) && all(is.finite(c(x, y))) && max(c(0, abs(x - y))) <= tol
  raises <- function(expr) tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  rng_before <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  config <- aft_public_config()
  control <- aft_cv_v012_controls(config, .1)
  check("configuration_locked_budgets", control$fista$max_iter == 30000 &&
          control$polish$max_newton_updates == 10 && control$fista$kkt_tol == 1e-7)

  # Direct double-loop oracle independent of the blocked core implementation.
  direct <- function(Y, delta, design, x, h, lambda) {
    W <- aft_explicit_design(design); n <- nrow(W)
    residual <- Y - as.vector(W %*% x); loss <- 0; g <- numeric(length(x))
    H <- matrix(0, length(x), length(x))
    for (i in seq_len(n)) for (j in seq_len(n)) if (i != j && delta[i] != 0) {
      r <- residual[j] - residual[i]; w <- W[i, ] - W[j, ]; N <- n * (n - 1)
      loss <- loss + delta[i] * (r * pnorm(r / h) + h * dnorm(r / h)) / N
      g <- g + delta[i] * pnorm(r / h) * w / N
      H <- H + delta[i] * dnorm(r / h) * tcrossprod(w) / (N * h)
    }
    ib <- which(seq_along(x) > design$q); kk <- abs(g)
    nz <- ib[x[ib] != 0]; zz <- ib[x[ib] == 0]
    kk[nz] <- abs(g[nz] + lambda * sign(x[nz]))
    kk[zz] <- pmax(abs(g[zz]) - lambda, 0)
    list(loss = loss, objective = loss + lambda * sum(abs(x[ib])),
         gradient = g, H = H, kkt = kk)
  }
  W <- rbind(c(-2, 0, 1), c(-1, 1, -1), c(0, -1, 2),
             c(1, 2, 0), c(2, -2, -2), c(0, 0, 0))
  design <- aft_design_matrix(W, q = 1)
  Y <- c(-.9, -.2, .5, .8, 1.6, -.4); delta <- c(1, 1, 0, 1, 0, 1)
  x <- c(.12, -.08, 0); h <- .63; lambda <- .01
  oracle <- direct(Y, delta, design, x, h, lambda)
  for (block in c(1, 2, 64)) {
    core <- aft_loss(Y, delta, design, x, h, block_size = block)
    check(paste0("direct_pair_loss_block_", block), close(core$loss, oracle$loss))
    check(paste0("direct_pair_gradient_block_", block), close(core$gradient, oracle$gradient))
  }
  independent <- phase2c_endpoint_arithmetic(Y, delta, design, x, h, lambda)
  check("independent_endpoint_objective", close(independent$objective, oracle$objective))
  check("independent_endpoint_full_KKT", close(independent$kkt$coordinates, oracle$kkt))
  epsilon <- 1e-6
  fd <- vapply(seq_along(x), function(j) {
    d <- numeric(length(x)); d[j] <- epsilon
    (aft_loss(Y, delta, design, x + d, h)$loss - aft_loss(Y, delta, design, x - d, h)$loss) / (2 * epsilon)
  }, numeric(1))
  check("gradient_finite_difference", close(fd, oracle$gradient, 1e-8))
  hs <- aft_active_face_hessian(Y, delta, design, x, h, seq_along(x))
  check("Hessian_direct_pairs", close(as.vector(hs$H), as.vector(oracle$H)))
  check("Hessian_dual_paths", isTRUE(hs$double_path_pass))
  check("partial_prox_preserves_unpenalized_step",
        close(aft_partial_prox(c(2, -.5, 0), c(.3, -.4, .1), .5, .2, 1), c(1.85, -.2, 0)))
  check("KKT_includes_nuisance_active_and_zero",
        close(aft_kkt(c(2, -.5, 0), c(.3, -.4, .1), 1, .2)$coordinates, c(.3, .6, 0)))
  segment <- aft_remainder_bound(c(0, 2), c(1, 0),
    aft_design_matrix(matrix(c(-1, 1), 2, 1), 0), 0, 2, .5)
  check("remainder_uses_whole_residual_segment", close(as.numeric(segment), dnorm(0) * 16 / (2 * 2 * .5)))

  # Binary64 references were generated independently at 100-digit precision.
  ref <- read.csv(file.path(aft_project_root(), "tests", "fixtures", "tail_reference.csv"),
                  colClasses = c(expected_100digit = "character", expected_binary64_le = "character"))
  expected <- vapply(ref$expected_binary64_le, function(hex) {
    starts <- seq(1, 15, 2)
    readBin(as.raw(strtoi(substring(hex, starts, starts + 1), 16L)), "double", n = 1, size = 8, endian = "little")
  }, numeric(1))
  # The reference file contains more than one h; evaluate each row separately.
  actual <- vapply(seq_len(nrow(ref)), function(i) aft_smooth(-ref$z[i] * ref$h[i], ref$h[i]), numeric(1))
  normal <- expected >= .Machine$double.xmin
  check("tail_normal_reference", all(abs(actual[normal] - expected[normal]) / expected[normal] < 7e-12))
  check("tail_subnormal_reference", all(abs(actual[!normal] - expected[!normal]) <= 2 * .Machine$double.xmin * .Machine$double.eps))
  check("finite_log_budget_keeps_positive_mass", .ap_bound(-1000) > 0)

  # Four-row analytic partial-penalty fit: one real FISTA call and one real polish call.
  X <- matrix(c(-1, 1, -1, 1), 4, 1, dimnames = list(NULL, "x1"))
  Z <- matrix(c(-1, -1, 1, 1), 4, 1, dimnames = list(NULL, "z1"))
  D <- aft_preprocess(X, Z, K = 0)
  ys <- as.vector(Z + .4 * X); ds <- rep(1, 4)
  star <- c(1, .4 - .3 * qnorm(.575)); old <- star + .01
  small <- config; small$fista$max_iter <- 2000L
  fit <- aft_solve_cv_v012(ys, ds, D, old, .6, .1, small)
  check("fixed_FISTA_success", isTRUE(fit$success))
  check("fixed_FISTA_analytic_solution", close(fit$standard$coefficients, star, 1e-6))
  check("fixed_FISTA_budget_and_no_restart", fit$fista$iterations > 0 && fit$fista$iterations <= 2000 &&
          identical(fit$fista$tolerances$restart, FALSE) && all(diff(fit$fista$trace$L) >= 0))
  check("fixed_FISTA_full_KKT", fit$standard$kkt$residual_with_underflow_budget <= 1e-7)
  check("fixed_FISTA_adapter_once", fit$counts$fista_calls == 1 && fit$counts$polish_calls == 0)
  polish <- aft_polish_active_face(ys, ds, D, old, .6, .1, control$polish)
  check("real_polish_success", isTRUE(polish$success) && polish$iterations > 0 && polish$iterations <= 10)
  check("real_polish_analytic_solution", close(polish$coefficients, star, 1e-7))
  values <- vapply(polish$points, function(p) p$objective, numeric(1))
  check("real_polish_objective_descent", all(diff(values) < 0))
  check("real_polish_sign_preservation", all(vapply(polish$points, function(p) p$coefficients[2] > 0, logical(1))))
  check("global_certificate_stays_false", identical(fit$global_optimum_certified, FALSE) &&
          identical(polish$global_optimum_certified, FALSE))

  # Mock scheduling uses real endpoint arithmetic. It deliberately supplies stale
  # auxiliary fields so the test detects reuse of the FISTA beta/b_pred/scores.
  raw <- function(coefficients, status, iterations, D, Y = ys, delta = ds) {
    core <- aft_loss(Y, delta, D, coefficients, .6)
    objective <- core$loss + .1 * sum(abs(coefficients[seq_along(coefficients) > D$q]))
    list(coefficients = coefficients, design = D, h = .6, lambda = .1,
         gradient = core$gradient, numerics = core$numerics, objective = objective,
         kkt = .cv012_adapter_kkt(coefficients, core$gradient, D$q, .1, .cv012_adapter_budget(core$numerics)$value),
         status = status, success = identical(status, "stationary_to_tolerance"),
         iterations = iterations, global_optimum_certified = FALSE,
         beta = 999, b_pred = 999, fitted_scores = rep(999, D$n))
  }
  mock_fista <- function(Y, delta, design, h, lambda, initial, control, checkpoint, ...) {
    for (i in seq(26000L, 30000L, 1000L)) checkpoint(list(iterations = i, coefficients = old, state = list(x = old)))
    raw(old, "iteration_limit", 30000L, D = design)
  }
  mock_polish <- function(...) raw(star, "stationary_to_tolerance", 1L, D = D)
  rebuilt <- aft_solve_cv_v012(ys, ds, D, old, .6, .1, config,
                               fista_fn = mock_fista, polish_fn = mock_polish)
  check("adapter_mock_one_conditional_polish", isTRUE(rebuilt$success) && rebuilt$counts$polish_calls == 1)
  check("adapter_keeps_original_FISTA_coefficients", close(rebuilt$fista$coefficients, old, 0))
  check("adapter_final_coefficients_beta_bpred", close(rebuilt$standard$coefficients, star) &&
          close(rebuilt$standard$beta, star[2]) && close(rebuilt$standard$b_pred, star[2]))
  at <- direct(ys, ds, D, star, .6, .1)
  check("adapter_final_gradient_and_predictions", close(rebuilt$standard$gradient, at$gradient) &&
          close(rebuilt$standard$fitted_scores, as.vector(aft_explicit_design(D) %*% star)))
  mock_failure <- function(...) raw(c(0, 0), "line_search_failure", 0L, D = D)
  no_polish <- function(...) stop("Unexpected polisher invocation")
  rejected <- aft_solve_cv_v012(ys, ds, D, c(0, 0), .6, .1, config,
                                fista_fn = mock_failure, polish_fn = no_polish)
  check("adapter_failure_does_not_retry", identical(rejected$status, "line_search_failure") && rejected$counts$polish_calls == 0)
  four_checkpoints <- function(Y, delta, design, h, lambda, initial, control, checkpoint, ...) {
    for (i in seq(27000L, 30000L, 1000L)) checkpoint(list(iterations = i, coefficients = old, state = list(x = old)))
    raw(old, "iteration_limit", 30000L, D = design)
  }
  gated <- aft_solve_cv_v012(ys, ds, D, old, .6, .1, config,
                             fista_fn = four_checkpoints, polish_fn = no_polish)
  check("adapter_requires_five_distinct_checkpoints", identical(gated$status, "iteration_limit_polish_gate_rejected") &&
          gated$counts$polish_calls == 0 && !gated$gate$enough_checkpoints)
  corrupt_gradient <- function(...) {
    value <- raw(old, "line_search_failure", 0L, D = D)
    value$gradient[1] <- value$gradient[1] + 1e-5
    value
  }
  corrupt <- aft_solve_cv_v012(ys, ds, D, old, .6, .1, config,
                               fista_fn = corrupt_gradient, polish_fn = no_polish)
  check("adapter_rejects_stale_gradient", identical(corrupt$status, "endpoint_formula_discrepancy") &&
          isTRUE(corrupt$blocked) && corrupt$counts$polish_calls == 0)

  # Nonempty factor block: mocked polishing changes gamma and deployment slope.
  factorX <- cbind(c(-2, -2, 2, 2), c(-1, 1, -1, 1))
  factorZ <- matrix(c(-1, 1, 1, -1), 4, 1)
  factorD <- aft_preprocess(factorX, factorZ, K = 1L)
  factorstar <- c(.8, .5, 0, 0); factorold <- c(.81, .51, 0, 0)
  factorY <- as.vector(aft_explicit_design(factorD) %*% factorstar)
  factor_fista <- function(Y, delta, design, h, lambda, initial, control, checkpoint, ...) {
    for (i in seq(26000L, 30000L, 1000L)) checkpoint(list(iterations = i, coefficients = factorold, state = list(x = factorold)))
    raw(factorold, "iteration_limit", 30000L, D = design, Y = Y, delta = delta)
  }
  factor_polish <- function(Y, delta, design, ...) raw(factorstar, "stationary_to_tolerance", 1L, D = design, Y = Y, delta = delta)
  factor_result <- aft_solve_cv_v012(factorY, ds, factorD, factorold, .6, .1, config,
                                     fista_fn = factor_fista, polish_fn = factor_polish)
  new_bp <- aft_prediction_coefficients(factorD, factorstar)
  old_bp <- aft_prediction_coefficients(factorD, factorold)
  check("FA_adapter_polish_final_gamma_and_bpred", isTRUE(factor_result$success) &&
          factor_result$counts$polish_calls == 1 && close(factor_result$standard$gamma, .5) &&
          close(factor_result$standard$b_pred, new_bp) && max(abs(new_bp - old_bp)) > 1e-3)
  check("FA_adapter_prediction_final_coefficients", close(
    aft_predict_map(factorD, factor_result$standard$coefficients, factorX, factorZ), factorY))

  # A held-out perturbation must leave training PCA and frozen schema unchanged.
  XF <- rbind(c(1, 3, -1), c(2, -2, 1), c(4, 0, 3), c(-1, 1, -2),
              c(0, 4, 2), c(3, -3, 0), c(50, -90, 70), c(-60, 80, -40))
  colnames(XF) <- sprintf("gene_%04d", 1:3)
  ZF <- matrix(c(2, 0, -1, 3, 1, -2, 100, -100), 8, 1, dimnames = list(NULL, "clin_01"))
  yf <- c(-2, -1, 0, 1, 2, 3, 9, -9); df <- c(1, 1, 0, 1, 1, 0, 1, 1)
  prepared <- p3a_prepare_split(XF, ZF, yf, df, 1:8, 1:6, 7:8, 1)
  changed <- XF; changed[7:8, ] <- changed[7:8, ] * 100
  changed_prepared <- p3a_prepare_split(changed, ZF, yf, df, 1:8, 1:6, 7:8, 1)
  check("held_out_data_do_not_change_PCA", identical(prepared$design, changed_prepared$design))
  xf <- c(.2, .3, .4, 0, -.2); Df <- prepared$design
  saved <- serialize(prepared, NULL)
  prediction <- p3a_predict_frozen(prepared, xf, XF[7:8, , drop = FALSE], ZF[7:8, , drop = FALSE])
  bp <- xf[3:5] + as.vector(Df$Bhat %*% ((xf[2] - as.vector(crossprod(Df$Bhat, xf[3:5]))) / diag(Df$Lambdahat)))
  manual <- as.vector(sweep(ZF[7:8, , drop = FALSE], 2, Df$z_mean) %*% xf[1] +
                      sweep(XF[7:8, , drop = FALSE], 2, Df$x_mean) %*% bp)
  check("FA_frozen_prediction_formula", close(prediction, manual) && identical(saved, serialize(prepared, NULL)))
  check("FA_schema_column_reorder_rejected", raises(p3a_predict_frozen(prepared, xf, XF[7:8, 3:1], ZF[7:8, , drop = FALSE])))
  rr <- raw(xf, "iteration_limit", 1L, D = Df, Y = yf[1:6], delta = df[1:6])
  endpoint <- aft_cv_v012_endpoint(yf[1:6], df[1:6], Df, xf, .6, .1, rr, config)
  check("FA_final_beta_bpred_reconstructed", close(endpoint$standard$beta, xf[3:5]) &&
          close(endpoint$standard$b_pred, bp) && isTRUE(endpoint$verification$pass))

  # Validation oracle and pair-weighted aggregation with unequal fold sizes.
  val_direct <- function(y, d, s) {
    total <- 0
    for (i in seq_along(y)) for (j in seq_along(y)) if (i != j)
      total <- total + d[i] * (max(y[j] - y[i] + s[i] - s[j], 0) - max(y[j] - y[i], 0))
    total
  }
  folds <- list(list(y = c(-1, .5, 2, 3), d = c(1, 0, 1, 1), s = c(.2, -.1, .6, 1)),
                list(y = c(-2, 0, 1), d = c(0, 1, 1), s = c(-.4, .3, .2)))
  scores <- lapply(folds, function(f) p3a_validation_scalar(f$y, f$d, f$s))
  for (i in seq_along(folds)) {
    f <- folds[[i]]; blocked <- p3a_validation_blocked(f$y, f$d, f$s, block_size = 2)
    check(paste0("validation_direct_pairs_", i), close(scores[[i]]$numerator, val_direct(f$y, f$d, f$s)))
    check(paste0("validation_fixed_sum_order_", i), identical(scores[[i]]$numerator, blocked$numerator))
  }
  aggregate <- p3a_aggregate_folds(scores, required_folds = 2)
  check("validation_pair_weighted_aggregation", aggregate$denominator == 18 &&
          close(aggregate$CV, sum(vapply(scores, `[[`, numeric(1), "numerator")) / 18))
  failure <- list(status = "solver_failed", success = FALSE)
  check("failed_fold_disqualifies_candidate", !p3a_aggregate_folds(c(scores, list(failure)))$eligible)
  records <- lapply(c(-1, -1 + 1.5e-10, -.5, 0), function(v) list(status = "eligible", CV = v))
  selection <- p3a_select_lambda(records, c(.15, .06, .025, .01))
  check("selection_tie_prefers_larger_lambda", selection$selected_lambda == .15 && selection$boundary_side == "upper")
  check("selection_all_failed", p3a_select_lambda(rep(list(failure), 4), c(.15, .06, .025, .01))$status == "tuning_failed")
  state <- p3a_warm_start_new(c(0, 0), "fixed-scope")
  failed <- list(status = "iteration_limit", success = FALSE, coefficients = c(999, 999))
  check("failed_fit_does_not_advance_warm_start", identical(p3a_warm_start_update(state, failed, "fixed-scope")$coefficients, c(0, 0)))
  check("warm_start_scope_enforced", raises(p3a_warm_start_update(state, failed, "other-scope")))
  check("RNG_state_unchanged", identical(rng_before,
    if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL))
  invisible(list(real_fista_calls = 1L, real_polish_calls = 1L,
                 fista_accepted_updates = fit$fista$iterations,
                 polish_accepted_updates = polish$iterations,
                 research_data_read = FALSE, random_cohorts_generated = FALSE))
}
