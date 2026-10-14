# Fixed-objective nonaccelerated partial proximal gradient, Algorithm 1.
aft_control <- function(lambda, control = list()) {
  defaults <- list(kkt_tol = min(1e-8, lambda * 1e-6), max_iter = 30000L,
                   max_backtracks = 80L, t_max = 1, zeta = .5, omega = .5,
                   block_size = 64L, lp_tolerance = 1e-9, rank_tolerance = 1e-10,
                   trace_vectors = FALSE)
  if (!is.list(control) || (length(control) && is.null(names(control))) ||
      any(!names(control) %in% names(defaults))) stop('Unknown or unnamed solver control')
  out <- utils::modifyList(defaults, control)
  positive <- c('kkt_tol', 't_max', 'lp_tolerance', 'rank_tolerance')
  for (nm in positive) if (length(out[[nm]]) != 1L || !is.finite(out[[nm]]) || out[[nm]] <= 0)
    stop(paste('Invalid control:', nm))
  if (out$kkt_tol > lambda / 100) stop('kkt_tol must be <= lambda/100; tighten the numerical tolerance')
  for (nm in c('max_iter', 'max_backtracks', 'block_size'))
    if (length(out[[nm]]) != 1L || !is.finite(out[[nm]]) || out[[nm]] != as.integer(out[[nm]]) ||
        out[[nm]] < if (nm == 'block_size') 1 else 0) stop(paste('Invalid control:', nm))
  for (nm in c('zeta', 'omega')) if (length(out[[nm]]) != 1L || !is.finite(out[[nm]]) || out[[nm]] <= 0 || out[[nm]] >= 1)
    stop(paste('Invalid control:', nm))
  if (length(out$trace_vectors) != 1L || !is.logical(out$trace_vectors) || is.na(out$trace_vectors)) stop('Invalid trace_vectors')
  out
}

aft_nuisance_matrix <- function(design) {
  if (!design$q) return(matrix(numeric(), design$n, 0))
  if (!is.null(design$Zc) && !is.null(design$Fhat)) return(cbind(design$Zc, design$Fhat))
  if (!is.null(design$W)) return(design$W[, seq_len(design$q), drop = FALSE])
  vapply(seq_len(design$q), function(j) {
    v <- numeric(design$P); v[j] <- 1; aft_forward(design, v)
  }, numeric(design$n))
}

aft_solve <- function(Y, delta, design, h, lambda, initial = NULL, control = list()) {
  if (!inherits(design, 'aft_design')) stop('design must be produced by aft_preprocess or aft_design_matrix')
  if (!is.numeric(Y) || length(Y) != design$n || any(!is.finite(Y)) ||
      length(delta) != design$n || any(!is.finite(delta)) || any(!delta %in% c(0, 1))) stop('Invalid Y or event coding')
  if (length(h) != 1L || !is.finite(h) || h <= 0 ||
      length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) stop('h and lambda must be finite and positive')
  ctl <- aft_control(lambda, control)
  x <- if (is.null(initial)) numeric(design$P) else initial
  if (!is.numeric(x) || length(x) != design$P || any(!is.finite(x))) stop('Invalid initial coefficients')
  q <- design$q; ib <- if (q < design$P) seq.int(q + 1L, design$P) else integer()
  diagnosis <- aft_nuisance_diagnostic(aft_nuisance_matrix(design), delta, ctl$lp_tolerance, ctl$rank_tolerance)
  counts <- list(objective = 0L, gradient = 0L, remainder = 0L, backtracks = 0L, trials = 0L)
  trace <- list(); steps <- 0L; last_step <- NA_real_; current <- NULL
  aggregate_numerics <- list()
  evaluate <- function(at) {
    counts$objective <<- counts$objective + 1L; counts$gradient <<- counts$gradient + 1L
    ans <- aft_loss(Y, delta, design, at, h, ctl$block_size)
    if (!is.finite(ans$loss) || ans$loss < 0 || any(!is.finite(ans$gradient))) stop('Nonfinite/negative objective or nonfinite gradient')
    lost <- ans$numerics$lost_gradient_log_upper_bound
    if (!is.null(lost) && lost > log(ctl$kkt_tol) - log(10)) {
      aggregate_numerics[[length(aggregate_numerics) + 1L]] <<- list(unreliable_gradient = ans$numerics)
      stop(structure(list(message = 'Underflow error bound exceeds one tenth of KKT tolerance', call = NULL),
                     class = c('aft_uncertain_error','error','condition')))
    }
    ans
  }
  finish <- function(status, message = NULL) {
    available <- !is.null(current)
    kk <- if (available) tryCatch(aft_kkt(x, current$gradient, q, lambda), error = function(e) NULL) else NULL
    if (!is.null(kk)) {
      lb <- current$numerics$lost_gradient_log_upper_bound
      kk$lost_gradient_upper_bound <- if (is.null(lb)) 0 else exp(lb)
      kk$residual_with_underflow_budget <- kk$residual + kk$lost_gradient_upper_bound
    }
    gm <- if (available) tryCatch(aft_gradient_mapping(x, current$gradient, q, lambda,
                                             if (is.finite(last_step)) last_step else ctl$t_max), error = function(e) NULL) else NULL
    d <- if (!is.null(design$d)) design$d else q
    K <- if (!is.null(design$K)) design$K else 0L
    theta <- x[seq_len(d)]; gamma <- if (K) x[d + seq_len(K)] else numeric()
    bpred <- if (!is.null(design$Xc)) tryCatch(aft_prediction_coefficients(design, x), error = function(e) NULL) else NULL
    fit <- list(status = status, success = identical(status, 'stationary_to_tolerance'),
                global_optimum_certified = FALSE, message = message,
                coefficients = x, theta = theta, gamma = gamma, beta = x[ib], b_pred = bpred,
                design = design, lambda = lambda, h = h, K = K, q = q,
                objective = if (available) current$loss + lambda * sum(abs(x[ib])) else NA_real_,
                smooth_loss = if (available) current$loss else NA_real_,
                kkt = kk, gradient_mapping = gm, gradient = if (available) current$gradient else NULL,
                accepted_step = last_step, iterations = steps, evaluations = counts,
                event_count = sum(delta), diagnosis = diagnosis, tolerances = ctl,
                numerics = if (available) current$numerics else NULL,
                numerical_events = aggregate_numerics, trace = trace,
                fitted_scores = if (available) aft_forward(design, x) else NULL,
                coefficient_scale = 'original molecular scale; only training centering and PCA',
                line_search = 'nonnegative full-segment Taylor bound; no additive acceptance slack')
    class(fit) <- 'aft_fit'; fit
  }
  if (diagnosis$status %in% c('no_events', 'strict_escape_numeric', 'diagnostic_inconclusive'))
    return(finish(diagnosis$status, diagnosis$reason))
  answer <- tryCatch({
    current <- evaluate(x)
    repeat {
      kk <- aft_kkt(x, current$gradient, q, lambda)
      log_lost <- current$numerics$lost_gradient_log_upper_bound
      lost_budget <- if (is.null(log_lost)) 0 else exp(log_lost)
      if (kk$residual + lost_budget <= ctl$kkt_tol) break
      if (steps >= ctl$max_iter) return(finish('iteration_limit', 'Maximum accepted iterations reached'))
      accepted <- FALSE; t <- ctl$t_max
      for (bt in seq.int(0L, ctl$max_backtracks)) {
        if (!is.finite(t) || t <= 0) return(finish('line_search_failure', 'Step underflow'))
        counts$trials <- counts$trials + 1L
        y <- aft_partial_prox(x, current$gradient, t, lambda, q)
        change <- y - x; square <- sum(change^2)
        if (!is.finite(square)) stop('Nonfinite step norm')
        if (square == 0) return(finish('numerical_stagnation', 'Step vanishes in floating point before full KKT tolerance'))
        counts$remainder <- counts$remainder + 1L
        bound <- aft_remainder_bound(Y, delta, design, x, change, h, ctl$block_size)
        margin <- ctl$omega * square / (2 * t)
        if (!is.finite(bound) || !is.finite(margin) || bound < 0) stop('Nonfinite or negative line-search certificate')
        if (margin == 0) return(finish('line_search_failure', 'Quadratic acceptance margin underflowed'))
        if (bound <= margin) { accepted <- TRUE; break }
        if (bt < ctl$max_backtracks) { t <- t * ctl$zeta; counts$backtracks <- counts$backtracks + 1L }
      }
      if (!accepted) return(finish('line_search_failure', 'No certified step within backtracking limit'))
      next_value <- evaluate(y)
      Q <- current$loss + lambda * sum(abs(x[ib]))
      Qy <- next_value$loss + lambda * sum(abs(y[ib]))
      if (!all(is.finite(c(Q, Qy)))) stop('Nonfinite penalized objective')
      # This diagnostic allowance NEVER enters the acceptance test.
      rounding <- 32 * .Machine$double.eps * max(1, abs(Q), abs(Qy))
      if (Qy > Q + rounding) return(finish('numerical_failure', 'Objective increase contradicts certified step beyond rounding allowance'))
      row <- list(iteration = steps + 1L, objective_before = Q, objective_after = Qy,
                  step = t, backtracks = bt, positive_remainder_bound = as.numeric(bound),
                  remainder_numerics = attr(bound, 'numerics'),
                  quadratic_margin = margin, step_squared = square,
                  majorant_defect = next_value$loss - current$loss - sum(current$gradient * change) - margin,
                  descent_rhs = (2 - ctl$omega) * square / (2 * t),
                  diagnostic_rounding = rounding,
                  kkt_after = aft_kkt(y, next_value$gradient, q, lambda)$residual)
      if (ctl$trace_vectors) {
        row$x <- x; row$y <- y; row$gradient <- current$gradient; row$gradient_y <- next_value$gradient
      }
      trace[[length(trace) + 1L]] <- row
      aggregate_numerics[[length(aggregate_numerics) + 1L]] <- list(loss = next_value$numerics, remainder = attr(bound, 'numerics'))
      x <- y; current <- next_value; steps <- steps + 1L; last_step <- t
    }
    finish('stationary_to_tolerance')
  }, error = function(e) finish(if (inherits(e, 'aft_uncertain_error')) 'numerical_uncertainty' else 'numerical_failure', conditionMessage(e)))
  answer
}

aft_fit <- function(Y, delta, X, Z = NULL, K = 0, h, lambda, initial = NULL, control = list()) {
  design <- aft_preprocess(X, Z, K)
  aft_solve(Y, delta, design, h, lambda, initial, control)
}

predict.aft_fit <- function(object, Xnew, Znew = NULL, ...) {
  if (!isTRUE(object$success)) stop('Prediction requires a stationary_to_tolerance fit')
  aft_predict_map(object$design, object$coefficients, Xnew, Znew)
}
