# Phase2c v0.10 experimental non-restarted FISTA. The frozen R/ core is reused.
# The acceptance certificate covers the full segment from the extrapolated y.
# This file neither replaces aft_solve nor supplies a resume/retry mechanism.

aft_fista_control <- function(lambda, control = list()) {
  defaults <- list(kkt_tol = 1e-7, max_iter = 30000L, max_backtracks = 80L,
    L_initial = 1, L_growth = 2, omega = .5, block_size = 64L,
    lp_tolerance = 1e-9, rank_tolerance = 1e-10, restart = FALSE,
    L_may_decrease = FALSE, full_KKT_every_accepted_point = TRUE,
    return_last_accepted_point_on_failure = TRUE,
    trace_schema = paste0('one row per accepted update; no per-iteration ',
                          'coefficients; numeric columns plus exceptional numerical-event records'),
    checkpoint_every_updates = 1000L, automatic_resume = FALSE)
  if (!is.list(control) || (length(control) && is.null(names(control))) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control)) ||
      any(!names(control) %in% names(defaults))) stop('Unknown or unnamed FISTA control')
  out <- utils::modifyList(defaults, control)
  for (nm in c('kkt_tol', 'L_initial', 'lp_tolerance', 'rank_tolerance'))
    if (!is.numeric(out[[nm]]) || length(out[[nm]]) != 1L ||
        !is.finite(out[[nm]]) || out[[nm]] <= 0) stop(paste('Invalid control:', nm))
  if (out$kkt_tol > lambda / 100) stop('kkt_tol must be <= lambda/100')
  for (nm in c('max_iter', 'max_backtracks', 'block_size', 'checkpoint_every_updates'))
    if (!is.numeric(out[[nm]]) || length(out[[nm]]) != 1L ||
        !is.finite(out[[nm]]) || out[[nm]] != floor(out[[nm]]) ||
        out[[nm]] < if (nm %in% c('max_iter', 'max_backtracks')) 0 else 1)
      stop(paste('Invalid control:', nm))
  if (out$max_iter > 30000L || out$max_backtracks > 80L)
    stop('Phase2c accepted-update/backtracking budget exceeded')
  if (!identical(as.numeric(out$L_initial), 1) ||
      !identical(as.numeric(out$L_growth), 2) ||
      !identical(as.numeric(out$omega), .5))
    stop('Phase2c requires L_initial=1, L_growth=2, omega=.5')
  for (nm in c('restart', 'L_may_decrease', 'automatic_resume'))
    if (!identical(out[[nm]], FALSE)) stop(paste('Phase2c requires FALSE:', nm))
  for (nm in c('full_KKT_every_accepted_point', 'return_last_accepted_point_on_failure'))
    if (!identical(out[[nm]], TRUE)) stop(paste('Phase2c requires TRUE:', nm))
  if (!identical(out$trace_schema, defaults$trace_schema)) stop('Unknown trace schema')
  if (out$checkpoint_every_updates != 1000L) stop('Phase2c checkpoint interval is 1000 updates')
  out
}

aft_solve_fista_experimental <- function(Y, delta, design, h, lambda,
    initial = NULL, control = list(), checkpoint = NULL, observer = NULL) {
  if (!inherits(design, 'aft_design')) stop('design must be an aft_design')
  if (!is.numeric(Y) || length(Y) != design$n || any(!is.finite(Y)) ||
      !is.numeric(delta) || length(delta) != design$n || any(!is.finite(delta)) ||
      any(!delta %in% c(0, 1))) stop('Invalid Y or event coding')
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0 ||
      !is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0)
    stop('h and lambda must be finite and positive')
  ctl <- aft_fista_control(lambda, control)
  x <- if (is.null(initial)) numeric(design$P) else initial
  if (!is.numeric(x) || length(x) != design$P || any(!is.finite(x)))
    stop('Invalid initial coefficients')
  if ((!is.null(checkpoint) && !is.function(checkpoint)) ||
      (!is.null(observer) && !is.function(observer))) stop('Callbacks must be functions or NULL')
  q <- design$q
  if (length(q) != 1L || !is.finite(q) || q != floor(q) || q < 0 || q > design$P)
    stop('Invalid nuisance block dimension')
  ib <- if (q < design$P) seq.int(q + 1L, design$P) else integer()
  initial_coefficients <- x
  initial_arithmetic <- NULL
  steps <- 0L; tau <- 1; L_previous <- ctl$L_initial; y <- x
  last_step <- NA_real_; current <- NULL; current_kkt <- NULL; current_gm <- NULL
  diagnosis <- NULL; working <- list(stage = 'initialization')
  accepted_state <- list(x = x, x_prev = x, y = y, tau = tau, L = L_previous)
  events <- list()
  counts <- list(objective = 0L, gradient = 0L, remainder = 0L, backtracks = 0L,
    trials = 0L, initial_gradient = 0L, y_gradient = 0L, z_gradient = 0L,
    y_cache_hits = 0L, kkt = 0L, gradient_mapping = 0L)
  loss_numeric_names <- c('gap_underflow', 'positive_gap_absorbed', 'tail_J_underflow',
    'cdf_underflow', 'cdf_rounded_one', 'cdf_normalization_underflow',
    'subnormal_gradient_weight_count', 'lost_gradient_log_upper_bound',
    'asymptotic_count', 'max_tail_relative_remainder')
  trace_names <- c('iteration', 'tau', 'tau_next', 'L', 'step', 'objective',
    'objective_before', 'objective_y', 'smooth_loss', 'smooth_loss_y', 'kkt',
    'kkt_with_budget', 'lost_gradient_upper_bound', 'y_lost_gradient_upper_bound',
    'gradient_mapping_inf', 'gradient_mapping_2', 'step_norm',
    'accepted_change_norm', 'step_squared', 'backtracks', 'remainder_bound',
    'quadratic_margin', 'majorant_defect', 'objective_increase',
    paste0('y_', loss_numeric_names), paste0('z_', loss_numeric_names),
    'remainder_term_underflow', 'interval_certificate', paste0('eval_', names(counts)))
  scalar_trace <- matrix(NA_real_, nrow = ctl$max_iter, ncol = length(trace_names),
                         dimnames = list(NULL, trace_names))

  add_event <- function(type, details = NULL) {
    events[[length(events) + 1L]] <<- list(type = type, accepted_updates = steps,
      stage = working$stage, details = details)
  }
  fail <- function(message, status = 'numerical_failure') {
    stop(structure(list(message = message, call = NULL, status = status),
                   class = c('aft_fista_error', 'error', 'condition')))
  }
  metadata_value <- function(meta, name) {
    value <- meta[[name]]
    if (is.null(value)) return(if (name == 'lost_gradient_log_upper_bound') -Inf else 0)
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        (name != 'lost_gradient_log_upper_bound' && !is.finite(value)))
      fail(paste('Invalid numerical metadata:', name))
    as.numeric(value)
  }
  budget <- function(ev) {
    loglost <- metadata_value(ev$numerics, 'lost_gradient_log_upper_bound')
    if (loglost > log(ctl$kkt_tol) - log(10))
      fail('Underflow error bound exceeds one tenth of KKT tolerance', 'numerical_uncertainty')
    exp(loglost)
  }
  evaluate <- function(at, role) {
    working$evaluation_point <<- at
    working$evaluation_role <<- role
    counts$objective <<- counts$objective + 1L
    counts$gradient <<- counts$gradient + 1L
    counts[[paste0(role, '_gradient')]] <<- counts[[paste0(role, '_gradient')]] + 1L
    ans <- aft_loss(Y, delta, design, at, h, ctl$block_size)
    working$last_evaluation <<- list(point = at, role = role, arithmetic = ans)
    if (length(ans$loss) != 1L || !is.finite(ans$loss) || ans$loss < 0 ||
        length(ans$gradient) != design$P || any(!is.finite(ans$gradient)))
      fail('Nonfinite/negative loss or nonfinite/mismatched gradient')
    ans
  }
  objective_at <- function(at, ev) {
    val <- ev$loss + lambda * sum(abs(at[ib]))
    if (!is.finite(val)) fail('Nonfinite penalized objective')
    val
  }
  kkt_at <- function(at, ev, check_budget = TRUE) {
    counts$kkt <<- counts$kkt + 1L
    kk <- aft_kkt(at, ev$gradient, q, lambda)
    lb <- metadata_value(ev$numerics, 'lost_gradient_log_upper_bound')
    kk$lost_gradient_upper_bound <- if (check_budget) budget(ev) else exp(lb)
    kk$residual_with_underflow_budget <- kk$residual + kk$lost_gradient_upper_bound
    kk
  }
  mapping_at <- function(at, ev, step) {
    counts$gradient_mapping <<- counts$gradient_mapping + 1L
    gm <- aft_gradient_mapping(at, ev$gradient, q, lambda, step)
    if (any(!is.finite(c(gm$vector, gm$norm_inf, gm$norm_2))))
      fail('Nonfinite gradient mapping')
    gm
  }
  finish <- function(status, message = NULL) {
    available <- !is.null(current)
    # Return diagnostics are always attached to x, the last completed accepted
    # point (or the explicitly recorded zero-displacement stationary candidate).
    kk <- current_kkt
    if (available && is.null(kk)) kk <- tryCatch(kkt_at(x, current, FALSE), error = function(e) NULL)
    gm <- current_gm
    if (available && is.null(gm)) gm <- tryCatch(mapping_at(x, current,
      if (is.finite(last_step)) last_step else 1 / ctl$L_initial), error = function(e) NULL)
    d <- if (!is.null(design$d)) design$d else q
    K <- if (!is.null(design$K)) design$K else 0L
    bpred <- if (!is.null(design$Xc))
      tryCatch(aft_prediction_coefficients(design, x), error = function(e) NULL) else NULL
    trace <- as.data.frame(scalar_trace[seq_len(steps), , drop = FALSE])
    fit <- list(solver_id = 'phase2c_fista_v0.10_nonrestarted_nondecreasing_L',
      status = status, success = identical(status, 'stationary_to_tolerance'),
      global_optimum_certified = FALSE, message = message, coefficients = x,
      theta = x[seq_len(d)], gamma = if (K) x[d + seq_len(K)] else numeric(),
      beta = x[ib], b_pred = bpred, design = design, h = h, lambda = lambda, K = K, q = q,
      objective = if (available) tryCatch(objective_at(x, current), error = function(e) NA_real_) else NA_real_,
      smooth_loss = if (available) current$loss else NA_real_,
      gradient = if (available) current$gradient else NULL, kkt = kk, gradient_mapping = gm,
      accepted_step = last_step, iterations = steps, evaluations = counts,
      numerics = if (available) current$numerics else NULL, numerical_events = events,
      trace = trace, diagnosis = diagnosis, event_count = sum(delta), tolerances = ctl,
      initial_coefficients = initial_coefficients, initial_arithmetic = initial_arithmetic,
      state = accepted_state, checkpoint = accepted_state, last_attempt = working,
      fitted_scores = if (available) tryCatch(aft_forward(design, x), error = function(e) NULL) else NULL,
      coefficient_scale = 'original molecular scale; saved design unchanged',
      line_search = 'nonnegative full-segment Taylor bound at extrapolated y; no additive acceptance slack',
      numerical_metadata_schema = list(loss_columns = loss_numeric_names,
        loss_constants = if (available) current$numerics[setdiff(names(current$numerics), loss_numeric_names)] else NULL,
        remainder_interval_certificate = FALSE,
        accepted_metadata_storage = 'every accepted update in y_/z_ numeric trace columns',
        exceptional_metadata_storage = 'numerical_events and last_attempt'),
      automatic_resume = FALSE)
    class(fit) <- c('aft_fit', 'aft_fista_experimental_fit')
    fit
  }

  tryCatch(withCallingHandlers({
    working$stage <- 'initial_arithmetic'
    current <- evaluate(x, 'initial')
    initial_arithmetic <- current
    # Retain completed initial arithmetic even when its budget is inadequate.
    current_kkt <- kkt_at(x, current, FALSE)
    working$stage <- 'initial_nuisance_diagnostic'
    diagnosis <- aft_nuisance_diagnostic(aft_nuisance_matrix(design), delta,
                                         ctl$lp_tolerance, ctl$rank_tolerance)
    accepted_state$gradient_x <- current$gradient
    accepted_state$numerics_x <- current$numerics
    if (diagnosis$status %in% c('no_events', 'strict_escape_numeric', 'diagnostic_inconclusive'))
      return(finish(diagnosis$status, diagnosis$reason))
    budget(current)
    current_gm <- mapping_at(x, current, 1 / ctl$L_initial)
    objective_at(x, current)
    if (current_kkt$residual_with_underflow_budget <= ctl$kkt_tol)
      return(finish('stationary_to_tolerance'))
    if (steps >= ctl$max_iter) return(finish('iteration_limit', 'Maximum accepted updates reached'))

    repeat {
      working <- list(stage = 'extrapolated_arithmetic', iteration = steps + 1L,
        x_prev = x, gradient_x_prev = current$gradient, y = y, tau = tau,
        L_previous = L_previous)
      if (!is.finite(tau) || tau <= 0 || any(!is.finite(y)))
        fail('Nonfinite tau or extrapolated coefficients')
      # Exact equality is required before reusing a cached point's gradient.
      if (identical(y, x)) {
        at_y <- current
        counts$y_cache_hits <- counts$y_cache_hits + 1L
      } else at_y <- evaluate(y, 'y')
      working$gradient_y <- at_y$gradient
      working$numerics_y <- at_y$numerics
      y_budget <- budget(at_y)
      Qy <- objective_at(y, at_y)
      L <- L_previous; accepted <- FALSE
      for (bt in seq.int(0L, ctl$max_backtracks)) {
        working$stage <- 'line_search'
        working$L <- L; working$backtracks <- bt
        if (!is.finite(L) || L <= 0) fail('Nonfinite inverse step', 'line_search_failure')
        step <- 1 / L
        if (!is.finite(step) || step <= 0) fail('Step underflow', 'line_search_failure')
        counts$trials <- counts$trials + 1L
        z <- aft_partial_prox(y, at_y$gradient, step, lambda, q)
        change <- z - y
        working$z <- z; working$change <- change
        if (any(!is.finite(change))) fail('Nonfinite candidate displacement')
        square <- sum(change^2)
        working$step_squared <- square
        if (!is.finite(square)) fail('Nonfinite step norm')
        if (all(change == 0)) {
          # The candidate is independently evaluated even if z equals a cache.
          working$stage <- 'zero_displacement_candidate'
          at_z <- evaluate(z, 'z')
          z_kkt <- kkt_at(z, at_z)
          working$gradient_z <- at_z$gradient
          working$numerics_z <- at_z$numerics
          working$kkt_z <- z_kkt
          if (z_kkt$residual_with_underflow_budget <= ctl$kkt_tol) {
            current_gm <- mapping_at(z, at_z, step)
            accepted_state <- list(x = z, x_prev = x, y = y, z = z,
              gradient_x_prev = current$gradient, gradient_y = at_y$gradient,
              gradient_z = at_z$gradient, gradient_x = at_z$gradient,
              numerics_x = at_z$numerics, tau = tau, L = L,
              terminal_candidate = 'independently stationary zero displacement',
              accepted_update_count_unchanged = TRUE)
            x <- z; current <- at_z; current_kkt <- z_kkt
            add_event('zero_displacement_stationary_candidate', working)
            return(finish('stationary_to_tolerance'))
          }
          fail('Zero displacement before full KKT tolerance', 'numerical_stagnation')
        }
        if (square == 0) fail('Nonzero displacement squared norm underflow', 'line_search_failure')
        counts$remainder <- counts$remainder + 1L
        bound <- aft_remainder_bound(Y, delta, design, x = y, d = change,
                                     h = h, block_size = ctl$block_size)
        remmeta <- attr(bound, 'numerics')
        margin <- ctl$omega * L * square / 2
        working$bound <- as.numeric(bound); working$margin <- margin
        working$remainder_numerics <- remmeta
        if (!is.finite(bound) || bound < 0 || !is.finite(margin))
          fail('Nonfinite or negative line-search certificate')
        if (margin == 0) fail('Quadratic acceptance margin underflow', 'line_search_failure')
        if (metadata_value(remmeta, 'term_underflow') > 0)
          add_event('remainder_term_underflow', list(iteration = steps + 1L,
            backtracks = bt, L = L, bound = as.numeric(bound), numerics = remmeta))
        if (bound <= margin) { accepted <- TRUE; break }
        if (bt < ctl$max_backtracks) {
          L <- L * ctl$L_growth
          counts$backtracks <- counts$backtracks + 1L
        }
      }
      if (!accepted) fail('No certified step within backtracking limit', 'line_search_failure')
      working$stage <- 'accepted_candidate_arithmetic'
      at_z <- evaluate(z, 'z')
      working$gradient_z <- at_z$gradient
      working$numerics_z <- at_z$numerics
      z_kkt <- kkt_at(z, at_z)
      z_gm <- mapping_at(z, at_z, step)
      Qprev <- objective_at(x, current); Qz <- objective_at(z, at_z)
      # Qz may exceed Qprev: FISTA imposes no monotonicity or restart condition.
      tau_next <- (1 + sqrt(1 + 4 * tau^2)) / 2
      y_next <- z + ((tau - 1) / tau_next) * (z - x)
      accepted_change <- z - x
      row <- c(iteration = steps + 1L, tau = tau, tau_next = tau_next, L = L,
        step = step, objective = Qz, objective_before = Qprev, objective_y = Qy,
        smooth_loss = at_z$loss, smooth_loss_y = at_y$loss, kkt = z_kkt$residual,
        kkt_with_budget = z_kkt$residual_with_underflow_budget,
        lost_gradient_upper_bound = z_kkt$lost_gradient_upper_bound,
        y_lost_gradient_upper_bound = y_budget, gradient_mapping_inf = z_gm$norm_inf,
        gradient_mapping_2 = z_gm$norm_2, step_norm = sqrt(square),
        accepted_change_norm = sqrt(sum(accepted_change^2)), step_squared = square,
        backtracks = bt, remainder_bound = as.numeric(bound), quadratic_margin = margin,
        majorant_defect = at_z$loss - at_y$loss - sum(at_y$gradient * change) - margin,
        objective_increase = as.numeric(Qz > Qprev))
      for (nm in loss_numeric_names) {
        row[paste0('y_', nm)] <- metadata_value(at_y$numerics, nm)
        row[paste0('z_', nm)] <- metadata_value(at_z$numerics, nm)
      }
      row['remainder_term_underflow'] <- metadata_value(remmeta, 'term_underflow')
      row['interval_certificate'] <- 0
      row[paste0('eval_', names(counts))] <- unlist(counts, use.names = FALSE)
      steps <- steps + 1L
      scalar_trace[steps, names(row)] <- row
      accepted_state <- list(iteration = steps, x = z, x_prev = x, y = y, z = z,
        gradient_x_prev = current$gradient, gradient_y = at_y$gradient,
        gradient_z = at_z$gradient, gradient_x = at_z$gradient,
        numerics_y = at_y$numerics, numerics_z = at_z$numerics, numerics_x = at_z$numerics,
        tau = tau, tau_next = tau_next, y_next = y_next, L = L, step = step,
        bound = as.numeric(bound), margin = margin, remainder_numerics = remmeta,
        kkt = z_kkt, automatic_resume = FALSE)
      observation <- if (!is.null(observer)) list(iteration = steps, k = steps,
        x_prev = x, y = y, z = z, grad_y = at_y$gradient, grad_z = at_z$gradient,
        tau = tau, tau_next = tau_next, y_next = y_next, L = L, step = step,
        backtracks = bt, bound = as.numeric(bound), margin = margin,
        objective_prev = Qprev, objective_y = Qy, objective_z = Qz, kkt = z_kkt,
        y_numerics = at_y$numerics, z_numerics = at_z$numerics,
        remainder_numerics = remmeta) else NULL
      # Commit the full point/arithmetic atomically before any user callback.
      x <- z; current <- at_z; current_kkt <- z_kkt; current_gm <- z_gm
      last_step <- step
      working$stage <- 'accepted_point_committed'
      if (!is.null(checkpoint) && steps %% ctl$checkpoint_every_updates == 0L) {
        working$stage <- 'checkpoint_callback'
        checkpoint(finish('checkpoint', 'Evidence snapshot; automatic resume disabled'))
      }
      if (!is.null(observer)) {
        working$stage <- 'observer_callback'
        observer(observation)
      }
      if (current_kkt$residual_with_underflow_budget <= ctl$kkt_tol)
        return(finish('stationary_to_tolerance'))
      if (steps >= ctl$max_iter)
        return(finish('iteration_limit', 'Maximum accepted updates reached'))
      working$stage <- 'momentum_update'
      if (!is.finite(tau_next) || any(!is.finite(y_next))) fail('Nonfinite FISTA extrapolation')
      y <- y_next; tau <- tau_next; L_previous <- L
    }
  }, warning = function(w) {
    add_event('warning', list(message = conditionMessage(w), call = deparse(conditionCall(w))))
    invokeRestart('muffleWarning')
  }), error = function(e) {
    status <- if (inherits(e, 'aft_fista_error')) e$status else 'numerical_failure'
    add_event('failure', list(status = status, message = conditionMessage(e),
                             class = class(e), last_attempt = working))
    finish(status, conditionMessage(e))
  })
}
