# Section 5 eq:comp-prox, eq:comp-kkt, eq:comp-mapping, eq:comp-escape-lp.
aft_max_abs <- function(x) if (length(x)) max(abs(x)) else 0

aft_partial_prox <- function(x, gradient, step, lambda, q) {
  if (!is.numeric(x) || !is.numeric(gradient) || length(x) != length(gradient) ||
      any(!is.finite(c(x, gradient))) || length(q) != 1L || !is.finite(q) ||
      q != as.integer(q) || q < 0 || q > length(x) ||
      length(step) != 1L || !is.finite(step) || step <= 0 ||
      length(lambda) != 1L || !is.finite(lambda) || lambda <= 0)
    stop('Invalid partial proximal input')
  y <- x - step * gradient
  if (any(!is.finite(y)) || !is.finite(step * lambda)) stop('Nonfinite proximal candidate')
  b <- seq.int(q + 1L, length(x) + 1L)
  b <- b[b <= length(x)]
  y[b] <- sign(y[b]) * pmax(abs(y[b]) - step * lambda, 0)
  y
}

aft_kkt <- function(coef, gradient, q, lambda) {
  if (!is.numeric(coef) || !is.numeric(gradient) || length(coef) != length(gradient) ||
      any(!is.finite(c(coef, gradient))) || length(q) != 1L || !is.finite(q) ||
      q != floor(q) || q < 0 || q > length(coef) ||
      length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) stop('Invalid KKT input')
  ia <- seq_len(q); ib <- if (q < length(coef)) seq.int(q + 1L, length(coef)) else integer()
  rho <- pmax(abs(gradient[ib]) - lambda, 0)
  nz <- coef[ib] != 0
  rho[nz] <- abs(gradient[ib][nz] + lambda * sign(coef[ib][nz]))
  coordinates <- c(abs(gradient[ia]), rho)
  residual <- if (length(coordinates)) max(coordinates) else 0
  list(residual = residual, normalized = residual / max(lambda, aft_max_abs(gradient)),
       nuisance_gradient = gradient[ia], beta_residual = rho,
       coordinates = coordinates, exact_nonzero = nz,
       all_coordinates_checked = TRUE)
}

aft_gradient_mapping <- function(coef, gradient, q, lambda, step) {
  v <- (coef - aft_partial_prox(coef, gradient, step, lambda, q)) / step
  list(vector = v, norm_inf = aft_max_abs(v), norm_2 = sqrt(sum(v^2)), step = step)
}

aft_nuisance_diagnostic <- function(A, delta, lp_tolerance = 1e-9,
                                    rank_tolerance = 1e-10) {
  if (!is.matrix(A) || !is.numeric(A) || nrow(A) < 2L || any(!is.finite(A)) ||
      length(delta) != nrow(A) || any(!is.finite(delta)) || any(!delta %in% c(0, 1)) ||
      any(!is.finite(c(lp_tolerance, rank_tolerance))) ||
      lp_tolerance <= 0 || rank_tolerance <= 0)
    stop('Invalid nuisance diagnostic input')
  n <- nrow(A); q <- ncol(A); Eidx <- which(delta == 1); Cidx <- which(delta == 0)
  singular <- tryCatch({
    Ac <- if (q) sweep(A, 2, colMeans(A)) else A
    if (any(!is.finite(Ac))) stop('overflow when centering nuisance design')
    val <- if (q) svd(Ac, nu = 0, nv = 0)$d else numeric()
    if (any(!is.finite(val))) stop('nonfinite nuisance singular values')
    val
  }, error = identity)
  if (inherits(singular, 'error')) return(list(status = 'diagnostic_inconclusive',
    event_count = length(Eidx), flat_dimension = NA_integer_, flat_flag = NA,
    rank_uncertain = TRUE, lp_tolerance = lp_tolerance, rank_tolerance = rank_tolerance,
    lp_status = NA_integer_, reason = conditionMessage(singular), mathematical_certificate = FALSE))
  cutoff <- rank_tolerance * if (length(singular)) max(singular) else 0
  rank <- sum(singular > cutoff)
  ambiguous <- cutoff > 0 && any(singular > cutoff / 10 & singular < cutoff * 10)
  ans <- list(status = NULL, event_count = length(Eidx), flat_dimension = q - rank,
              flat_flag = q > rank, singular_values = singular, rank_cutoff = cutoff,
              rank_uncertain = ambiguous, lp_tolerance = lp_tolerance,
              rank_tolerance = rank_tolerance, lp_status = NA_integer_,
              lp_solver = 'lpSolve', mathematical_certificate = FALSE,
              direction = numeric(), max_event_residual = NA_real_,
              min_censored_margin = NA_real_, sum_margin = NA_real_)
  if (!length(Eidx)) { ans$status <- 'no_events'; return(ans) }
  if (!q || !length(Cidx)) {
    ans$status <- 'no_strict_escape'; ans$reason <- if (!q) 'empty nuisance block' else 'all patients are events'
    return(ans)
  }
  if (ambiguous) {
    ans$status <- 'diagnostic_inconclusive'; ans$reason <- 'nuisance rank near numerical cutoff'
    return(ans)
  }
  if (!requireNamespace('lpSolve', quietly = TRUE)) {
    ans$status <- 'diagnostic_inconclusive'; ans$reason <- 'lpSolve unavailable'; return(ans)
  }
  dif <- sweep(A, 2, A[Eidx[1], ], '-')
  if (any(!is.finite(dif))) {
    ans$status <- 'diagnostic_inconclusive'; ans$reason <- 'overflow in nuisance differences'; return(ans)
  }
  col_scale <- apply(abs(dif), 2, max)
  col_scale[col_scale == 0] <- 1
  dif_scaled <- sweep(dif, 2, col_scale, '/')
  E <- dif_scaled[Eidx[-1L], , drop = FALSE]
  C <- dif_scaled[Cidx, , drop = FALSE]
  constraints <- rbind(E, C, colSums(C))
  directions <- c(rep('=', nrow(E)), rep('>=', nrow(C) + 1L))
  rhs <- c(rep(0, nrow(E) + nrow(C)), 1)
  # v = v_plus-v_minus makes each original coordinate unrestricted.
  lp <- tryCatch(lpSolve::lp('min', rep(1, 2L * q), cbind(constraints, -constraints),
                             directions, rhs), error = identity)
  ans$column_scale <- col_scale
  if (inherits(lp, 'error')) {
    ans$status <- 'diagnostic_inconclusive'; ans$reason <- conditionMessage(lp); return(ans)
  }
  ans$lp_status <- lp$status
  if (lp$status == 2L) {
    ans$status <- 'no_strict_escape_numeric'; ans$reason <- 'LP reported infeasible in floating point'
    return(ans)
  }
  if (lp$status != 0L || any(!is.finite(lp$solution))) {
    ans$status <- 'diagnostic_inconclusive'; ans$reason <- 'LP neither feasible nor infeasible'; return(ans)
  }
  v <- (lp$solution[seq_len(q)] - lp$solution[q + seq_len(q)]) / col_scale
  Em <- dif[Eidx[-1L], , drop = FALSE] %*% v
  Cm <- as.vector(dif[Cidx, , drop = FALSE] %*% v)
  ans$direction <- as.vector(v)
  ans$max_event_residual <- aft_max_abs(Em)
  ans$min_censored_margin <- min(Cm)
  ans$sum_margin <- sum(Cm)
  valid <- all(is.finite(c(v, Em, Cm, ans$sum_margin))) &&
    ans$max_event_residual <= lp_tolerance &&
    ans$min_censored_margin >= -lp_tolerance && ans$sum_margin >= 1 - lp_tolerance
  ans$status <- if (valid) 'strict_escape_numeric' else 'diagnostic_inconclusive'
  ans$reason <- if (valid) 'LP witness checked in original coordinates' else 'LP witness failed original-coordinate residual checks'
  ans
}
