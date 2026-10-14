# Fixed Gaussian-smoothed, event-anchored Gehan loss. Original coordinates.
# Section 3: smoothed-gehan/loss-derivatives; Section 5: comp-gradient,
# comp-positive-certificate, comp-stable-loss; Appendix E: comp-derivatives/tail.

.aft_numeric_fault <- function(message) {
  stop(structure(list(message = message, call = NULL),
                 class = c("aft_numeric_error", "error", "condition")))
}

.aft_require_finite <- function(x, what) {
  if (any(!is.finite(x))) .aft_numeric_fault(paste("nonfinite", what))
  invisible(x)
}

.aft_h <- function(h) {
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0)
    stop("h must be finite and positive", call. = FALSE)
}

# For details=TRUE, the values retain r's dimensions and the metadata record
# underflow of positive mathematical quantities. The asymptotic remainder is
# an exact-arithmetic truncation bound, not a floating-point interval bound.
aft_smooth <- function(r, h, details = FALSE) {
  .aft_h(h)
  if (!is.numeric(r)) stop("r must be numeric", call. = FALSE)
  .aft_require_finite(r, "residual difference")
  shape <- dim(r)
  x <- as.numeric(r)
  z <- abs(x) / h
  .aft_require_finite(z, "scaled residual")
  central <- z <= 12
  J <- rep(NA_real_, length(x))
  remainder <- numeric(length(x))
  nterms <- integer(length(x))
  J[z == 0] <- 1
  use <- central & z > 0
  # -expm1 avoids subtraction after the normal-tail ratio has been computed.
  zz <- z[use]
  J[use] <- -expm1(log(zz) + pnorm(-zz, log.p = TRUE) -
                      dnorm(zz, log = TRUE))
  for (i in which(!central)) {
    inv <- (1 / z[i])^2
    # When inv itself underflows, work directly on the leading log term.
    if (inv == 0) {
      J[i] <- 0
      nterms[i] <- 1L
      next
    }
    term <- inv
    total <- term
    for (k in seq_len(200L)) {
      nxt <- -term * (2 * k + 1) * inv
      if (abs(nxt) <= 2e-16 * abs(total)) break
      if (abs(nxt) >= abs(term))
        .aft_numeric_fault("tail expansion did not attain its remainder tolerance")
      total <- total + nxt
      term <- nxt
      if (k == 200L) .aft_numeric_fault("tail expansion iteration limit")
    }
    J[i] <- total
    remainder[i] <- abs(nxt)
    nterms[i] <- k
  }
  if (any(J < 0 | !is.finite(J)) || any(J[central] <= 0))
    .aft_numeric_fault("nonpositive or nonfinite tail correction")
  logJ <- log(J)
  tinyJ <- !central & J == 0
  logJ[tinyJ] <- -2 * log(z[tinyJ])
  loggap <- log(h) + dnorm(z, log = TRUE) + logJ
  gap <- exp(loggap)
  value <- pmax(x, 0) + gap
  .aft_require_finite(value, "smoothed loss (overflow)")
  if (any(value < 0)) .aft_numeric_fault("negative smoothed loss")
  if (!is.null(shape)) {
    dim(value) <- shape
    dim(gap) <- shape
    dim(J) <- shape
    dim(remainder) <- shape
    dim(nterms) <- shape
  }
  if (!details) return(value)
  list(value = value, gap = gap, J = J, tail_remainder = remainder,
       tail_terms = nterms,
       numerics = list(gap_underflow = sum(gap == 0),
                       positive_gap_absorbed = sum(x > 0 & as.numeric(gap) > 0 &
                                                    as.numeric(value) == x),
                       tail_J_underflow = sum(tinyJ),
                       asymptotic_count = sum(!central),
                       max_tail_relative_remainder = if (any(!central & J > 0))
                         max(remainder[!central & J > 0] / J[!central & J > 0]) else 0,
                       central_switch = 12,
                       interval_certificate = FALSE))
}

.aft_loss_validate <- function(Y, delta, design, coef, h, block_size) {
  .aft_h(h)
  if (!is.list(design) || length(design$n) != 1L || design$n < 2 ||
      length(design$P) != 1L || design$P < 0)
    stop("invalid design dimensions", call. = FALSE)
  if (!is.numeric(Y) || length(Y) != design$n ||
      !is.numeric(delta) || length(delta) != design$n ||
      any(!is.finite(delta)) || any(!delta %in% c(0, 1)))
    stop("invalid sample dimensions or event code", call. = FALSE)
  if (!is.numeric(coef) || length(coef) != design$P)
    stop("invalid coefficient dimensions", call. = FALSE)
  if (!is.numeric(block_size) || length(block_size) != 1L ||
      !is.finite(block_size) || block_size < 1 || block_size != floor(block_size))
    stop("block_size must be a positive integer", call. = FALSE)
  .aft_require_finite(Y, "outcome")
  .aft_require_finite(coef, "coefficient")
  prediction <- as.numeric(aft_forward(design, coef))
  if (length(prediction) != length(Y)) stop("invalid forward product dimensions")
  .aft_require_finite(prediction, "predictor")
  residual <- Y - prediction
  .aft_require_finite(residual, "residual")
  residual
}

.aft_blocks <- function(delta, block_size) {
  ii <- which(delta == 1)
  if (!length(ii)) return(list())
  split(ii, ceiling(seq_along(ii) / block_size))
}

.aft_kernel_metadata <- function() {
  list(gap_underflow = 0L, positive_gap_absorbed = 0L,
       tail_J_underflow = 0L, cdf_underflow = 0L, cdf_rounded_one = 0L,
       cdf_normalization_underflow = 0L, subnormal_gradient_weight_count = 0L,
       lost_gradient_log_upper_bound = -Inf,
       lost_gradient_bound_kind = "Gaussian Mills and design max-norm upper bound; floating-point evaluation",
       asymptotic_count = 0L, max_tail_relative_remainder = 0,
       interval_certificate = FALSE,
       summation = "R sum/rowSums/colSums; normalized per pair before accumulation")
}

.aft_logsumexp <- function(x) {
  if (!length(x)) return(-Inf)
  largest <- max(x)
  if (largest == -Inf) return(-Inf)
  largest + log(sum(exp(x-largest)))
}

# Column scans keep scratch storage O(n+P), even for a large explicit W.
.aft_matrix_maxabs <- function(x) {
  ans <- 0
  for (j in seq_len(ncol(x))) ans <- max(ans,max(abs(x[,j])))
  ans
}

.aft_design_log_maxabs <- function(design) {
  if (design$P == 0) return(-Inf)
  if (design$kind == "explicit") return(log(.aft_matrix_maxabs(design$W)))
  if (design$kind != "pca") .aft_numeric_fault("unknown design for underflow bound")
  xlog <- log(.aft_matrix_maxabs(design$Xc))
  zlog <- log(.aft_matrix_maxabs(design$Zc))
  flog <- log(.aft_matrix_maxabs(design$Fhat))
  # |U_ij| <= max|Xc| + K max|F| max|B|, in the original scales.
  blog <- log(.aft_matrix_maxabs(design$Bhat))
  residual_log <- .aft_logsumexp(c(xlog,log(design$K)+flog+blog))
  max(zlog,flog,residual_log)
}

# Bound the infinity norm of score error from unrepresentable/subnormal pair
# weights before multiplying by potentially large W. In exact arithmetic,
# Phi(-z) <= phi(z)/z (z>0), and ||Wi-Wj||inf <= 2 max|W|. For nonzero
# subnormals we conservatively bound exact + computed contributions; this also
# covers their poor relative precision. Logs avoid underflow of the bound itself.
.aft_lost_score_log_bound <- function(scaled_residual, computed_weight, N,
                                      design_log_max) {
  if (design_log_max == -Inf || !length(scaled_residual)) return(-Inf)
  # Capping and reducing z weakens the bound; it prevents squared overflow and
  # leaves headroom for the floating division used to form the residual ratio.
  z <- pmax(0,pmin(-scaled_residual,1e100)*(1-8*.Machine$double.eps))
  log_probability_upper <- numeric(length(z)) # universal Phi <= 1
  tail <- z > 1
  log_probability_upper[tail] <- -.5*z[tail]^2-.5*log(2*pi)-log(z[tail])
  log_mass <- .aft_logsumexp(c(log_probability_upper-log(N),log(computed_weight)))
  bound <- log(2)+design_log_max+log_mass
  # Numerical evaluation of an analytic bound, not directed interval arithmetic.
  bound + 64*.Machine$double.eps*(1+abs(bound))
}

# O(block_size*n+n+P) scratch, without any pair-by-feature matrix.
aft_loss <- function(Y, delta, design, coef, h, block_size = 64,
                     gradient = TRUE) {
  e <- .aft_loss_validate(Y, delta, design, coef, h, block_size)
  n <- length(Y)
  N <- as.double(n) * (n - 1)
  cvec <- numeric(n)
  total <- 0
  metadata <- .aft_kernel_metadata()
  design_log_max <- NULL
  for (ii in .aft_blocks(delta, block_size)) {
    rr <- outer(e[ii], e, function(a, b) b - a)
    .aft_require_finite(rr, "pair residual")
    diagonal <- cbind(seq_along(ii), ii)
    scalar <- aft_smooth(rr, h, details = TRUE)
    vals <- scalar$value
    vals[diagonal] <- 0
    total <- total + sum(vals / N)
    keep <- matrix(TRUE, length(ii), n)
    keep[diagonal] <- FALSE
    metadata$gap_underflow <- metadata$gap_underflow + sum(scalar$gap[keep] == 0)
    metadata$positive_gap_absorbed <- metadata$positive_gap_absorbed +
      sum(rr[keep] > 0 & scalar$gap[keep] > 0 & scalar$value[keep] == rr[keep])
    metadata$tail_J_underflow <- metadata$tail_J_underflow + sum(scalar$J[keep] == 0)
    metadata$asymptotic_count <- metadata$asymptotic_count + sum(scalar$tail_terms[keep] > 0)
    metadata$max_tail_relative_remainder <- max(metadata$max_tail_relative_remainder,
                                               scalar$numerics$max_tail_relative_remainder)
    if (gradient) {
      weights <- pnorm(rr / h)
      metadata$cdf_underflow <- metadata$cdf_underflow + sum(weights[keep] == 0)
      metadata$cdf_rounded_one <- metadata$cdf_rounded_one + sum(weights[keep] == 1)
      weights[diagonal] <- 0
      normalized <- weights / N
      metadata$cdf_normalization_underflow <- metadata$cdf_normalization_underflow +
        sum(keep & weights > 0 & normalized == 0)
      affected <- keep & normalized < .Machine$double.xmin
      if (any(affected)) {
        if (is.null(design_log_max)) design_log_max <- .aft_design_log_maxabs(design)
        metadata$subnormal_gradient_weight_count <-
          metadata$subnormal_gradient_weight_count + sum(affected)
        contribution_bound <- .aft_lost_score_log_bound(rr[affected]/h,
                                    normalized[affected],N,design_log_max)
        metadata$lost_gradient_log_upper_bound <- .aft_logsumexp(c(
          metadata$lost_gradient_log_upper_bound,contribution_bound))
      }
      weights <- normalized
      cvec[ii] <- cvec[ii] + rowSums(weights)
      cvec <- cvec - colSums(weights)
    }
  }
  grad <- if (gradient) as.numeric(aft_adjoint(design, cvec)) else NULL
  .aft_require_finite(total, "total loss (overflow)")
  if (gradient) .aft_require_finite(grad, "gradient")
  list(loss = total, gradient = grad, numerics = metadata,
       event_count = sum(delta), pair_count = sum(delta) * (n - 1),
       normalization = N, score_weights = if (gradient) cvec else NULL)
}

# Product exp(log_weight)*v, preserving signs and avoiding premature density
# underflow or overflow. Metadata counts true nonzero products lost to underflow.
.aft_log_product <- function(log_weight, v) {
  out <- numeric(length(v))
  nz <- v != 0
  out[nz] <- sign(v[nz]) * exp(log_weight[nz] + log(abs(v[nz])))
  .aft_require_finite(out, "weighted derivative (overflow)")
  attr(out, "underflow") <- sum(nz & out == 0)
  out
}

# Small-sample reference only: explicit ordered-pair loop and full Hessian.
aft_loss_explicit <- function(Y, delta, W, coef, h) {
  if (!is.matrix(W) || !is.numeric(W)) stop("W must be a numeric matrix")
  .aft_require_finite(W, "design")
  design <- list(n = nrow(W), P = ncol(W))
  # Validate independently of the implicit design representation.
  if (nrow(W) < 2 || length(Y) != nrow(W) || length(delta) != nrow(W) ||
      length(coef) != ncol(W) || any(!is.finite(delta)) || any(!delta %in% c(0, 1)))
    stop("invalid explicit dimensions or event code")
  .aft_h(h)
  .aft_require_finite(Y, "outcome")
  .aft_require_finite(coef, "coefficient")
  e <- Y - as.numeric(W %*% coef)
  .aft_require_finite(e, "residual")
  n <- nrow(W)
  P <- ncol(W)
  N <- as.double(n) * (n - 1)
  f <- 0
  g <- numeric(P)
  H <- matrix(0, P, P)
  underflow <- 0L
  for (i in which(delta == 1)) for (j in seq_len(n)) {
    if (i == j) next
    w <- W[i, ] - W[j, ]
    r <- e[j] - e[i]
    f <- f + aft_smooth(r, h) / N
    z <- r / h
    g <- g + (pnorm(z) / N) * w
    if (P) {
      log_weight <- dnorm(z, log = TRUE) - log(h) - log(N)
      # Form products in log space so tiny density and large differences do
      # not produce spurious zero or infinity before multiplication.
      signs <- outer(sign(w), sign(w))
      logs <- outer(log(abs(w)), log(abs(w)), "+")
      add <- signs * exp(log_weight + logs)
      add[signs == 0] <- 0
      .aft_require_finite(add, "Hessian term (overflow)")
      underflow <- underflow + sum(signs != 0 & add == 0)
      H <- H + add
    }
  }
  .aft_require_finite(c(f, g, H), "explicit loss/derivatives")
  list(loss = f, gradient = g, hessian = H,
       numerics = list(hessian_term_underflow = underflow))
}

aft_hvp <- function(Y, delta, design, coef, h, v, block_size = 64) {
  e <- .aft_loss_validate(Y, delta, design, coef, h, block_size)
  if (!is.numeric(v) || length(v) != design$P) stop("invalid HVP direction")
  .aft_require_finite(v, "HVP direction")
  advance <- as.numeric(aft_forward(design, v))
  .aft_require_finite(advance, "HVP predictor")
  n <- length(Y)
  N <- as.double(n) * (n - 1)
  b <- numeric(n)
  underflow <- 0L
  for (ii in .aft_blocks(delta, block_size)) {
    rr <- outer(e[ii], e, function(a, b) b - a)
    dd <- outer(advance[ii], advance, "-")
    .aft_require_finite(c(rr, dd), "HVP pair difference")
    z <- rr / h
    .aft_require_finite(z, "scaled residual")
    log_weight <- dnorm(z, log = TRUE) - log(h) - log(N)
    R <- .aft_log_product(as.numeric(log_weight), as.numeric(dd))
    underflow <- underflow + attr(R, "underflow")
    dim(R) <- dim(rr)
    R[cbind(seq_along(ii), ii)] <- 0
    b[ii] <- b[ii] + rowSums(R)
    b <- b - colSums(R)
  }
  ans <- as.numeric(aft_adjoint(design, b))
  .aft_require_finite(ans, "Hessian-vector product")
  attr(ans, "numerics") <- list(product_underflow = underflow)
  ans
}

# d denotes y-x. This upper bound covers the entire segment, not one endpoint.
# Scalar result carries attr(result, "numerics") for reporting, including
# underflow. It is a floating-point evaluation, not interval certification.
aft_remainder_bound <- function(Y, delta, design, x, d, h, block_size = 64) {
  e <- .aft_loss_validate(Y, delta, design, x, h, block_size)
  if (!is.numeric(d) || length(d) != design$P) stop("invalid increment")
  .aft_require_finite(d, "increment")
  advance <- as.numeric(aft_forward(design, d))
  .aft_require_finite(advance, "increment predictor")
  n <- length(Y)
  N <- as.double(n) * (n - 1)
  total <- 0
  underflow <- 0L
  for (ii in .aft_blocks(delta, block_size)) {
    rr <- outer(e[ii], e, function(a, b) b - a)
    dd <- outer(advance[ii], advance, "-")
    next_r <- rr + dd
    .aft_require_finite(c(rr, dd, next_r), "line-segment pair difference")
    minimum <- pmin(abs(rr), abs(next_r))
    crossing <- (rr <= 0 & next_r >= 0) | (rr >= 0 & next_r <= 0)
    minimum[crossing] <- 0
    z <- minimum / h
    .aft_require_finite(z, "scaled segment residual")
    log_terms <- dnorm(z, log = TRUE) - log(h) - log(2 * N) + 2 * log(abs(dd))
    terms <- exp(log_terms)
    terms[dd == 0] <- 0
    terms[cbind(seq_along(ii), ii)] <- 0
    .aft_require_finite(terms, "remainder bound term (overflow)")
    underflow <- underflow + sum(dd != 0 & terms == 0)
    total <- total + sum(terms)
  }
  .aft_require_finite(total, "remainder bound (overflow)")
  attr(total, "numerics") <- list(term_underflow = underflow,
                                   interval_certificate = FALSE)
  total
}
