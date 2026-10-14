# Fixed training transformations from Sections 2--3 and (comp-implicit-W).
# No residual rescaling, outcome inputs, or changes of coefficient coordinates.

.aft_preprocess_stop <- function(code, message, diagnostics = NULL) {
  stop(structure(list(message = message, call = NULL, code = code,
                      diagnostics = diagnostics),
                 class = c("aft_preprocess_error", "error", "condition")))
}

.aft_finite_matrix <- function(x, name, n = NULL) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x))
    .aft_preprocess_stop("invalid_input", paste(name, "must be a numeric matrix"))
  if (!is.null(n) && nrow(x) != n)
    .aft_preprocess_stop("invalid_input", paste(name, "has the wrong row count"))
  if (any(!is.finite(x)))
    .aft_preprocess_stop("nonfinite_input", paste(name, "contains a non-finite value"))
  storage.mode(x) <- "double"
  x
}

.aft_finite_vector <- function(x, size, name) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || length(x) != size)
    .aft_preprocess_stop("invalid_input", paste(name, "must be a numeric vector of length", size))
  if (any(!is.finite(x)))
    .aft_preprocess_stop("nonfinite_input", paste(name, "contains a non-finite value"))
  as.numeric(x)
}

.aft_indices <- function(start, size) if (size) seq.int(start, length.out = size) else integer(0)

# Canonical basis of a tied left singular subspace. Project coordinate axes in
# row order, twice orthogonalize, and make the first largest entry positive.
# This depends only on the subspace projector (and the specified row order).
.aft_canonical_subspace <- function(V) {
  n <- nrow(V)
  k <- ncol(V)
  out <- matrix(0, n, k)
  found <- 0L
  tolerance <- 128 * .Machine$double.eps * max(n, k)
  for (j in seq_len(n)) {
    v <- as.vector(V %*% V[j, ])
    if (found) {
      Q <- out[, seq_len(found), drop = FALSE]
      v <- v - as.vector(Q %*% crossprod(Q, v))
      v <- v - as.vector(Q %*% crossprod(Q, v))
    }
    length_v <- sqrt(sum(v * v))
    if (length_v > tolerance) {
      v <- v / length_v
      # Treat last-bit ties as ties so an exact mathematical tie does not pick
      # a different sign merely because one equal entry rounded upward.
      largest <- max(abs(v))
      anchor <- which(abs(v) >= largest - 16 * .Machine$double.eps * max(n, k) * largest)[1L]
      if (v[anchor] < 0) v <- -v
      found <- found + 1L
      out[, found] <- v
      if (found == k) break
    }
  }
  if (found != k)
    .aft_preprocess_stop("rank_uncertain", "Cannot resolve a deterministic basis of the tied PCA subspace")
  out
}

# X and Z retain their supplied coordinate scales. K is prespecified.
aft_preprocess <- function(X, Z = NULL, K = 0) {
  X <- .aft_finite_matrix(X, "X")
  n <- nrow(X)
  p <- ncol(X)
  if (n < 2L) .aft_preprocess_stop("invalid_input", "n must be at least 2")
  if (is.null(Z)) Z <- matrix(numeric(0), n, 0L)
  Z <- .aft_finite_matrix(Z, "Z", n)
  d <- ncol(Z)
  if (length(K) != 1L || !is.numeric(K) || !is.finite(K) || K < 0 || K != floor(K))
    .aft_preprocess_stop("invalid_input", "K must be one nonnegative integer")
  if (K > min(n - 1L, p))
    .aft_preprocess_stop("rank_deficient", "K exceeds the maximum possible centered rank")
  K <- as.integer(K)
  x_mean <- colMeans(X)
  z_mean <- colMeans(Z)
  Xc <- sweep(X, 2L, x_mean, "-")
  Zc <- sweep(Z, 2L, z_mean, "-")
  if (any(!is.finite(Xc)) || any(!is.finite(Zc)))
    .aft_preprocess_stop("nonfinite_computation", "Centering overflowed")
  Fhat <- matrix(numeric(0), n, 0L)
  Bhat <- matrix(numeric(0), p, 0L)
  Lambdahat <- matrix(numeric(0), 0L, 0L)
  numeric_rank <- NA_integer_
  rank_tolerance <- 0
  tie_tolerance <- 0
  selected_eigenvalues <- numeric(0)
  tie_groups <- list()
  pca_diagnostics <- list(status = "not_applicable", tolerance = 1e-8,
                          reason = "K=0; no factor mapping")
  if (K) {
    # Scaling here protects the decomposition only; it does not scale the model.
    magnitude <- max(abs(Xc))
    if (magnitude == 0)
      .aft_preprocess_stop("rank_deficient", "Centered X has zero rank")
    decomposition <- tryCatch(svd(Xc / magnitude, nu = min(n, p), nv = min(n, p)),
      error = function(e) .aft_preprocess_stop("pca_failure", conditionMessage(e)))
    singular <- decomposition$d
    rank_tolerance <- 64 * .Machine$double.eps * max(n, p) * singular[1L]
    numeric_rank <- sum(singular > rank_tolerance)
    if (K > numeric_rank)
      .aft_preprocess_stop("rank_uncertain", "Requested PCA direction is zero or cannot be distinguished from numerical rank error")
    vectors <- decomposition$u
    right <- decomposition$v
    # Never rotate unequal computed singular values. In particular, a global
    # absolute tolerance must not merge weak directions of different scales.
    # An exactly equal computed group uses the full subspace before selecting K;
    # the right vectors receive the same rotation, preserving the SVD relation.
    first <- 1L
    while (first <= numeric_rank) {
      last <- first
      while (last < numeric_rank && singular[last + 1L] == singular[first])
        last <- last + 1L
      ids <- seq.int(first, last)
      if (length(ids) > 1L) {
        old <- vectors[, ids, drop = FALSE]
        canonical <- .aft_canonical_subspace(old)
        rotation <- crossprod(old, canonical)
        vectors[, ids] <- canonical
        right[, ids] <- right[, ids, drop = FALSE] %*% rotation
        tie_groups[[length(tie_groups) + 1L]] <- ids
      } else {
        v <- vectors[, first]
        largest <- max(abs(v))
        anchor <- which(abs(v) >= largest -
          16 * .Machine$double.eps * max(n, p) * largest)[1L]
        if (v[anchor] < 0) {
          vectors[, first] <- -v
          right[, first] <- -right[, first]
        }
      }
      first <- last + 1L
    }
    Fhat <- sqrt(n) * vectors[, seq_len(K), drop = FALSE]
    # Xc = U D V' implies Xc' Fhat/n = V D/sqrt(n). The latter avoids
    # cancelling strong columns while constructing a weak loading direction.
    loading_scales <- singular[seq_len(K)] * (magnitude / sqrt(n))
    Bhat <- sweep(right[, seq_len(K), drop = FALSE], 2L, loading_scales, "*")
    selected_eigenvalues <- loading_scales * loading_scales
    if (any(!is.finite(Bhat)) || any(!is.finite(selected_eigenvalues)) || any(selected_eigenvalues <= 0))
      .aft_preprocess_stop("nonfinite_computation", "PCA loading/eigenvalue computation overflowed or underflowed")
    Lambdahat <- diag(selected_eigenvalues, nrow = K, ncol = K)
    # Check every selected direction on its own scale. A small absolute Gram
    # defect alone cannot rule out leakage outside a selected K-dimensional
    # subspace. Both multiplication orders used by the public frozen map are
    # checked, along with the loading identity and normalized Gram matrix.
    column_max <- function(M) apply(abs(M), 2L, max)
    raw_loadings <- crossprod(Xc, Fhat) / n
    loading_relative <- column_max(sweep(raw_loadings - Bhat, 2L, loading_scales, "/"))
    projected <- Xc %*% Bhat
    factor_replay <- sweep(projected, 2L, selected_eigenvalues, "/")
    frozen_factors <- Xc %*% sweep(Bhat, 2L, selected_eigenvalues, "/")
    direction_scale <- 1 + column_max(Fhat)
    spectral_relative <- column_max(sweep(projected - sweep(Fhat, 2L, selected_eigenvalues, "*"),
                                          2L, selected_eigenvalues, "/")) / direction_scale
    factor_replay_relative <- column_max(factor_replay - Fhat) / direction_scale
    prediction_replay_relative <- column_max(frozen_factors - Fhat) / direction_scale
    # Divide each loading column before taking the Gram to avoid lambda_i *
    # lambda_j overflow/underflow in a square-root normalization.
    normalized_gram <- crossprod(sweep(Bhat, 2L, loading_scales, "/")) - diag(K)
    normalized_gram_defect <- max(abs(normalized_gram))
    factor_normalization_defect <- max(abs(crossprod(Fhat) / n - diag(K)))
    all_errors <- c(loading_relative, spectral_relative, factor_replay_relative,
                    prediction_replay_relative, normalized_gram_defect,
                    factor_normalization_defect)
    pca_diagnostics <- list(status = "checked", tolerance = 1e-8,
      scaled_singular_values = singular, exact_computed_tie_groups = tie_groups,
      loading_relative_by_direction = loading_relative,
      spectral_relative_by_direction = spectral_relative,
      factor_replay_relative_by_direction = factor_replay_relative,
      prediction_replay_relative_by_direction = prediction_replay_relative,
      factor_replay_absolute_by_direction = column_max(factor_replay - Fhat),
      prediction_replay_absolute_by_direction = column_max(frozen_factors - Fhat),
      normalized_gram = normalized_gram,
      normalized_gram_defect = normalized_gram_defect,
      factor_normalization_defect = factor_normalization_defect,
      selected_eigenvalues = selected_eigenvalues,
      loading_construction = "V D/sqrt(n), equivalent to Xc'Fhat/n",
      interval_certificate = FALSE)
    if (any(!is.finite(all_errors)) || any(all_errors > pca_diagnostics$tolerance)) {
      pca_diagnostics$status <- "numerical_uncertainty"
      .aft_preprocess_stop("pca_uncertain",
        "Selected PCA directions do not satisfy the directional loading and frozen-map accuracy checks",
        pca_diagnostics)
    }
  }
  structure(list(kind = "pca", n = n, p = p, d = d, K = K, q = d + K,
    P = p + d + K, Xc = Xc, Zc = Zc, Fhat = Fhat, Bhat = Bhat,
    Lambdahat = Lambdahat, x_mean = x_mean, z_mean = z_mean,
    numeric_rank = numeric_rank, rank_tolerance = rank_tolerance,
    tie_tolerance = tie_tolerance, selected_eigenvalues = selected_eigenvalues,
    pca_diagnostics = pca_diagnostics,
    scale = rep(1, p), pca_convention = "SVD; retain unequal singular directions; exact computed ties rotate both singular bases; largest-entry-positive signs"),
    class = "aft_design")
}

# A low-level fixture/solver interface; W is used exactly as supplied.
aft_design_matrix <- function(W, q = 0) {
  W <- .aft_finite_matrix(W, "W")
  if (nrow(W) < 2L) .aft_preprocess_stop("invalid_input", "n must be at least 2")
  if (length(q) != 1L || !is.numeric(q) || !is.finite(q) || q < 0 || q != floor(q) || q > ncol(W))
    .aft_preprocess_stop("invalid_input", "q must be an integer between 0 and ncol(W)")
  q <- as.integer(q)
  structure(list(kind = "explicit", W = W, n = nrow(W), P = ncol(W),
                 p = ncol(W) - q, q = q, d = q, K = 0L), class = "aft_design")
}

.aft_check_design <- function(design) {
  if (!inherits(design, "aft_design") || !design$kind %in% c("explicit", "pca"))
    .aft_preprocess_stop("invalid_input", "design must be an aft_design object")
}

aft_forward <- function(design, v) {
  .aft_check_design(design)
  v <- .aft_finite_vector(v, design$P, "v")
  if (design$kind == "explicit") result <- as.vector(design$W %*% v) else {
    theta <- v[.aft_indices(1L, design$d)]
    gamma <- v[.aft_indices(design$d + 1L, design$K)]
    beta <- v[.aft_indices(design$q + 1L, design$p)]
    result <- as.vector(design$Zc %*% theta + design$Xc %*% beta +
      design$Fhat %*% (gamma - crossprod(design$Bhat, beta)))
  }
  if (any(!is.finite(result))) .aft_preprocess_stop("nonfinite_computation", "Design forward product overflowed")
  result
}

aft_adjoint <- function(design, c) {
  .aft_check_design(design)
  c <- .aft_finite_vector(c, design$n, "c")
  if (design$kind == "explicit") result <- as.vector(crossprod(design$W, c)) else {
    factor_product <- crossprod(design$Fhat, c)
    result <- as.vector(c(crossprod(design$Zc, c), factor_product,
      crossprod(design$Xc, c) - design$Bhat %*% factor_product))
  }
  if (any(!is.finite(result))) .aft_preprocess_stop("nonfinite_computation", "Design adjoint product overflowed")
  result
}

aft_explicit_design <- function(design) {
  .aft_check_design(design)
  if (design$kind == "explicit") return(design$W)
  cbind(design$Zc, design$Fhat, design$Xc - tcrossprod(design$Fhat, design$Bhat))
}

# This dense deployment slope is distinct from the sparse beta estimand.
aft_prediction_coefficients <- function(design, coef) {
  .aft_check_design(design)
  if (design$kind != "pca")
    .aft_preprocess_stop("invalid_input", "Prediction mapping requires a training preprocessing object")
  coef <- .aft_finite_vector(coef, design$P, "coef")
  beta <- coef[.aft_indices(design$q + 1L, design$p)]
  gamma <- coef[.aft_indices(design$d + 1L, design$K)]
  result <- beta
  if (design$K) result <- beta + as.vector(design$Bhat %*%
    ((gamma - as.vector(crossprod(design$Bhat, beta))) / diag(design$Lambdahat)))
  if (any(!is.finite(result))) .aft_preprocess_stop("nonfinite_computation", "Prediction coefficients overflowed")
  result
}

aft_predict_map <- function(design, coef, Xnew, Znew = NULL) {
  slope <- aft_prediction_coefficients(design, coef)
  Xnew <- .aft_finite_matrix(Xnew, "Xnew")
  if (ncol(Xnew) != design$p)
    .aft_preprocess_stop("invalid_input", "Xnew has the wrong column count")
  if (is.null(Znew)) {
    if (design$d) .aft_preprocess_stop("invalid_input", "Znew is required for a model with clinical covariates")
    Znew <- matrix(numeric(0), nrow(Xnew), 0L)
  }
  Znew <- .aft_finite_matrix(Znew, "Znew", nrow(Xnew))
  if (ncol(Znew) != design$d)
    .aft_preprocess_stop("invalid_input", "Znew has the wrong column count")
  result <- as.vector(sweep(Xnew, 2L, design$x_mean, "-") %*% slope +
    sweep(Znew, 2L, design$z_mean, "-") %*% coef[.aft_indices(1L, design$d)])
  if (any(!is.finite(result))) .aft_preprocess_stop("nonfinite_computation", "Prediction mapping overflowed")
  result
}
