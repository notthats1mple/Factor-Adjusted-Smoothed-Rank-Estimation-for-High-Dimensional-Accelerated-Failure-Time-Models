# Independent Phase2b endpoint arithmetic. This file calls no core or pilot
# function, performs no estimation/preprocessing, and never accesses RNG state.
# It materializes the saved design, then accumulates ordered pairs directly.
.phase2b_endpoint_stop <- function(message) {
  stop(structure(list(message=message,call=NULL),
    class=c('phase2b_endpoint_error','error','condition')))
}
.phase2b_endpoint_finite <- function(x) is.numeric(x) && all(is.finite(x))
.phase2b_endpoint_logadd <- function(x,y) {
  if (is.infinite(x) && x<0) return(y)
  if (is.infinite(y) && y<0) return(x)
  max(x,y)+log1p(exp(-abs(x-y)))
}
.phase2b_endpoint_bound <- function(log_bound) {
  if (identical(log_bound,-Inf)) return(0)
  ans <- exp(log_bound)
  # A positive mathematical bound that rounds to zero must remain disclosed.
  if (is.finite(log_bound) && ans==0) ans <- .Machine$double.xmin*.Machine$double.eps
  ans
}
.phase2b_endpoint_design <- function(design) {
  if (!is.list(design)) .phase2b_endpoint_stop('Saved design must be a list')
  if (identical(design$kind,'explicit')) {
    W <- design$W
    q <- design$q
  } else if (identical(design$kind,'pca')) {
    matrices <- design[c('Zc','Fhat','Xc','Bhat')]
    if (length(matrices)!=4L || !all(vapply(matrices,function(x)
        is.matrix(x) && .phase2b_endpoint_finite(x),logical(1))))
      .phase2b_endpoint_stop('Saved PCA design contains missing or nonfinite matrices')
    n <- nrow(design$Xc)
    p <- ncol(design$Xc)
    K <- ncol(design$Fhat)
    d <- ncol(design$Zc)
    if (nrow(design$Zc)!=n || nrow(design$Fhat)!=n ||
        nrow(design$Bhat)!=p || ncol(design$Bhat)!=K ||
        !identical(as.numeric(c(design$n,design$p,design$d,design$K,design$q,design$P)),
                   as.numeric(c(n,p,d,K,d+K,d+K+p))))
      .phase2b_endpoint_stop('Saved PCA design dimensions or coefficient partition disagree')
    # This is a materialization of saved factors/loadings, never a fresh PCA.
    U <- design$Xc - design$Fhat %*% t(design$Bhat)
    W <- cbind(design$Zc,design$Fhat,U)
    q <- d+K
  } else .phase2b_endpoint_stop('Unsupported saved design kind')
  if (!is.matrix(W) || !.phase2b_endpoint_finite(W) ||
      !is.numeric(q) || length(q)!=1L || !is.finite(q) || q!=as.integer(q) ||
      q<0 || q>ncol(W)) .phase2b_endpoint_stop('Invalid explicit design or nuisance partition')
  list(W=W,q=as.integer(q))
}

phase2b_endpoint_arithmetic <- function(Y,delta,design,coefficients,h,lambda) {
  explicit <- .phase2b_endpoint_design(design)
  W <- explicit$W
  n <- nrow(W)
  P <- ncol(W)
  q <- explicit$q
  if (n<2L || !.phase2b_endpoint_finite(Y) || length(Y)!=n ||
      !.phase2b_endpoint_finite(delta) || length(delta)!=n || any(!delta %in% c(0,1)) ||
      !.phase2b_endpoint_finite(coefficients) || length(coefficients)!=P ||
      !.phase2b_endpoint_finite(h) || length(h)!=1L || h<=0 ||
      !.phase2b_endpoint_finite(lambda) || length(lambda)!=1L || lambda<0)
    .phase2b_endpoint_stop('Invalid response, event vector, coefficients, h, or lambda')
  fitted <- as.numeric(W %*% coefficients)
  residual <- Y-fitted
  if (any(!is.finite(residual))) .phase2b_endpoint_stop('Nonfinite saved-design fitted score or residual')
  normalization <- as.double(n)*(n-1)
  log_normalization <- log(normalization)
  log_smallest_positive <- log(.Machine$double.xmin)+log(.Machine$double.eps)
  loss <- 0
  gradient <- numeric(P)
  tail_log_bound <- -Inf
  gradient_log_bound <- -Inf
  omitted_loss_pairs <- 0L
  omitted_cdf_pairs <- 0L
  rounded_one_pairs <- 0L
  rounded_one_log_bound <- -Inf
  pairs <- 0L
  for (i in seq_len(n)) {
    if (delta[[i]]==0) next
    for (j in seq_len(n)) {
      if (i==j) next
      pairs <- pairs+1L
      r <- residual[[j]]-residual[[i]]
      if (!is.finite(r)) .phase2b_endpoint_stop('Ordered residual difference overflow')
      z <- r/h
      log_phi <- stats::dnorm(z,log=TRUE)
      probability <- stats::pnorm(z)
      difference <- W[i,]-W[j,]
      if (any(!is.finite(difference))) .phase2b_endpoint_stop('Ordered design difference overflow')
      # For r=-h*a<0, 0 <= smooth(r,h)=h*(phi(a)-a*Phi(-a)) <= h*phi(a).
      # Discard only z<-8 and retain the explicit nonnegative absolute budget.
      # The locked tolerance is not increased by this budget in the comparison.
      if (z < -8) {
        omitted_loss_pairs <- omitted_loss_pairs+1L
        pair_tail_log_bound <- if (is.infinite(z)) log_smallest_positive else
          log(h)+log_phi-log_normalization
        # If finite r/h overflows, its Gaussian tail times any finite h is
        # smaller than the least positive double; retain that conservative cap.
        tail_log_bound <- .phase2b_endpoint_logadd(tail_log_bound,pair_tail_log_bound)
        pair_loss <- 0
      } else if (is.infinite(z) && z>0) {
        pair_loss <- r/normalization
      } else {
        pair_loss <- (r*probability+h*exp(log_phi))/normalization
        if (!is.finite(pair_loss) || pair_loss<0)
          .phase2b_endpoint_stop('Central smoothed loss became negative or nonfinite')
      }
      loss <- loss+pair_loss
      # Multiply the CDF by each coordinate difference before normalization.
      # Explicit pair-vector accumulation is separate from the core row-score
      # adjoint route. Extreme underflow is conservatively bounded below.
      weighted <- difference*probability/normalization
      if (any(!is.finite(weighted))) .phase2b_endpoint_stop('Pair gradient overflow')
      gradient <- gradient+weighted
      dropped <- difference!=0 & weighted==0
      if (any(dropped)) {
        omitted_cdf_pairs <- omitted_cdf_pairs+1L
        # Mills upper bound Phi(z)<=phi(z)/(-z), with the max coordinate
        # difference; sum of these bounds controls the gradient infinity norm.
        pair_gradient_log_bound <- log(max(abs(difference[dropped])))-log_normalization
        if (is.finite(z) && z<0)
          pair_gradient_log_bound <- pair_gradient_log_bound+min(0,log_phi-log(-z))
        else if (is.infinite(z) && z<0) pair_gradient_log_bound <- log_smallest_positive
        gradient_log_bound <- .phase2b_endpoint_logadd(gradient_log_bound,pair_gradient_log_bound)
      }
      if (probability==1 && z>0 && length(difference) && any(difference!=0)) {
        rounded_one_pairs <- rounded_one_pairs+1L
        pair_one_log_bound <- if (is.infinite(z)) log_smallest_positive else
          log(max(abs(difference)))+log_phi-log(z)-log_normalization
        rounded_one_log_bound <- .phase2b_endpoint_logadd(rounded_one_log_bound,
          pair_one_log_bound)
      }
    }
  }
  beta_indices <- if (q<P) seq.int(q+1L,P) else integer()
  beta <- coefficients[beta_indices]
  objective <- loss+lambda*sum(abs(beta))
  if (any(!is.finite(c(loss,objective,gradient))))
    .phase2b_endpoint_stop('Nonfinite accumulated loss, objective, or gradient')
  coordinate_residuals <- abs(gradient)
  gb <- gradient[beta_indices]
  active <- beta!=0
  beta_residuals <- numeric(length(beta))
  beta_residuals[active] <- abs(gb[active]+lambda*sign(beta[active]))
  beta_residuals[!active] <- pmax(abs(gb[!active])-lambda,0)
  coordinate_residuals[beta_indices] <- beta_residuals
  total_gradient_log_bound <- .phase2b_endpoint_logadd(gradient_log_bound,rounded_one_log_bound)
  list(loss=loss,objective=objective,gradient=gradient,fitted_scores=fitted,
    kkt=list(residual=max(c(0,coordinate_residuals)),coordinates=coordinate_residuals,
      nuisance_gradient=gradient[seq_len(q)],beta_residual=beta_residuals,exact_nonzero=active),
    dimensions=list(n=n,P=P,q=q,p=P-q),
    numerics=list(route='Explicit saved-design matrix, sequential event-weighted ordered pairs',
      loss_tail_cutoff=-8,loss_tail_bound=.phase2b_endpoint_bound(tail_log_bound),
      loss_tail_log_bound=tail_log_bound,omitted_loss_pairs=omitted_loss_pairs,
      gradient_omission_bound=.phase2b_endpoint_bound(gradient_log_bound),
      gradient_omission_log_bound=gradient_log_bound,omitted_cdf_pairs=omitted_cdf_pairs,
      cdf_rounded_one_bound=.phase2b_endpoint_bound(rounded_one_log_bound),
      cdf_rounded_one_log_bound=rounded_one_log_bound,cdf_rounded_one_pairs=rounded_one_pairs,
      gradient_tail_bound=.phase2b_endpoint_bound(total_gradient_log_bound),
      gradient_tail_log_bound=total_gradient_log_bound,event_ordered_pairs=pairs,
      normalization=normalization,interval_certificate=FALSE,
      bounds_scope='Mathematical Gaussian tail bounds evaluated in floating point; excludes complete roundoff enclosure',
      no_core_calls=TRUE,no_solver_calls=TRUE,no_rng_calls=TRUE))
}

phase2b_endpoint <- function(fit,Y,delta,config,reference=NULL,reference_source=NULL) {
  # Errors are retained as structured failures so the runner can halt the
  # affected continuation while saving the exact endpoint evidence.
  tryCatch({
    tol <- config$numeric_reference_tolerances
    gt <- tol$gradient_max_abs
    oa <- tol$objective_abs_plus_rel$atol
    or <- tol$objective_abs_plus_rel$rtol
    if (!identical(as.numeric(c(gt,oa,or)),c(1e-10,1e-10,1e-12)))
      .phase2b_endpoint_stop('Comparison tolerances must equal the locked v0.9 values')
    if (!is.null(reference) && (!is.character(reference_source) ||
        length(reference_source)!=1L || !nzchar(reference_source)))
      .phase2b_endpoint_stop('Provided reference requires its actual source path')
    independent <- phase2b_endpoint_arithmetic(Y,delta,fit$design,fit$coefficients,fit$h,fit$lambda)
    checks <- list()
    compare <- function(name,actual,expected,limit,budget=0) {
      valid <- .phase2b_endpoint_finite(actual) && .phase2b_endpoint_finite(expected) &&
        length(actual)==length(expected)
      difference <- if (valid) max(c(0,abs(actual-expected))) else NA_real_
      passed <- valid && is.finite(limit) && is.finite(budget) && difference+budget<=limit
      checks[[length(checks)+1L]] <<- list(name=name,passed=isTRUE(passed),
        max_abs_difference=difference,absolute_tolerance=limit,tail_budget=budget,
        difference_plus_tail_budget=if(valid) difference+budget else NA_real_)
    }
    objective_limit <- oa+or*abs(fit$objective)
    compare('saved_core_full_gradient',independent$gradient,fit$gradient,gt,
      independent$numerics$gradient_tail_bound)
    compare('saved_core_objective',independent$objective,fit$objective,objective_limit,
      independent$numerics$loss_tail_bound)
    compare('saved_core_full_KKT_scalar',independent$kkt$residual,fit$kkt$residual,gt,
      independent$numerics$gradient_tail_bound)
    compare('saved_core_full_KKT_coordinates',independent$kkt$coordinates,fit$kkt$coordinates,gt,
      independent$numerics$gradient_tail_bound)
    if (!is.null(reference)) {
      compare('provided_reference_objective',independent$objective,reference$independent_objective,
        oa+or*abs(reference$independent_objective),independent$numerics$loss_tail_bound)
      compare('provided_reference_KKT_scalar',independent$kkt$residual,reference$independent_kkt,gt,
        independent$numerics$gradient_tail_bound)
      reference_record <- list(status='compared',source=reference_source,record=reference,
        full_gradient_comparison='not_available_in_provided_reference',
        coefficient_provenance='Reference used decimal JSON; this recomputation uses authoritative RDS coefficients')
    } else reference_record <- list(status='not_requested',full_gradient_comparison='not_requested')
    pass <- all(vapply(checks,function(x)isTRUE(x$passed),logical(1)))
    list(status=if(pass)'matched' else 'discrepancy',success=pass,pass=pass,
      independent=independent,comparison=list(pass=pass,checks=checks,
        thresholds=tol,budget_policy='Tail bounds consume the locked tolerance; no tolerance increase'),
      saved_core=list(objective=fit$objective,gradient=fit$gradient,kkt=fit$kkt,
        status=fit$status,success=fit$success,coefficient_source='supplied fit RDS object'),
      reference=reference_record,global_optimum_certified=FALSE)
  },error=function(e)list(status='arithmetic_failure',success=FALSE,pass=FALSE,
    independent=NULL,comparison=list(pass=FALSE,checks=list()),
    error=list(message=conditionMessage(e),class=class(e)),global_optimum_certified=FALSE))
}
