# Phase2c read-only endpoint comparisons. Source the preserved independent
# ordered-pair implementation into a private environment; no solver is loaded.
phase2c_arithmetic_load <- function(project) {
  e <- new.env(parent=baseenv())
  # The old independent helper uses stats:: explicitly and no project core.
  sys.source(file.path(project,'diagnostics/phase2b/endpoint.R'),envir=e)
  .phase2c_arithmetic_state$implementation <- e$phase2b_endpoint_arithmetic
  invisible(TRUE)
}
.phase2c_arithmetic_state <- new.env(parent=emptyenv())

phase2c_endpoint_arithmetic <- function(Y,delta,design,coefficients,h,lambda) {
  if (!exists('implementation',.phase2c_arithmetic_state,inherits=FALSE))
    stop('Call phase2c_arithmetic_load(project) before endpoint arithmetic')
  .phase2c_arithmetic_state$implementation(Y,delta,design,coefficients,h,lambda)
}

.phase2c_budget <- function(numerics) {
  z <- numerics$lost_gradient_log_upper_bound
  if (is.null(z)) return(0)
  if (identical(z,-Inf)) return(0)
  if (!is.numeric(z) || length(z)!=1L || is.na(z)) return(Inf)
  b <- exp(z)
  if (is.finite(z) && b==0) b <- .Machine$double.xmin*.Machine$double.eps
  b
}

.phase2c_comparison <- function(name,actual,expected,tolerance,budget=0) {
  valid <- is.numeric(actual) && is.numeric(expected) &&
    length(actual)==length(expected) && all(is.finite(c(actual,expected)))
  difference <- if (valid) max(c(0,abs(actual-expected))) else NA_real_
  list(name=name,passed=isTRUE(valid && is.finite(budget) && budget>=0 &&
    difference+budget<=tolerance),max_abs_difference=difference,
    absolute_tolerance=tolerance,tail_budget=budget,
    difference_plus_tail_budget=if(valid) difference+budget else NA_real_)
}

# input is the exact saved Y/delta/design/initial/h/lambda record. fit can be
# a core arithmetic record or the new solver return, including failure returns.
phase2c_endpoint_check <- function(input,fit,config,reference=NULL) {
  tryCatch({
    tol <- config$endpoint_crosscheck_tolerances
    gt <- tol$gradient_and_KKT_atol
    oa <- tol$objective_atol
    or <- tol$objective_rtol
    if (!identical(as.numeric(c(gt,oa,or)),c(1e-10,1e-10,1e-12)))
      stop('Endpoint thresholds differ from locked v0.10')
    independent <- phase2c_endpoint_arithmetic(input$Y,input$delta,input$design,
      fit$coefficients,input$h,input$lambda)
    checks <- list()
    compare <- function(name,a,b,t,budget=0)
      checks[[length(checks)+1L]] <<- .phase2c_comparison(name,a,b,t,budget)
    gb <- independent$numerics$gradient_tail_bound+.phase2c_budget(fit$numerics)
    lb <- independent$numerics$loss_tail_bound
    compare('core_full_gradient',independent$gradient,fit$gradient,gt,gb)
    compare('core_objective',independent$objective,fit$objective,oa+or*abs(fit$objective),lb)
    compare('core_full_KKT_scalar',independent$kkt$residual,fit$kkt$residual,gt,gb)
    compare('core_full_KKT_coordinates',independent$kkt$coordinates,fit$kkt$coordinates,gt,gb)
    if (!is.null(reference)) {
      rv <- function(x) as.numeric(unlist(x,use.names=FALSE))
      compare('locked_reference_objective',independent$objective,reference$objective,
        oa+or*abs(reference$objective),lb+reference$negative_tail_bound)
      compare('locked_reference_loss',independent$loss,reference$loss,
        oa+or*abs(reference$loss),lb+reference$negative_tail_bound)
      compare('locked_reference_full_gradient',independent$gradient,rv(reference$gradient),gt,
        independent$numerics$gradient_tail_bound)
      compare('locked_reference_full_KKT_scalar',independent$kkt$residual,reference$KKT,gt,
        independent$numerics$gradient_tail_bound)
      compare('locked_reference_full_KKT_coordinates',independent$kkt$coordinates,
        rv(reference$KKT_coordinates),gt,independent$numerics$gradient_tail_bound)
    }
    pass <- all(vapply(checks,function(x)isTRUE(x$passed),logical(1)))
    list(status=if(pass)'matched' else 'discrepancy',pass=pass,success=pass,
      independent=independent,comparison=list(checks=checks,thresholds=tol,
        budget_policy='Tail bounds consume locked thresholds; no additive relaxation'),
      core=list(coefficients=fit$coefficients,objective=fit$objective,gradient=fit$gradient,
        kkt=fit$kkt,numerics=fit$numerics),global_optimum_certified=FALSE)
  },error=function(e)list(status='arithmetic_failure',pass=FALSE,success=FALSE,
    error=list(message=conditionMessage(e),class=class(e)),global_optimum_certified=FALSE))
}
