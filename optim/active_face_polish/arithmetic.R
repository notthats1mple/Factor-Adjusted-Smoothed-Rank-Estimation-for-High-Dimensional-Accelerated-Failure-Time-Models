# Pure arithmetic for the separately authorized v0.11 active-face experiment.
# Requires preserved R/preprocess.R and R/loss.R; no solver is sourced here.
.ap_stop <- function(status,message,details=NULL) {
  stop(structure(list(message=message,call=NULL,status=status,details=details),
    class=c('aft_active_face_error','error','condition')))
}
.ap_bound <- function(log_bound) {
  if(identical(log_bound,-Inf))return(0)
  if(!is.numeric(log_bound)||length(log_bound)!=1L||is.na(log_bound)||log_bound==Inf)
    .ap_stop('numerical_uncertainty','Invalid or missing logarithmic numerical bound')
  b<-exp(log_bound)
  if(!is.finite(b)) .ap_stop('numerical_uncertainty','Numerical bound overflowed')
  if(is.finite(log_bound)&&b==0)b<-.Machine$double.xmin*.Machine$double.eps
  b
}
.ap_logadd <- function(a,b) {
  if(identical(a,-Inf))return(b)
  if(identical(b,-Inf))return(a)
  max(a,b)+log1p(exp(-abs(a-b)))
}
.ap_gradient_budget <- function(metadata) {
  necessary<-c('lost_gradient_log_upper_bound','interval_certificate',
    'cdf_underflow','cdf_normalization_underflow','subnormal_gradient_weight_count')
  if(!is.list(metadata)||!all(necessary%in%names(metadata)))
    .ap_stop('numerical_uncertainty','Required full-gradient numerical metadata are missing',metadata)
  if(!is.logical(metadata$interval_certificate)||length(metadata$interval_certificate)!=1L||
     is.na(metadata$interval_certificate))
    .ap_stop('numerical_uncertainty','Invalid gradient interval-certificate metadata',metadata)
  for(nm in necessary[3:5]) {
    v<-metadata[[nm]]
    if(!is.numeric(v)||length(v)!=1L||!is.finite(v)||v<0||v!=floor(v))
      .ap_stop('numerical_uncertainty',paste('Invalid gradient metadata',nm),metadata)
  }
  .ap_bound(metadata$lost_gradient_log_upper_bound)
}
.ap_kkt <- function(x,g,q,lambda,E) {
  P<-length(x);beta<-if(q<P)seq.int(q+1L,P)else integer()
  r<-abs(g)
  nonzero<-beta[x[beta]!=0];zero<-beta[x[beta]==0]
  r[nonzero]<-abs(g[nonzero]+lambda*sign(x[nonzero]))
  r[zero]<-pmax(abs(g[zero])-lambda,0)
  k<-max(c(0,r))
  list(residual=k,coordinates=r,residual_with_underflow_budget=k+E,
    lost_gradient_upper_bound=E,original_coordinate_scale=TRUE)
}
.ap_saved_W <- function(design) {
  .aft_check_design(design)
  W<-if(design$kind=='explicit')design$W else
    cbind(design$Zc,design$Fhat,design$Xc-design$Fhat%*%t(design$Bhat))
  if(!is.matrix(W)||!is.numeric(W)||any(!is.finite(W))||
     nrow(W)!=design$n||ncol(W)!=design$P)
    .ap_stop('invalid_input','Saved design cannot be materialized without nonfinite values')
  W
}

# Two independently accumulated Hessians at an unchanged point. The direct
# route forms all pair-by-coordinate products in log space. The second route
# uses a weighted graph Laplacian. No penalty, ridge, or inverse is added.
aft_active_face_hessian <- function(Y,delta,design,x,h,A) {
  W<-.ap_saved_W(design);Wa<-W[,A,drop=FALSE]
  n<-nrow(W);m<-length(A);N<-as.double(n)*(n-1)
  e<-Y-as.numeric(W%*%x)
  if(any(!is.finite(e))) .ap_stop('numerical_failure','Nonfinite Hessian residual')
  gap<-outer(e,e,function(a,b)(b-a)/h)
  if(any(!is.finite(gap))) .ap_stop('numerical_failure','Nonfinite scaled Hessian residual')
  log_weight<-matrix(stats::dnorm(as.vector(gap),log=TRUE),n,n)-log(N)-log(h)
  weights<-exp(log_weight)*delta
  diag(weights)<-0
  adjacency<-weights+t(weights)
  laplacian<-diag(rowSums(adjacency),n)-adjacency
  H_laplacian<-crossprod(Wa,laplacian%*%Wa)
  H_pair<-matrix(0,m,m)
  pair_term_underflow<-0L;laplacian_weight_underflow<-0L
  pair_lost_log<--Inf;laplacian_lost_log<--Inf
  for(i in which(delta==1))for(j in seq_len(n)) {
    if(i==j)next
    w<-Wa[i,]-Wa[j,]
    signs<-outer(sign(w),sign(w))
    logs<-log_weight[i,j]+outer(log(abs(w)),log(abs(w)),'+')
    add<-signs*exp(logs);add[signs==0]<-0
    if(any(!is.finite(add))) .ap_stop('numerical_failure','Hessian pair term overflow')
    H_pair<-H_pair+add
    dropped<-signs!=0&add==0
    if(any(dropped)) {
      pair_term_underflow<-pair_term_underflow+sum(dropped)
      # A max-entry bound per pair, summed across omitted pairs.
      b<-max(logs[dropped]);b<-b+64*.Machine$double.eps*(1+abs(b))
      pair_lost_log<-.ap_logadd(pair_lost_log,b)
    }
    if(weights[i,j]==0&&any(w!=0)) {
      laplacian_weight_underflow<-laplacian_weight_underflow+1L
      b<-log_weight[i,j]+2*log(max(abs(w)))
      b<-b+64*.Machine$double.eps*(1+abs(b))
      laplacian_lost_log<-.ap_logadd(laplacian_lost_log,b)
    }
  }
  if(any(!is.finite(c(H_pair,H_laplacian)))) .ap_stop('numerical_failure','Nonfinite accumulated Hessian')
  pair_budget<-.ap_bound(pair_lost_log);laplacian_budget<-.ap_bound(laplacian_lost_log)
  difference<-max(c(0,abs(H_pair-H_laplacian)))
  tolerance<-1e-10+1e-12*max(c(0,abs(H_pair),abs(H_laplacian)))
  list(H=(H_pair+t(H_pair))/2,H_pair=H_pair,H_laplacian=H_laplacian,
    active_indices=A,raw_pair_symmetry_error=max(c(0,abs(H_pair-t(H_pair)))),
    raw_laplacian_symmetry_error=max(c(0,abs(H_laplacian-t(H_laplacian)))),
    double_path_max_abs_difference=difference,double_path_tolerance=tolerance,
    double_path_difference_plus_budgets=difference+pair_budget+laplacian_budget,
    double_path_pass=difference+pair_budget+laplacian_budget<=tolerance,
    numerics=list(pair_term_underflow=pair_term_underflow,
      laplacian_weight_underflow=laplacian_weight_underflow,
      pair_lost_entry_log_upper_bound=pair_lost_log,pair_lost_entry_upper_bound=pair_budget,
      laplacian_lost_entry_log_upper_bound=laplacian_lost_log,laplacian_lost_entry_upper_bound=laplacian_budget,
      core_vs_explicit_predictor_max_abs_difference=max(c(0,abs(as.numeric(aft_forward(design,x))-as.numeric(W%*%x)))),
      interval_certificate=FALSE))
}

# Audit only terms lost by the preserved whole-segment remainder calculation.
# Its returned value remains the principal bound. The omitted positive mass is
# added conservatively, never treated as an exact zero or an acceptance slack.
aft_active_face_remainder <- function(Y,delta,design,x,d,h,block_size=64) {
  original<-aft_remainder_bound(Y,delta,design,x,d,h,block_size)
  metadata<-attr(original,'numerics')
  if(!is.list(metadata)||is.null(metadata$term_underflow)||
     !is.numeric(metadata$term_underflow)||length(metadata$term_underflow)!=1L||
     !is.finite(metadata$term_underflow)||metadata$term_underflow<0||
     metadata$term_underflow!=floor(metadata$term_underflow)||
     !is.logical(metadata$interval_certificate)||length(metadata$interval_certificate)!=1L||
     is.na(metadata$interval_certificate))
    .ap_stop('numerical_uncertainty','Required remainder metadata are missing',metadata)
  lost_log<--Inf;observed<-0L
  if(metadata$term_underflow>0) {
    e<-.aft_loss_validate(Y,delta,design,x,h,block_size)
    advance<-as.numeric(aft_forward(design,d));n<-length(Y);N<-as.double(n)*(n-1)
    for(ii in .aft_blocks(delta,block_size)) {
      rr<-outer(e[ii],e,function(a,b)b-a)
      dd<-outer(advance[ii],advance,'-')
      next_r<-rr+dd
      minimum<-pmin(abs(rr),abs(next_r))
      minimum[(rr<=0&next_r>=0)|(rr>=0&next_r<=0)]<-0
      z<-minimum/h
      logs<-stats::dnorm(z,log=TRUE)-log(h)-log(2*N)+2*log(abs(dd))
      terms<-exp(logs);terms[dd==0]<-0
      terms[cbind(seq_along(ii),ii)]<-0
      dropped<-dd!=0&terms==0
      observed<-observed+sum(dropped)
      if(any(dropped)) {
        # Reducing z weakens the Gaussian upper bound. Small multiplicative
        # guards also cover floating evaluation of the logarithmic factors;
        # this is explicitly not directed interval arithmetic.
        safe_z<-pmin(z[dropped],1e100)*(1-8*.Machine$double.eps)
        bounds<--.5*safe_z^2-.5*log(2*pi)-log(h)-log(2*N)+2*log(abs(dd[dropped]))
        bounds<-bounds+64*.Machine$double.eps*(1+abs(bounds))
        if(any(!is.finite(bounds)))
          .ap_stop('numerical_uncertainty','Cannot control lost remainder terms',metadata)
        for(b in bounds)lost_log<-.ap_logadd(lost_log,b)
      }
    }
    if(observed!=metadata$term_underflow)
      .ap_stop('numerical_uncertainty','Remainder underflow audit count disagrees with original metadata',
        list(original=metadata,audited_count=observed))
  }
  lost<-.ap_bound(lost_log);bound<-as.numeric(original)+lost
  if(!is.finite(bound)||bound<0) .ap_stop('numerical_uncertainty','Invalid controlled remainder bound')
  list(bound=bound,original_bound=as.numeric(original),lost_upper_bound=lost,
    lost_log_upper_bound=lost_log,original_numerics=metadata,
    numerics=list(term_underflow=metadata$term_underflow,audited_underflow_count=observed,
      lost_upper_bound=lost,lost_log_upper_bound=lost_log,
      underflow_controlled=TRUE,interval_certificate=FALSE,
      bound_kind='Preserved whole-segment remainder plus explicit positive omitted-mass bound'))
}
