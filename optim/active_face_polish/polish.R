# Bounded local Newton polishing on one exact, fixed active sign face.
# This function never calls the default solver, FISTA, PCA, or random routines.
aft_active_face_control <- function(control=list()) {
  defaults<-list(kkt_tol=1e-7,max_newton_updates=10L,max_backtrack_halvings=30L,
    armijo_constant=1e-4,backtrack_factor=.5,fraction_to_sign_boundary=.99,
    hessian_rank_factor=100,linear_system_relative_backward_error=1e-10,
    gradient_budget_fraction=.1,independent_endpoint_atol=1e-10,
    independent_objective_rtol=1e-12,block_size=64L,
    full_coordinate_KKT_every_accepted_point=TRUE,active_set_fixed_within_polish=TRUE,
    momentum=FALSE,extra_ridge=FALSE,automatic_retry=FALSE,automatic_fallback=FALSE)
  if(!is.list(control)||is.null(names(control))&&length(control)||
     anyDuplicated(names(control))||any(!names(control)%in%names(defaults)))
    stop('Invalid or unknown active-face control fields')
  out<-modifyList(defaults,control,keep.null=TRUE)
  integers<-c('max_newton_updates','max_backtrack_halvings','block_size')
  for(nm in integers) {
    z<-out[[nm]]
    if(!is.numeric(z)||length(z)!=1L||!is.finite(z)||z!=floor(z)||z<0)
      stop(paste('Invalid integer control',nm))
    out[[nm]]<-as.integer(z)
  }
  if(out$max_newton_updates>20L||out$max_backtrack_halvings>30L||out$block_size<1L)
    stop('Control exceeds bounded experiment limits')
  if(!is.numeric(out$kkt_tol)||length(out$kkt_tol)!=1L||!is.finite(out$kkt_tol)||out$kkt_tol<=0)
    stop('kkt_tol must be finite and positive')
  fixed<-setdiff(names(defaults),c('kkt_tol',integers))
  for(nm in fixed)if(!isTRUE(all.equal(out[[nm]],defaults[[nm]],check.attributes=FALSE,tolerance=0)))
    stop(paste('The fixed mathematical contract does not permit changing',nm))
  out
}

aft_polish_active_face <- function(Y,delta,design,initial,h,lambda,control=list(),checkpoint=NULL) {
  ctl<-aft_active_face_control(control)
  if(!is.null(checkpoint)&&!is.function(checkpoint))stop('checkpoint must be NULL or a function')
  .aft_loss_validate(Y,delta,design,initial,h,ctl$block_size)
  if(!is.numeric(lambda)||length(lambda)!=1L||!is.finite(lambda)||lambda<0)
    stop('lambda must be finite and nonnegative')
  P<-design$P;q<-design$q
  if(!is.numeric(q)||length(q)!=1L||!is.finite(q)||q!=floor(q)||q<0||q>P)
    stop('Invalid saved nuisance partition')
  beta<-if(q<P)seq.int(q+1L,P)else integer()
  nuisance<-if(q)seq_len(q)else integer()
  active_beta<-beta[initial[beta]!=0];Z<-beta[initial[beta]==0]
  A<-sort(c(nuisance,active_beta));original_sign<-sign(initial[active_beta])
  x<-as.numeric(initial);current<-NULL;initial_point<-NULL
  points<-list();accepted_points<-list();directions<-list();trials<-list();warnings<-list()
  updates<-0L;last_attempt<-NULL;failure<-NULL
  counts<-list(objective=0L,gradient=0L,kkt=0L,initial_gradient=0L,candidate_gradient=0L,
    hessian=0L,hessian_pair=0L,hessian_laplacian=0L,hessian_eigen=0L,
    cholesky=0L,linear_solve=0L,remainder=0L,remainder_underflow_audit=0L,
    trials=0L,backtrack_halvings=0L,accepted_updates=0L,checkpoints=0L)
  increment<-function(nm,n=1L)counts[[nm]]<<-counts[[nm]]+n
  make_point<-function(at,iteration,role) {
    increment('objective');increment('gradient');increment('kkt')
    increment(if(role=='initial')'initial_gradient'else'candidate_gradient')
    ev<-aft_loss(Y,delta,design,at,h,ctl$block_size,gradient=TRUE)
    last_attempt<<-list(stage='point_arithmetic',role=role,iteration=iteration,
      coefficients=at,raw_arithmetic=ev)
    E<-.ap_gradient_budget(ev$numerics)
    kk<-.ap_kkt(at,ev$gradient,q,lambda,E)
    obj<-ev$loss+lambda*sum(abs(at[beta]))
    if(any(!is.finite(c(obj,ev$gradient,kk$coordinates))))
      .ap_stop('numerical_failure','Nonfinite complete objective, gradient, or KKT')
    list(iteration=iteration,role=role,coefficients=at,objective=obj,smooth_loss=ev$loss,
      gradient=ev$gradient,kkt=kk,numerics=ev$numerics,gradient_budget=E,
      inactive_applicable=length(Z)>0L,
      inactive_margin=if(length(Z))min(lambda-abs(ev$gradient[Z]))-E else NA_real_,
      inactive_max_abs_gradient=if(length(Z))max(abs(ev$gradient[Z]))else NA_real_,
      inactive_condition=if(length(Z))max(abs(ev$gradient[Z]))+E<=lambda else TRUE,
      active_sign_preserved=all(sign(at[active_beta])==original_sign)&all(at[active_beta]!=0),
      initial_zero_coordinates_preserved=all(at[Z]==0),evaluations=counts)
  }
  budget_ok<-function(point) {
    if(point$gradient_budget>ctl$gradient_budget_fraction*ctl$kkt_tol)
      .ap_stop('numerical_uncertainty','Full-gradient error budget exceeds kkt_tol/10',point$numerics)
    invisible(TRUE)
  }
  finish<-function(status,success=FALSE,message=NULL) {
    available<-!is.null(current)
    structure(list(solver_id='active_face_newton_polish_v0.11_fixed_sign_face',
      status=status,success=success,global_optimum_certified=FALSE,message=message,
      coefficients=if(available)current$coefficients else x,
      initial_coefficients=as.numeric(initial),
      objective=if(available)current$objective else NA_real_,
      smooth_loss=if(available)current$smooth_loss else NA_real_,
      gradient=if(available)current$gradient else NULL,kkt=if(available)current$kkt else NULL,
      numerics=if(available)current$numerics else NULL,
      gradient_budget=if(available)current$gradient_budget else NA_real_,
      iterations=updates,accepted_updates=updates,evaluations=counts,
      initial=initial_point,final=current,points=points,accepted_points=accepted_points,
      directions=directions,trials=trials,last_attempt=last_attempt,failure=failure,warnings=warnings,
      active_indices=A,inactive_indices=Z,active_beta_indices=active_beta,
      original_active_beta_sign=original_sign,design=design,h=h,lambda=lambda,q=q,
      control=ctl,automatic_retry=FALSE,automatic_fallback=FALSE,
      interval_certificate=FALSE,
      fixed_face='All initial nonzero beta and all nuisance coordinates; exact initial beta zeros fixed.',
      coefficient_scale='Original saved design and original penalty scale'),
      class=c('aft_active_face_polish_fit','list'))
  }
  tryCatch(withCallingHandlers({
    current<-make_point(x,0L,'initial');initial_point<-current;points<-list(current)
    budget_ok(current)
    if(sum(delta)==0)return(finish('no_events',FALSE,'No event-anchored loss terms'))
    if(current$kkt$residual_with_underflow_budget<=ctl$kkt_tol)
      return(finish('stationary_to_tolerance',TRUE,'Initial point already satisfies full KKT plus budget'))
    if(length(Z)&&!(current$inactive_margin>0))
      return(finish('active_face_not_ready',FALSE,'Initial exact-zero beta conditions do not have strictly positive margin'))
    for(iteration in seq_len(ctl$max_newton_updates)) {
      direction_index<-length(directions)+1L
      dir<-list(iteration=iteration,point_before=current,active_indices=A,
        stage='Hessian evaluation',accepted=FALSE)
      directions[[direction_index]]<-dir;last_attempt<-dir
      increment('hessian');increment('hessian_pair');increment('hessian_laplacian')
      hs<-aft_active_face_hessian(Y,delta,design,x,h,A)
      dir$hessian<-hs;dir$H<-hs$H
      directions[[direction_index]]<-dir;last_attempt<-dir
      if(!isTRUE(hs$double_path_pass)) .ap_stop('hessian_crosscheck_failure','Independent Hessian routes disagree',hs)
      increment('hessian_eigen')
      eigenvalues<-eigen(hs$H,symmetric=TRUE,only.values=TRUE)$values
      spectral_tolerance<-ctl$hessian_rank_factor*.Machine$double.eps*length(A)*max(abs(eigenvalues))
      dir$eigenvalues<-eigenvalues;dir$spectral_tolerance<-spectral_tolerance
      dir$numerical_rank<-sum(abs(eigenvalues)>spectral_tolerance)
      dir$stage<-'Hessian positive-definiteness check'
      directions[[direction_index]]<-dir;last_attempt<-dir
      if(!length(eigenvalues)||min(eigenvalues)<=spectral_tolerance)
        .ap_stop('active_hessian_not_positive_definite','Active Hessian minimum eigenvalue is not above the fixed rank threshold',dir)
      increment('cholesky')
      R<-tryCatch(chol(hs$H),error=function(e)
        .ap_stop('active_hessian_not_positive_definite',paste('Cholesky failed:',conditionMessage(e)),dir))
      sA<-current$gradient[A]
      penalized_positions<-match(active_beta,A)
      sA[penalized_positions]<-sA[penalized_positions]+lambda*original_sign
      dir$cholesky<-R;dir$signed_residual<-sA;dir$stage<-'Cholesky linear solve'
      directions[[direction_index]]<-dir;last_attempt<-dir
      increment('linear_solve')
      dA<-backsolve(R,forwardsolve(t(R),-sA))
      d<-numeric(P);d[A]<-as.numeric(dA)
      residual<-as.numeric(hs$H%*%dA+sA)
      numerator<-max(c(0,abs(residual)))
      denominator<-max(rowSums(abs(hs$H)))*max(c(0,abs(dA)))+max(c(0,abs(sA)))
      backward_error<-if(denominator==0)if(numerator==0)0 else Inf else numerator/denominator
      v<-sum(sA*dA);direction_l1<-sum(abs(d));robust_v<-v+current$gradient_budget*direction_l1
      dir$cholesky<-R;dir$signed_residual<-sA;dir$d<-d;dir$d_A<-as.numeric(dA)
      dir$linear_system_residual<-residual;dir$relative_backward_error<-backward_error
      dir$linear_system_backward_error_numerator<-numerator
      dir$linear_system_backward_error_denominator<-denominator
      dir$v<-v;dir$predicted_linear_decrease<--v
      dir$gradient_budget<-current$gradient_budget;dir$direction_l1<-direction_l1
      dir$robust_directional_derivative<-robust_v;dir$stage<-'direction certification'
      directions[[direction_index]]<-dir;last_attempt<-dir
      if(any(!is.finite(c(d,residual,v,robust_v,backward_error))))
        .ap_stop('numerical_failure','Nonfinite Cholesky direction or certificate',dir)
      if(backward_error>ctl$linear_system_relative_backward_error)
        .ap_stop('linear_system_inaccurate','Relative backward error exceeds 1e-10',dir)
      if(!(v<0&&robust_v<0)) .ap_stop('non_descent_direction','Direction lacks strictly negative budget-aware derivative',dir)
      toward_boundary<-active_beta[d[active_beta]!=0&sign(x[active_beta])!=sign(d[active_beta])]
      boundary<-if(length(toward_boundary))min(-x[toward_boundary]/d[toward_boundary])else Inf
      initial_step<-min(1,ctl$fraction_to_sign_boundary*boundary)
      dir$sign_boundary<-boundary;dir$initial_step<-initial_step
      directions[[direction_index]]<-dir;last_attempt<-dir
      if(!is.finite(initial_step)||initial_step<=0)
        .ap_stop('face_or_line_search_failure','No positive representable sign-preserving step',dir)
      accepted<-FALSE;trial_indices<-integer()
      for(halvings in seq.int(0L,ctl$max_backtrack_halvings)) {
        if(halvings>0L)increment('backtrack_halvings')
        t<-initial_step*ctl$backtrack_factor^halvings
        increment('trials')
        trial_index<-length(trials)+1L;trial_indices<-c(trial_indices,trial_index)
        candidate<-x+t*d
        trial<-list(iteration=iteration,direction_index=direction_index,trial_index=trial_index,
          step=t,halvings=halvings,candidate=candidate,accepted=FALSE,
          sign_preserved=all(sign(candidate[active_beta])==original_sign)&all(candidate[active_beta]!=0),
          zeros_preserved=all(candidate[Z]==0),reason='checking candidate',
          objective_difference=NA_real_,objective_difference_available=FALSE)
        trials[[trial_index]]<-trial;last_attempt<-trial
        if(!all(is.finite(candidate))||t<=0) {
          trial$reason<-'nonfinite candidate or nonpositive step'
          trials[[trial_index]]<-trial;last_attempt<-trial;next
        }
        if(!trial$sign_preserved||!trial$zeros_preserved) {
          trial$reason<-'original signs or exact-zero coordinates changed'
          trials[[trial_index]]<-trial;last_attempt<-trial;next
        }
        if(identical(as.numeric(candidate),as.numeric(x))) {
          trial$reason<-'candidate rounded to the current point'
          trials[[trial_index]]<-trial;last_attempt<-trial;next
        }
        increment('remainder')
        rem<-aft_active_face_remainder(Y,delta,design,x,t*d,h,ctl$block_size)
        if(rem$original_numerics$term_underflow>0)increment('remainder_underflow_audit')
        lhs<-t*v+rem$bound+t*current$gradient_budget*direction_l1
        rhs<-ctl$armijo_constant*t*v
        trial$remainder<-rem$bound;trial$remainder_numerics<-rem$numerics
        trial$remainder_original<-rem$original_bound
        trial$remainder_original_numerics<-rem$original_numerics
        trial$armijo_lhs<-lhs;trial$armijo_rhs<-rhs;trial$armijo_pass<-is.finite(lhs)&&lhs<=rhs
        trial$reason<-'whole-segment Armijo check'
        trials[[trial_index]]<-trial;last_attempt<-trial
        if(!trial$armijo_pass) {
          trial$reason<-'whole-segment Armijo sufficient condition failed; candidate arithmetic not evaluated'
          trials[[trial_index]]<-trial;last_attempt<-trial;next
        }
        candidate_point<-make_point(candidate,iteration,'candidate')
        trial$candidate_point<-candidate_point
        trial$objective_difference<-candidate_point$objective-current$objective
        trial$objective_difference_available<-TRUE
        trial$inactive_condition<-candidate_point$inactive_condition
        trial$reason<-'complete candidate arithmetic and inactive-coordinate check'
        trials[[trial_index]]<-trial;last_attempt<-trial
        budget_ok(candidate_point)
        if(!candidate_point$inactive_condition) {
          trial$reason<-'candidate exact-zero beta gradient condition failed'
          trials[[trial_index]]<-trial;last_attempt<-trial;next
        }
        accepted<-TRUE;updates<-updates+1L;increment('accepted_updates')
        candidate_point$role<-'accepted';candidate_point$iteration<-updates
        candidate_point$step<-t;candidate_point$halvings<-halvings
        candidate_point$objective_difference<-trial$objective_difference
        candidate_point$evaluations<-counts
        x<-candidate;current<-candidate_point
        points[[length(points)+1L]]<-current;accepted_points[[length(accepted_points)+1L]]<-current
        trial$candidate_point<-current;trial$accepted<-TRUE;trial$reason<-'accepted'
        trials[[trial_index]]<-trial;last_attempt<-trial
        dir$accepted<-TRUE;dir$accepted_trial_index<-trial_index;dir$trial_indices<-trial_indices
        dir$accepted_step<-t;dir$stage<-'accepted';directions[[direction_index]]<-dir
        if(!is.null(checkpoint)) {
          increment('checkpoints')
          tryCatch(checkpoint(list(iteration=updates,accepted_update=updates,point=current,
            direction=dir,trials=trials[trial_indices],evaluations=counts)),error=function(e)
              .ap_stop('checkpoint_failure',paste('Accepted point checkpoint failed:',conditionMessage(e))))
        }
        break
      }
      if(!accepted) {
        dir$trial_indices<-trial_indices;dir$stage<-'line search exhausted'
        directions[[direction_index]]<-dir
        return(finish('face_or_line_search_failure',FALSE,'No accepted candidate after the allowed halvings'))
      }
      if(current$kkt$residual_with_underflow_budget<=ctl$kkt_tol)
        return(finish('stationary_to_tolerance',TRUE,'Accepted point satisfies complete original-coordinate KKT plus budget'))
    }
    finish('polish_iteration_limit',FALSE,'Newton update budget exhausted on the fixed sign face')
  },warning=function(w) {
    warnings[[length(warnings)+1L]]<<-list(message=conditionMessage(w),class=class(w))
    invokeRestart('muffleWarning')
  }),error=function(e) {
    failure<<-list(message=conditionMessage(e),class=class(e),
      details=if(inherits(e,'aft_active_face_error'))e$details else NULL,last_attempt=last_attempt)
    status<-if(inherits(e,'aft_active_face_error'))e$status else if(inherits(e,'aft_numeric_error'))'numerical_failure' else 'numerical_failure'
    finish(status,FALSE,conditionMessage(e))
  })
}
