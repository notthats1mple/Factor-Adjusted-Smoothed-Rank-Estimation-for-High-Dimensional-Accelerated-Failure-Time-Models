# v0.12 scheduling and final-field adapter. This file contains no optimizer
# iteration formula. The frozen FISTA and fixed-face Newton functions are called
# at most once each, and every final numerical field is rebuilt at final x.

.cv012_adapter_error <- function(code,message,details=NULL) {
  stop(structure(list(message=message,call=NULL,code=code,details=details),
    class=c('aft_cv_v012_adapter_error','error','condition')))
}
.cv012_adapter_budget <- function(numerics) {
  if(!is.list(numerics)||is.null(numerics$lost_gradient_log_upper_bound))
    .cv012_adapter_error('numerical_metadata_missing','Missing original full-gradient logarithmic error budget')
  z<-numerics$lost_gradient_log_upper_bound
  if(!is.numeric(z)||length(z)!=1L||is.na(z)||z==Inf)
    .cv012_adapter_error('numerical_metadata_invalid','Invalid original full-gradient logarithmic error budget',numerics)
  E<-if(z== -Inf)0 else exp(z)
  if(!is.finite(E)) .cv012_adapter_error('numerical_metadata_invalid','Full-gradient budget overflow',numerics)
  if(is.finite(z)&&E==0)E<-.Machine$double.xmin*.Machine$double.eps
  list(value=E,log_value=as.numeric(z),underflow_preserved=is.finite(z)&&exp(z)==0)
}
.cv012_adapter_kkt <- function(coef,gradient,q,lambda,budget) {
  beta<-if(q<length(coef))seq.int(q+1L,length(coef))else integer()
  r<-abs(gradient);active<-beta[coef[beta]!=0];zero<-beta[coef[beta]==0]
  r[active]<-abs(gradient[active]+lambda*sign(coef[active]))
  r[zero]<-pmax(abs(gradient[zero])-lambda,0)
  k<-max(c(0,r))
  list(residual=k,coordinates=r,residual_with_underflow_budget=k+budget,
    lost_gradient_upper_bound=budget,exact_nonzero=coef[beta]!=0,
    nuisance_gradient=gradient[seq_len(q)],beta_residual=r[beta])
}
.cv012_adapter_native_block <- function(raw,solver) {
  if(raw$status%in%c('hessian_crosscheck_failure','checkpoint_failure'))
    return(list(reason='Frozen solver reported a formula cross-check or persistence failure',
      native_status=raw$status,message=raw$message))
  if(!identical(raw$status,'numerical_failure'))return(NULL)
  classes<-character();messages<-raw$message
  if(solver=='polish'&&!is.null(raw$failure)) {
    classes<-as.character(raw$failure$class)
    messages<-c(messages,raw$failure$message)
  }
  if(solver=='fista'&&length(raw$numerical_events))for(ev in raw$numerical_events) {
    if(identical(ev$type,'failure')) {
      classes<-c(classes,as.character(ev$details$class))
      messages<-c(messages,ev$details$message)
    }
  }
  contradiction<-length(messages)>0L&&any(grepl('contradict|disagree|inconsisten|cross.?check',messages,ignore.case=TRUE))
  controlled_classes<-c('aft_fista_error','aft_numeric_error','aft_active_face_error')
  recognized<-length(classes)>0L&&any(classes%in%controlled_classes)
  if(contradiction||!recognized)return(list(
    reason=if(contradiction)'Frozen solver reported a formula contradiction'else
      'Frozen solver numerical_failure contains an unknown or unrecorded exception class',
    native_status=raw$status,condition_classes=classes,messages=messages))
  NULL
}
aft_cv_v012_controls <- function(config,lambda) {
  if(!is.list(config)||!identical(as.character(config$protocol_version),'0.12'))
    stop('Expected v0.12 configuration')
  fc<-config$fista;pc<-config$polish
  if(!identical(fc[['function']],'aft_solve_fista_experimental')||
     !identical(pc[['function']],'aft_polish_active_face'))stop('Unexpected frozen solver identifiers')
  fc[['function']]<-NULL;pc[['function']]<-NULL
  # The preserved normalizers reject unknown controls. The trace-schema and
  # return-last-accepted-point defaults remain the frozen v0.10 values.
  fc<-aft_fista_control(lambda,fc);pc<-aft_active_face_control(pc)
  a<-config$candidate_acceptance;g<-config$polish_gate
  eq<-function(x,y)isTRUE(all.equal(x,y,check.attributes=FALSE,tolerance=0))
  stopifnot(eq(fc$kkt_tol,1e-7),eq(pc$kkt_tol,1e-7),pc$max_newton_updates<=10,
    eq(a$KKT_plus_budget_max,1e-7),eq(a$post_return_independent_arithmetic_atol,1e-10),
    eq(a$post_return_objective_rtol,1e-12),identical(a$all_original_coordinates,TRUE),
    identical(a$global_optimum_certified,FALSE),eq(g$max_calls_per_penalized_fit,1),
    identical(g$required_fista_status,'iteration_limit'),eq(g$required_saved_checkpoints,5),
    eq(g$gradient_budget_fraction,.1),identical(g$exact_beta_support_and_sign_stable,TRUE),
    identical(g$strict_inactive_margin,TRUE),identical(g$early_trigger,FALSE),
    identical(g$retry_after_failure,FALSE),identical(g$resume_fista_after_polish,FALSE))
  list(fista=fc,polish=pc,gate=g,acceptance=a)
}

# Deterministic, independently checkable endpoint reconstruction. The saved
# solver gradient is checked but never reused as the new standard gradient.
aft_cv_v012_endpoint <- function(Y,delta,design,coefficients,h,lambda,source_fit,config,on_evaluation=NULL) {
  P<-design$P;q<-design$q
  if(!is.numeric(coefficients)||length(coefficients)!=P||any(!is.finite(coefficients)))
    .cv012_adapter_error('invalid_final_coefficients','Final coefficient vector is invalid')
  if(!is.list(source_fit)||!identical(source_fit$design,design)||
     !identical(as.numeric(source_fit$h),as.numeric(h))||
     !identical(as.numeric(source_fit$lambda),as.numeric(lambda)))
    .cv012_adapter_error('source_design_mismatch','Returned solver design, h or lambda differs from the requested problem')
  source_budget<-.cv012_adapter_budget(source_fit$numerics)
  if(!is.null(on_evaluation))on_evaluation('core_objective_gradient')
  core<-aft_loss(Y,delta,design,coefficients,h,config$fista$block_size,gradient=TRUE)
  core_budget<-.cv012_adapter_budget(core$numerics)
  E<-max(source_budget$value,core_budget$value)
  Elog<-max(source_budget$log_value,core_budget$log_value)
  if(!is.null(on_evaluation))on_evaluation('independent_objective_gradient')
  independent<-phase2b_endpoint_arithmetic(Y,delta,design,coefficients,h,lambda)
  beta_index<-if(q<P)seq.int(q+1L,P)else integer()
  objective<-core$loss+lambda*sum(abs(coefficients[beta_index]))
  kk<-.cv012_adapter_kkt(coefficients,core$gradient,q,lambda,E)
  a<-config$candidate_acceptance;atol<-a$post_return_independent_arithmetic_atol
  rtol<-a$post_return_objective_rtol
  comparisons<-list()
  compare<-function(name,actual,reference,limit,budget=0) {
    valid<-is.numeric(actual)&&is.numeric(reference)&&length(actual)==length(reference)&&
      all(is.finite(c(actual,reference)))&&is.numeric(budget)&&length(budget)==1L&&is.finite(budget)&&budget>=0
    difference<-if(valid)max(c(0,abs(actual-reference)))else NA_real_
    comparisons[[name]]<<-list(pass=isTRUE(valid&&difference+budget<=limit),
      valid=valid,max_absolute_difference=difference,comparison_limit=limit,budget=budget,
      difference_plus_budget=if(valid)difference+budget else NA_real_,
      formula_discrepancy=!valid||isTRUE(difference>limit),
      budget_consumes_tolerance=TRUE)
  }
  gb<-independent$numerics$gradient_tail_bound
  lb<-independent$numerics$loss_tail_bound
  compare('rebuilt_core_vs_independent_gradient',core$gradient,independent$gradient,atol,E+gb)
  compare('rebuilt_core_vs_independent_KKT',kk$coordinates,independent$kkt$coordinates,atol,E+gb)
  compare('rebuilt_core_vs_independent_KKT_scalar',kk$residual,independent$kkt$residual,atol,E+gb)
  compare('rebuilt_core_vs_independent_objective',objective,independent$objective,atol+rtol*abs(objective),lb)
  compare('returned_gradient_vs_independent',source_fit$gradient,independent$gradient,atol,source_budget$value+gb)
  compare('returned_KKT_coordinates_vs_independent',source_fit$kkt$coordinates,independent$kkt$coordinates,atol,source_budget$value+gb)
  compare('returned_KKT_scalar_vs_independent',source_fit$kkt$residual,independent$kkt$residual,atol,source_budget$value+gb)
  compare('returned_objective_vs_independent',source_fit$objective,independent$objective,
    atol+rtol*abs(objective),lb)
  d<-if(!is.null(design$d))design$d else q
  K<-if(!is.null(design$K))design$K else 0L
  if(d+K!=q) .cv012_adapter_error('design_partition_mismatch','Nuisance block does not equal d+K')
  theta<-coefficients[seq_len(d)]
  gamma<-if(K)coefficients[d+seq_len(K)]else numeric()
  beta<-coefficients[beta_index]
  bpred<-if(identical(design$kind,'pca')) {
    if(!is.null(on_evaluation))on_evaluation('prediction_coefficients')
    aft_prediction_coefficients(design,coefficients)
  } else beta
  if(!is.null(on_evaluation))on_evaluation('training_scores')
  scores<-as.numeric(aft_forward(design,coefficients))
  compare('rebuilt_training_scores_vs_independent',scores,independent$fitted_scores,atol)
  if(any(!is.finite(c(objective,core$gradient,theta,gamma,beta,bpred,scores))))
    .cv012_adapter_error('nonfinite_standard_fields','Reconstructed standard fields are nonfinite')
  # Retain both original records. The standard error budget takes their maximum
  # so a representable log budget is never silently replaced by a rounded zero.
  numerics<-core$numerics
  numerics$lost_gradient_log_upper_bound<-Elog
  numerics$adapter_source_lost_gradient_log_upper_bound<-source_budget$log_value
  numerics$adapter_recomputed_lost_gradient_log_upper_bound<-core_budget$log_value
  numerics$adapter_gradient_budget<-E
  numerics$adapter_finite_log_underflow_preserved<-source_budget$underflow_preserved||core_budget$underflow_preserved
  pass<-all(vapply(comparisons,function(x)isTRUE(x$pass),logical(1)))
  discrepancy<-any(vapply(comparisons,function(x)isTRUE(x$formula_discrepancy),logical(1)))
  core_stationary<-kk$residual_with_underflow_budget<=a$KKT_plus_budget_max
  independent_stationary<-independent$kkt$residual+gb<=a$KKT_plus_budget_max
  list(standard=list(coefficients=as.numeric(coefficients),theta=as.numeric(theta),
    gamma=as.numeric(gamma),beta=as.numeric(beta),b_pred=as.numeric(bpred),
    fitted_scores=scores,gradient=core$gradient,kkt=kk,objective=objective,smooth_loss=core$loss,
    numerics=numerics,gradient_budget=E,gradient_budget_log=Elog,
    design=design,h=h,lambda=lambda,K=K,q=q,
    exact_nonzero_beta=sum(beta!=0),global_optimum_certified=FALSE,
    original_solver_status=source_fit$status,original_fit_status=source_fit$status,
    coefficients_source='Final returned coefficient vector; all derived fields reconstructed',
    prediction_interface=if(design$kind=='pca')'aft_predict_map(design, coefficients, Xvalid, Zvalid)'else
      'Explicit fixture design: Wnew %*% coefficients'),
    verification=list(pass=pass,formula_discrepancy=discrepancy,
      core_stationary=core_stationary,independent_stationary=independent_stationary,
      comparisons=comparisons,independent=independent,
      core_recomputed=core,source_budget=source_budget,core_budget=core_budget,
      final_tracked_budget=E,final_tracked_budget_log=Elog,
      budgets_consume_fixed_tolerances=TRUE,global_optimum_certified=FALSE))
}

aft_solve_cv_v012 <- function(Y,delta,design,initial,h,lambda,config,
    checkpoint=NULL,event=NULL,fista_fn=aft_solve_fista_experimental,
    polish_fn=aft_polish_active_face) {
  ctl<-aft_cv_v012_controls(config,lambda)
  if(!is.function(fista_fn)||!is.function(polish_fn)||
     !is.null(checkpoint)&&!is.function(checkpoint)||!is.null(event)&&!is.function(event))
    stop('Solver and callback arguments must be functions or NULL where permitted')
  for(nm in c('aft_loss','aft_forward','aft_prediction_coefficients','phase2b_endpoint_arithmetic'))
    if(!exists(nm,mode='function'))stop(paste('Load preserved pure arithmetic before the adapter:',nm))
  .aft_loss_validate(Y,delta,design,initial,h,ctl$fista$block_size)
  fista<-NULL;polish<-NULL;standard<-NULL;verification<-list()
  snapshots<-list();checkpoint_fault<-NULL
  gate<-list(evaluated=FALSE,eligible=FALSE,reason='FISTA has not returned')
  counts<-list(fista_calls=0L,polish_calls=0L,fista_checkpoint_callbacks=0L,
    distinct_fista_checkpoints=0L,polish_checkpoint_callbacks=0L,
    endpoint_verifications=0L,recomputed_core_objective=0L,recomputed_core_gradient=0L,
    independent_pair_objective=0L,independent_pair_gradient=0L,
    reconstructed_training_scores=0L,prediction_coefficient_maps=0L,polish_gate_evaluations=0L)
  timings<-list(fista_seconds=0,polish_seconds=0,endpoint_verification_seconds=0,gate_seconds=0)
  stage<-'prepared';native_status<-NULL
  started<-proc.time()[['elapsed']]
  emit<-function(payload)if(!is.null(event))event(payload)else invisible(NULL)
  result<-function(status,success=FALSE,blocked=FALSE,error=NULL) {
    if(!is.null(standard)) {
      standard$status<-status;standard$success<-success
      standard$fista_status<-if(!is.null(fista))fista$status else NULL
      standard$polish_status<-if(!is.null(polish))polish$status else NULL
      standard$original_fista_status<-if(!is.null(fista))fista$status else NULL
      standard$original_polish_status<-if(!is.null(polish))polish$status else NULL
      standard$solver_id<-'phase3a_v0.12_frozen_fista_optional_fixed_face_polish'
      standard$fista_iterations<-if(!is.null(fista))fista$iterations else 0L
      standard$polish_iterations<-if(!is.null(polish))polish$iterations else 0L
      standard$iterations<-standard$fista_iterations+standard$polish_iterations
      standard$blocked<-blocked
    }
    structure(list(status=status,success=success,blocked=blocked,error=error,stage=stage,
      fista=fista,polish=polish,standard=standard,gate=gate,verification=verification,
      checkpoints=snapshots,counts=counts,
      native_evaluations=list(fista=if(!is.null(fista))fista$evaluations else NULL,
        polish=if(!is.null(polish))polish$evaluations else NULL),
      timings=c(timings,list(adapter_elapsed_seconds=proc.time()[['elapsed']]-started)),
      controls=ctl,global_optimum_certified=FALSE,automatic_retry=FALSE),
      class=c('aft_cv_v012_result','list'))
  }
  retain_snapshot<-function(record,origin) {
    k<-record$iterations;x<-record$coefficients
    if(!is.numeric(k)||length(k)!=1L||!is.finite(k)||k<=0||k!=floor(k)||
       k%%ctl$fista$checkpoint_every_updates!=0||!is.numeric(x)||length(x)!=design$P||any(!is.finite(x)))
      .cv012_adapter_error('invalid_accepted_checkpoint','Checkpoint does not contain a valid periodic accepted x',record)
    if(!is.null(record$state$x)&&!identical(as.numeric(record$state$x),as.numeric(x)))
      .cv012_adapter_error('checkpoint_point_mismatch','Checkpoint coefficients disagree with its accepted x')
    key<-as.character(as.integer(k))
    if(!is.null(snapshots[[key]])) {
      if(!identical(as.numeric(snapshots[[key]]$coefficients),as.numeric(x)))
        .cv012_adapter_error('duplicate_checkpoint_mismatch','One accepted-update index has inconsistent coefficients')
      snapshots[[key]]$origins<<-unique(c(snapshots[[key]]$origins,origin))
    } else snapshots[[key]]<<-list(iteration=as.integer(k),coefficients=as.numeric(x),origins=origin)
    counts$distinct_fista_checkpoints<<-length(snapshots)
    invisible(TRUE)
  }
  fista_checkpoint<-function(record) {
    counts$fista_checkpoint_callbacks<<-counts$fista_checkpoint_callbacks+1L
    tryCatch({
      if(!is.null(checkpoint))checkpoint(list(solver='fista',iteration=record$iterations,record=record))
      retain_snapshot(record,'fista_checkpoint_callback')
    },error=function(e){checkpoint_fault<<-list(message=conditionMessage(e),class=class(e));stop(e)})
  }
  polish_checkpoint<-function(record) {
    counts$polish_checkpoint_callbacks<<-counts$polish_checkpoint_callbacks+1L
    tryCatch({if(!is.null(checkpoint))checkpoint(list(solver='polish',iteration=record$iteration,record=record))},
      error=function(e){checkpoint_fault<<-list(message=conditionMessage(e),class=class(e));stop(e)})
  }
  call_saved<-function(solver,fn,args) {
    public_args<-args;public_args$checkpoint<-NULL
    emit(list(type='solver_call_started',solver=solver,input=public_args,control=args$control))
    counts[[paste0(solver,'_calls')]]<<-counts[[paste0(solver,'_calls')]]+1L
    t0<-proc.time()[['elapsed']]
    value<-tryCatch(do.call(fn,args),error=function(e)structure(list(message=conditionMessage(e),class=class(e)),
      class='cv012_solver_raised'),interrupt=function(e)structure(list(message='Execution interrupted',class=class(e)),
      class='cv012_solver_raised'))
    elapsed<-proc.time()[['elapsed']]-t0
    # The raw returned object is handed to durable storage before validation,
    # field reconstruction, gate evaluation, or any conditional second call.
    if(inherits(value,'cv012_solver_raised')) {
      emit(list(type='solver_raised',solver=solver,error=value,elapsed_seconds=elapsed))
      .cv012_adapter_error('solver_raised',paste(solver,'raised an exception'),value)
    }
    emit(list(type='solver_returned',solver=solver,result=value,elapsed_seconds=elapsed))
    timings[[paste0(solver,'_seconds')]]<<-elapsed
    value
  }
  evaluate_endpoint<-function(raw,label) {
    counts$endpoint_verifications<<-counts$endpoint_verifications+1L
    actual_evaluation<-function(kind) {
      if(kind=='core_objective_gradient') {
        counts$recomputed_core_objective<<-counts$recomputed_core_objective+1L
        counts$recomputed_core_gradient<<-counts$recomputed_core_gradient+1L
      } else if(kind=='independent_objective_gradient') {
        counts$independent_pair_objective<<-counts$independent_pair_objective+1L
        counts$independent_pair_gradient<<-counts$independent_pair_gradient+1L
      } else if(kind=='prediction_coefficients')counts$prediction_coefficient_maps<<-counts$prediction_coefficient_maps+1L
      else if(kind=='training_scores')counts$reconstructed_training_scores<<-counts$reconstructed_training_scores+1L
      else .cv012_adapter_error('unknown_evaluation_counter','Unknown pure arithmetic counter')
    }
    t0<-proc.time()[['elapsed']]
    ep<-tryCatch(aft_cv_v012_endpoint(Y,delta,design,raw$coefficients,h,lambda,raw,config,actual_evaluation),
      finally={timings$endpoint_verification_seconds<<-timings$endpoint_verification_seconds+proc.time()[['elapsed']]-t0})
    verification[[label]]<<-ep$verification
    standard<<-ep$standard
    ep$verification
  }
  check_status<-function(raw,allowed,label) {
    if(!is.list(raw)||!is.character(raw$status)||length(raw$status)!=1L||!raw$status%in%allowed)
      .cv012_adapter_error('unknown_solver_status',paste(label,'returned an unknown status'),raw$status)
    if(!is.logical(raw$success)||length(raw$success)!=1L||is.na(raw$success)||
       !identical(raw$success,identical(raw$status,'stationary_to_tolerance')))
      .cv012_adapter_error('solver_status_contradiction',paste(label,'status and success disagree'))
    if(!identical(raw$global_optimum_certified,FALSE))
      .cv012_adapter_error('unexpected_optimum_certificate','Frozen solver must retain global_optimum_certified=FALSE')
  }
  tryCatch({
    stage<-'fista_call'
    fista<-call_saved('fista',fista_fn,list(Y=Y,delta=delta,design=design,h=h,lambda=lambda,
      initial=initial,control=ctl$fista,checkpoint=fista_checkpoint))
    check_status(fista,c('stationary_to_tolerance','iteration_limit','no_events','strict_escape_numeric',
      'diagnostic_inconclusive','numerical_uncertainty','numerical_failure','numerical_stagnation','line_search_failure'),'FISTA')
    native_status<-fista$status
    native_block<-.cv012_adapter_native_block(fista,'fista')
    if(!is.null(checkpoint_fault)).cv012_adapter_error('checkpoint_failure','Checkpoint persistence or identity failed',checkpoint_fault)
    if(!is.numeric(fista$iterations)||length(fista$iterations)!=1L||!is.finite(fista$iterations)||
       fista$iterations<0||fista$iterations>ctl$fista$max_iter||fista$iterations!=floor(fista$iterations))
      .cv012_adapter_error('fista_budget_contradiction','FISTA accepted-update count is invalid')
    # A terminal record exactly on the periodic grid is a saved accepted point;
    # if it duplicates the final callback, equality is checked and it counts once.
    if(fista$iterations>0&&fista$iterations%%ctl$fista$checkpoint_every_updates==0)
      retain_snapshot(fista,'saved_terminal_fista_return')
    if(length(snapshots)&&any(vapply(snapshots,function(x)x$iteration,numeric(1))>fista$iterations))
      .cv012_adapter_error('checkpoint_after_terminal','A checkpoint lies after the returned accepted update')
    stage<-'fista_endpoint_verification'
    fv<-evaluate_endpoint(fista,'fista')
    if(fv$formula_discrepancy)
      .cv012_adapter_error('endpoint_formula_discrepancy','FISTA endpoint failed independent arithmetic comparison',fv)
    if(!is.null(native_block)) {
      gate<-list(evaluated=FALSE,eligible=FALSE,reason='Native solver contradiction or unknown exception stops this execution')
      return(result(fista$status,FALSE,TRUE,native_block))
    }
    if(identical(fista$status,'stationary_to_tolerance')) {
      gate<-list(evaluated=FALSE,eligible=FALSE,reason='FISTA already returned stationary; polishing not called')
      if(!fv$pass||!fv$core_stationary||!fv$independent_stationary)
        return(result('final_numerical_verification_failed',FALSE,FALSE))
      stage<-'complete';return(result('stationary_to_tolerance',TRUE,FALSE))
    }
    if(!identical(fista$status,'iteration_limit')) {
      gate<-list(evaluated=FALSE,eligible=FALSE,reason=paste('FISTA status is',fista$status))
      stage<-'complete';return(result(fista$status,FALSE,FALSE))
    }
    stage<-'polish_gate';counts$polish_gate_evaluations<-counts$polish_gate_evaluations+1L
    gate_t0<-proc.time()[['elapsed']]
    ordered<-snapshots[order(vapply(snapshots,function(x)x$iteration,numeric(1)))]
    chosen<-tail(ordered,ctl$gate$required_saved_checkpoints)
    beta<-if(design$q<design$P)seq.int(design$q+1L,design$P)else integer()
    enough<-length(chosen)==ctl$gate$required_saved_checkpoints
    final_support<-as.numeric(fista$coefficients)[beta]!=0
    final_sign<-sign(as.numeric(fista$coefficients)[beta])
    stable_support<-enough&&all(vapply(chosen,function(x)identical(x$coefficients[beta]!=0,final_support),logical(1)))
    stable_sign<-enough&&all(vapply(chosen,function(x)identical(sign(x$coefficients[beta]),final_sign),logical(1)))
    zero<-beta[!final_support];E<-standard$gradient_budget
    margin<-if(length(zero))min(lambda-abs(standard$gradient[zero]))-E else NA_real_
    budget_limit<-ctl$fista$kkt_tol*ctl$gate$gradient_budget_fraction
    budget_pass<-is.finite(E)&&E<=budget_limit
    inactive_pass<-!length(zero)||isTRUE(margin>0)
    eligible<-enough&&stable_support&&stable_sign&&budget_pass&&inactive_pass
    reasons<-c(if(!enough)'insufficient_distinct_accepted_checkpoints',
      if(enough&&!stable_support)'exact_beta_support_changed',
      if(enough&&!stable_sign)'exact_beta_sign_changed',
      if(!budget_pass)'terminal_gradient_budget_exceeds_gate',
      if(!inactive_pass)'inactive_beta_lacks_strict_margin')
    gate<-list(evaluated=TRUE,eligible=eligible,reason=if(eligible)'all_fixed_terminal_gate_conditions_passed'else paste(reasons,collapse='; '),
      fista_status=fista$status,required_saved_checkpoints=ctl$gate$required_saved_checkpoints,
      available_distinct_saved_checkpoints=length(ordered),last_five=chosen,
      accepted_update_indices=vapply(chosen,function(x)x$iteration,numeric(1)),
      enough_checkpoints=enough,exact_support_stable=stable_support,exact_sign_stable=stable_sign,
      terminal_gradient_budget=E,terminal_gradient_budget_log=standard$gradient_budget_log,
      gradient_budget_limit=budget_limit,gradient_budget_pass=budget_pass,
      exact_zero_beta_indices=zero,inactive_applicable=length(zero)>0L,
      inactive_margin=margin,inactive_condition_pass=inactive_pass,
      terminal_coefficients=fista$coefficients,terminal_gradient=standard$gradient,
      max_polish_calls=1L,early_trigger=FALSE)
    timings$gate_seconds<-proc.time()[['elapsed']]-gate_t0
    if(!eligible) {stage<-'complete';return(result('iteration_limit_polish_gate_rejected',FALSE,FALSE))}
    stage<-'polish_call'
    polish<-call_saved('polish',polish_fn,list(Y=Y,delta=delta,design=design,
      initial=fista$coefficients,h=h,lambda=lambda,control=ctl$polish,checkpoint=polish_checkpoint))
    check_status(polish,c('stationary_to_tolerance','polish_iteration_limit','no_events','active_face_not_ready',
      'numerical_uncertainty','numerical_failure','hessian_crosscheck_failure','active_hessian_not_positive_definite',
      'linear_system_inaccurate','non_descent_direction','face_or_line_search_failure','checkpoint_failure','invalid_input'),'polisher')
    native_status<-polish$status
    native_block<-.cv012_adapter_native_block(polish,'polish')
    if(!is.null(checkpoint_fault)).cv012_adapter_error('checkpoint_failure','Polisher checkpoint persistence failed',checkpoint_fault)
    if(counts$polish_calls!=1L||!is.numeric(polish$iterations)||length(polish$iterations)!=1L||
       !is.finite(polish$iterations)||polish$iterations<0||polish$iterations>ctl$polish$max_newton_updates||
       polish$iterations!=floor(polish$iterations))
      .cv012_adapter_error('polish_budget_contradiction','Conditional polish exceeded the one-call/update budget')
    stage<-'polish_endpoint_verification'
    pv<-evaluate_endpoint(polish,'polish')
    if(pv$formula_discrepancy)
      .cv012_adapter_error('endpoint_formula_discrepancy','Polished endpoint failed independent arithmetic comparison',pv)
    if(!is.null(native_block))return(result(polish$status,FALSE,TRUE,native_block))
    if(!identical(polish$status,'stationary_to_tolerance')) {stage<-'complete';return(result(polish$status,FALSE,FALSE))}
    if(!pv$pass||!pv$core_stationary||!pv$independent_stationary)
      return(result('final_numerical_verification_failed',FALSE,FALSE))
    stage<-'complete';result('stationary_to_tolerance',TRUE,FALSE)
  },error=function(e) {
    code<-if(inherits(e,'aft_cv_v012_adapter_error'))e$code else 'adapter_unexpected_exception'
    result(code,FALSE,TRUE,list(message=conditionMessage(e),class=class(e),
      native_status=native_status,details=if(inherits(e,'aft_cv_v012_adapter_error'))e$details else NULL))
  })
}
