# Probability updating and bootstrap calibration assessment helpers.
safe_logit <- function(p) qlogis(pmin(pmax(p,1e-6),1-1e-6))
calibration_fit <- function(y,p,method=c("slope","intercept")) {
  method<-match.arg(method)
  if(length(unique(y))!=2L || anyNA(y) || any(!is.finite(p))) return(NULL)
  lp<-safe_logit(p)
  fit<-tryCatch(suppressWarnings(if(method=="intercept") glm(y~1+offset(lp),family=binomial()) else glm(y~lp,family=binomial())),error=function(e)NULL)
  if(is.null(fit)||!isTRUE(fit$converged)||any(!is.finite(coef(fit)))||isTRUE(fit$boundary)) return(NULL)
  fit
}
cal_metrics <- function(y,p) {
  mi<-calibration_fit(y,p,"intercept"); ms<-calibration_fit(y,p,"slope")
  prev<-mean(y);ref<-prev*(1-prev)
  c(alpha_int=if(is.null(mi))NA_real_ else unname(coef(mi)[1]),
    alpha=if(is.null(ms))NA_real_ else unname(coef(ms)[1]),beta=if(is.null(ms))NA_real_ else unname(coef(ms)[2]),
    Brier=mean((p-y)^2),prevalence_reference_Brier=ref,Brier_skill=if(ref>0)1-mean((p-y)^2)/ref else NA_real_)
}
update_probabilities <- function(y,p,method) {
  fit<-calibration_fit(y,p,method)
  if(is.null(fit)) return(list(p=rep(NA_real_,length(y)),alpha=NA_real_,beta=NA_real_))
  alpha<-unname(coef(fit)[1]); beta<-if(method=="intercept")1 else unname(coef(fit)[2])
  list(p=plogis(alpha+beta*safe_logit(p)),alpha=alpha,beta=beta)
}
calibration_bootstrap <- function(y,p,B,analysis_id,mode=c("original","apparent_intercept","apparent_slope","heldout_fixed"),cores=1L) {
  mode<-match.arg(mode)
  evaluate<-function(ii) {
    yy<-y[ii];pp<-p[ii];up<-list(p=pp,alpha=NA_real_,beta=NA_real_)
    if(startsWith(mode,"apparent")) up<-update_probabilities(yy,pp,if(mode=="apparent_intercept")"intercept" else "slope")
    c(cal_metrics(yy,up$p),updating_alpha=up$alpha,updating_beta=up$beta)
  }
  point<-evaluate(seq_along(y))
  draws<-deterministic_rows(B,analysis_id,function(i)evaluate(sample.int(length(y),length(y),replace=TRUE)),cores)
  colnames(draws)<-names(point)
  summary<-map_dfr(seq_along(point),function(j) {
    v<-draws[,j];ok<-is.finite(v);q<-if(any(ok))quantile(v[ok],c(.025,.975),names=FALSE)else c(NA,NA)
    applicable<-!startsWith(names(point)[j],"updating_")||startsWith(mode,"apparent")
    tibble(analysis_id=analysis_id,mode=mode,metric=names(point)[j],estimate=point[j],lo=q[1],hi=q[2],
      attempted=if(applicable)B else 0L,valid=if(applicable)sum(ok)else 0L,invalid=if(applicable)sum(!ok)else 0L,
      extreme=NA_integer_,extreme_applicable=FALSE,method_review_required=applicable&&any(!ok),
      estimand=if(startsWith(mode,"apparent"))"apparent assessment with update refitted inside each participant bootstrap" else if(mode=="heldout_fixed")"conditional assessment of held-out predictions; training fits held fixed" else "assessment of inherited probabilities")
  })
  list(summary=summary,draws=draws)
}
