# Source/cohort/held-out assessment uses n_boot_cohort configured draws.
# Cohort apparent assessment holds the full-cohort update fixed.
rows<-list();heldout_constant<-list()
add_cal<-function(x,p,scope,method,training,mode="heldout_fixed") {
  b<-calibration_bootstrap(x$EXP_LABEL,p,n_boot_cohort,paste(scope,method,sep="."),mode,n_cores)
  z<-b$summary%>%mutate(scope=scope,method=method,training_cohorts=training,N=nrow(x),events=sum(x$EXP_LABEL),prevalence=mean(x$EXP_LABEL))
  resampling_registry[[paste(scope,method)]]<<-z
  z
}
source_x<-filter(d,cohort=="PRONIA")
source_cal<-add_cal(source_x,source_x$Probs_platt,"PRONIA model-origin","Original imported probability","historical frozen calibrator","original")
write.csv(source_cal,file.path(tables_dir,"source_calibration.csv"),row.names=FALSE)
for(co in sort(unique(external$cohort))) {
  te<-filter(external,cohort==co);tr<-filter(external,cohort!=co)
  rows[[paste0(co,".original")]]<-add_cal(te,te$Probs_platt,co,"Original","historical frozen calibrator","original")
  p_app<-update_probabilities(te$EXP_LABEL,te$Probs_platt,"slope")$p
  rows[[paste0(co,".apparent")]]<-add_cal(te,p_app,co,"Cohort apparent update; conditional fixed-update assessment",co)
  train_fit<-calibration_fit(tr$EXP_LABEL,tr$Probs_platt,"slope")
  heldout<-if(is.null(train_fit))rep(NA_real_,nrow(te))else plogis(coef(train_fit)[1]+coef(train_fit)[2]*safe_logit(te$Probs_platt))
  rows[[paste0(co,".heldout")]]<-add_cal(te,heldout,paste0("IECV ",co),"Held-out update assessment",paste(sort(unique(tr$cohort)),collapse="+"))
  ref<-mean(te$EXP_LABEL)*(1-mean(te$EXP_LABEL));bs<-mean((mean(tr$EXP_LABEL)-te$EXP_LABEL)^2)
  heldout_constant[[co]]<-tibble(evaluation_cohort=co,training_cohorts=paste(sort(unique(tr$cohort)),collapse="+"),training_prevalence=mean(tr$EXP_LABEL),N=nrow(te),events=sum(te$EXP_LABEL),Brier=bs,evaluation_prevalence_reference_Brier=ref,Brier_skill=1-bs/ref)
}
write.csv(bind_rows(cal_boot,bind_rows(rows)),file.path(tables_dir,"calibration_audit.csv"),row.names=FALSE)
write.csv(bind_rows(heldout_constant),file.path(tables_dir,"iecv_training_prevalence_constant.csv"),row.names=FALSE)
