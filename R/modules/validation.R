# Deterministic fixed-score validation permutations.
permute_labels <- function(y, strata = NULL) {
  out <- integer(length(y)); blocks <- if(is.null(strata)) list(all=seq_along(y)) else split(seq_along(y),strata)
  for (ii in blocks) { k <- sum(y[ii] == 1L); if(k) out[ii[sample.int(length(ii),k,replace=FALSE)]] <- 1L }
  out
}
validation_nulls <- list(); validation_summaries <- list(); validation_rows <- list()
perf_one <- function(x, scope, blocked = FALSE) {
  strata <- if(blocked) x$cohort else NULL
  auc_obs <- auc_fast(x$Mean_Score,x$EXP_LABEL); bacc_obs <- bacc_fast(x$PRED_LABEL,x$EXP_LABEL)
  auc_null <- deterministic_rows(n_perm,paste0("validation.auc.",scope),function(i){yy<-permute_labels(x$EXP_LABEL,strata);auc_fast(x$Mean_Score,yy)},cores=n_cores)[,1]
  bacc_null <- deterministic_rows(n_perm,paste0("validation.bacc.",scope),function(i){yy<-permute_labels(x$EXP_LABEL,strata);bacc_fast(x$PRED_LABEL,yy)},cores=n_cores)[,1]
  validation_nulls[[scope]] <<- list(auc=auc_null,bacc=bacc_null)
  validation_summaries[[paste0(scope,".auc")]] <<- null_summary(auc_null,auc_obs,paste0("validation.auc.",scope),"AUC","absolute_reference",.5)
  validation_summaries[[paste0(scope,".bacc")]] <<- null_summary(bacc_null,bacc_obs,paste0("validation.bacc.",scope),"balanced accuracy","upper")
  roc <- pROC::roc(x$EXP_LABEL,x$Mean_Score,direction="<",quiet=TRUE); ci<-as.numeric(pROC::ci.auc(roc,method="delong"))
  tp<-sum(x$PRED_LABEL==1L&x$EXP_LABEL==1L);fp<-sum(x$PRED_LABEL==1L&x$EXP_LABEL==0L);tn<-sum(x$PRED_LABEL==0L&x$EXP_LABEL==0L);fn<-sum(x$PRED_LABEL==0L&x$EXP_LABEL==1L)
  tibble(scope=scope,N=nrow(x),events=sum(x$EXP_LABEL),non_events=nrow(x)-sum(x$EXP_LABEL),AUC=auc_obs,AUC_lo=ci[1],AUC_hi=ci[3],
         AUC_p_two_sided=validation_summaries[[paste0(scope,".auc")]]$p_add_one,
         TP=tp,FP=fp,TN=tn,FN=fn,Sensitivity=tp/(tp+fn),Specificity=tn/(tn+fp),BACC=bacc_obs,
         BACC_p_upper=validation_summaries[[paste0(scope,".bacc")]]$p_add_one,
         PPV=tp/(tp+fp),NPV=tn/(tn+fn),LR_pos=(tp/(tp+fn))/(fp/(tn+fp)),
         LR_neg=(fn/(tp+fn))/(tn/(tn+fp)),Brier=mean((x$Probs_platt-x$EXP_LABEL)^2))
}
external <- filter(d,cohort!="PRONIA")
validation_rows[[1]] <- perf_one(filter(d,cohort=="PRONIA"),"PRONIA model-origin")
validation_rows[[2]] <- perf_one(external,"External pooled",TRUE)
for (co in sort(unique(d$cohort))) validation_rows[[length(validation_rows)+1L]] <- perf_one(filter(d,cohort==co),co)
primary_validation <- bind_rows(validation_rows)
write.csv(primary_validation,file.path(tables_dir,"primary_validation.csv"),row.names=FALSE)
sex_performance <- map_dfr(c(1,2),function(s) {
  x<-filter(d,sex==s);label<-c("Male","Female")[[s]]
  perf_one(x,paste0("Sex ",label),TRUE)%>%mutate(recorded_sex=label,interpretation="descriptive pooled discovery/external; not a fairness test")
})
write.csv(sex_performance,file.path(tables_dir,"sex_stratified_performance.csv"),row.names=FALSE)
saveRDS(validation_nulls,file.path(private_dir,"validation_permutation_nulls.rds"))

# Fisher conditional fixed-margin simulation, using the probability-ordering
# statistic and floating-point tie rule in stats::fisher.test. Each table has an indexed stream.
categorical_tests <- map_dfr(c("sex","EXP_LABEL"),function(v) {
  tab<-table(d$cohort,d[[v]],useNA="no");id<-paste0("baseline.fisher.",v)
  observed<-sum(lfactorial(tab))
  values<-deterministic_rows(n_fisher,id,function(i)sum(lfactorial(r2dtable(1,rowSums(tab),colSums(tab))[[1]])),n_cores)[,1]
  ns<-null_summary(values,observed/(1+64*.Machine$double.eps),id,"negative log conditional table probability up to fixed-margin constant","upper")
  ns$observed<-observed;ns$tie_rule<-"stats::fisher.test 1+64*machine_epsilon"
  validation_summaries[[id]]<<-ns
  tibble(variable=v,test="Fisher conditional fixed-margin simulation",N=sum(tab),p=ns$p_add_one,replicates=n_fisher,valid=ns$valid,invalid=ns$invalid,exceedances=ns$exceedances)
})%>%mutate(q_BH_categorical_2=p.adjust(p,"BH",n=2))
write.csv(categorical_tests,file.path(tables_dir,"baseline_group_tests_categorical_seeded.csv"),row.names=FALSE)

# External random-effects synthesis and forest data.
meta_cohorts <- external %>% group_by(cohort) %>% group_modify(~{
  roc<-pROC::roc(.x$EXP_LABEL,.x$Mean_Score,direction="<",quiet=TRUE);auc<-as.numeric(pROC::auc(roc));se<-sqrt(as.numeric(pROC::var(roc,method="delong")));ci<-as.numeric(pROC::ci.auc(roc,method="delong"))
  tibble(N=nrow(.x),events=sum(.x$EXP_LABEL),AUC=auc,AUC_lo=ci[1],AUC_hi=ci[3],DeLong_SE=se,logit_AUC=qlogis(auc),logit_variance=(se/(auc*(1-auc)))^2)
}) %>% ungroup()
meta_fit <- metafor::rma(yi=logit_AUC,vi=logit_variance,method="REML",data=meta_cohorts)
meta_summary <- tibble(cohort="Random-effects summary",N=sum(meta_cohorts$N),events=sum(meta_cohorts$events),AUC=plogis(meta_fit$b[1]),AUC_lo=plogis(meta_fit$ci.lb),AUC_hi=plogis(meta_fit$ci.ub),
                       DeLong_SE=NA_real_,logit_AUC=as.numeric(meta_fit$b[1]),logit_variance=NA_real_,tau_squared=meta_fit$tau2,I_squared_percent=meta_fit$I2,Q=meta_fit$QE,Q_p=meta_fit$QEp)
forest_data <- bind_rows(meta_cohorts %>% mutate(tau_squared=NA_real_,I_squared_percent=NA_real_,Q=NA_real_,Q_p=NA_real_),meta_summary)
write.csv(forest_data,file.path(figdata_dir,"external_auc_random_effects_forest_data.csv"),row.names=FALSE)
write.csv(meta_summary,file.path(tables_dir,"external_auc_meta_summary.csv"),row.names=FALSE)
