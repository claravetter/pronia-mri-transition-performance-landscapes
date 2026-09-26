# Canonical analysis-set standardization and authoritative unpenalized GLMs.
symptom_vars<-c("Mean_Score","COGDIS_score","SIPS_Positiv_Gesamt","SIPS_Negativ_Gesamt")
symptom<-d%>%filter(if_all(all_of(c("EXP_LABEL",symptom_vars,"cohort")),~!is.na(.x)&if(is.numeric(.x))is.finite(.x)else TRUE))
functioning<-symptom%>%filter(is.finite(Psychosoz_aequiv));check_analysis_counts(symptom,functioning,cfg,mode,tables_dir)
sc_sym<-scale_fit(symptom,symptom_vars,"SCALE_M1_CC");dm<-sc_sym$data%>%mutate(y=EXP_LABEL,cohort_f=relevel(factor(cohort),"NAPLS-3"))
f_m1<-y~Mean_Score_z*(COGDIS_score_z+SIPS_Positiv_Gesamt_z+SIPS_Negativ_Gesamt_z)+cohort_f
m1<-glm(f_m1,data=dm,family=binomial())
fvars<-c(symptom_vars,"Psychosoz_aequiv");sc_f<-scale_fit(functioning,fvars,"SCALE_F1_CC");df<-sc_f$data%>%mutate(y=EXP_LABEL,cohort_f=relevel(factor(cohort),"NAPLS-3"))
f_f1<-y~Mean_Score_z*(COGDIS_score_z+SIPS_Positiv_Gesamt_z+SIPS_Negativ_Gesamt_z+Psychosoz_aequiv_z)+cohort_f
f1<-glm(f_f1,data=df,family=binomial())
sc_match<-list(data=sc_f$data,params=filter(sc_f$params,variable%in%symptom_vars));dmatch<-sc_match$data%>%mutate(y=EXP_LABEL,cohort_f=relevel(factor(cohort),"NAPLS-3"));m1_match<-glm(f_m1,data=dmatch,family=binomial())
tidy_fit<-function(fit,model_id,scaling_id){s<-coef(summary(fit));terms<-rownames(s);out<-tibble(model_id=model_id,scaling_id=scaling_id,reference_cohort="NAPLS-3",N=nobs(fit),events=sum(model.response(model.frame(fit))),formula=paste(deparse(formula(fit)),collapse=" "),term=terms,estimate=s[,1],SE=s[,2],statistic=s[,3],p=s[,4],ci_lo=s[,1]-qnorm(.975)*s[,2],ci_hi=s[,1]+qnorm(.975)*s[,2],is_score_interaction=str_detect(terms,"Mean_Score_z:"),q_BH_within_model_interactions=NA_real_);ii<-which(out$is_score_interaction);out$q_BH_within_model_interactions[ii]<-p.adjust(out$p[ii],"BH");out}
glm_coefficients<-bind_rows(tidy_fit(m1,"GLM_M1_CANONICAL","SCALE_M1_CC"),tidy_fit(f1,"GLM_F1_CANONICAL","SCALE_F1_CC"),tidy_fit(m1_match,"GLM_M1_MATCHED","SCALE_F1_CC"))
write.csv(glm_coefficients,file.path(tables_dir,"canonical_glm_coefficients.csv"),row.names=FALSE)
scaling_parameters<-bind_rows(sc_sym$params,sc_f$params,sc_match$params);write.csv(scaling_parameters,file.path(tables_dir,"canonical_scaling_parameters.csv"),row.names=FALSE)
saveRDS(list(GLM_M1_CANONICAL=m1,GLM_F1_CANONICAL=f1,GLM_M1_MATCHED=m1_match),file.path(private_dir,"canonical_glm_models.rds"))

# Explicit all-available-data scaling comparator and fitted-value/log-likelihood invariance.
allavail_params<-map_dfr(symptom_vars,function(v)tibble(variable=v,mean=mean(d[[v]],na.rm=TRUE),SD=sd(d[[v]],na.rm=TRUE)))
allavail_dm<-symptom
for(v in symptom_vars){p<-filter(allavail_params,variable==v);allavail_dm[[paste0(v,"_z")]]<-(allavail_dm[[v]]-p$mean)/p$SD}
allavail_dm<-allavail_dm%>%mutate(y=EXP_LABEL,cohort_f=relevel(factor(cohort),"NAPLS-3"));allavail_fit<-glm(f_m1,data=allavail_dm,family=binomial())
allavail_coef<-tidy_fit(allavail_fit,"GLM_M1_ALL_AVAILABLE_SCALE_COMPARATOR","SCALE_ALL_AVAILABLE")
write.csv(allavail_coef,file.path(tables_dir,"all_available_scaling_comparator_coefficients.csv"),row.names=FALSE)
scaling_invariance<-tibble(canonical_model_id="GLM_M1_CANONICAL",comparison_model_id="GLM_M1_ALL_AVAILABLE_SCALE_COMPARATOR",same_analysis_rows=identical(as.character(dm$Cases),as.character(allavail_dm$Cases)),
  max_abs_fitted_probability_difference=max(abs(fitted(m1)-fitted(allavail_fit))),logLik_canonical=as.numeric(logLik(m1)),logLik_comparator=as.numeric(logLik(allavail_fit)),logLik_difference=as.numeric(logLik(m1)-logLik(allavail_fit)),
  coefficient_units_changed=TRUE,fit_invariant=max(abs(fitted(m1)-fitted(allavail_fit)))<1e-10)
write.csv(scaling_invariance,file.path(tables_dir,"glm_scaling_invariance.csv"),row.names=FALSE);stopifnot(scaling_invariance$fit_invariant)

# Canonical linear moderation curve, derived from the saved fit and NAPLS-3 reference.
nd<-expand_grid(Mean_Score_z=seq(quantile(dm$Mean_Score_z,.025),quantile(dm$Mean_Score_z,.975),length.out=160),SIPS_Negativ_Gesamt_z=c(-1,0,1),COGDIS_score_z=0,SIPS_Positiv_Gesamt_z=0,cohort_f=factor("NAPLS-3",levels=levels(dm$cohort_f)))
pp<-predict(m1,nd,type="link",se.fit=TRUE);nd<-nd%>%mutate(probability=plogis(pp$fit),lo=plogis(pp$fit-qnorm(.975)*pp$se.fit),hi=plogis(pp$fit+qnorm(.975)*pp$se.fit),SIPS_N=factor(SIPS_Negativ_Gesamt_z,levels=c(-1,0,1),labels=c("Low (-1 SD)","Mean","High (+1 SD)")),model_id="GLM_M1_CANONICAL",reference_cohort="NAPLS-3")
write.csv(nd,file.path(figdata_dir,"canonical_linear_moderation_curve.csv"),row.names=FALSE)

# Symptom (M1/M2/M2b/M3) and functioning (F1/F2/F2b/F3) model hierarchies.
# M3/F3 extend demographic-main-effects models, not demographic-interaction models.
hierarchical_fits<-list(GLM_M1_CANONICAL=m1,GLM_F1_CANONICAL=f1,GLM_M1_MATCHED=m1_match)
hierarchy<-list();hierarchy_diagnostics<-list();hierarchy_coefs<-list(glm_coefficients);slopes<-list()
for(prefix in c("M","F")) {
  raw<-if(prefix=="M")symptom else functioning
  vars<-c(symptom_vars,if(prefix=="F")"Psychosoz_aequiv","age","BrainAGE_SBC_corr")
  if(!all(eligible(raw,c(vars,"sex"))))stop("Hierarchical covariate availability changes declared rows; reconcile before comparing nested fits")
  sid<-if(prefix=="M")"SCALE_M1_CC"else"SCALE_F1_CC"
  scaled<-scale_fit(raw,vars,sid)
  hdata<-scaled$data%>%mutate(y=EXP_LABEL,cohort_f=relevel(factor(cohort),"NAPLS-3"),sex_f=factor(sex,levels=c(1,2)))
  base<-if(prefix=="M")m1 else f1
  fits<-list(`1`=base)
  fits[["2"]]<-glm(update(formula(base),.~.+age_z+sex_f),data=hdata,family=binomial())
  fits[["2b"]]<-glm(update(formula(fits[["2"]]),.~.+Mean_Score_z:age_z+Mean_Score_z:sex_f),data=hdata,family=binomial())
  fits[["3"]]<-glm(update(formula(fits[["2"]]),.~.+BrainAGE_SBC_corr_z+Mean_Score_z:BrainAGE_SBC_corr_z),data=hdata,family=binomial())
  for(name in names(fits)) {
    id<-paste0("GLM_",prefix,name,if(name=="1")"_CANONICAL"else"")
    fit<-fits[[name]];hierarchical_fits[[id]]<-fit
    if(name!="1")hierarchy_coefs[[id]]<-tidy_fit(fit,id,sid)
    hierarchy_diagnostics[[id]]<-tibble(model_id=id,N=nobs(fit),events=sum(model.response(model.frame(fit))),AIC=AIC(fit),BIC=BIC(fit),logLik=as.numeric(logLik(fit)),converged=isTRUE(fit$converged),rank=fit$rank,coefficient_count=length(coef(fit)))
  }
  for(comparison in list(c("1","2"),c("2","2b"),c("2","3"))) {
    small<-fits[[comparison[1]]];large<-fits[[comparison[2]]]
    stopifnot(identical(rownames(model.frame(small)),rownames(model.frame(large))))
    lr<-anova(small,large,test="LRT")
    hierarchy[[paste0(prefix,paste(comparison,collapse="_"))]]<-tibble(hierarchy=prefix,smaller=paste0(prefix,comparison[1]),larger=paste0(prefix,comparison[2]),N=nobs(large),deviance_change=lr$Deviance[2],df=lr$Df[2],p=lr$`Pr(>Chi)`[2],same_rows=TRUE)
  }
  scaling_parameters<-bind_rows(scaling_parameters,filter(scaled$params,variable%in%c("age","BrainAGE_SBC_corr")))
}
for(id in names(hierarchical_fits))for(level in c(-1,0,1)) {
  fit<-hierarchical_fits[[id]]
  contrast<-setNames(rep(0,length(coef(fit))),names(coef(fit)));contrast["Mean_Score_z"]<-1;contrast["Mean_Score_z:SIPS_Negativ_Gesamt_z"]<-level
  estimate<-sum(contrast*coef(fit));se<-sqrt(as.numeric(t(contrast)%*%vcov(fit)%*%contrast))
  slopes[[paste(id,level)]]<-tibble(model_id=id,moderator="SIPS_Negativ_Gesamt_z",level_SD=level,log_odds_slope=estimate,SE=se,OR=exp(estimate),OR_lo=exp(estimate-qnorm(.975)*se),OR_hi=exp(estimate+qnorm(.975)*se),p=2*pnorm(-abs(estimate/se)))
}
glm_coefficients<-bind_rows(hierarchy_coefs)
write.csv(glm_coefficients,file.path(tables_dir,"canonical_glm_coefficients.csv"),row.names=FALSE)
write.csv(scaling_parameters,file.path(tables_dir,"canonical_scaling_parameters.csv"),row.names=FALSE)
write.csv(bind_rows(hierarchy),file.path(tables_dir,"glm_hierarchy.csv"),row.names=FALSE)
write.csv(bind_rows(hierarchy_diagnostics),file.path(tables_dir,"glm_hierarchy_diagnostics.csv"),row.names=FALSE)
write.csv(bind_rows(slopes),file.path(tables_dir,"glm_simple_slopes.csv"),row.names=FALSE)
saveRDS(hierarchical_fits,file.path(private_dir,"canonical_glm_models.rds"))
saveRDS(list(symptom=dm$Cases,functioning=df$Cases,matched=dmatch$Cases),file.path(private_dir,"model_row_ids.rds"))
