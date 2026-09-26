# Conditional second-stage LOCO: frozen MRI input, train-only transforms, paired held-out assessment.
loco_rows<-list();loco_scales<-list();loco_invariance<-list();loco_fits<-list()
for(include_functioning in c(FALSE,TRUE)) {
  variant<-if(include_functioning)"with_functioning"else"without_functioning"
  vars<-c("COGDIS_score","SIPS_Positiv_Gesamt","SIPS_Negativ_Gesamt","age",if(include_functioning)"Psychosoz_aequiv")
  z<-d[eligible(d,c(vars,"sex","cohort")),]
  for(co in sort(unique(z$cohort))) {
    tr<-filter(z,cohort!=co);te<-filter(z,cohort==co)
    scales<-scale_fit(tr,c(vars,"Mean_Score"),paste("LOCO",variant,co,sep="."))
    train<-scales$data;test<-te
    for(v in c(vars,"Mean_Score")){p<-filter(scales$params,variable==v);test[[paste0(v,"_z")]]<-(te[[v]]-p$mean)/p$SD}
    train$sex_f<-factor(train$sex,levels=c(1,2));test$sex_f<-factor(test$sex,levels=c(1,2))
    loco_scales[[paste(variant,co)]]<-mutate(scales$params,heldout_cohort=co,training_cohorts=paste(sort(unique(tr$cohort)),collapse="+"))
    result<-tibble(variant=variant,heldout_cohort=co,y=te$EXP_LABEL,row_id=te$Cases)
    for(nm in c("clinical","clinical_plus_mri")) {
      terms<-c(paste0(vars,"_z"),"sex_f",if(nm=="clinical_plus_mri")"Mean_Score_z")
      fit<-glm(reformulate(terms,"EXP_LABEL"),data=train,family=binomial())
      if(!isTRUE(fit$converged)||any(!is.finite(coef(fit))))stop("LOCO fit failed: ",variant," ",co," ",nm)
      result[[nm]]<-predict(fit,test,type="response")
      # Raw-unit model spans the same unpenalized design; no held-out transform is fitted.
      raw_fit<-glm(reformulate(c(vars,"sex_f",if(nm=="clinical_plus_mri")"Mean_Score"),"EXP_LABEL"),data=train,family=binomial())
      loco_invariance[[paste(variant,co,nm)]]<-tibble(variant=variant,heldout_cohort=co,model=nm,max_abs_prediction_difference=max(abs(result[[nm]]-predict(raw_fit,test,type="response"))))
      loco_fits[[paste(variant,co,nm)]]<-fit
    }
    loco_rows[[paste(variant,co)]]<-result
  }
}
loco_pred<-bind_rows(loco_rows)
write.csv(bind_rows(loco_scales),file.path(tables_dir,"loco_training_scaling.csv"),row.names=FALSE)
write.csv(bind_rows(loco_invariance),file.path(tables_dir,"loco_scaling_invariance.csv"),row.names=FALSE)
saveRDS(list(predictions=loco_pred,fits=loco_fits),file.path(private_dir,"loco_fits_predictions.rds"))
loco_metrics<-function(g) {
  a<-auc_fast(g$clinical,g$y);b<-auc_fast(g$clinical_plus_mri,g$y);c<-mean((g$clinical-g$y)^2);e<-mean((g$clinical_plus_mri-g$y)^2)
  c(clinical_AUC=a,clinical_plus_mri_AUC=b,clinical_Brier=c,clinical_plus_mri_Brier=e,AUC_difference=b-a,Brier_difference=e-c)
}
loco_paired<-list();loco_summary<-list()
for(variant in unique(loco_pred$variant)) {
  g<-filter(loco_pred,.data$variant==!!variant);point<-loco_metrics(g)
  blocks<-split(seq_len(nrow(g)),g$heldout_cohort)
  draws<-deterministic_rows(n_boot_loco,paste0("loco.paired.",variant),function(i)loco_metrics(g[unlist(lapply(blocks,sample_indices,replace=TRUE)),]),n_cores)
  colnames(draws)<-names(point)
  z<-map_dfr(seq_along(point),function(j){ok<-is.finite(draws[,j]);q<-if(any(ok))quantile(draws[ok,j],c(.025,.975))else c(NA,NA);tibble(analysis_id=paste0("loco.paired.",variant),variant=variant,N=nrow(g),events=sum(g$y),metric=names(point)[j],estimate=point[j],lo=q[1],hi=q[2],attempted=n_boot_loco,valid=sum(ok),invalid=sum(!ok),extreme=NA_integer_,extreme_applicable=FALSE,method_review_required=any(!ok),estimand="paired cohort-stratified conditional bootstrap of fixed held-out second-stage predictions")})
  loco_paired[[variant]]<-z;resampling_registry[[paste0("loco.",variant)]]<-z
  loco_summary[[variant]]<-bind_rows(tibble(variant=variant,model="clinical",N=nrow(g),events=sum(g$y),AUC=point[1],Brier=point[3]),tibble(variant=variant,model="clinical_plus_mri",N=nrow(g),events=sum(g$y),AUC=point[2],Brier=point[4]))
}
write.csv(bind_rows(loco_paired),file.path(tables_dir,"loco_paired_comparison.csv"),row.names=FALSE)
write.csv(bind_rows(loco_summary),file.path(tables_dir,"loco_incremental_value.csv"),row.names=FALSE)
