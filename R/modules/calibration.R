# Pooled external participant bootstrap; every apparent update is refitted in its resample.
cal_methods<-c(Original="original",`Apparent intercept-only`="apparent_intercept",`Apparent intercept+slope`="apparent_slope")
cal_points<-list();cal_rows<-list()
for(label in names(cal_methods)) {
  method<-cal_methods[[label]]
  p<-if(method=="original")external$Probs_platt else update_probabilities(external$EXP_LABEL,external$Probs_platt,if(method=="apparent_intercept")"intercept"else"slope")$p
  cal_points[[label]]<-tibble(y=external$EXP_LABEL,p=p)%>%mutate(bin=ntile(p,10))%>%group_by(bin)%>%summarise(N=n(),events=sum(y),mean_predicted=mean(p),observed_fraction=mean(y),.groups="drop")%>%mutate(method=label)
  result<-calibration_bootstrap(external$EXP_LABEL,external$Probs_platt,n_boot_cal,paste0("calibration.pooled.",method),method,n_cores)
  cal_rows[[label]]<-result$summary%>%mutate(method=label,scope="External pooled")
  saveRDS(result$draws,file.path(private_dir,paste0("calibration_pooled_",method,".rds")))
}
cal_boot<-bind_rows(cal_rows)
write.csv(bind_rows(cal_points),file.path(figdata_dir,"external_original_vs_apparent_calibration_curve.csv"),row.names=FALSE)
# Each row records the configured number of attempted, valid and invalid draws.
write.csv(cal_boot,file.path(tables_dir,"external_calibration_original_vs_apparent_bootstrap.csv"),row.names=FALSE)
resampling_registry[["pooled_calibration"]]<-cal_boot
