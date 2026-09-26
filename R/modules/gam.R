# Canonical GAM and prespecified diagnostic variants use the same analysis-set scaling.
gam_formula<-y~s(Mean_Score_z,k=10)+cohort_f+s(COGDIS_score_z,k=10)+ti(Mean_Score_z,COGDIS_score_z,k=c(10,10))+s(SIPS_Positiv_Gesamt_z,k=10)+ti(Mean_Score_z,SIPS_Positiv_Gesamt_z,k=c(10,10))+s(SIPS_Negativ_Gesamt_z,k=10)+ti(Mean_Score_z,SIPS_Negativ_Gesamt_z,k=c(10,10))
fit_gam<-function(form,select=TRUE,data=dm)mgcv::gam(form,data=data,family=binomial(),method="REML",select=select)
gam_models<-list(
  GAM_CANONICAL_K10_SELECT=fit_gam(gam_formula,TRUE),
  GAM_DIAGNOSTIC_K6_SELECT=fit_gam(update(gam_formula,.~s(Mean_Score_z,k=6)+cohort_f+s(COGDIS_score_z,k=6)+ti(Mean_Score_z,COGDIS_score_z,k=c(6,6))+s(SIPS_Positiv_Gesamt_z,k=6)+ti(Mean_Score_z,SIPS_Positiv_Gesamt_z,k=c(6,6))+s(SIPS_Negativ_Gesamt_z,k=6)+ti(Mean_Score_z,SIPS_Negativ_Gesamt_z,k=c(6,6))),TRUE),
  GAM_DIAGNOSTIC_K10_NOSELECT=fit_gam(gam_formula,FALSE),
  GAM_DIAGNOSTIC_CR_K10_SELECT=fit_gam(update(gam_formula,.~s(Mean_Score_z,k=10,bs="cr")+cohort_f+s(COGDIS_score_z,k=10,bs="cr")+ti(Mean_Score_z,COGDIS_score_z,k=c(10,10),bs=c("cr","cr"))+s(SIPS_Positiv_Gesamt_z,k=10,bs="cr")+ti(Mean_Score_z,SIPS_Positiv_Gesamt_z,k=c(10,10),bs=c("cr","cr"))+s(SIPS_Negativ_Gesamt_z,k=10,bs="cr")+ti(Mean_Score_z,SIPS_Negativ_Gesamt_z,k=c(10,10),bs=c("cr","cr"))),TRUE)
)
allavail_gam <- fit_gam(gam_formula, TRUE, allavail_dm)
gam_term_table <- function(gm, id, scaling_id) {
  st <- as.data.frame(summary(gm)$s.table); st$term <- rownames(st)
  pcol <- grep("p", names(st), ignore.case = TRUE, value = TRUE)[1]
  statcol <- grep("Chi|F", names(st), value = TRUE)[1]
  out <- as_tibble(st) %>% transmute(model_id=id,scaling_id=scaling_id,reference_cohort="NAPLS-3",N=nobs(gm),events=sum(model.response(model.frame(gm))),term,edf,reference_df=Ref.df,statistic=.data[[statcol]],p=.data[[pcol]],interaction=str_detect(term,"^ti\\("),q_BH_3_tensor=NA_real_)
  ii <- which(out$interaction); out$q_BH_3_tensor[ii] <- p.adjust(out$p[ii], "BH"); out
}
gam_terms<-imap_dfr(gam_models,function(gm,id)gam_term_table(gm,id,"SCALE_M1_CC"))
write.csv(gam_terms,file.path(tables_dir,"gam_bounded_diagnostic_sensitivity.csv"),row.names=FALSE)
gam_scaling_comparison <- bind_rows(
  gam_term_table(gam_models[["GAM_CANONICAL_K10_SELECT"]], "GAM_CANONICAL_K10_SELECT", "SCALE_M1_CC"),
  gam_term_table(allavail_gam, "GAM_ALL_AVAILABLE_SCALE_COMPARATOR", "SCALE_ALL_AVAILABLE")
)
write.csv(gam_scaling_comparison,file.path(tables_dir,"gam_scaling_comparison.csv"),row.names=FALSE)
gam_scaling_fit_comparison <- tibble(
  canonical_model_id="GAM_CANONICAL_K10_SELECT",
  comparison_model_id="GAM_ALL_AVAILABLE_SCALE_COMPARATOR",
  same_analysis_rows=identical(as.character(dm$Cases),as.character(allavail_dm$Cases)),
  max_abs_fitted_probability_difference=max(abs(fitted(gam_models[["GAM_CANONICAL_K10_SELECT"]])-fitted(allavail_gam))),
  logLik_canonical=as.numeric(logLik(gam_models[["GAM_CANONICAL_K10_SELECT"]])),
  logLik_comparator=as.numeric(logLik(allavail_gam)),
  logLik_difference=as.numeric(logLik(gam_models[["GAM_CANONICAL_K10_SELECT"]])-logLik(allavail_gam)),
  fit_invariant=FALSE,
  note="Penalized GAM scaling comparison is reported empirically; invariance is not assumed."
)
gam_scaling_fit_comparison$probability_tolerance <- 1e-10
gam_scaling_fit_comparison$logLik_tolerance <- 1e-10
gam_scaling_fit_comparison$fit_invariant <- with(gam_scaling_fit_comparison, same_analysis_rows & max_abs_fitted_probability_difference <= probability_tolerance & abs(logLik_difference) <= logLik_tolerance)
gam_scaling_fit_comparison$empirical_comparison_status <- ifelse(gam_scaling_fit_comparison$fit_invariant,"NUMERICALLY EQUIVALENT WITHIN DECLARED TOLERANCES","DIFFERENT BEYOND DECLARED TOLERANCES")
gam_scaling_fit_comparison$note <- "Empirical comparison on the same analysis rows only; not a universal invariance claim for penalized models."
write.csv(gam_scaling_fit_comparison,file.path(tables_dir,"gam_scaling_fit_comparison.csv"),row.names=FALSE)
gam_fit_diagnostics<-imap_dfr(gam_models,function(gm,id)tibble(model_id=id,converged=isTRUE(gm$converged),rank=gm$rank,coefficients=length(coef(gm)),rank_deficient=gm$rank<length(coef(gm)),REML_score=gm$gcv.ubre,logLik=as.numeric(logLik(gm)),deviance_explained=summary(gm)$dev.expl))
write.csv(gam_fit_diagnostics,file.path(tables_dir,"gam_fit_diagnostics.csv"),row.names=FALSE)
saveRDS(c(gam_models,list(GAM_ALL_AVAILABLE_SCALE_COMPARATOR=allavail_gam)),file.path(private_dir,"canonical_and_diagnostic_gam_models.rds"))

gm<-gam_models[["GAM_CANONICAL_K10_SELECT"]]
record_kcheck<-function(fit,id) {
  raw<-deterministic_rows(n_kcheck,id,function(i)as.numeric(t(mgcv::k.check(fit,n.rep=1))),n_cores)
  term_names<-vapply(fit$smooth,function(z)z$label,character(1));nterms<-length(term_names)
  first<-matrix(raw[1,],nrow=nterms,byrow=TRUE,dimnames=list(term_names,c("k'","edf","k-index","p-value")))
  for(j in seq_len(nterms)) {
    pdraw<-raw[,4*j];ok<-is.finite(pdraw)
    first[j,4]<-if(any(ok))mean(pdraw[ok])else NA_real_
    resampling_registry[[paste(id,j)]]<<-tibble(analysis_id=id,metric=term_names[j],attempted=n_kcheck,valid=sum(ok),invalid=sum(!ok),extreme=if(any(ok))sum(pdraw[ok])else NA_real_,extreme_applicable=TRUE,diagnostic_p=first[j,4],method_review_required=any(!ok),estimand="mgcv k-index diagnostic empirical fraction; indexed one-replicate calls; not a primary add-one permutation test")
  }
  as.data.frame(first)
}
kc<-record_kcheck(gm,"gam.kcheck.canonical");kc$term<-rownames(kc);write.csv(kc,file.path(tables_dir,"gam_k_diagnostics_canonical.csv"),row.names=FALSE)
cv<-as.data.frame(mgcv::concurvity(gm,full=TRUE));cv$measure<-rownames(cv);write.csv(cv,file.path(tables_dir,"gam_concurvity_canonical.csv"),row.names=FALSE)

# Row-order/tie sensitivity of neighbour-based k diagnostics.
orders<-list(original=seq_len(nrow(dm)),reverse=rev(seq_len(nrow(dm))))
set_replicate_seed("gam.roworder.shuffle",1);orders$deterministic_shuffle<-sample.int(nrow(dm))
gam_row_order<-imap_dfr(orders,function(ord,label){g<-fit_gam(gam_formula,TRUE,dm[ord,]);k<-record_kcheck(g,paste0("gam.kcheck.",label));k$term<-rownames(k);as_tibble(k)%>%mutate(order=label,logLik=as.numeric(logLik(g)),converged=isTRUE(g$converged),.before=1)})
write.csv(gam_row_order,file.path(tables_dir,"gam_k_diagnostic_row_order_sensitivity.csv"),row.names=FALSE)

# Simulation-based Bernoulli residual checks, conditional on the canonical fit.
gam_p<-fitted(gm);obs_brier<-mean((dm$y-gam_p)^2);obs_events<-sum(dm$y);obs_pearson<-sum((dm$y-gam_p)^2/(gam_p*(1-gam_p)))
sim_diag<-deterministic_rows(n_gam_sim,"gam.bernoulli.simulation",function(i){ys<-rbinom(length(gam_p),1,gam_p);c(events=sum(ys),Brier=mean((ys-gam_p)^2),Pearson=sum((ys-gam_p)^2/(gam_p*(1-gam_p))))},cores=n_cores)
colnames(sim_diag)<-c("events","Brier","Pearson")
gam_sim_summary<-tibble(metric=c("events","Brier","Pearson"),observed=c(obs_events,obs_brier,obs_pearson),simulated_mean=colMeans(sim_diag),simulated_sd=apply(sim_diag,2,sd),simulated_q025=apply(sim_diag,2,quantile,.025),simulated_q975=apply(sim_diag,2,quantile,.975),two_sided_simulation_p=c((1+sum(abs(sim_diag[,1]-mean(sim_diag[,1]))>=abs(obs_events-mean(sim_diag[,1]))))/(n_gam_sim+1),(1+sum(abs(sim_diag[,2]-mean(sim_diag[,2]))>=abs(obs_brier-mean(sim_diag[,2]))))/(n_gam_sim+1),(1+sum(abs(sim_diag[,3]-mean(sim_diag[,3]))>=abs(obs_pearson-mean(sim_diag[,3]))))/(n_gam_sim+1)),replicates=n_gam_sim)
write.csv(gam_sim_summary,file.path(tables_dir,"gam_bernoulli_simulation_diagnostics.csv"),row.names=FALSE);saveRDS(sim_diag,file.path(private_dir,"gam_bernoulli_simulation_nulls.rds"))

# GAM uncertainty slices with aggregate support counts, no participant points.
gx<-seq(quantile(dm$Mean_Score_z,.025),quantile(dm$Mean_Score_z,.975),length.out=140);gy<-c(-1,0,1)
gnd<-expand_grid(Mean_Score_z=gx,SIPS_Negativ_Gesamt_z=gy,COGDIS_score_z=0,SIPS_Positiv_Gesamt_z=0,cohort_f=factor("NAPLS-3",levels=levels(dm$cohort_f)))
gp<-predict(gm,gnd,type="link",se.fit=TRUE);gnd<-gnd%>%mutate(probability=plogis(gp$fit),lo=plogis(gp$fit-qnorm(.975)*gp$se.fit),hi=plogis(gp$fit+qnorm(.975)*gp$se.fit),support_count=map2_int(Mean_Score_z,SIPS_Negativ_Gesamt_z,~sum(abs(dm$Mean_Score_z-.x)<=.25&abs(dm$SIPS_Negativ_Gesamt_z-.y)<=.25)),supported_region=support_count>=10,SIPS_N=factor(SIPS_Negativ_Gesamt_z,levels=c(-1,0,1),labels=c("Low (-1 SD)","Mean","High (+1 SD)")),model_id="GAM_CANONICAL_K10_SELECT",reference_cohort="NAPLS-3")
write.csv(gnd,file.path(figdata_dir,"gam_sips_n_uncertainty_slices.csv"),row.names=FALSE)

# GAM support segmentation: no line or ribbon can span an unsupported grid point.
seg <- gnd %>% arrange(SIPS_N,Mean_Score_z) %>% group_by(SIPS_N) %>%
  mutate(grid_order=row_number(),segment_number=segment_support(supported_region),segment_id=paste0(as.integer(SIPS_N),"_",segment_number)) %>%
  group_by(SIPS_N,segment_id) %>% mutate(segment_size=n(),singleton_segment=n()==1L) %>% ungroup()
segment_summary <- seg %>% group_by(SIPS_N,segment_id,supported_region) %>% summarise(segment_size=n(),x_min=min(Mean_Score_z),x_max=max(Mean_Score_z),singleton=first(singleton_segment),.groups="drop")
gap_counts <- seg %>% group_by(SIPS_N) %>% summarise(internal_unsupported_points=sum(!supported_region & Mean_Score_z>min(Mean_Score_z[supported_region]) & Mean_Score_z<max(Mean_Score_z[supported_region])),.groups="drop")
stopifnot(all(seg %>% group_by(SIPS_N,segment_id) %>% summarise(n_status=n_distinct(supported_region),.groups="drop") %>% pull(n_status)==1L))
write.csv(seg,file.path(figdata_dir,"gam_sips_n_uncertainty_segmented.csv"),row.names=FALSE)
write.csv(left_join(segment_summary,gap_counts,by="SIPS_N"),file.path(tables_dir,"gam_support_segmentation_validation.csv"),row.names=FALSE)

resampling_registry[["gam.bernoulli"]]<-gam_sim_summary%>%transmute(analysis_id="gam.bernoulli.simulation",metric,attempted=replicates,valid=replicates,invalid=0L,extreme=round(two_sided_simulation_p*(replicates+1)-1),p_add_one=two_sided_simulation_p,extreme_applicable=TRUE,method_review_required=FALSE,estimand="conditional Bernoulli simulation around canonical fitted probabilities")
