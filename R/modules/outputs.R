# Population ledger, coverage, uncertainty and correction registries.
symptom_vars<-c("Mean_Score","COGDIS_score","SIPS_Positiv_Gesamt","SIPS_Negativ_Gesamt")
symptom<-d%>%filter(if_all(all_of(c("EXP_LABEL",symptom_vars,"cohort")),~!is.na(.x)&if(is.numeric(.x))is.finite(.x)else TRUE))
functioning<-symptom%>%filter(is.finite(Psychosoz_aequiv))
population_rows[[length(population_rows)+1L]]<-add_population("POP_GLM_SYMPTOM","Canonical symptom-complete moderation",symptom,required="MRI risk score; COGDIS; SIPS-P; SIPS-N; outcome; cohort")
population_rows[[length(population_rows)+1L]]<-add_population("POP_GLM_FUNCTIONING","Canonical functioning-complete moderation",functioning,required="symptom model variables plus functioning")
analysis_populations<-bind_rows(population_rows);write.csv(analysis_populations,file.path(tables_dir,"analysis_populations.csv"),row.names=FALSE)

# Correction families with actual model IDs and sizes.
glm_coef <- glm_coefficients
glm_registry <- glm_coef %>% filter(is_score_interaction) %>% count(model_id,name="family_size") %>% transmute(family_id=paste0(model_id,"_INTERACTIONS"),scope=paste0("score interaction terms within ",model_id),family_size,method="BH",status="complete",interpretation="canonical model-specific family")
registry <- bind_rows(
  tribble(~family_id,~scope,~family_size,~method,~status,~interpretation,
          "PRIMARY_1D_6","six primary 1D Dmax tests",6L,"BH","complete","primary family",
          "PRIMARY_2D_6","six primary 2D Dmax tests",6L,"BH","complete","primary family",
          "SENS_1D_WIDTH_STEP_BY_SETTING","six moderators at each of four settings",6L,"BH","complete","four separate sensitivity families",
          "SENS_SELECTED_SURFACE_WIDTH_STEP_4","four settings for selected SIPS-N × BrainAGE surface",4L,"BH","complete","does not correct surface selection",
          "S2_ORIGINAL_PEARSON_21","21 moderator pairs",21L,"BH","complete","original S2 Pearson family",
          "S2_ORIGINAL_SPEARMAN_21","21 moderator pairs",21L,"BH","complete","separate original S2 Spearman family",
          "S3_MRI_SCORE_PEARSON_7","MRI risk score with seven covariates",7L,"BH","complete","original S3 counterpart",
          "S3_MRI_SCORE_SPEARMAN_7","MRI risk score with seven covariates",7L,"BH","complete","separate sensitivity family",
          "SCORE_COMPONENT_4","MRI risk score with four components",4L,"BH","complete","score/component family",
          "SUBJECT_COMPONENT_MODERATOR_24","six moderators × four components",24L,"BH","complete","subject-level family",
          "CATEGORICAL_SEEDED_2","two seeded categorical comparisons",2L,"BH","complete","categorical family",
          "GAM_TENSOR_3","three tensor interactions within each GAM specification",3L,"BH","complete","model-specific; nonlinear result remains exploratory",
          "COMPONENT_CONTEXT_24","four components × six summaries within each context",24L,"BH","complete","one family for each of 12 contexts"),
  glm_registry)
write.csv(registry,file.path(tables_dir,"correction_registry.csv"),row.names=FALSE)

# Six-panel figure data from this run's configured local-bootstrap intervals.
one_d <- bind_rows(one_rows) %>% arrange(moderator,window_id) %>% group_by(moderator) %>% mutate(segment_id=segment_support(supported)) %>% ungroup()
global <- primary_1d
labels <- c(age="Age",brainage_sbc="BrainAGE",cogdis="COGDIS",sips_p="SIPS-P",sips_n="SIPS-N",psychosoz="Functioning")
one_d <- one_d %>% left_join(global %>% select(moderator,within_cohort_p,within_cohort_q_BH_6),by="moderator") %>%
  mutate(display_moderator=factor(unname(labels[moderator]),levels=unname(labels)),panel_note=sprintf("p = %.3f; BH q = %.3f",within_cohort_p,within_cohort_q_BH_6))
stopifnot(all(one_d$bootstrap_valid[one_d$supported]==n_boot_local),all(one_d$bootstrap_invalid==0L))
write.csv(one_d,file.path(figdata_dir,"primary_1d_window_details_figure_ready.csv"),row.names=FALSE)

null_summary_table<-bind_rows(validation_summaries);write.csv(null_summary_table,file.path(tables_dir,"permutation_null_summaries_core.csv"),row.names=FALSE)
saveRDS(landscape_nulls,file.path(private_dir,"landscape_permutation_nulls.rds"));saveRDS(surface_nulls,file.path(private_dir,"surface_permutation_nulls.rds"))

runtime<-tibble(field=c("run_id","R_version","platform","RNGkind","normal_kind","sample_kind","permutations","cores","p_value_rule","seed_scheme","mgcv_version","pROC_version","metafor_version"),value=c(run_id,R.version.string,R.version$platform,RNG_KIND,RNG_NORMAL_KIND,RNG_SAMPLE_KIND,n_perm,n_cores,"add-one using one plus valid replicates","deterministic seed keyed to analysis ID and replicate; nextRNGStream indexed replicates; mc.set.seed=FALSE",as.character(packageVersion("mgcv")),as.character(packageVersion("pROC")),as.character(packageVersion("metafor"))))
write.csv(runtime,file.path(tables_dir,"runtime_rng_manifest.csv"),row.names=FALSE)
message("Core analysis completed: ",run_dir)

# Unified uncertainty/validity registry: all attempts and exclusions remain inspectable.
local_registry<-bind_rows(
  bind_rows(one_rows)%>%transmute(analysis_id=paste0("local1d.",moderator,".",window_id),metric="AUC",attempted=bootstrap_valid+bootstrap_invalid,valid=bootstrap_valid,invalid=bootstrap_invalid),
  all_surface_rows%>%transmute(analysis_id=paste0("local2d.",pair,".",cell_i,".",cell_j),metric="AUC",attempted=bootstrap_valid+bootstrap_invalid,valid=bootstrap_valid,invalid=bootstrap_invalid))%>%mutate(extreme=NA_integer_,extreme_applicable=FALSE,method_review_required=invalid>0,estimand="stratified percentile AUC bootstrap; supported cells only")
permutation_registry<-bind_rows(bind_rows(validation_summaries),bind_rows(null_summaries))%>%rename(extreme=exceedances)%>%mutate(extreme_applicable=TRUE,estimand="declared permutation statistic conditional on finite null draws")
registry_all<-bind_rows(permutation_registry,local_registry,bind_rows(resampling_registry))
write.csv(registry_all,file.path(tables_dir,"resampling_validity_registry.csv"),row.names=FALSE)
uncertainty<-tribble(~family,~method,~resamples,~conditioning,
 "validation","DeLong AUC intervals",0,"fixed imported MRI score",
 "local_AUC","class-stratified participant percentile bootstrap",n_boot_local,"same eligible cases as support counts",
 "pooled_calibration","participant bootstrap; refit apparent update inside replicate",n_boot_cal,"external participants",
 "source_cohort_heldout_calibration","participant bootstrap assessment",n_boot_cohort,"historical probabilities or fixed update/held-out predictions",
 "LOCO","paired cohort-stratified participant bootstrap",n_boot_loco,"frozen upstream MRI and second-stage fold fits",
 "GLM","two-sided Wald, normal pointwise link intervals",0,"canonical analysis-set scaling",
 "GAM","mgcv approximate smooth tests and pointwise link intervals",0,"canonical and bounded sensitivity specifications")
write.csv(uncertainty,file.path(tables_dir,"uncertainty_methods_registry.csv"),row.names=FALSE)
streams<-map_dfr(ls(stream_cache),function(k){z<-get(k,stream_cache);tibble(analysis_key=k,replicates=length(z),initial_state=paste(z[[1]],collapse=":"),last_state=paste(tail(z,1)[[1]],collapse=":"),algorithm="L'Ecuyer-CMRG / nextRNGStream / Inversion / Rejection")})
if(anyDuplicated(streams$initial_state))stop("Analysis RNG seed collision; method review required")
write.csv(streams,file.path(tables_dir,"rng_stream_registry.csv"),row.names=FALSE)
# Explicitly record method issues; a completed computation is not empirical validation.
issues<-registry_all%>%filter(method_review_required%in%TRUE)
write.csv(issues,file.path(tables_dir,"method_review_required.csv"),row.names=FALSE)
