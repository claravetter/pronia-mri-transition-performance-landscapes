read_study_config <- function(path) {
  cfg<-yaml::read_yaml(path)
  # Scientific choices are fixed by the accepted methods. Reject unsupported changes rather than silently ignore them.
  stopifnot(cfg$reference_cohort=="NAPLS-3",cfg$landscape$primary_width==.20,cfg$landscape$primary_step==.05,
    cfg$landscape$minimum_events==7,cfg$landscape$minimum_non_events==7,cfg$landscape$sensitivity_support==10,
    identical(as.numeric(cfg$landscape$sensitivity_widths),c(.15,.20,.25,.30)),cfg$landscape$sensitivity_step_divisor==4,
    cfg$landscape$interval=="left_closed_right_open_final_right_closed",cfg$landscape$tolerance==1e-12,
    cfg$components$contexts==12,cfg$components$comparisons==288,cfg$components$BH_family_size==24,
    identical(unlist(cfg$components$mapping),c("COMP01","COMP02","COMP04","COMP06")),
    cfg$calibration$clip_epsilon==1e-6,cfg$calibration$pooled_apparent=="refit_each_participant_bootstrap",
    cfg$calibration$cohort_apparent=="conditional_assessment_of_fixed_full_cohort_update",
    cfg$calibration$heldout=="conditional_assessment_of_fixed_predictions",
    cfg$components$selected_statistic=="maximum_absolute_rho_among_six_summaries",
    cfg$permutations$rng=="L_Ecuyer_CMRG",cfg$permutations$replicate_stream=="nextRNGStream",
    cfg$permutations$analysis_stream=="SHA256_analysis_id_and_seed_first_28_bits",
    cfg$permutations$validity=="conditional_on_finite_draws_with_method_review_if_any_invalid")
  wanted<-c(permutations=5000,local_auc_bootstrap=2000,pooled_calibration_bootstrap=2000,source_cohort_heldout_bootstrap=500,loco_bootstrap=2000,fisher_simulations=20000,gam_kcheck=400,gam_simulations=1000)
  for(k in names(wanted)) if(cfg$full[[k]]!=wanted[[k]])stop("Configured full-mode resample setting changed: ",k)
  for(mode in c("full","synthetic"))for(k in names(wanted))if(!is.finite(cfg[[mode]][[k]])||cfg[[mode]][[k]]<2||cfg[[mode]][[k]]!=as.integer(cfg[[mode]][[k]]))stop("Invalid resample count: ",k)
  cfg
}
check_analysis_counts<-function(symptom,functioning,cfg,mode,tables_dir) {
  a<-tibble(population=c("symptom_complete","functioning_complete"),N_actual=c(nrow(symptom),nrow(functioning)),events_actual=c(sum(symptom$EXP_LABEL),sum(functioning$EXP_LABEL)))
  a$N_expected<-vapply(a$population,function(k)cfg$populations[[k]]$N,numeric(1));a$events_expected<-vapply(a$population,function(k)cfg$populations[[k]]$events,numeric(1))
  a$status<-if(mode=="synthetic")"SYNTHETIC_NOT_EMPIRICAL_COMPARISON" else ifelse(a$N_actual==a$N_expected&a$events_actual==a$events_expected,"MATCH","INPUT_RECONCILIATION_REQUIRED")
  write.csv(a,file.path(tables_dir,"benchmark_population_comparison.csv"),row.names=FALSE)
  if(mode=="full"&&any(a$status!="MATCH"))stop("Complete-case population counts do not match the configuration; see benchmark_population_comparison.csv")
}
synthetic_data<-function(cfg) {
  set_replicate_seed("synthetic_data",1,cfg$synthetic$seed)
  n<-cfg$synthetic$participants_per_cohort*5L
  cohort<-rep(names(cfg$cohorts),each=cfg$synthetic$participants_per_cohort)
  latent<-runif(n)
  covariate<-function()pmin(1,pmax(0,.85*latent+.15*runif(n)))
  age<-16+24*covariate();brain<-(-8+16*covariate())
  score<-rnorm(n)+.3*latent;y<-rbinom(n,1,plogis(-.6+.5*score+.3*latent))
  d<-tibble(Cases=sprintf("SYNTHETIC_%05d",seq_len(n)),cohort=cohort,EXP_LABEL=y,PRED_LABEL=as.integer(score>0),Mean_Score=score,Probs_platt=plogis(score),age=age,sex=sample(c(1,2),n,replace=TRUE),Studygroup="CHR",COGDIS_score=20*covariate(),SIPS_Positiv_Gesamt=30*covariate(),SIPS_Negativ_Gesamt=25*covariate(),Psychosoz_aequiv=40+60*covariate(),BrainAGE_corr=brain,BrainAGE_SBC_corr=brain,PredictedAge_SBC_corr=age+brain)
  for(j in 1:4)d[[paste0("Comp",j)]]<-rnorm(n)+.2*score
  d$Psychosoz_aequiv[seq(3,n,5)]<-NA_real_
  # Missingness exercises population selection without reproducing study records.
  d$COGDIS_score[seq(7,n,23)]<-NA_real_
  d
}
