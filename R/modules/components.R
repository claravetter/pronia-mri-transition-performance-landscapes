comp_cols <- paste0("Comp", 1:4)
stopifnot(all(vapply(d[comp_cols], function(z) all(is.finite(z)), logical(1))))

mapping <- tibble(
  exported_component_identifier=c("component 01","component 02","component 04","component 06"),
  analysis_column=comp_cols,
  manuscript_display_label=paste("Comp",1:4),
  anatomical_map=c("Component_comp01_CVRabsGT3_AND_SBCGT1p3.nii","Component_comp02_CVRabsGT3_AND_SBCGT1p3.nii","Component_comp04_CVRabsGT3_AND_SBCGT1p3.nii","Component_comp06_CVRabsGT3_AND_SBCGT1p3.nii"),
  map_verified_present=file.exists(file.path(map_root,anatomical_map)),
  mapping_basis="NeuroMiner export order and component-map filenames"
)
# Maps are fixed upstream inputs for anatomical rendering, not required for score associations.
write.csv(mapping,file.path(tables_dir,"component_identifier_mapping.csv"),row.names=FALSE)

import_audit <- bind_rows(
  d %>% group_by(cohort) %>% summarise(participants=n(),transitions=sum(EXP_LABEL),across(all_of(comp_cols),~sum(is.finite(.x)),.names="{.col}_available"),.groups="drop"),
  d %>% summarise(cohort="ALL",participants=n(),transitions=sum(EXP_LABEL),across(all_of(comp_cols),~sum(is.finite(.x)),.names="{.col}_available"))
) %>% mutate(unique_participant_cohort_keys=!anyDuplicated(d[c("cohort","Cases")]),visit_field_available=any(str_detect(names(d),regex("visit|timepoint",ignore_case=TRUE))),
             note=if_else(cohort=="ALL","No separate visit field was required or found in the fixed one-row-per-participant analysis export.",""))
write.csv(import_audit,file.path(tables_dir,"component_import_matching_audit.csv"),row.names=FALSE)

mods <- c(age="age",brainage_sbc="BrainAGE_SBC_corr",cogdis="COGDIS_score",sips_p="SIPS_Positiv_Gesamt",sips_n="SIPS_Negativ_Gesamt",psychosoz="Psychosoz_aequiv")
surface_vars <- c(age="age",brainage_sbc="BrainAGE_SBC_corr",cogdis="COGDIS_score",sips_n="SIPS_Negativ_Gesamt")
pairs <- combn(names(surface_vars),2,simplify=FALSE)
metric_names <- unlist(lapply(tolower(comp_cols),function(cc)paste0(cc,"_",c("mad_all","mad_neg","mad_pos","mean_pos","mean_neg","delta_mean"))))

summarize_association <- function(observed, null, context_id, dimension, moderator=NA_character_, pair=NA_character_) {
  map_dfr(seq_along(observed), function(j) {
    ns <- null_summary(null[,j],observed[[j]],paste0("component.",context_id,".",names(observed)[j]),"Spearman rho","absolute_zero")
    tibble(dimension=dimension,moderator=moderator,pair=pair,component=sub("^(comp[1-4])_.*$","\\1",names(observed)[j]),metric=sub("^comp[1-4]_","",names(observed)[j]),rho=unname(observed[[j]]),
           raw_p=ns$p_add_one,effective_permutations=ns$valid,invalid_permutations=ns$invalid,exceedances=ns$exceedances,
           null_mean=ns$null_mean,null_sd=ns$null_sd,null_q025=ns$null_q025,null_q50=ns$null_q50,null_q975=ns$null_q975)
  }) %>% mutate(q_BH_table_S5_24=p.adjust(raw_p,"BH"),correction_family=paste0("24 component-summary tests within ",context_id))
}

association_rows <- list(); null_summaries <- list(); window_rows <- list(); private_nulls <- list()
for (mi in seq_along(mods)) {
  lab <- names(mods)[mi]; v <- mods[[mi]]
  x <- d %>% filter(is.finite(.data[[v]]),is.finite(Mean_Score),!is.na(EXP_LABEL))
  g <- window_grid(x[[v]],.20,.05)
  message("1D component associations: ",lab)
  obs <- component_1d(x,v,x[[v]],g)
  null <- deterministic_rows(n_perm,paste0("component1d.within.",lab),function(b){mv<-permute_within(x[[v]],x$cohort);component_rho_1d_fast(x,mv,g)},cores=n_cores)
  colnames(null) <- names(obs$rho)
  context <- paste0("1d.",lab)
  a <- summarize_association(obs$rho,null,context,"1D",moderator=lab)
  association_rows[[context]] <- a
  null_summaries[[context]] <- map_dfr(seq_along(obs$rho),function(j) null_summary(null[,j],unname(obs$rho[[j]]),paste0("component.",context,".",names(obs$rho)[j]),"Spearman rho","absolute_zero",0))
  window_rows[[context]] <- obs$rows %>% mutate(dimension="1D",moderator=lab,pair=NA_character_,analysis_N=nrow(x),analysis_events=sum(x$EXP_LABEL),.before=1)
  private_nulls[[context]] <- null
}

for (pi in seq_along(pairs)) {
  labs <- pairs[[pi]]; v1 <- surface_vars[[labs[1]]]; v2 <- surface_vars[[labs[2]]]; pair <- paste(labs,collapse="_x_")
  x <- d %>% filter(is.finite(.data[[v1]]),is.finite(.data[[v2]]),is.finite(Mean_Score),!is.na(EXP_LABEL))
  g1 <- window_grid(x[[v1]],.20,.05); g2 <- window_grid(x[[v2]],.20,.05)
  message("2D component associations: ",pair)
  obs <- component_2d(x,x[[v1]],x[[v2]],g1,g2)
  null <- deterministic_rows(n_perm,paste0("component2d.within.",pair),function(b){idx<-seq_len(nrow(x));for(co in unique(x$cohort)){ii<-which(x$cohort==co);idx[ii]<-sample_indices(ii)};component_rho_2d_fast(x,x[[v1]][idx],x[[v2]][idx],g1,g2)},cores=n_cores)
  colnames(null) <- names(obs$rho)
  context <- paste0("2d.",pair)
  a <- summarize_association(obs$rho,null,context,"2D",pair=pair)
  association_rows[[context]] <- a
  null_summaries[[context]] <- map_dfr(seq_along(obs$rho),function(j) null_summary(null[,j],unname(obs$rho[[j]]),paste0("component.",context,".",names(obs$rho)[j]),"Spearman rho","absolute_zero",0))
  window_rows[[context]] <- obs$rows %>% mutate(dimension="2D",moderator=NA_character_,pair=pair,analysis_N=nrow(x),analysis_events=sum(x$EXP_LABEL),.before=1)
  private_nulls[[context]] <- null
}

table_s5 <- bind_rows(association_rows)
stopifnot(nrow(table_s5)==288L,all(table_s5$effective_permutations+table_s5$invalid_permutations==n_perm))
write.csv(table_s5,file.path(tables_dir,"table_S5_component_landscape_associations.csv"),row.names=FALSE)
write.csv(bind_rows(null_summaries),file.path(tables_dir,"permutation_null_summaries_components.csv"),row.names=FALSE)
all_component_rows <- bind_rows(window_rows)
write.csv(all_component_rows,file.path(figdata_dir,"component_window_and_cell_summaries.csv"),row.names=FALSE)
saveRDS(private_nulls,file.path(private_dir,"component_permutation_nulls.rds"))

# Make the maximum-|rho| display rule and the selected statistic explicit for every context/component cell.
metric_labels <- c(mad_all="Overall MAD",mad_neg="Non-transition MAD",mad_pos="Transition MAD",mean_pos="Transition mean",mean_neg="Non-transition mean",delta_mean="Mean difference")
context_labels <- c(age="Age",brainage_sbc="BrainAGE",cogdis="COGDIS",sips_p="SIPS-P",sips_n="SIPS-N",psychosoz="Functioning")
selected_key <- table_s5 %>% group_by(dimension,moderator,pair,component) %>% slice_max(abs(rho),n=1,with_ties=FALSE) %>% ungroup() %>% mutate(selection_rule="largest absolute observed Spearman rho among six component summaries; inferential p and q use the full permutation and 24-test context family",display_metric=unname(metric_labels[metric]),display_component=recode(component,comp1="Comp 1",comp2="Comp 2",comp3="Comp 3",comp4="Comp 4"),display_context=if_else(dimension=="1D",unname(context_labels[moderator]),map_chr(str_split(pair,"_x_"),~paste(unname(context_labels[.x]),collapse=" × "))),display_label=paste0(display_metric,"; rho=",sprintf("%.2f",rho)))
write.csv(selected_key,file.path(tables_dir,"component_selected_statistic_key.csv"),row.names=FALSE)
subject_cor <- imap_dfr(mods,function(v,lab)map_dfr(comp_cols,function(cc){ok<-is.finite(d[[v]])&is.finite(d[[cc]]);ct<-suppressWarnings(cor.test(d[[v]][ok],d[[cc]][ok],method="spearman",exact=FALSE));tibble(moderator=lab,component=tolower(cc),N=sum(ok),rho=unname(ct$estimate),p=ct$p.value)})) %>% mutate(q_BH_24=p.adjust(p,"BH"))
score_cor <- map_dfr(comp_cols,function(cc){ct<-suppressWarnings(cor.test(d$Mean_Score,d[[cc]],method="spearman",exact=FALSE));tibble(component=tolower(cc),N=nrow(d),rho=unname(ct$estimate),p=ct$p.value)}) %>% mutate(q_BH_4=p.adjust(p,"BH"))
write.csv(subject_cor,file.path(tables_dir,"subject_level_component_moderator_correlations.csv"),row.names=FALSE)
write.csv(score_cor,file.path(tables_dir,"score_component_correlations.csv"),row.names=FALSE)

runtime <- tibble(field=c("run_id","R_version","platform","RNGkind","normal_kind","sample_kind","permutations","cores","p_value_rule","seed_scheme","component_mapping"),value=c(run_id,R.version.string,R.version$platform,RNG_KIND,RNG_NORMAL_KIND,RNG_SAMPLE_KIND,n_perm,n_cores,"two-sided absolute-rho add-one using valid replicates","deterministic seed keyed to analysis ID and replicate; nextRNGStream indexed replicates; mc.set.seed=FALSE","NeuroMiner components 01/02/04/06 map to display Comp 1/2/3/4"))
write.csv(runtime,file.path(tables_dir,"component_runtime_rng_manifest.csv"),row.names=FALSE)
