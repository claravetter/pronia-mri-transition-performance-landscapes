# Primary and joint width/step landscape analyses.
widths <- unlist(cfg$landscape$sensitivity_widths); one_rows<-list();one_sens<-list();landscape_nulls<-list()
for (mi in seq_along(mods)) {
  lab<-names(mods)[mi];v<-mods[[mi]];x<-d%>%filter(is.finite(.data[[v]]),is.finite(Mean_Score),!is.na(EXP_LABEL))
  message("1D landscape: ",lab)
  grids<-setNames(map(widths,~window_grid(x[[v]],.x,.x/4)),as.character(widths))
  obs_all<-map(widths,function(w)window_stats_1d(x,v,w,w/4,7,7,if(w==.20)n_boot_local else 0L,paste0("local1d.",lab,".",w)))
  within<-deterministic_rows(n_perm,paste0("landscape1d.within.",lab),function(b){mv<-permute_within(x[[v]],x$cohort);vapply(grids,function(g)window_dmax(x,v,mv,g),numeric(1))},cores=n_cores)
  global<-deterministic_rows(n_perm,paste0("landscape1d.global.",lab),function(b){mv<-sample(x[[v]],length(x[[v]]),replace=FALSE);vapply(grids,function(g)window_dmax(x,v,mv,g),numeric(1))},cores=n_cores)
  landscape_nulls[[paste0("1d_",lab)]]<-list(within=within,global=global)
  for(wi in seq_along(widths)){
    o<-obs_all[[wi]];good<-filter(o$rows,supported,is.finite(AUC));represented_idx<-which(Reduce(`|`,map(seq_len(nrow(good)),function(i)in_window(x[[v]],good,i)),init=rep(FALSE,nrow(x))))
    s1<-null_summary(within[,wi],o$Dmax,paste0("landscape1d.within.",lab,".w",widths[wi]),"Dmax","upper")
    sg<-null_summary(global[,wi],o$Dmax,paste0("landscape1d.global.",lab,".w",widths[wi]),"Dmax","upper")
    validation_summaries[[s1$analysis_id]]<-s1;validation_summaries[[sg$analysis_id]]<-sg
    one_sens[[length(one_sens)+1L]]<-tibble(moderator=lab,analysis_N=nrow(x),analysis_events=sum(x$EXP_LABEL),width_fraction=widths[wi],step_fraction=widths[wi]/4,
      final_start=tail(grids[[wi]]$start,1),final_end=tail(grids[[wi]]$end,1),observed_max=max(x[[v]]),final_reaches_max=abs(tail(grids[[wi]]$end,1)-max(x[[v]]))<1e-10,
      n_candidate_windows=nrow(o$rows),n_eligible_windows=nrow(good),represented_N=length(represented_idx),represented_events=sum(x$EXP_LABEL[represented_idx]),
      Dmax=o$Dmax,Dmean=o$Dmean,trajectory_rho=o$rho,within_cohort_p=s1$p_add_one,unrestricted_p=sg$p_add_one,valid_permutations=s1$valid,invalid_permutations=s1$invalid)
    if(widths[wi]==.20){one_rows[[lab]]<-o$rows%>%mutate(moderator=lab,analysis_N=nrow(x),analysis_events=sum(x$EXP_LABEL));population_rows[[length(population_rows)+1L]]<-add_population(paste0("POP_1D_",toupper(lab)),paste("One-dimensional",lab),x,x[represented_idx,],paste("MRI risk score; outcome;",lab))}
  }
}
one_sensitivity<-bind_rows(one_sens)%>%group_by(width_fraction,step_fraction)%>%mutate(within_cohort_q_BH_6=p.adjust(within_cohort_p,"BH"),unrestricted_q_BH_6=p.adjust(unrestricted_p,"BH"))%>%ungroup()
primary_1d<-one_sensitivity%>%filter(width_fraction==.20)
write.csv(primary_1d,file.path(tables_dir,"global_1d_landscape_tests.csv"),row.names=FALSE)
write.csv(one_sensitivity,file.path(tables_dir,"joint_width_step_sensitivity_1d.csv"),row.names=FALSE)
write.csv(bind_rows(one_rows),file.path(figdata_dir,"primary_1d_window_details.csv"),row.names=FALSE)

# Six primary surfaces and selected-pair joint width/step sensitivity.
surface_rows<-list();surface_tests<-list();surface_nulls<-list()
for(pi in seq_along(pairs)){
  labs<-pairs[[pi]];v1<-surface_vars[[labs[1]]];v2<-surface_vars[[labs[2]]];pair<-paste(labs,collapse="_x_")
  x<-d%>%filter(is.finite(.data[[v1]]),is.finite(.data[[v2]]),is.finite(Mean_Score),!is.na(EXP_LABEL))
  message("2D landscape: ",pair)
  obs<-surface_stats(x,v1,v2,.20,.05,7,7,n_boot_local,paste0("local2d.",pair));g1<-obs$g1;g2<-obs$g2
  within<-deterministic_rows(n_perm,paste0("landscape2d.within.",pair),function(b){idx<-seq_len(nrow(x));for(co in unique(x$cohort)){ii<-which(x$cohort==co);idx[ii]<-sample_indices(ii)};surface_dmax(x,x[[v1]][idx],x[[v2]][idx],g1,g2)},cores=n_cores)[,1]
  global<-deterministic_rows(n_perm,paste0("landscape2d.global.",pair),function(b){idx<-sample.int(nrow(x));surface_dmax(x,x[[v1]][idx],x[[v2]][idx],g1,g2)},cores=n_cores)[,1]
  surface_nulls[[pair]]<-list(within=within,global=global);sw<-null_summary(within,obs$Dmax,paste0("landscape2d.within.",pair),"Dmax","upper");sg<-null_summary(global,obs$Dmax,paste0("landscape2d.global.",pair),"Dmax","upper")
  validation_summaries[[sw$analysis_id]]<-sw;validation_summaries[[sg$analysis_id]]<-sg
  good<-filter(obs$rows,supported,is.finite(AUC));represented<-Reduce(`|`,map(seq_len(nrow(good)),function(k)in_window(x[[v1]],good%>%transmute(start=start_1,end=end_1,right_closed=right_closed_1),k)&in_window(x[[v2]],good%>%transmute(start=start_2,end=end_2,right_closed=right_closed_2),k)),init=rep(FALSE,nrow(x)))
  surface_tests[[pair]]<-tibble(pair=pair,moderator_1=labs[1],moderator_2=labs[2],analysis_N=nrow(x),analysis_events=sum(x$EXP_LABEL),n_eligible_cells=nrow(good),represented_N=sum(represented),represented_events=sum(x$EXP_LABEL[represented]),Dmax=obs$Dmax,Dmean=obs$Dmean,within_cohort_p=sw$p_add_one,unrestricted_p=sg$p_add_one,valid_permutations=sw$valid,invalid_permutations=sw$invalid)
  surface_rows[[pair]]<-obs$rows%>%mutate(pair=pair,moderator_1=labs[1],moderator_2=labs[2],analysis_N=nrow(x),analysis_events=sum(x$EXP_LABEL))
  population_rows[[length(population_rows)+1L]]<-add_population(paste0("POP_2D_",toupper(labs[1]),"_X_",toupper(labs[2])),paste("Two-dimensional",pair),x,x[represented,],paste("MRI risk score; outcome;",labs[1],labs[2]))
}
primary_2d<-bind_rows(surface_tests)%>%mutate(within_cohort_q_BH_6=p.adjust(within_cohort_p,"BH"),unrestricted_q_BH_6=p.adjust(unrestricted_p,"BH"))
write.csv(primary_2d,file.path(tables_dir,"global_2d_surface_tests.csv"),row.names=FALSE)
all_surface_rows<-bind_rows(surface_rows);write.csv(all_surface_rows,file.path(figdata_dir,"primary_2d_cell_details.csv"),row.names=FALSE)

# Highlighted surface width/step and support/population checks.
v1<-"SIPS_Negativ_Gesamt";v2<-"BrainAGE_SBC_corr";x<-d%>%filter(is.finite(.data[[v1]]),is.finite(.data[[v2]]),is.finite(Mean_Score),!is.na(EXP_LABEL))
grids1<-map(widths,~window_grid(x[[v1]],.x,.x/4));grids2<-map(widths,~window_grid(x[[v2]],.x,.x/4));obsw<-map(widths,~surface_stats(x,v1,v2,.x,.x/4,7,7,0L))
nullw<-deterministic_rows(n_perm,"surface_width_sensitivity.sips_n_x_brainage_sbc",function(b){idx<-seq_len(nrow(x));for(co in unique(x$cohort)){ii<-which(x$cohort==co);idx[ii]<-sample_indices(ii)};vapply(seq_along(widths),function(k)surface_dmax(x,x[[v1]][idx],x[[v2]][idx],grids1[[k]],grids2[[k]]),numeric(1))},cores=n_cores)
saveRDS(nullw,file.path(private_dir,"selected_surface_width_step_permutation_nulls.rds"))
surface_width<-map_dfr(seq_along(widths),function(k){ns<-null_summary(nullw[,k],obsw[[k]]$Dmax,paste0("surface_width.sips_n_x_brainage_sbc.",widths[k]),"Dmax","upper");validation_summaries[[ns$analysis_id]]<<-ns;tibble(pair="sips_n_x_brainage_sbc",width_fraction=widths[k],step_fraction=widths[k]/4,final_end_1=tail(grids1[[k]]$end,1),max_1=max(x[[v1]]),final_end_2=tail(grids2[[k]]$end,1),max_2=max(x[[v2]]),both_final_intervals_reach_max=abs(tail(grids1[[k]]$end,1)-max(x[[v1]]))<1e-10&&abs(tail(grids2[[k]]$end,1)-max(x[[v2]]))<1e-10,Dmax=obsw[[k]]$Dmax,p=ns$p_add_one,valid=ns$valid,invalid=ns$invalid)})%>%mutate(q_BH_4_selected_pair_width_step=p.adjust(p,"BH"),family_scope="four joint width/step settings for preselected SIPS-N x BrainAGE pair; does not correct surface selection")
write.csv(surface_width,file.path(tables_dir,"joint_width_step_sensitivity_selected_surface.csv"),row.names=FALSE)
surface_support<-bind_rows(
  surface_stats(x,v1,v2,.20,.05,10,10,0L)$rows%>%filter(supported)%>%summarise(population="pooled",support="10/10",N_range=nrow(x),events_range=sum(x$EXP_LABEL),eligible_cells=n(),Dmax=diff(range(AUC)),AUC_min=min(AUC),AUC_max=max(AUC)),
  {z<-filter(x,cohort!="PRONIA");surface_stats(z,v1,v2,.20,.05,7,7,0L)$rows%>%filter(supported)%>%summarise(population="external_only",support="7/7",N_range=nrow(z),events_range=sum(z$EXP_LABEL),eligible_cells=n(),Dmax=diff(range(AUC)),AUC_min=min(AUC),AUC_max=max(AUC))},
  {z<-filter(x,cohort!="PRONIA");surface_stats(z,v1,v2,.20,.05,10,10,0L)$rows%>%filter(supported)%>%summarise(population="external_only",support="10/10",N_range=nrow(z),events_range=sum(z$EXP_LABEL),eligible_cells=n(),Dmax=if(n()>=3)diff(range(AUC))else NA_real_,AUC_min=if(n())min(AUC)else NA_real_,AUC_max=if(n())max(AUC)else NA_real_)}
)
write.csv(surface_support,file.path(tables_dir,"highlighted_surface_population_support_sensitivity.csv"),row.names=FALSE)

# Complete descriptive population and 10/10 coverage for all six 1D moderators.
sets<-list(pooled=d,external_only=filter(d,cohort!="PRONIA"),chr_only_excluding_pronia_rop=chr_sensitivity_population(d),newly_evaluated_external=filter(d,cohort%in%c("NAPLS-3","MUC-FRUE")))
for(co in sort(unique(d$cohort)))sets[[paste0("leave_out_",co)]]<-filter(d,cohort!=co)
population_sensitivity<-imap_dfr(sets,function(dat,pop)imap_dfr(mods,function(v,lab){z<-dat%>%filter(is.finite(.data[[v]]),is.finite(Mean_Score),!is.na(EXP_LABEL));o<-window_stats_1d(z,v,.20,.05,7,7,0L);good<-filter(o$rows,supported);tibble(population=pop,moderator=lab,N=nrow(z),events=sum(z$EXP_LABEL),eligible_windows=nrow(good),Dmax=o$Dmax,AUC_min=if(nrow(good))min(good$AUC)else NA_real_,AUC_max=if(nrow(good))max(good$AUC)else NA_real_,final_interval_reaches_max=abs(tail(o$rows$end,1)-max(z[[v]]))<1e-10)}))
write.csv(population_sensitivity,file.path(tables_dir,"landscape_population_sensitivity_1d_all_six.csv"),row.names=FALSE)
support_sensitivity<-imap_dfr(mods,function(v,lab){z<-d%>%filter(is.finite(.data[[v]]),is.finite(Mean_Score),!is.na(EXP_LABEL));map_dfr(c(7L,10L),function(k){o<-window_stats_1d(z,v,.20,.05,k,k,0L);tibble(moderator=lab,min_events=k,min_non_events=k,analysis_N=nrow(z),analysis_events=sum(z$EXP_LABEL),eligible_windows=sum(o$rows$supported),Dmax=o$Dmax)})})
write.csv(support_sensitivity,file.path(tables_dir,"window_support_sensitivity_1d_all_six.csv"),row.names=FALSE)

# Figure 2 fidelity: include every supported AUC on a 0-1 scale and assert no supported cell is censored.
fig2data<-all_surface_rows%>%filter(pair=="brainage_sbc_x_sips_n")
stopifnot(nrow(fig2data)>0L)
stopifnot(all(fig2data$AUC[fig2data$supported]>=0&fig2data$AUC[fig2data$supported]<=1),all(!fig2data$supported | is.finite(fig2data$AUC)))
write.csv(tibble(check_id="figure2_supported_auc_scale",supported_cells=sum(fig2data$supported),supported_below_0_3=sum(fig2data$supported&fig2data$AUC<.3),all_supported_inside_scale=all(fig2data$AUC[fig2data$supported]>=0&fig2data$AUC[fig2data$supported]<=1),scale_min=0,scale_max=1,missing_colour_reserved_for_unsupported=TRUE),file.path(tables_dir,"figure2_fidelity_assertion.csv"),row.names=FALSE)
