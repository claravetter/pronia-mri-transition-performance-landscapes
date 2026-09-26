summarise_num <- function(x) {
  ok <- is.finite(x)
  if (!any(ok)) {
    return(c(N = 0, mean = NA_real_, SD = NA_real_, median = NA_real_,
             IQR = NA_real_, min = NA_real_, max = NA_real_))
  }
  z <- x[ok]
  c(N = length(z), mean = mean(z), SD = sd(z), median = median(z),
    IQR = IQR(z), min = min(z), max = max(z))
}
baseline_vars <- c(age = "age", mri_score = "Mean_Score", cogdis = "COGDIS_score", sips_p = "SIPS_Positiv_Gesamt",
                   sips_n = "SIPS_Negativ_Gesamt", functioning = "Psychosoz_aequiv", brainage = "BrainAGE_SBC_corr")
baseline <- d %>% group_by(cohort) %>% group_modify(~{
  z <- .x
  vals <- unlist(lapply(baseline_vars, function(v) setNames(summarise_num(z[[v]]), paste(names(summarise_num(z[[v]])), v, sep = "_"))))
  bind_cols(tibble(participants = nrow(z), transitions = sum(z$EXP_LABEL), non_transitions = sum(z$EXP_LABEL == 0L),
                   female_N = sum(z$sex == 2, na.rm = TRUE), male_N = sum(z$sex == 1, na.rm = TRUE)), as_tibble_row(vals))
}) %>% ungroup()
write.csv(baseline, file.path(tables_dir, "baseline_descriptive_by_cohort.csv"), row.names = FALSE)

group_tests <- map_dfr(names(baseline_vars), function(label) {
  v <- baseline_vars[[label]]
  x <- d %>% filter(is.finite(.data[[v]]))
  kt <- kruskal.test(x[[v]] ~ x$cohort)
  tibble(variable = label, test = "Kruskal-Wallis", N = nrow(x), statistic = unname(kt$statistic), df = unname(kt$parameter), p = kt$p.value)
}) %>% mutate(q_BH_continuous_7 = p.adjust(p, "BH"))
write.csv(group_tests, file.path(tables_dir, "baseline_group_tests_continuous.csv"), row.names = FALSE)


# Missingness and BrainAGE clinical transport diagnostics.
audit_vars<-c("EXP_LABEL","Mean_Score","Probs_platt","PRED_LABEL","age","sex","COGDIS_score","SIPS_Positiv_Gesamt","SIPS_Negativ_Gesamt","Psychosoz_aequiv","BrainAGE_corr","BrainAGE_SBC_corr","Comp1","Comp2","Comp3","Comp4")
miss<-d%>%group_by(cohort)%>%summarise(across(all_of(audit_vars),list(N_available=~sum(!is.na(.x)),N_missing=~sum(is.na(.x))),.names="{.col}__{.fn}"),.groups="drop")%>%pivot_longer(-cohort,names_to=c("variable","measure"),names_sep="__")%>%pivot_wider(names_from=measure,values_from=value)%>%mutate(missingness_type=case_when(cohort%in%c("FePsy","MUC-FRUE")&variable%in%c("COGDIS_score","SIPS_Positiv_Gesamt","SIPS_Negativ_Gesamt","Psychosoz_aequiv")~"structurally unavailable by design",TRUE~"incidental/unknown"))
write.csv(miss,file.path(tables_dir,"missingness_by_cohort.csv"),row.names=FALSE)

brainage<-map_dfr(sort(unique(d$cohort)),function(co){x<-filter(d,cohort==co,is.finite(age),is.finite(PredictedAge_SBC_corr));err<-x$PredictedAge_SBC_corr-x$age;fit<-lm(PredictedAge_SBC_corr~age,data=x);tibble(cohort=co,N=nrow(x),age_min=min(x$age),age_max=max(x$age),MAE=mean(abs(err)),RMSE=sqrt(mean(err^2)),mean_signed_error=mean(err),pearson=cor(x$PredictedAge_SBC_corr,x$age),spearman=cor(x$PredictedAge_SBC_corr,x$age,method="spearman"),predictive_R2=1-sum(err^2)/sum((x$age-mean(x$age))^2),calibration_intercept=coef(fit)[1],calibration_slope=coef(fit)[2],residual_age_association=cor(err,x$age))})
write.csv(brainage,file.path(tables_dir,"brainage_external_performance.csv"),row.names=FALSE)
brain_identity<-d%>%filter(is.finite(age),is.finite(PredictedAge_SBC_corr),is.finite(BrainAGE_SBC_corr))%>%group_by(cohort)%>%summarise(N=n(),max_abs_identity_error=max(abs(PredictedAge_SBC_corr-age-BrainAGE_SBC_corr)),mean_signed_error=mean(PredictedAge_SBC_corr-age),brainage_column_mean=mean(BrainAGE_SBC_corr),rmse=sqrt(mean((PredictedAge_SBC_corr-age)^2)),.groups="drop")
write.csv(brain_identity,file.path(tables_dir,"brainage_identity_checks.csv"),row.names=FALSE)
