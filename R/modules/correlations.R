# Prespecified S2 moderator correlations and S3 MRI-score correlations.
cor_families<-list(original_s2_s3=c(age="age",sex="sex",cogdis="COGDIS_score",sips_p="SIPS_Positiv_Gesamt",sips_n="SIPS_Negativ_Gesamt",functioning="Psychosoz_aequiv",brainage="BrainAGE_SBC_corr"))
correlations<-imap_dfr(cor_families,function(vars,family)map_dfr(c("pearson","spearman"),function(method){z<-map_dfr(combn(names(vars),2,simplify=FALSE),function(pair){v1<-vars[[pair[1]]];v2<-vars[[pair[2]]];x<-d%>%filter(is.finite(.data[[v1]]),is.finite(.data[[v2]]));ct<-suppressWarnings(cor.test(x[[v1]],x[[v2]],method=method,exact=FALSE));tibble(family=family,method=method,variable_1=pair[1],variable_2=pair[2],N=nrow(x),estimate=unname(ct$estimate),p=ct$p.value)});z%>%mutate(q_BH_21=p.adjust(p,"BH"))}))


# Separate 21 moderator-pair and seven MRI-score-pair families for each method.
s2 <- correlations %>% mutate(family="S2_ORIGINAL_MODERATOR_21")
stopifnot(nrow(filter(s2,method=="pearson"))==21L,nrow(filter(s2,method=="spearman"))==21L)
write.csv(s2,file.path(tables_dir,"correlation_family_S2_original_21.csv"),row.names=FALSE)
s3_vars <- c(age="age",sex="sex",BrainAGE="BrainAGE_SBC_corr",COGDIS="COGDIS_score",`SIPS-P`="SIPS_Positiv_Gesamt",`SIPS-N`="SIPS_Negativ_Gesamt",functioning="Psychosoz_aequiv")
s3 <- map_dfr(c("pearson","spearman"),function(method) map_dfr(names(s3_vars),function(label){
  variable <- s3_vars[[label]]; ok<-is.finite(d$Mean_Score)&is.finite(d[[variable]])
  ct<-suppressWarnings(cor.test(d$Mean_Score[ok],d[[variable]][ok],method=method,exact=FALSE))
  tibble(family="S3_MRI_RISK_SCORE_7",method=method,variable_1="MRI risk score",variable_2=label,N=sum(ok),estimate=unname(ct$estimate),p=ct$p.value)
})) %>% group_by(method) %>% mutate(q_BH_7=p.adjust(p,"BH"),correction_family=paste0("seven MRI-risk-score pairs; ",method)) %>% ungroup()
write.csv(s3,file.path(tables_dir,"correlation_family_S3_mri_risk_score_7.csv"),row.names=FALSE)
