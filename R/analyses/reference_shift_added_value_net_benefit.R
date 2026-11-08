#!/usr/bin/env Rscript
# PRONIA CHR-P reference, score distributions and per-cohort shift, landscape robustness to population
# composition, symptom regression with and without the MRI score per held-out cohort, and decision curves.
# Reuses the published pipeline helpers (analysis_figures_publication) so that
# populations, window definitions and resampling streams match the manuscript.
# Writes aggregate tables and figures only; no participant-level output.
#
# Usage: Rscript R/analyses/reference_shift_added_value_net_benefit.R INPUT.csv OUTPUT_DIR PIPELINE_REPO [cores]

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3) stop("Usage: reference_shift_added_value_net_benefit.R INPUT.csv OUTPUT_DIR PIPELINE_REPO [cores]")
input_path <- normalizePath(args[1]); out_dir <- args[2]; repo <- normalizePath(args[3])
n_cores <- if(length(args)>=4) as.integer(args[4]) else 1L
source(file.path(repo,"R/lib/common.R")); source(file.path(repo,"R/lib/calibration.R"))
suppressPackageStartupMessages(library(patchwork))
options(pronia.seed=20260913L)
n_perm <- 5000L; n_boot <- 2000L; n_boot_cal <- 500L

tables_dir <- file.path(out_dir,"tables"); figures_dir <- file.path(out_dir,"figures")
for(p in c(tables_dir,figures_dir)) dir.create(p,recursive=TRUE,showWarnings=FALSE)

# Canonical population: same duplicate resolution as the manuscript (correction record via PRONIA_CORRECTION_JSON).
raw <- read.csv(input_path,check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(Cases="character"))
cpath <- Sys.getenv("PRONIA_CORRECTION_JSON")  # authorized duplicate-resolution record, kept outside the repository
correction <- if(nzchar(cpath)) jsonlite::read_json(cpath,simplifyVector=TRUE) else NULL
expected <- list(PRONIA=list(N=329,events=26),`NAPLS-3`=list(N=458,events=51),ZInEP=list(N=157,events=23),
                 FePsy=list(N=37,events=16),`MUC-FRUE`=list(N=33,events=16))
d <- canonicalize_data(input_path,NULL,correction,expected)
external <- filter(d,cohort!="PRONIA")
cohort_levels <- c("PRONIA","NAPLS-3","ZInEP","FePsy","MUC-FRUE")
cohort_status <- c(PRONIA="Discovery",`NAPLS-3`="New",ZInEP="Previously used",FePsy="Previously used",`MUC-FRUE`="New")
mods <- c(Age="age",BrainAGE="BrainAGE_SBC_corr",COGDIS="COGDIS_score",`SIPS-P`="SIPS_Positiv_Gesamt",
          `SIPS-N`="SIPS_Negativ_Gesamt",Functioning="Psychosoz_aequiv")
theme_set(theme_minimal(base_size=10)+theme(panel.grid.minor=element_blank(),strip.text=element_text(face="bold")))
save_fig <- function(p,name,w,h) {ggsave(file.path(figures_dir,paste0(name,".pdf")),p,width=w,height=h);ggsave(file.path(figures_dir,paste0(name,".png")),p,width=w,height=h,dpi=300)}
delong <- function(score,y) {roc<-pROC::roc(y,score,direction="<",quiet=TRUE);as.numeric(pROC::ci.auc(roc,method="delong"))}

# ---------------------------------------------------------------------------
# PRONIA CHR-P-only reference
# ---------------------------------------------------------------------------
reference_metrics <- function(x,scope) {
  ci <- delong(x$Mean_Score,x$EXP_LABEL)
  null <- deterministic_rows(n_perm,paste0("a1.auc.",scope),function(i)auc_fast(x$Mean_Score,sample(x$EXP_LABEL)),n_cores)[,1]
  ns <- null_summary(null,ci[2],paste0("a1.auc.",scope),"AUC","absolute_reference",.5)
  tp<-sum(x$PRED_LABEL==1&x$EXP_LABEL==1);fp<-sum(x$PRED_LABEL==1&x$EXP_LABEL==0)
  tn<-sum(x$PRED_LABEL==0&x$EXP_LABEL==0);fn<-sum(x$PRED_LABEL==0&x$EXP_LABEL==1)
  cal <- calibration_bootstrap(x$EXP_LABEL,x$Probs_platt,n_boot_cal,paste0("a1.cal.",scope),"original",n_cores)$summary
  get <- function(m,col) cal[[col]][cal$metric==m]
  tibble(scope=scope,N=nrow(x),events=sum(x$EXP_LABEL),AUC=ci[2],AUC_lo=ci[1],AUC_hi=ci[3],AUC_p=ns$p_add_one,
         Sensitivity=tp/(tp+fn),Specificity=tn/(tn+fp),BACC=bacc_fast(x$PRED_LABEL,x$EXP_LABEL),
         PPV=tp/(tp+fp),NPV=tn/(tn+fn),Brier=get("Brier","estimate"),
         prevalence_reference_Brier=get("prevalence_reference_Brier","estimate"),
         calibration_in_the_large=get("alpha_int","estimate"),
         calibration_slope=get("beta","estimate"),slope_lo=get("beta","lo"),slope_hi=get("beta","hi"))
}
pronia <- filter(d,cohort=="PRONIA")
a1 <- bind_rows(reference_metrics(pronia,"PRONIA all (CHR-P + non-CHR group)"),
                reference_metrics(filter(pronia,Studygroup=="CHR"),"PRONIA CHR-P only"),
                reference_metrics(filter(pronia,Studygroup!="CHR"),"PRONIA non-CHR group (coded ROP)"))
write.csv(a1,file.path(tables_dir,"A1_pronia_reference_by_group.csv"),row.names=FALSE)
message("A1 done")

# ---------------------------------------------------------------------------
# MRI score distributions by cohort; what a per-cohort score shift can explain
# ---------------------------------------------------------------------------
d$cohort_f <- factor(d$cohort,levels=cohort_levels)
d$outcome <- factor(ifelse(d$EXP_LABEL==1,"Transition","No transition"),levels=c("No transition","Transition"))
a4_dist <- d %>% group_by(cohort=cohort_f,outcome) %>%
  summarise(N=n(),mean=mean(Mean_Score),SD=sd(Mean_Score),median=median(Mean_Score),
            q25=quantile(Mean_Score,.25),q75=quantile(Mean_Score,.75),
            mean_original_probability=mean(Probs_platt),.groups="drop")
ref_nt <- filter(d,cohort=="PRONIA",EXP_LABEL==0)$Mean_Score
a4_shift <- d %>% group_by(cohort=cohort_f) %>% group_modify(~{
  nt<-.x$Mean_Score[.x$EXP_LABEL==0];tr<-.x$Mean_Score[.x$EXP_LABEL==1]
  pooled_sd<-sqrt((var(nt)+var(ref_nt))/2)
  tibble(N=nrow(.x),events=length(tr),observed_rate=mean(.x$EXP_LABEL),mean_original_probability=mean(.x$Probs_platt),
         non_transition_shift_vs_PRONIA_SMD=(mean(nt)-mean(ref_nt))/pooled_sd,
         non_transition_SD_ratio_vs_PRONIA=sd(nt)/sd(ref_nt),
         within_cohort_separation_SMD=(mean(tr)-mean(nt))/sqrt((var(tr)+var(nt))/2),
         within_cohort_AUC=auc_fast(.x$Mean_Score,.x$EXP_LABEL))
}) %>% ungroup()
# Outcome-blind per-cohort standardization changes only between-cohort comparisons:
# within-cohort AUC is invariant to any monotone per-cohort transformation of the score.
ext <- external %>% group_by(cohort) %>%
  mutate(score_centered=Mean_Score-median(Mean_Score),score_z=(Mean_Score-mean(Mean_Score))/sd(Mean_Score)) %>% ungroup()
a4_pooled <- map_dfr(c(`Original score`="Mean_Score",`Per-cohort median-centred`="score_centered",`Per-cohort z-scored`="score_z"),
  function(v){ci<-delong(ext[[v]],ext$EXP_LABEL);tibble(AUC=ci[2],AUC_lo=ci[1],AUC_hi=ci[3])},.id="external_pooled_score")
write.csv(a4_dist,file.path(tables_dir,"A4_score_distribution_by_cohort_outcome.csv"),row.names=FALSE)
write.csv(a4_shift,file.path(tables_dir,"A4_score_shift_by_cohort.csv"),row.names=FALSE)
write.csv(a4_pooled,file.path(tables_dir,"A4_pooled_auc_after_per_cohort_standardization.csv"),row.names=FALSE)
p_a4 <- ggplot(d,aes(x=cohort_f,y=Mean_Score,fill=outcome)) +
  geom_violin(position=position_dodge(.8),width=.75,alpha=.55,colour=NA,scale="width") +
  geom_boxplot(position=position_dodge(.8),width=.18,outlier.size=.4,alpha=.9) +
  scale_fill_manual(values=c(`No transition`="#8DA0CB",Transition="#E78A3C"),name=NULL) +
  scale_x_discrete(labels=function(x)paste0(x,"\n(",cohort_status[x],")")) +
  labs(x=NULL,y="MRI risk score (SVM decision value)") + theme(legend.position="top")
save_fig(p_a4,"A4_score_distribution_by_cohort",7,4)
message("A4 done")

# ---------------------------------------------------------------------------
# Landscape robustness to population composition (Dmax, 20%/5%, 7/7)
# ---------------------------------------------------------------------------
populations <- c(list(Pooled=d,`External only`=external,`Newly evaluated`=filter(d,cohort%in%c("NAPLS-3","MUC-FRUE")),
                      `CHR-P only`=chr_sensitivity_population(d)),
                 setNames(lapply(cohort_levels,function(co)filter(d,cohort!=co)),paste("Omit",cohort_levels)))
a5 <- imap_dfr(populations,function(dat,pop)imap_dfr(mods,function(v,lab){
  z<-dat[eligible(dat,v),];o<-window_stats_1d(z,v,.20,.05,7L,7L,0L);good<-filter(o$rows,supported)
  tibble(population=pop,moderator=lab,N=nrow(z),events=sum(z$EXP_LABEL),
         cohorts_with_data=paste(intersect(cohort_levels,unique(z$cohort)),collapse="+"),
         eligible_windows=nrow(good),Dmax=o$Dmax,AUC_min=if(nrow(good))min(good$AUC)else NA_real_,AUC_max=if(nrow(good))max(good$AUC)else NA_real_)}))
write.csv(a5,file.path(tables_dir,"A5_landscape_population_robustness.csv"),row.names=FALSE)
pop_order <- names(populations)
a5_plot <- a5 %>% mutate(population=factor(population,levels=rev(pop_order)),moderator=factor(moderator,levels=names(mods)),
  type=case_when(population=="Pooled"~"Primary",grepl("^Omit",population)~"One cohort omitted",TRUE~"Restricted population"))
p_a5 <- ggplot(filter(a5_plot,is.finite(Dmax)),aes(x=Dmax,y=population,colour=type)) +
  geom_vline(data=filter(a5_plot,population=="Pooled"),aes(xintercept=Dmax),linetype=2,colour="grey50") +
  geom_point(aes(size=events)) + facet_wrap(~moderator,nrow=1) +
  scale_colour_manual(values=c(Primary="black",`Restricted population`="#1B9E77",`One cohort omitted`="#7570B3"),name=NULL) +
  scale_size_area(max_size=3.5,name="Transitions") + scale_x_continuous(limits=c(0,.45),breaks=c(0,.2,.4)) +
  labs(x="Range of local AUC across windows (Dmax)",y=NULL) + theme(legend.position="bottom")
save_fig(p_a5,"A5_landscape_population_robustness",10,3.8)
message("A5 done")

# ---------------------------------------------------------------------------
# Clinical-only vs clinical + MRI score, per held-out cohort (LOCO)
# Same specification as the manuscript (modules/loco.R).
# ---------------------------------------------------------------------------
loco_pred <- list()
for(include_functioning in c(FALSE,TRUE)) {
  variant <- if(include_functioning)"with_functioning"else"without_functioning"
  vars <- c("COGDIS_score","SIPS_Positiv_Gesamt","SIPS_Negativ_Gesamt","age",if(include_functioning)"Psychosoz_aequiv")
  z <- d[eligible(d,c(vars,"sex","cohort")),]
  for(co in sort(unique(z$cohort))) {
    tr<-filter(z,cohort!=co);te<-filter(z,cohort==co)
    sc<-scale_fit(tr,c(vars,"Mean_Score"),paste("A6",variant,co,sep="."));train<-sc$data;test<-te
    for(v in c(vars,"Mean_Score")){p<-filter(sc$params,variable==v);test[[paste0(v,"_z")]]<-(te[[v]]-p$mean)/p$SD}
    train$sex_f<-factor(train$sex,levels=c(1,2));test$sex_f<-factor(test$sex,levels=c(1,2))
    res<-tibble(variant=variant,heldout_cohort=co,y=te$EXP_LABEL,Mean_Score=te$Mean_Score)
    for(nm in c("clinical","clinical_plus_mri")) {
      fit<-glm(reformulate(c(paste0(vars,"_z"),"sex_f",if(nm=="clinical_plus_mri")"Mean_Score_z"),"EXP_LABEL"),data=train,family=binomial())
      res[[nm]]<-predict(fit,test,type="response")
    }
    loco_pred[[paste(variant,co)]]<-res
  }
}
loco_pred <- bind_rows(loco_pred)
loco_metric <- function(g) c(clinical_AUC=auc_fast(g$clinical,g$y),clinical_plus_mri_AUC=auc_fast(g$clinical_plus_mri,g$y),
                             AUC_difference=auc_fast(g$clinical_plus_mri,g$y)-auc_fast(g$clinical,g$y),
                             clinical_Brier=mean((g$clinical-g$y)^2),clinical_plus_mri_Brier=mean((g$clinical_plus_mri-g$y)^2))
a6 <- loco_pred %>% group_by(variant) %>% group_modify(~{
  groups <- c(split(.x,.x$heldout_cohort),list(`Pooled held-out`=.x))
  imap_dfr(groups,function(g,co){
    pt<-loco_metric(g)
    blocks<-split(seq_len(nrow(g)),paste(g$heldout_cohort,g$y))  # resample within cohort and outcome
    draws<-deterministic_rows(n_boot,paste("a6",.y$variant,co,sep="."),function(i){m<-loco_metric(g[unlist(lapply(blocks,sample_indices,replace=TRUE)),])
      c(m[["AUC_difference"]],m[["clinical_plus_mri_Brier"]]-m[["clinical_Brier"]])},n_cores)
    paired<-if(co!="Pooled held-out")pROC::roc.test(pROC::roc(g$y,g$clinical,direction="<",quiet=TRUE),pROC::roc(g$y,g$clinical_plus_mri,direction="<",quiet=TRUE),paired=TRUE,method="delong")$p.value else NA_real_
    as_tibble_row(c(pt,AUC_difference_lo=quantile(draws[,1],.025,na.rm=TRUE,names=FALSE),AUC_difference_hi=quantile(draws[,1],.975,na.rm=TRUE,names=FALSE),
                    Brier_difference=pt[["clinical_plus_mri_Brier"]]-pt[["clinical_Brier"]],
                    Brier_difference_lo=quantile(draws[,2],.025,na.rm=TRUE,names=FALSE),Brier_difference_hi=quantile(draws[,2],.975,na.rm=TRUE,names=FALSE),
                    DeLong_paired_p=paired,MRI_score_alone_AUC=auc_fast(g$Mean_Score,g$y))) %>%
      mutate(heldout_cohort=co,N=nrow(g),events=sum(g$y),.before=1)
  })
}) %>% ungroup()
write.csv(a6,file.path(tables_dir,"A6_loco_clinical_vs_clinical_plus_mri_by_cohort.csv"),row.names=FALSE)
a6_plot <- a6 %>% mutate(variant=recode(variant,without_functioning="Symptoms, age, sex",with_functioning="Plus functioning"),
                         label=paste0(heldout_cohort,"\n",N,"/",events))
p_a6 <- ggplot(a6_plot,aes(x=AUC_difference,y=fct_rev(factor(label,levels=unique(label))))) +
  geom_vline(xintercept=0,colour="grey50") + geom_errorbarh(aes(xmin=AUC_difference_lo,xmax=AUC_difference_hi),height=.2) +
  geom_point(aes(shape=heldout_cohort=="Pooled held-out"),size=2.4) + scale_shape_manual(values=c(16,18),guide="none") +
  facet_wrap(~variant,scales="free_y") +
  labs(x="Held-out AUC gain from adding the MRI score (95% CI)",y="Held-out cohort (N/transitions)")
save_fig(p_a6,"A6_loco_mri_gain_by_cohort",7,3.6)
message("A6 done")

# ---------------------------------------------------------------------------
# Decision curves (net benefit)
# ---------------------------------------------------------------------------
thresholds <- seq(.05,.50,by=.01)
net_benefit <- function(p,y,t) {n<-length(y);sapply(t,function(tt){pos<-p>=tt;sum(pos&y==1)/n-sum(pos&y==0)/n*tt/(1-tt)})}
treat_all <- function(y,t) mean(y)-(1-mean(y))*t/(1-t)
# Held-out probability update: slope/intercept fitted in the other external cohorts.
ext$p_heldout_update <- NA_real_
for(co in unique(ext$cohort)) {
  fit<-calibration_fit(ext$EXP_LABEL[ext$cohort!=co],ext$Probs_platt[ext$cohort!=co],"slope")
  ext$p_heldout_update[ext$cohort==co]<-plogis(coef(fit)[1]+coef(fit)[2]*safe_logit(ext$Probs_platt[ext$cohort==co]))
}
loco_sym <- filter(loco_pred,variant=="without_functioning")
dca_sets <- list(
  `External cohorts (n=685)`=list(y=ext$EXP_LABEL,strata=ext$cohort,models=list(`MRI model, original probabilities`=ext$Probs_platt,`MRI model, held-out recalibration`=ext$p_heldout_update)),
  `Symptom-complete cohorts, leave-one-cohort-out (n=658)`=list(y=loco_sym$y,strata=loco_sym$heldout_cohort,models=list(`Clinical only`=loco_sym$clinical,`Clinical + MRI score`=loco_sym$clinical_plus_mri)))
a3_curves <- imap_dfr(dca_sets,function(s,set)bind_rows(
  tibble(set=set,strategy="Treat all",threshold=thresholds,net_benefit=treat_all(s$y,thresholds)),
  tibble(set=set,strategy="Treat none",threshold=thresholds,net_benefit=0),
  imap_dfr(s$models,function(p,m)tibble(set=set,strategy=m,threshold=thresholds,net_benefit=net_benefit(p,s$y,thresholds)))))
key_t <- c(.10,.15,.20,.30)
a3_key <- imap_dfr(dca_sets,function(s,set)imap_dfr(s$models,function(p,m){
  stat<-function(ii){yy<-s$y[ii];pp<-p[ii];net_benefit(pp,yy,key_t)-pmax(sapply(key_t,function(t)treat_all(yy,t)),0)}
  pt<-stat(seq_along(s$y));blocks<-split(seq_along(s$y),paste(s$strata,s$y))
  draws<-deterministic_rows(n_boot,paste("a3",set,m,sep="."),function(i)stat(unlist(lapply(blocks,sample_indices,replace=TRUE))),n_cores)
  tibble(set=set,model=m,threshold=key_t,net_benefit=net_benefit(p,s$y,key_t),
         best_default=pmax(sapply(key_t,function(t)treat_all(s$y,t)),0),
         difference_vs_best_default=pt,lo=apply(draws,2,quantile,.025),hi=apply(draws,2,quantile,.975))
}))
write.csv(a3_curves,file.path(tables_dir,"A3_decision_curves.csv"),row.names=FALSE)
write.csv(a3_key,file.path(tables_dir,"A3_net_benefit_key_thresholds.csv"),row.names=FALSE)
p_a3 <- ggplot(a3_curves,aes(x=threshold,y=net_benefit,colour=strategy,linetype=strategy%in%c("Treat all","Treat none"))) +
  geom_line(linewidth=.7) + facet_wrap(~set) + coord_cartesian(ylim=c(-.03,.16)) +
  scale_linetype_manual(values=c(`FALSE`=1,`TRUE`=2),guide="none") +
  scale_colour_manual(values=c(`Treat all`="grey40",`Treat none`="grey70",`MRI model, original probabilities`="#D95F02",
    `MRI model, held-out recalibration`="#7570B3",`Clinical only`="#1B9E77",`Clinical + MRI score`="#E7298A"),name=NULL) +
  labs(x="Threshold probability",y="Net benefit") + theme(legend.position="bottom") + guides(colour=guide_legend(nrow=2))
save_fig(p_a3,"A3_decision_curves",9,4.2)
message("A3 done")

writeLines(c(paste("input_sha256:",digest::digest(file=input_path,algo="sha256")),paste("run:",format(Sys.time())),
             capture.output(sessionInfo())),file.path(out_dir,"run_info.txt"))
