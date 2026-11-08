#!/usr/bin/env Rscript
# Can measured case-mix differences explain the external loss of discrimination?
# Following the case-mix step of Debray et al. (J Clin Epidemiol 2015;68:279-89):
#  1. Relatedness: a membership model (PRONIA CHR-P vs external cohort) on the measured
#     characteristics; its AUC shows how different the case mix is.
#  2. Case-mix-standardized reference: PRONIA CHR-P participants are reweighted by the
#     odds of external membership, so that their characteristics *within each outcome
#     group* match the external cohort. The weighted AUC estimates the discovery
#     performance expected under the external cohort's measured case mix.
#     If it stays near 0.73 while the external AUC is ~0.53, measured case mix does not
#     explain the loss.
# Covariate sets: age + sex (all cohorts); plus SIPS-P, SIPS-N and COGDIS (NAPLS-3,
# ZInEP). Weights are trimmed at the 1st/99th percentile.
#
# Usage: Rscript R/analyses/casemix_standardization.R INPUT.csv OUTPUT_DIR PIPELINE_REPO

args <- commandArgs(trailingOnly=TRUE)
input_path <- normalizePath(args[1]); out_dir <- args[2]; repo <- normalizePath(args[3])
source(file.path(repo,"R/lib/common.R"))
options(pronia.seed=20260913L)
n_boot <- 2000L
dir.create(file.path(out_dir,"tables"),recursive=TRUE,showWarnings=FALSE)

raw <- read.csv(input_path,check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(Cases="character"))
cpath <- Sys.getenv("PRONIA_CORRECTION_JSON")  # authorized duplicate-resolution record, kept outside the repository
correction <- if(nzchar(cpath)) jsonlite::read_json(cpath,simplifyVector=TRUE) else NULL
d <- canonicalize_data(input_path,NULL,correction,NULL)
ref <- filter(d,cohort=="PRONIA",Studygroup=="CHR")

weighted_auc <- function(s,y,w) {
  o<-order(s);s<-s[o];y<-y[o];w<-w[o]
  w1<-ifelse(y==1,w,0);w0<-ifelse(y==0,w,0)
  # for each transition, weight of non-transitions with lower score (+half ties)
  below<-cumsum(w0)-w0; tie<-ave(w0,s,FUN=sum); below_strict<-ave(below,s,FUN=min)
  sum(w1*(below_strict+.5*tie))/(sum(w1)*sum(w0))
}
reweight <- function(r,e,vars) {
  z <- bind_rows(mutate(r,member=0),mutate(e,member=1))
  f <- reformulate(c(paste0("(",paste(vars,collapse="+"),")*EXP_LABEL")),"member")
  fit <- suppressWarnings(glm(f,data=z,family=binomial()))
  p <- predict(fit,r,type="response"); w <- p/(1-p)
  w <- pmin(pmax(w,quantile(w,.01)),quantile(w,.99))
  # rescale within outcome group so each group keeps its PRONIA size
  w <- ave(w,r$EXP_LABEL,FUN=function(v)v/mean(v))
  ess <- function(v) sum(v)^2/sum(v^2)
  list(w=w,membership_AUC=auc_fast(predict(fit,z),z$member),ess=ess(w),ess_events=ess(w[r$EXP_LABEL==1]))
}
sets <- list(`Age, sex`=c("age","sex"),
             `Age, sex, SIPS-P, SIPS-N, COGDIS`=c("age","sex","SIPS_Positiv_Gesamt","SIPS_Negativ_Gesamt","COGDIS_score"))
rows <- list()
for(co in c("NAPLS-3","ZInEP","FePsy","MUC-FRUE")) for(sn in names(sets)) {
  vars <- sets[[sn]]
  e <- filter(d,cohort==co); e <- e[eligible(e,vars),]; r <- ref[eligible(ref,vars),]
  if(nrow(e)<20||sum(e$EXP_LABEL)<5) next
  rw <- reweight(r,e,vars)
  stat <- function(ii,jj){w<-reweight(r[ii,],e[jj,],vars)$w;weighted_auc(r$Mean_Score[ii],r$EXP_LABEL[ii],w)}
  rb <- split(seq_len(nrow(r)),r$EXP_LABEL); eb <- split(seq_len(nrow(e)),e$EXP_LABEL)
  draws <- deterministic_rows(n_boot,paste("casemix",co,sn,sep="."),function(i)
    stat(unlist(lapply(rb,sample_indices,replace=TRUE)),unlist(lapply(eb,sample_indices,replace=TRUE))),8L)[,1]
  rows[[paste(co,sn)]] <- tibble(external_cohort=co,covariates=sn,external_N=nrow(e),external_events=sum(e$EXP_LABEL),
    reference_N=nrow(r),reference_events=sum(r$EXP_LABEL),membership_AUC=rw$membership_AUC,
    effective_reference_N=rw$ess,effective_reference_events=rw$ess_events,reference_AUC_unweighted=auc_fast(r$Mean_Score,r$EXP_LABEL),
    reference_AUC_casemix_standardized=weighted_auc(r$Mean_Score,r$EXP_LABEL,rw$w),
    lo=quantile(draws,.025,na.rm=TRUE,names=FALSE),hi=quantile(draws,.975,na.rm=TRUE,names=FALSE),
    observed_external_AUC=auc_fast(e$Mean_Score,e$EXP_LABEL))
}
res <- bind_rows(rows)
write.csv(res,file.path(out_dir,"tables","C1_casemix_standardized_reference.csv"),row.names=FALSE)
print(as.data.frame(res %>% select(external_cohort,covariates,external_N,external_events,membership_AUC,effective_reference_N,
  reference_AUC_casemix_standardized,lo,hi,observed_external_AUC)),digits=3)
