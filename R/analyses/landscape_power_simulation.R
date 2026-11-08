#!/usr/bin/env Rscript
# Operating characteristics of the Dmax landscape test in the study's own design.
# Moderator values, cohort membership and outcomes are the observed ones; only the
# MRI score is simulated, with a known local AUC along the moderator:
#   score | no transition ~ N(0,1);  score | transition ~ N(sqrt(2)*qnorm(AUC(u)), 1),
# where u is the moderator's position in its observed range (the window scale).
# The global test is the manuscript's: Dmax over 20%/5% windows with 7/7 support,
# within-cohort moderator permutation. Permuting moderator values within cohort is
# equivalent to permuting (outcome, score) rows within cohort against fixed values,
# which lets window membership be computed once.
# A score x moderator logistic interaction test (with cohort effects) is recorded
# on the same data for comparison.
#
# Usage: Rscript R/analyses/landscape_power_simulation.R INPUT.csv OUTPUT_DIR PIPELINE_REPO [cores] [reps] [perms]

args <- commandArgs(trailingOnly=TRUE)
input_path <- normalizePath(args[1]); out_dir <- args[2]; repo <- normalizePath(args[3])
n_cores <- if(length(args)>=4) as.integer(args[4]) else 1L
n_rep <- if(length(args)>=5) as.integer(args[5]) else 1000L
n_perm <- if(length(args)>=6) as.integer(args[6]) else 499L
source(file.path(repo,"R/lib/common.R"))
options(pronia.seed=20260913L)
tables_dir <- file.path(out_dir,"tables"); figures_dir <- file.path(out_dir,"figures")
for(p in c(tables_dir,figures_dir)) dir.create(p,recursive=TRUE,showWarnings=FALSE)

raw <- read.csv(input_path,check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(Cases="character"))
cpath <- Sys.getenv("PRONIA_CORRECTION_JSON")  # authorized duplicate-resolution record, kept outside the repository
correction <- if(nzchar(cpath)) jsonlite::read_json(cpath,simplifyVector=TRUE) else NULL
d <- canonicalize_data(input_path,NULL,correction,NULL)

designs <- list(`SIPS-N`="SIPS_Negativ_Gesamt",Age="age")

run_design <- function(label,v) {
  x <- d[eligible(d,v),]; m <- x[[v]]; y <- x$EXP_LABEL; co <- x$cohort
  grid <- window_grid(m,.20,.05); u <- (m-min(m))/diff(range(m))
  members <- lapply(seq_len(nrow(grid)),function(i)which(in_window(m,grid,i)))
  blocks <- split(seq_along(y),co)
  # Observed trajectory, smoothed and taken as the truth: an optimistic scenario,
  # because the observed range also contains sampling noise.
  obs <- window_stats_1d(x,v,.20,.05,7L,7L,0L)$rows %>% filter(supported)
  obs_u <- (obs$center-min(m))/diff(range(m))
  obs_fit <- loess(AUC~obs_u,data=data.frame(AUC=obs$AUC,obs_u=obs_u),span=.75)
  observed_curve <- function(u) pmin(pmax(predict(obs_fit,data.frame(obs_u=pmin(pmax(u,min(obs_u)),max(obs_u)))),.40),.90)
  scenarios <- list(
    `Null, AUC 0.60 everywhere`=function(u)rep(.60,length(u)),
    `Linear, 0.50 to 0.70`=function(u).50+.20*u,
    `Linear, 0.45 to 0.80`=function(u).45+.35*u,
    `Step at midpoint, 0.55 to 0.75`=function(u)ifelse(u<.5,.55,.75),
    `Observed trajectory taken as truth`=observed_curve)
  # Same estimator as auc_fast without its per-call input validation (hot loop).
  auc_raw <- function(ss,yy,k1,k0){r<-rank(ss);(sum(r[yy==1L])-k1*(k1+1)/2)/(k1*k0)}
  dmax_rows <- function(yy,ss) {
    a <- numeric(0)
    for(ii in members){k1<-sum(yy[ii]);k0<-length(ii)-k1;if(k1>=7&&k0>=7)a<-c(a,auc_raw(ss[ii],yy[ii],k1,k0))}
    if(length(a)>=3) diff(range(a)) else NA_real_
  }
  imap_dfr(scenarios,function(f,scen){
    auc_u <- f(u); delta <- sqrt(2)*qnorm(auc_u)
    sims <- deterministic_rows(n_rep,paste("a2",label,scen,sep="."),function(r){
      s <- rnorm(length(y),mean=ifelse(y==1,delta,0))
      obs_d <- dmax_rows(y,s)
      null <- numeric(n_perm)
      for(b in seq_len(n_perm)){perm<-unlist(lapply(blocks,function(ii)ii[sample.int(length(ii))]),use.names=FALSE)
        idx<-integer(length(y));idx[unlist(blocks,use.names=FALSE)]<-perm;null[b]<-dmax_rows(y[idx],s[idx])}
      p_dmax <- (1+sum(null>=obs_d,na.rm=TRUE))/(1+sum(is.finite(null)))
      fit <- glm(y~scale(s)*scale(m)+factor(co),family=binomial())
      p_int <- summary(fit)$coefficients["scale(s):scale(m)","Pr(>|z|)"]
      c(Dmax=obs_d,p_dmax=p_dmax,p_interaction=p_int,pooled_AUC=auc_fast(s,y))
    },n_cores)
    tibble(design=label,N=length(y),events=sum(y),scenario=scen,
           true_AUC_min=min(auc_u),true_AUC_max=max(auc_u),true_AUC_range=max(auc_u)-min(auc_u),
           replicates=n_rep,permutations=n_perm,
           rejection_Dmax=mean(sims[,"p_dmax"]<=.05,na.rm=TRUE),
           rejection_Dmax_mcse=sqrt(mean(sims[,"p_dmax"]<=.05,na.rm=TRUE)*(1-mean(sims[,"p_dmax"]<=.05,na.rm=TRUE))/n_rep),
           rejection_interaction=mean(sims[,"p_interaction"]<=.05,na.rm=TRUE),
           mean_Dmax=mean(sims[,"Dmax"],na.rm=TRUE),Dmax_q95=quantile(sims[,"Dmax"],.95,na.rm=TRUE,names=FALSE),
           mean_pooled_AUC=mean(sims[,"pooled_AUC"]))
  })
}
res <- imap_dfr(designs,function(v,label){message("A2 design: ",label);run_design(label,v)})
write.csv(res,file.path(tables_dir,"A2_dmax_power_simulation.csv"),row.names=FALSE)
writeLines(c(paste("input_sha256:",digest::digest(file=input_path,algo="sha256")),paste("reps:",n_rep,"perms:",n_perm),
             paste("run:",format(Sys.time()))),file.path(out_dir,"A2_run_info.txt"))
message("A2 done")
