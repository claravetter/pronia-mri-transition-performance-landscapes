#!/usr/bin/env Rscript
# PRONIA participants with CHR-P: clinical profile of participants with and without transition,
# premorbid adjustment, and timing of documented transitions.
#   Profile    age, sex, SIPS positive and negative scores, COGDIS and MRI score by outcome in PRONIA
#              CHR-P participants with SIPS negative scores; standardized differences (pooled SD);
#              SIPS negative scores of participants with transition, NAPLS-3 vs PRONIA (Welch test).
#   Premorbid  mean of the Premorbid Adjustment Scale childhood and early-adolescence items;
#              outcome comparison (Wilcoxon), MRI-score AUC at or below vs above the median, and the
#              MRI score x premorbid adjustment logistic interaction (both standardized).
#   Timing     onset of psychosis documented at a study visit (earliest SIPS POPS date after
#              baseline); AUC counting only documented transitions, after excluding the others,
#              and for the later or undocumented transitions against nontransitions.
# Inputs: the authorized participant table (with the PRONIA participant number in `PatID`), the
# PRONIA baseline PAS export (columns PSN, PAS_A*_T0, PAS_B*_T0) and the PRONIA visit export with
# examination and POPS dates. Writes aggregate tables only.
#
# Usage: Rscript R/analyses/pronia_profile_premorbid_onset_timing.R INPUT.csv PAS.csv VISITS.csv OUTPUT_DIR PIPELINE_REPO
# Set PRONIA_CORRECTION_JSON to the duplicate-resolution record if the input contains it.

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<5) stop("Usage: pronia_profile_premorbid_onset_timing.R INPUT.csv PAS.csv VISITS.csv OUTPUT_DIR PIPELINE_REPO")
input_path <- normalizePath(args[1]); pas_path <- normalizePath(args[2]); visit_path <- normalizePath(args[3])
out_dir <- args[4]; repo <- normalizePath(args[5])
source(file.path(repo,"R/lib/common.R"))
suppressMessages({library(data.table); library(pROC)})
tables_dir <- file.path(out_dir,"tables"); dir.create(tables_dir,recursive=TRUE,showWarnings=FALSE)

cpath <- Sys.getenv("PRONIA_CORRECTION_JSON")  # authorized duplicate-resolution record, kept outside the repository
correction <- if(nzchar(cpath)) jsonlite::read_json(cpath,simplifyVector=TRUE) else NULL
d <- as.data.table(canonicalize_data(input_path,NULL,correction,NULL))
if(!"PatID" %in% names(d)) stop("The input needs the PRONIA participant number in column PatID")
pc <- d[cohort=="PRONIA" & Studygroup=="CHR"]

auc_ci <- function(y,s) {
  r <- roc(y,s,quiet=TRUE,direction="<"); ci <- ci.auc(r,method="delong")
  data.table(N=length(y),transitions=sum(y),AUC=as.numeric(auc(r)),lo=ci[1],hi=ci[3])
}
smd <- function(x,v) {
  a <- na.omit(x[EXP_LABEL==1][[v]]); b <- na.omit(x[EXP_LABEL==0][[v]])
  (mean(a)-mean(b))/sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2))
}

# ---- Profile -------------------------------------------------------------------------------
x <- pc[!is.na(SIPS_Negativ_Gesamt)]
vars <- c(age="age",SIPS_P="SIPS_Positiv_Gesamt",SIPS_N="SIPS_Negativ_Gesamt",COGDIS="COGDIS_score",MRI_score="Mean_Score")
profile <- rbindlist(lapply(1:0,function(y){
  z <- x[EXP_LABEL==y]
  data.table(outcome=y,N=nrow(z),female=sum(z$sex==2),
             variable=names(vars),mean=sapply(vars,function(v)mean(z[[v]],na.rm=TRUE)),
             sd=sapply(vars,function(v)sd(z[[v]],na.rm=TRUE)))
}))
profile_smd <- data.table(variable=names(vars),standardized_difference=sapply(vars,function(v)smd(x,v)))
a <- d[cohort=="NAPLS-3" & EXP_LABEL==1 & !is.na(SIPS_Negativ_Gesamt)]$SIPS_Negativ_Gesamt
b <- x[EXP_LABEL==1]$SIPS_Negativ_Gesamt; w <- t.test(a,b)
sipsn_transitions <- data.table(comparison="SIPS-N in participants with transition, NAPLS-3 minus PRONIA",
                                difference=mean(a)-mean(b),lo=w$conf.int[1],hi=w$conf.int[2],p_welch=w$p.value)
fwrite(profile,file.path(tables_dir,"pronia_chr_profile_by_outcome.csv"))
fwrite(profile_smd,file.path(tables_dir,"pronia_chr_profile_standardized_differences.csv"))
fwrite(sipsn_transitions,file.path(tables_dir,"sipsn_transitions_napls3_vs_pronia.csv"))

# ---- Premorbid adjustment ----------------------------------------------------------------------
p <- fread(pas_path,colClasses=c(PSN="character"))
items <- grep("^PAS_[AB][0-9]_T0$",names(p),value=TRUE)
p[,pas:=rowMeans(.SD),.SDcols=items]
m <- merge(pc[,.(PatID=as.character(PatID),EXP_LABEL,Mean_Score)],p[,.(PatID=PSN,pas)],by="PatID")[!is.na(pas)]
med <- median(m$pas)
m[,`:=`(zs=scale(Mean_Score)[,1],zp=scale(pas)[,1])]
co <- summary(glm(EXP_LABEL~zs*zp,data=m,family=binomial))$coefficients["zs:zp",]
premorbid <- rbind(
  data.table(analysis="outcome comparison",N=nrow(m),transitions=sum(m$EXP_LABEL),
             mean_transition=mean(m[EXP_LABEL==1]$pas),sd_transition=sd(m[EXP_LABEL==1]$pas),
             mean_nontransition=mean(m[EXP_LABEL==0]$pas),sd_nontransition=sd(m[EXP_LABEL==0]$pas),
             p=wilcox.test(pas~EXP_LABEL,m)$p.value,median_pas=med),
  data.table(analysis="score x premorbid interaction",N=nrow(m),transitions=sum(m$EXP_LABEL),
             beta=co[1],se=co[2],lo=co[1]-1.96*co[2],hi=co[1]+1.96*co[2],p=co[4]),fill=TRUE)
premorbid_auc <- rbind(cbind(group="at or below median",auc_ci(m[pas<=med]$EXP_LABEL,m[pas<=med]$Mean_Score)),
                       cbind(group="above median",auc_ci(m[pas>med]$EXP_LABEL,m[pas>med]$Mean_Score)))
fwrite(premorbid,file.path(tables_dir,"pronia_chr_premorbid_adjustment.csv"))
fwrite(premorbid_auc,file.path(tables_dir,"pronia_chr_auc_by_premorbid_adjustment.csv"))

# ---- Timing of documented transitions ------------------------------------------------------------
cols <- c("PSN","ALL_QUEST_Examination_Date_T0","CASE_DATE_T0","SIPS_CHR_01_POPS_05_Date_T1","SIPS_CHR_01_POPS_05_Date_T2",
          "CHRC_SIPS_01_POPS_04_Date_T1","CHRC_SIPS_01_POPS_04_Date_T2")
v <- fread(visit_path,select=cols,colClasses="character")
pd <- function(z) as.Date(z,"%d-%m-%Y")
t <- merge(pc[,.(PSN=as.character(PatID),T=EXP_LABEL,s=Mean_Score)],v,by="PSN")
t[,t0:=pd(ALL_QUEST_Examination_Date_T0)]; t[is.na(t0),t0:=as.Date(substr(CASE_DATE_T0,1,10))]
t[,mon:=as.numeric(pmin(pd(SIPS_CHR_01_POPS_05_Date_T1),pd(SIPS_CHR_01_POPS_05_Date_T2),
                        pd(CHRC_SIPS_01_POPS_04_Date_T1),pd(CHRC_SIPS_01_POPS_04_Date_T2),na.rm=TRUE)-t0)/30.44]
t[,Tdoc:=as.integer(T==1 & !is.na(mon))]
u <- t[!(T==1 & is.na(mon))]; w2 <- t[Tdoc==0]
timing <- rbind(cbind(analysis="all transitions",auc_ci(t$T,t$s)),
                cbind(analysis="documented onset only; others nontransition",auc_ci(t$Tdoc,t$s)),
                cbind(analysis="documented onset only; others excluded",auc_ci(u$T,u$s)),
                cbind(analysis="later or undocumented transitions vs nontransitions",auc_ci(w2$T,w2$s)))
onset <- data.table(documented=sum(t$Tdoc),median_months=median(t[Tdoc==1]$mon),
                    q25=quantile(t[Tdoc==1]$mon,.25),q75=quantile(t[Tdoc==1]$mon,.75),max_months=max(t[Tdoc==1]$mon))
fwrite(timing,file.path(tables_dir,"pronia_chr_auc_by_onset_documentation.csv"))
fwrite(onset,file.path(tables_dir,"pronia_chr_documented_onset_months.csv"))

print(profile_smd); print(sipsn_transitions); print(premorbid); print(premorbid_auc); print(timing); print(onset)
