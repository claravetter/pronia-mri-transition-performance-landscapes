#!/usr/bin/env Rscript
# Why does the MRI score not separate transitioners in NAPLS-3?
# Runs whichever checks the supplied NAPLS-3 covariate file supports:
#   Horizon     months_to_transition + months_followup -> AUC at 12/18/24-month horizons
#               (cases: transition by h; controls: no transition by h and followed >= h)
#               and AUC of early vs late transitioners against all non-transitioners.
#   Categorical every column listed in CATEGORICAL that is present (e.g. subjecttype,
#               chr_subtype, antipsychotic_baseline, antipsychotic_followup, site):
#               AUC with DeLong CI per level.
#   Continuous  every column listed in CONTINUOUS that is present (e.g. cat12_qc,
#               gfs_social, gfr_role, parental_education, cannabis, trauma):
#               AUC by within-NAPLS tertile and the score x variable logistic interaction.
# The covariate file is keyed by `username` (e.g. 01_S0001) or by `Cases`.
#
# Usage: Rscript R/analyses/napls3_subgroups_and_image_quality.R INPUT.csv NAPLS_COVARIATES.csv OUTPUT_DIR PIPELINE_REPO

args <- commandArgs(trailingOnly=TRUE)
input_path <- normalizePath(args[1]); cov_path <- normalizePath(args[2]); out_dir <- args[3]; repo <- normalizePath(args[4])
source(file.path(repo,"R/lib/common.R"))
dir.create(file.path(out_dir,"tables"),recursive=TRUE,showWarnings=FALSE)
CATEGORICAL <- c("subjecttype","chr_subtype","antipsychotic_baseline","antipsychotic_followup","site")
CONTINUOUS <- c("cat12_qc","global_gmv","gfs_social","gfr_role","parental_education","cannabis","trauma",
                "aps_duration_months","sips_g","calgary_depression")

raw <- read.csv(input_path,check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(Cases="character"))
x <- raw %>% filter(cohort%in%c("NAPLES","NAPLS-3")) %>%
  transmute(Cases,EXP_LABEL=as.integer(EXP_LABEL),Mean_Score,
            username=paste0(formatC(as.integer(sub("_.*$","",Cases)),width=2,flag="0"),"_",sapply(strsplit(Cases,"_"),`[`,2)))
cov <- read.csv(cov_path,check.names=FALSE,stringsAsFactors=FALSE,colClasses="character")
names(cov) <- tolower(names(cov)); names(cov)[names(cov)=="cases"] <- "Cases"
names(cov)[names(cov)=="cat12qc"] <- "cat12_qc"
key <- if("username"%in%names(cov)) "username" else "Cases"
x <- left_join(x,cov,by=key)
matched <- sum(!is.na(x[[setdiff(intersect(c(CATEGORICAL,CONTINUOUS,"months_followup"),names(cov)),key)[1]]]))
message("NAPLS-3 participants: ",nrow(x),"; with covariate record: ",matched)

delong <- function(s,y) if(sum(y==1)>=3&&sum(y==0)>=3) as.numeric(pROC::ci.auc(pROC::roc(y,s,direction="<",quiet=TRUE),method="delong")) else c(NA,auc_fast(s,y),NA)
row <- function(check,group,s,y) {ci<-delong(s,y);tibble(check=check,group=group,N=length(y),events=sum(y),AUC=ci[2],lo=ci[1],hi=ci[3])}
out <- list()

if(all(c("months_to_transition","months_followup")%in%names(x))) {
  t <- as.numeric(x$months_to_transition); fu <- as.numeric(x$months_followup)
  for(h in c(12,18,24)) {
    case <- x$EXP_LABEL==1 & t<=h; ctrl <- (x$EXP_LABEL==0 & fu>=h) | (x$EXP_LABEL==1 & t>h)
    keep <- which((case|ctrl)%in%TRUE); out[[paste("h",h)]] <- row("Transition horizon",paste0("<= ",h," months"),x$Mean_Score[keep],as.integer(case[keep]))
  }
  med <- median(t[x$EXP_LABEL==1],na.rm=TRUE)
  for(g in c("early","late")) {
    sel <- x$EXP_LABEL==0 | (x$EXP_LABEL==1 & if(g=="early") t<=med else t>med)
    keep <- which(sel%in%TRUE); out[[g]] <- row("Time to transition",paste0(g," transitions (median split at ",round(med,1)," months) vs all non-transitions"),x$Mean_Score[keep],x$EXP_LABEL[keep])
  }
}
for(v in intersect(CATEGORICAL,names(x))) for(lv in sort(unique(na.omit(x[[v]])))) {
  keep <- which(x[[v]]==lv); out[[paste(v,lv)]] <- row(v,lv,x$Mean_Score[keep],x$EXP_LABEL[keep])
}
inter <- list()
for(v in intersect(CONTINUOUS,names(x))) {
  z <- suppressWarnings(as.numeric(x[[v]])); ok <- is.finite(z); if(sum(ok)<30) next
  tert <- cut(z,quantile(z[ok],c(0,1/3,2/3,1)),include.lowest=TRUE,labels=c("low","middle","high"))
  for(lv in levels(tert)){keep<-which(tert==lv);out[[paste(v,lv)]]<-row(v,paste(lv,"tertile"),x$Mean_Score[keep],x$EXP_LABEL[keep])}
  fit <- glm(EXP_LABEL~scale(Mean_Score)*scale(z),family=binomial(),data=data.frame(EXP_LABEL=x$EXP_LABEL[ok],Mean_Score=x$Mean_Score[ok],z=z[ok]))
  cf <- summary(fit)$coefficients["scale(Mean_Score):scale(z)",]
  inter[[v]] <- tibble(variable=v,N=sum(ok),events=sum(x$EXP_LABEL[ok]),interaction_beta=cf[1],SE=cf[2],p=cf[4])
}
res <- bind_rows(out); write.csv(res,file.path(out_dir,"tables","C2_napls_followup_checks.csv"),row.names=FALSE)
if(length(inter)) write.csv(bind_rows(inter),file.path(out_dir,"tables","C2_napls_followup_interactions.csv"),row.names=FALSE)
print(as.data.frame(res),digits=3); if(length(inter)) print(as.data.frame(bind_rows(inter)),digits=3)
