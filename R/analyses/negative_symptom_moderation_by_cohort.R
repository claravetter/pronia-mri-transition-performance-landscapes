#!/usr/bin/env Rscript
# Direction check for the SIPS-N pattern: within each cohort, MRI-score AUC below vs
# at/above the cohort median SIPS-N, and the score x SIPS-N interaction per cohort.
# Usage: Rscript R/analyses/negative_symptom_moderation_by_cohort.R INPUT.csv OUTPUT_DIR PIPELINE_REPO
args <- commandArgs(trailingOnly=TRUE)
source(file.path(normalizePath(args[3]),"R/lib/common.R"))
raw <- read.csv(args[1],check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(Cases="character"))
cpath <- Sys.getenv("PRONIA_CORRECTION_JSON")  # authorized duplicate-resolution record, kept outside the repository
correction <- if(nzchar(cpath)) jsonlite::read_json(cpath,simplifyVector=TRUE) else NULL
d <- canonicalize_data(normalizePath(args[1]),NULL,correction,NULL)
x <- d[eligible(d,"SIPS_Negativ_Gesamt"),]
res <- x %>% group_by(cohort) %>% group_modify(~{
  hi <- .x$SIPS_Negativ_Gesamt >= median(.x$SIPS_Negativ_Gesamt)
  fit <- glm(EXP_LABEL~scale(Mean_Score)*scale(SIPS_Negativ_Gesamt),family=binomial(),data=.x)
  co <- summary(fit)$coefficients["scale(Mean_Score):scale(SIPS_Negativ_Gesamt)",]
  tibble(N=nrow(.x),events=sum(.x$EXP_LABEL),median_SIPS_N=median(.x$SIPS_Negativ_Gesamt),
         low_N=sum(!hi),low_events=sum(.x$EXP_LABEL[!hi]),AUC_low=auc_fast(.x$Mean_Score[!hi],.x$EXP_LABEL[!hi]),
         high_N=sum(hi),high_events=sum(.x$EXP_LABEL[hi]),AUC_high=auc_fast(.x$Mean_Score[hi],.x$EXP_LABEL[hi]),
         interaction_beta=co[1],interaction_SE=co[2],interaction_p=co[4])
}) %>% ungroup()
dir.create(file.path(args[2],"tables"),showWarnings=FALSE,recursive=TRUE)
write.csv(res,file.path(args[2],"tables","A5b_sipsn_direction_by_cohort.csv"),row.names=FALSE)
print(as.data.frame(res),digits=3)
