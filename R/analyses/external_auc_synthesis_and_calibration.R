#!/usr/bin/env Rscript
# External discrimination and calibration by cohort:
#   Random-effects synthesis of external AUCs with Hartung-Knapp-Sidik-Jonkman CI and
#       95% prediction interval; never-used vs previously used comparison.
#   Calibration by external cohort after cohort-specific intercept updating
#       (slope with 95% CI, grouped observed vs predicted, loess curve data for eFigure 1).
#
# Usage: Rscript R/analyses/external_auc_synthesis_and_calibration.R INPUT.csv OUTPUT_DIR PIPELINE_REPO

args <- commandArgs(trailingOnly=TRUE)
input_path <- normalizePath(args[1]); out_dir <- args[2]; repo <- normalizePath(args[3])
source(file.path(repo,"R/lib/common.R"))
dir.create(file.path(out_dir,"tables"),recursive=TRUE,showWarnings=FALSE)

raw <- read.csv(input_path,check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(Cases="character"))
cpath <- Sys.getenv("PRONIA_CORRECTION_JSON")  # authorized duplicate-resolution record, kept outside the repository
correction <- if(nzchar(cpath)) jsonlite::read_json(cpath,simplifyVector=TRUE) else NULL
d <- canonicalize_data(input_path,NULL,correction,NULL)
external <- filter(d,cohort!="PRONIA")
used <- c(ZInEP="Previously used",FePsy="Previously used",`NAPLS-3`="Never used",`MUC-FRUE`="Never used")

# ---- Meta-analysis on the logit-AUC scale -----------------------------------------
meta <- external %>% group_by(cohort) %>% group_modify(~{
  r <- pROC::roc(.x$EXP_LABEL,.x$Mean_Score,direction="<",quiet=TRUE)
  a <- as.numeric(pROC::auc(r)); se <- sqrt(as.numeric(pROC::var(r,method="delong")))
  tibble(N=nrow(.x),events=sum(.x$EXP_LABEL),AUC=a,yi=qlogis(a),vi=(se/(a*(1-a)))^2)
}) %>% ungroup() %>% mutate(group=used[cohort])
fit_reml <- metafor::rma(yi=yi,vi=vi,data=meta,method="REML")
fit_hk <- metafor::rma(yi=yi,vi=vi,data=meta,method="REML",test="knha")
pred <- predict(fit_hk)
sub <- metafor::rma(yi=yi,vi=vi,mods=~group,data=meta,method="REML",test="knha")
# inverse-variance (fixed) pooling within each group and a z test of the difference
grp <- meta %>% group_by(group) %>% summarise(est=sum(yi/vi)/sum(1/vi),se=sqrt(1/sum(1/vi)),.groups="drop")
diff_logit <- diff(grp$est[match(c("Never used","Previously used"),grp$group)])
diff_se <- sqrt(sum(grp$se^2))
d1 <- bind_rows(
  tibble(analysis="Random effects, REML (Wald CI)",estimate=plogis(fit_reml$b[1]),lo=plogis(fit_reml$ci.lb),hi=plogis(fit_reml$ci.ub),
         I2=fit_reml$I2,tau2=fit_reml$tau2,p=NA_real_),
  tibble(analysis="Random effects, REML (Hartung-Knapp CI)",estimate=plogis(pred$pred),lo=plogis(pred$ci.lb),hi=plogis(pred$ci.ub),
         I2=fit_hk$I2,tau2=fit_hk$tau2,p=NA_real_),
  tibble(analysis="95% prediction interval (Hartung-Knapp)",estimate=plogis(pred$pred),lo=plogis(pred$pi.lb),hi=plogis(pred$pi.ub),
         I2=NA_real_,tau2=NA_real_,p=NA_real_),
  tibble(analysis=paste("Inverse-variance pooled,",grp$group),estimate=plogis(grp$est),lo=plogis(grp$est-1.96*grp$se),
         hi=plogis(grp$est+1.96*grp$se),I2=NA_real_,tau2=NA_real_,p=NA_real_),
  tibble(analysis="Previously used vs never used, z test (logit scale)",estimate=diff_logit,lo=diff_logit-1.96*diff_se,
         hi=diff_logit+1.96*diff_se,I2=NA_real_,tau2=NA_real_,p=2*pnorm(-abs(diff_logit/diff_se))),
  tibble(analysis="Previously used vs never used, meta-regression (Hartung-Knapp)",estimate=sub$b[2],lo=sub$ci.lb[2],hi=sub$ci.ub[2],
         I2=NA_real_,tau2=sub$tau2,p=sub$pval[2]))
write.csv(meta,file.path(out_dir,"tables","D1_meta_cohort_inputs.csv"),row.names=FALSE)
write.csv(d1,file.path(out_dir,"tables","D1_meta_analysis.csv"),row.names=FALSE)
print(as.data.frame(d1),digits=3)

# ---- Calibration by cohort after intercept updating --------------------------------
clip <- function(p) pmin(pmax(p,1e-6),1-1e-6)
d2 <- list(); pts <- list(); curves <- list()
for(co in names(used)) {
  x <- filter(external,cohort==co); lp <- qlogis(clip(x$Probs_platt)); y <- x$EXP_LABEL
  a <- coef(glm(y~offset(lp),family=binomial()))[1]          # calibration-in-the-large
  sl <- glm(y~lp,family=binomial()); ci <- suppressMessages(confint.default(sl))[2,]
  p_upd <- plogis(lp+a)                                       # intercept-updated probabilities
  k <- if(nrow(x)>=100) 5 else 3
  g <- tibble(p=p_upd,y=y) %>% mutate(bin=ntile(p,k)) %>% group_by(bin) %>%
    summarise(predicted=mean(p),observed=mean(y),n=n(),events=sum(y),.groups="drop") %>%
    mutate(lo=mapply(function(e,n)binom.test(e,n)$conf.int[1],events,n),hi=mapply(function(e,n)binom.test(e,n)$conf.int[2],events,n))
  grid <- seq(min(p_upd),max(p_upd),length.out=60)
  lo_fit <- suppressWarnings(loess(y~p,data=tibble(p=p_upd,y=y),span=1,degree=1))
  d2[[co]] <- tibble(cohort=co,group=used[[co]],N=nrow(x),events=sum(y),observed_proportion=mean(y),
                     mean_predicted_original=mean(x$Probs_platt),intercept_update=a,
                     slope=coef(sl)[2],slope_lo=ci[1],slope_hi=ci[2],mean_predicted_updated=mean(p_upd))
  pts[[co]] <- mutate(g,cohort=co)
  curves[[co]] <- tibble(cohort=co,p=grid,smooth=pmin(pmax(predict(lo_fit,tibble(p=grid)),0),1))
}
write.csv(bind_rows(d2),file.path(out_dir,"tables","D2_calibration_by_cohort.csv"),row.names=FALSE)
write.csv(bind_rows(pts),file.path(out_dir,"tables","D2_calibration_points.csv"),row.names=FALSE)
write.csv(bind_rows(curves),file.path(out_dir,"tables","D2_calibration_curves.csv"),row.names=FALSE)
print(as.data.frame(bind_rows(d2)),digits=3)
