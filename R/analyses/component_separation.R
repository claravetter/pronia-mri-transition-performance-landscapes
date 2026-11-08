#!/usr/bin/env Rscript
# Component-level separation by cohort.
# Question: is the loss of separation between participants with and without transition outside
# PRONIA shared by the 4 structural components of the MRI score, or concentrated in some of them?
#   - per cohort and variable (MRI score, Comp 1-4): separation as standardized mean difference
#     (transition minus nontransition, pooled SD; as in Figure 2A) and AUC (higher value = transition)
#   - loss vs PRONIA: SMD in PRONIA minus SMD in the cohort
#   - relative loss: loss of component k minus the mean loss of the 4 components; tested for NAPLS-3
#     and for the external summary (participant-weighted mean of the 4 external cohort SMDs);
#     BH across these 8 contrasts
# Bootstrap: 2000 replicates resampling within cohort and outcome, deterministic streams.
# Writes aggregate tables only; no participant-level output.
#
# Usage: Rscript R/analyses/component_separation.R INPUT.csv OUTPUT_DIR PIPELINE_REPO [cores]

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3) stop("Usage: component_separation.R INPUT.csv OUTPUT_DIR PIPELINE_REPO [cores]")
input_path <- normalizePath(args[1]); out_dir <- args[2]; repo <- normalizePath(args[3])
n_cores <- if(length(args)>=4) as.integer(args[4]) else 1L
source(file.path(repo,"R/lib/common.R"))
options(pronia.seed=20260913L)
n_boot <- 2000L
tables_dir <- file.path(out_dir,"tables"); dir.create(tables_dir,recursive=TRUE,showWarnings=FALSE)

# Canonical population: same duplicate resolution as the manuscript (correction record via PRONIA_CORRECTION_JSON).
raw <- read.csv(input_path,check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(Cases="character"))
cpath <- Sys.getenv("PRONIA_CORRECTION_JSON")  # authorized duplicate-resolution record, kept outside the repository
correction <- if(nzchar(cpath)) jsonlite::read_json(cpath,simplifyVector=TRUE) else NULL
expected <- list(PRONIA=list(N=329,events=26),`NAPLS-3`=list(N=458,events=51),ZInEP=list(N=157,events=23),
                 FePsy=list(N=37,events=16),`MUC-FRUE`=list(N=33,events=16))
d <- canonicalize_data(input_path,NULL,correction,expected)

cohort_levels <- c("PRONIA","NAPLS-3","ZInEP","FePsy","MUC-FRUE")
external_levels <- cohort_levels[-1]
vars <- c(`MRI score`="Mean_Score",`Comp 1`="Comp1",`Comp 2`="Comp2",`Comp 3`="Comp3",`Comp 4`="Comp4")
comps <- names(vars)[-1]
stopifnot(all(vapply(d[vars],function(z)all(is.finite(z)),logical(1))))
n_by_cohort <- table(d$cohort)[external_levels]

smd <- function(x,y) (mean(x[y==1])-mean(x[y==0]))/sqrt((var(x[y==1])+var(x[y==0]))/2)

# statistic vector for one data set: SMD and AUC per cohort and variable
stat_names <- as.vector(outer(paste(rep(cohort_levels,each=length(vars)),names(vars),sep="|"),c("SMD","AUC"),paste,sep="|"))
stats <- function(x) {
  smds <- aucs <- numeric(0)
  for(co in cohort_levels) {
    g <- x[x$cohort==co,]
    for(v in vars) { smds <- c(smds,smd(g[[v]],g$EXP_LABEL)); aucs <- c(aucs,auc_fast(g[[v]],g$EXP_LABEL)) }
  }
  setNames(c(smds,aucs),stat_names)
}
# derived contrasts from one statistic vector
derive <- function(s) {
  get <- function(co,v,m="SMD") s[[paste(co,v,m,sep="|")]]
  out <- c()
  ext_smd <- function(v) sum(vapply(external_levels,function(co)get(co,v),numeric(1))*n_by_cohort)/sum(n_by_cohort)
  for(v in names(vars)) out[paste("External summary",v,"SMD",sep="|")] <- ext_smd(v)
  for(target in c(external_levels,"External summary")) {
    tg <- function(v) if(target=="External summary") ext_smd(v) else get(target,v)
    loss <- vapply(names(vars),function(v)get("PRONIA",v)-tg(v),numeric(1))
    for(v in names(vars)) out[paste(target,v,"loss",sep="|")] <- loss[[v]]
    for(v in comps) out[paste(target,v,"relative_loss",sep="|")] <- loss[[v]]-mean(loss[comps])
  }
  out
}

obs <- stats(d); obs_d <- derive(obs)
blocks <- split(seq_len(nrow(d)),interaction(d$cohort,d$EXP_LABEL,drop=TRUE))
draws <- deterministic_rows(n_boot,"component_separation.boot",function(i){
  x <- d[unlist(lapply(blocks,sample_indices,replace=TRUE)),]
  s <- stats(x); c(s,derive(s))
},cores=n_cores)

summ <- function(name,est) {
  b <- draws[,name]; b <- b[is.finite(b)]
  p <- min(1,2*min(mean(b<=0),mean(b>=0)))
  q <- quantile(b,c(.025,.975),names=FALSE)
  parts <- strsplit(name,"|",fixed=TRUE)[[1]]
  tibble(cohort=parts[1],variable=parts[2],measure=parts[3],estimate=est,lo=q[1],hi=q[2],p_boot=p,valid_replicates=length(b))
}
all_est <- c(obs,obs_d)
res <- map_dfr(names(all_est),function(nm)summ(nm,all_est[[nm]]))
res <- res %>% mutate(N=as.integer(ifelse(cohort %in% cohort_levels,table(d$cohort)[cohort],sum(n_by_cohort))),
                      transitions=as.integer(ifelse(cohort %in% cohort_levels,tapply(d$EXP_LABEL,d$cohort,sum)[cohort],
                                                    sum(d$EXP_LABEL[d$cohort!="PRONIA"]))))
family <- res$measure=="relative_loss" & res$cohort %in% c("NAPLS-3","External summary")
res$q_BH_8 <- NA_real_; res$q_BH_8[family] <- p.adjust(res$p_boot[family],"BH")
write.csv(res,file.path(tables_dir,"C1_component_separation_by_cohort.csv"),row.names=FALSE)

# correlations among the components and with the MRI score, by cohort (aggregate)
cors <- map_dfr(cohort_levels,function(co){
  g <- d[d$cohort==co,]; m <- cor(g[vars],method="spearman"); dimnames(m) <- list(names(vars),names(vars))
  as.data.frame(as.table(m)) %>% filter(as.integer(Var1)<as.integer(Var2)) %>%
    transmute(cohort=co,variable_1=as.character(Var1),variable_2=as.character(Var2),spearman_rho=Freq)
})
write.csv(cors,file.path(tables_dir,"C1_component_correlations_by_cohort.csv"),row.names=FALSE)

writeLines(c(paste("input_sha256:",digest::digest(file=input_path,algo="sha256")),paste("run:",format(Sys.time())),
             paste("bootstrap replicates:",n_boot),capture.output(sessionInfo())),file.path(out_dir,"C1_run_info.txt"))
message("C1 done")
