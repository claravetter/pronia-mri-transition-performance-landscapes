#!/usr/bin/env Rscript
# Scanner-site check within multi-site cohorts.
# The MRI pipeline is linear from grey-matter maps to score, so a per-site location
# offset in the features moves every score at that site by a constant. Such offsets
# can lower a multi-site cohort's AUC only through comparisons between sites. This
# script compares the cohort's AUC with
#   (a) the within-site AUC: only transition/non-transition pairs from the same site,
#   (b) the AUC after outcome-blind per-site centring or z-scoring of the score,
# and describes site offsets among non-transitioners.
#
# NAPLS-3 site: taken from SITE_MAP.csv (columns Cases, site) when supplied; otherwise
# derived from the leading field of the NAPLS-3 participant key (values 1-9), which is
# assumed to be the site code and must be confirmed. PRONIA uses its recorded site.
#
# Usage: Rscript R/analyses/scanner_site_check.R INPUT.csv OUTPUT_DIR PIPELINE_REPO [SITE_MAP.csv]

args <- commandArgs(trailingOnly=TRUE)
input_path <- normalizePath(args[1]); out_dir <- args[2]; repo <- normalizePath(args[3])
site_map_path <- if(length(args)>=4) normalizePath(args[4]) else NULL
source(file.path(repo,"R/lib/common.R"))
options(pronia.seed=20260913L)
n_boot <- 2000L
tables_dir <- file.path(out_dir,"tables"); figures_dir <- file.path(out_dir,"figures")
for(p in c(tables_dir,figures_dir)) dir.create(p,recursive=TRUE,showWarnings=FALSE)

raw <- read.csv(input_path,check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(Cases="character"))
cpath <- Sys.getenv("PRONIA_CORRECTION_JSON")  # authorized duplicate-resolution record, kept outside the repository
correction <- if(nzchar(cpath)) jsonlite::read_json(cpath,simplifyVector=TRUE) else NULL
d <- canonicalize_data(input_path,NULL,correction,NULL)
site_lookup <- setNames(as.character(raw$site),paste(sub("^NAPLES$","NAPLS-3",raw$cohort),raw$Cases))
d$site <- unname(site_lookup[paste(d$cohort,d$Cases)])

napls <- filter(d,cohort=="NAPLS-3")
if(is.null(site_map_path)) {
  napls$site <- sub("_.*$","",napls$Cases); site_source <- "derived from leading field of participant key (to be confirmed)"
} else {
  map <- read.csv(site_map_path,stringsAsFactors=FALSE,colClasses="character")
  if(!all(c("Cases","site")%in%names(map))||anyDuplicated(map$Cases)) stop("SITE_MAP.csv needs unique Cases and a site column")
  napls$site <- map$site[match(napls$Cases,map$Cases)]
  if(anyNA(napls$site)) stop(sum(is.na(napls$site))," NAPLS-3 participants have no site in SITE_MAP.csv")
  site_source <- paste("site map:",basename(site_map_path))
}
pronia <- filter(d,cohort=="PRONIA")
if(anyNA(pronia$site)||any(pronia$site%in%c("","NA"))) stop("PRONIA site missing")
cohorts <- list(`NAPLS-3`=list(x=napls,source=site_source),PRONIA=list(x=pronia,source="recorded site (leave-one-site-out OOT scores)"))

# AUC restricted to pairs within the same site (stratified Mann-Whitney).
within_site_auc <- function(score,y,site) {
  num<-0;den<-0
  for(s in unique(site)){ii<-site==s;k1<-sum(y[ii]==1);k0<-sum(y[ii]==0)
    if(k1&&k0){r<-rank(score[ii]);num<-num+sum(r[y[ii]==1])-k1*(k1+1)/2;den<-den+k1*k0}}
  num/den
}
boot_ci <- function(x,stat,id) {
  blocks<-split(seq_len(nrow(x)),paste(x$site,x$EXP_LABEL))
  v<-deterministic_rows(n_boot,id,function(i)stat(x[unlist(lapply(blocks,sample_indices,replace=TRUE)),]),1L)[,1]
  quantile(v,c(.025,.975),na.rm=TRUE,names=FALSE)
}
centre <- function(x) x %>% group_by(site) %>% mutate(score_centred=Mean_Score-median(Mean_Score),
  score_z=(Mean_Score-mean(Mean_Score))/sd(Mean_Score)) %>% ungroup()

summary_rows <- list(); site_rows <- list()
for(co in names(cohorts)) {
  x <- centre(cohorts[[co]]$x); src <- cohorts[[co]]$source
  estimates <- list(
    `Cohort AUC (all pairs)`=function(z)auc_fast(z$Mean_Score,z$EXP_LABEL),
    `Within-site AUC (same-site pairs only)`=function(z)within_site_auc(z$Mean_Score,z$EXP_LABEL,z$site),
    `Per-site median-centred score`=function(z)auc_fast(z$score_centred,z$EXP_LABEL),
    `Per-site z-scored score`=function(z)auc_fast(z$score_z,z$EXP_LABEL))
  summary_rows[[co]] <- imap_dfr(estimates,function(f,lab){ci<-boot_ci(x,f,paste("site",co,lab,sep="."))
    tibble(cohort=co,site_source=src,sites=n_distinct(x$site),N=nrow(x),events=sum(x$EXP_LABEL),estimate=lab,AUC=f(x),lo=ci[1],hi=ci[2])})
  nt_all <- x$Mean_Score[x$EXP_LABEL==0]
  per_site <- x %>% group_by(site) %>% summarise(N=n(),events=sum(EXP_LABEL),transition_rate=mean(EXP_LABEL),
      non_transition_mean=mean(Mean_Score[EXP_LABEL==0]),non_transition_SD=sd(Mean_Score[EXP_LABEL==0]),
      AUC=if(sum(EXP_LABEL)>=3&&sum(1-EXP_LABEL)>=3)auc_fast(Mean_Score,EXP_LABEL)else NA_real_,.groups="drop") %>%
    mutate(cohort=co,non_transition_offset_SMD=(non_transition_mean-mean(nt_all))/sd(nt_all),
           non_transition_SD_ratio=non_transition_SD/sd(nt_all),.before=1)
  kw <- kruskal.test(Mean_Score~factor(site),data=filter(x,EXP_LABEL==0))
  offset_rate <- suppressWarnings(cor.test(per_site$non_transition_mean,per_site$transition_rate,method="spearman"))
  per_site$site_offset_KW_p <- kw$p.value
  per_site$offset_vs_rate_spearman <- unname(offset_rate$estimate); per_site$offset_vs_rate_p <- offset_rate$p.value
  site_rows[[co]] <- per_site
}
site_summary <- bind_rows(summary_rows); site_detail <- bind_rows(site_rows)
write.csv(site_summary,file.path(tables_dir,"B1_site_check_auc.csv"),row.names=FALSE)
write.csv(site_detail,file.path(tables_dir,"B1_site_check_by_site.csv"),row.names=FALSE)

theme_set(theme_minimal(base_size=10)+theme(panel.grid.minor=element_blank()))
plot_d <- bind_rows(lapply(names(cohorts),function(co)mutate(cohorts[[co]]$x,cohort=co))) %>%
  mutate(outcome=factor(ifelse(EXP_LABEL==1,"Transition","No transition"),levels=c("No transition","Transition")),
         site=factor(site,levels=unique(site[order(suppressWarnings(as.numeric(site)),site)])))
p <- ggplot(plot_d,aes(x=site,y=Mean_Score,fill=outcome)) + geom_boxplot(outlier.size=.4,width=.6,position=position_dodge(.7)) +
  facet_wrap(~cohort,scales="free_x") + scale_fill_manual(values=c(`No transition`="#8DA0CB",Transition="#E78A3C"),name=NULL) +
  labs(x="Site",y="MRI risk score") + theme(legend.position="top")
ggsave(file.path(figures_dir,"B1_site_score_distribution.pdf"),p,width=9,height=3.8)
ggsave(file.path(figures_dir,"B1_site_score_distribution.png"),p,width=9,height=3.8,dpi=300)
writeLines(c(paste("input_sha256:",digest::digest(file=input_path,algo="sha256")),paste("NAPLS-3 site source:",site_source),
             paste("run:",format(Sys.time()))),file.path(out_dir,"B1_run_info.txt"))
print(as.data.frame(site_summary),digits=3)
