#!/usr/bin/env Rscript
# Main Figures 1 and 2 and supplementary figures (no grid lines).
# Figure 1: (A) AUC by cohort, (B) calibration by cohort, (C) decision curves, (D) held-out AUC
# gain from adding the MRI score to the symptom model.
# Figure 2: (A) MRI score densities by cohort and outcome with the separation in SD,
# (B) observed external AUC vs the AUC expected in PRONIA CHR-P participants reweighted to
# each cohort's age and sex (C1_casemix_standardized_reference.csv), NAPLS-3 by PRONIA age range.
# Reads the canonical participant table and the decision-curve and added-value outputs of
# reference_shift_added_value_net_benefit.R.
#
# Usage: Rscript R/analyses/main_figures.R INPUT.csv OUTPUT_DIR PIPELINE_REPO

args <- commandArgs(trailingOnly=TRUE)
input_path <- normalizePath(args[1]); out_dir <- args[2]; repo <- normalizePath(args[3])
source(file.path(repo,"R/lib/common.R"))
suppressPackageStartupMessages(library(patchwork))
fig_dir <- file.path(out_dir,"figures"); dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)

raw <- read.csv(input_path,check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(Cases="character"))
cpath <- Sys.getenv("PRONIA_CORRECTION_JSON")  # authorized duplicate-resolution record, kept outside the repository
correction <- if(nzchar(cpath)) jsonlite::read_json(cpath,simplifyVector=TRUE) else NULL
d <- canonicalize_data(input_path,NULL,correction,NULL)
external <- filter(d,cohort!="PRONIA")

# Validated categorical slots (dataviz reference palette, light): blue, orange; ink for reference rows.
COL <- c(`Previously used`="#2a78d6",`Never used`="#eb6834",Reference="#3d3c39")
SHAPE <- c(`Previously used`=16,`Never used`=17,Reference=15)
INK2 <- "#52514e"; GRID <- "#e4e3df"
status <- c(PRONIA="Reference",ZInEP="Previously used",FePsy="Previously used",`NAPLS-3`="Never used",`MUC-FRUE`="Never used")
# Figures 1 and 2 show all external cohorts in one colour and marker (no used/never-used distinction);
# the status labels are kept only for the source CSVs and the supplementary random-effects eFigure.
EXT_COL <- COL[["Previously used"]]
COL_EXT <- c(`Previously used`=EXT_COL,`Never used`=EXT_COL,Reference=COL[["Reference"]])
SHAPE_EXT <- c(`Previously used`=16,`Never used`=16,Reference=15)
theme_figure <- theme_classic(base_size=8,base_family="Helvetica") +
  theme(axis.line=element_line(colour=INK2,linewidth=.3),axis.ticks=element_line(colour=INK2,linewidth=.3),
        axis.text=element_text(colour="#0b0b0b"),axis.title=element_text(colour="#0b0b0b"),
        panel.grid=element_blank(),
        legend.position="bottom",legend.key.size=unit(9,"pt"),legend.title=element_blank(),
        legend.margin=margin(0,0,0,0),plot.tag=element_text(face="bold",size=10),
        plot.margin=margin(4,6,4,4))
delong <- function(s,y) as.numeric(pROC::ci.auc(pROC::roc(y,s,direction="<",quiet=TRUE),method="delong"))

# ---- Figure 1A: AUC by cohort ---------------------------------------------------------
row <- function(label,x,grp) {ci<-delong(x$Mean_Score,x$EXP_LABEL);tibble(label=label,N=nrow(x),events=sum(x$EXP_LABEL),AUC=ci[2],lo=ci[1],hi=ci[3],group=grp)}
meta <- external %>% group_by(cohort) %>% group_modify(~{r<-pROC::roc(.x$EXP_LABEL,.x$Mean_Score,direction="<",quiet=TRUE);a<-as.numeric(pROC::auc(r));se<-sqrt(as.numeric(pROC::var(r,method="delong")));tibble(yi=qlogis(a),vi=(se/(a*(1-a)))^2)}) %>% ungroup()
re <- metafor::rma(yi=yi,vi=vi,data=meta,method="REML",test="knha")  # Hartung-Knapp CI, as in Table 2
forest <- bind_rows(
  row("PRONIA, all",filter(d,cohort=="PRONIA"),"Reference"),
  row("PRONIA, CHR-P only",filter(d,cohort=="PRONIA",Studygroup=="CHR"),"Reference"),
  row("ZInEP",filter(d,cohort=="ZInEP"),"Previously used"),
  row("FePsy",filter(d,cohort=="FePsy"),"Previously used"),
  row("NAPLS-3",filter(d,cohort=="NAPLS-3"),"Never used"),
  row("MUC-FRUE",filter(d,cohort=="MUC-FRUE"),"Never used"),
  row("External, pooled",external,"Reference"),
  tibble(label="External, random effects",N=nrow(external),events=sum(external$EXP_LABEL),AUC=plogis(re$b[1]),lo=plogis(re$ci.lb),hi=plogis(re$ci.ub),group="Reference")) %>%
  mutate(text=ifelse(label=="External, random effects",paste0(label," (4 cohorts)"),paste0(label," (",N,"/",events,")")),
         pos=rev(seq_len(n())) + ifelse(label%in%c("External, pooled","External, random effects"),0,.6))
published <- tibble(label=c("PRONIA, all","ZInEP","FePsy"),AUC=c(.73,.71,.71)) %>% left_join(select(forest,label,pos),by="label")
p1a <- ggplot(forest,aes(y=pos)) +
  geom_vline(xintercept=.5,colour=INK2,linewidth=.3,linetype=2) +
  geom_errorbar(aes(xmin=lo,xmax=hi,colour=group),width=0,linewidth=.6,orientation="y") +
  geom_point(aes(x=AUC,colour=group,shape=group),size=2.2) +
  geom_point(data=published,aes(x=AUC,y=pos+.32),shape=1,size=1.9,colour=INK2,stroke=.5) +
  scale_colour_manual(values=COL_EXT,guide="none") + scale_shape_manual(values=SHAPE_EXT,guide="none") +
  scale_y_continuous(breaks=forest$pos,labels=forest$text,expand=expansion(add=.6)) +
  scale_x_continuous(limits=c(.3,.9),breaks=seq(.3,.9,.1)) +
  labs(x="AUC (95% CI)",y=NULL) +
  theme_figure + theme(panel.grid.major.y=element_blank())

# ---- Figure 1B: mean predicted vs observed risk by cohort (dumbbell) ------------------
cal_cohort <- d %>% group_by(cohort) %>% summarise(predicted=mean(Probs_platt),observed=mean(EXP_LABEL),N=n(),.groups="drop") %>%
  mutate(group=status[cohort],cohort=factor(cohort,levels=rev(c("PRONIA","ZInEP","FePsy","NAPLS-3","MUC-FRUE"))))
cal_dec <- external %>% mutate(bin=ntile(Probs_platt,10)) %>% group_by(bin) %>%
  summarise(predicted=mean(Probs_platt),observed=mean(EXP_LABEL),.groups="drop")
cal_long <- cal_cohort %>% tidyr::pivot_longer(c(predicted,observed),names_to="quantity",values_to="value") %>%
  mutate(quantity=factor(ifelse(quantity=="predicted","Mean predicted","Observed"),levels=c("Mean predicted","Observed")))
p1b <- ggplot(cal_cohort,aes(y=cohort)) +
  geom_segment(aes(x=predicted,xend=observed,yend=cohort),colour="#b5b4ae",linewidth=.6) +
  geom_point(data=filter(cal_long,quantity=="Observed"),aes(x=value,shape=quantity),size=2,colour="#0b0b0b") +
  geom_point(data=filter(cal_long,quantity=="Mean predicted"),aes(x=value,shape=quantity),size=3,colour="#0b0b0b",fill=NA,stroke=.6) +
  scale_shape_manual(values=c(`Mean predicted`=23,Observed=16),breaks=c("Mean predicted","Observed")) +
  scale_x_continuous(limits=c(0,.6),breaks=seq(0,.6,.2),labels=function(x)paste0(x*100,"%")) +
  labs(x="Risk of transition",y=NULL) + guides(shape=guide_legend(nrow=2)) +
  theme_figure + theme(panel.grid.major.y=element_blank())

# ---- Figure 1C: decision curves (external cohorts) -----------------------------------
dca <- read.csv(file.path(out_dir,"tables","A3_decision_curves.csv")) %>% filter(grepl("^External",set)) %>%
  mutate(strategy=recode(strategy,`MRI model, original probabilities`="Original probabilities",
                         `MRI model, held-out recalibration`="Held-out recalibration"))
DCA_COL <- c(`Original probabilities`="#2a78d6",`Held-out recalibration`="#eb6834",`Treat all`="#3d3c39",`Treat none`="#8a8984")
DCA_LT <- c(`Original probabilities`="solid",`Held-out recalibration`="longdash",`Treat all`="dotted",`Treat none`="solid")
p1c <- ggplot(dca,aes(threshold,net_benefit,colour=strategy,linetype=strategy)) +
  geom_line(linewidth=.6) + scale_colour_manual(values=DCA_COL) + scale_linetype_manual(values=DCA_LT) +
  coord_cartesian(xlim=c(.05,.5),ylim=c(-.07,.12)) + scale_x_continuous(breaks=seq(.1,.5,.1),labels=function(x)paste0(x*100,"%")) +
  labs(x="Threshold probability",y="Net benefit") + guides(colour=guide_legend(nrow=2),linetype=guide_legend(nrow=2)) + theme_figure

# ---- Figure 1D: held-out gain from adding the MRI score ------------------------------
a6 <- read.csv(file.path(out_dir,"tables","A6_loco_clinical_vs_clinical_plus_mri_by_cohort.csv")) %>%
  filter(variant=="without_functioning") %>%
  mutate(heldout_cohort=recode(heldout_cohort,`Pooled held-out`="Pooled"),
         group=ifelse(heldout_cohort=="PRONIA","Reference",ifelse(heldout_cohort=="Pooled","Pooled",status[heldout_cohort])),
         text=sprintf("%s (%d/%d)",heldout_cohort,N,events),
         auc=sprintf("%.2f vs %.2f",clinical_AUC,clinical_plus_mri_AUC),
         pos=c(PRONIA=4.4,ZInEP=3.4,`NAPLS-3`=2.4,Pooled=1)[heldout_cohort])
p1d <- ggplot(a6,aes(y=pos)) +
  geom_vline(xintercept=0,colour=INK2,linewidth=.3,linetype=2) +
  geom_errorbar(aes(xmin=AUC_difference_lo,xmax=AUC_difference_hi,colour=group),width=0,linewidth=.6,orientation="y") +
  geom_point(aes(x=AUC_difference,colour=group,shape=group),size=2.2) +
  geom_text(aes(x=.27,label=auc),hjust=1,size=2.2,colour=INK2) +
  annotate("text",x=.27,y=5.05,label="AUC, symptoms vs\nsymptoms + MRI",hjust=1,size=2.2,colour="#0b0b0b",lineheight=.9) +
  scale_colour_manual(values=c(COL_EXT,Pooled="#0b0b0b"),guide="none") + scale_shape_manual(values=c(SHAPE_EXT,Pooled=18),guide="none") +
  scale_y_continuous(breaks=a6$pos,labels=a6$text,limits=c(.5,5.5)) +
  scale_x_continuous(limits=c(-.12,.27),breaks=seq(-.1,.15,.05)) +
  labs(x="AUC gain from adding the MRI score (95% CI)",y="Held-out cohort (N/transitions)") +
  theme_figure + theme(panel.grid.major.y=element_blank())

fig1 <- (p1a + p1b) / (free(p1c,type="label",side="l") + p1d) + plot_layout(heights=c(1,1)) + plot_annotation(tag_levels="A")
ggsave(file.path(fig_dir,"Figure1.pdf"),fig1,width=7.4,height=6.2,device=cairo_pdf)
ggsave(file.path(fig_dir,"Figure1.png"),fig1,width=7.4,height=6.2,dpi=600)

# ---- Figure 2: score densities, nontransition shift and separation ------------------------
cohort_levels <- c("PRONIA","ZInEP","FePsy","NAPLS-3","MUC-FRUE")
OUT_COL <- c(Nontransition="#8a8984",Transition="#eb6834")
COL2 <- c(COL,Development="#3d3c39"); SHAPE2 <- c(SHAPE,Development=15)
status2 <- c(status,PRONIA="Development")
smd_ci <- function(x1,x0) {
  n1 <- length(x1); n0 <- length(x0); s <- (mean(x1)-mean(x0))/sqrt((var(x1)+var(x0))/2)
  se <- sqrt((n1+n0)/(n1*n0)+s^2/(2*(n1+n0))); c(s,s-1.96*se,s+1.96*se)
}
ypos <- setNames(rev(seq_along(cohort_levels)),cohort_levels); HMAX <- .6; xr <- c(-1.45,1.45)
dens <- d %>% mutate(outcome=ifelse(EXP_LABEL==1,"Transition","Nontransition")) %>% group_by(cohort,outcome) %>%
  group_modify(~{k<-density(.x$Mean_Score,adjust=1.25,from=xr[1],to=xr[2],n=400);tibble(x=k$x,y=k$y)}) %>% ungroup() %>%
  group_by(cohort,outcome) %>% mutate(base=ypos[cohort],top=base+y/max(y)*HMAX,grp=paste(cohort,outcome)) %>% ungroup()
ref0 <- d$Mean_Score[d$cohort=="PRONIA" & d$EXP_LABEL==0]
stats <- d %>% group_by(cohort) %>% group_modify(~{
  sv <- smd_ci(.x$Mean_Score[.x$EXP_LABEL==1],.x$Mean_Score[.x$EXP_LABEL==0])
  hv <- smd_ci(.x$Mean_Score[.x$EXP_LABEL==0],ref0)
  tibble(N=nrow(.x),ev=sum(.x$EXP_LABEL),m0=mean(.x$Mean_Score[.x$EXP_LABEL==0]),m1=mean(.x$Mean_Score[.x$EXP_LABEL==1]),
         sep=sv[1],sep_lo=sv[2],sep_hi=sv[3],sh=hv[1],sh_lo=hv[2],sh_hi=hv[3])}) %>% ungroup() %>%
  mutate(y=ypos[cohort],grp=status2[cohort],
         lab=sprintf("%s (%d/%d)",cohort,N,ev),
         arrow_y=y+HMAX+.16,sep_txt=sub("-","−",sprintf("%.2f SD",sep)))
p2a <- ggplot() +
  geom_ribbon(data=filter(dens,outcome=="Nontransition"),aes(x,ymin=base,ymax=top,group=grp),fill=OUT_COL["Nontransition"],alpha=.35) +
  geom_line(data=filter(dens,outcome=="Nontransition"),aes(x,top,group=grp),colour=OUT_COL["Nontransition"],linewidth=.35) +
  geom_ribbon(data=filter(dens,outcome=="Transition"),aes(x,ymin=base,ymax=top,group=grp),fill=OUT_COL["Transition"],alpha=.1) +
  geom_line(data=filter(dens,outcome=="Transition"),aes(x,top,group=grp),colour=OUT_COL["Transition"],linewidth=.6) +
  geom_vline(xintercept=median(ref0),colour=INK2,linewidth=.3,linetype=2) +
  geom_segment(data=stats,aes(x=m0,xend=m0,y=y,yend=y+.14),colour=INK2,linewidth=.6) +
  geom_segment(data=stats,aes(x=m1,xend=m1,y=y,yend=y+.14),colour=OUT_COL["Transition"],linewidth=.6) +
  geom_segment(data=stats,aes(x=m0,xend=m1,y=arrow_y,yend=arrow_y),colour="#0b0b0b",linewidth=.4,arrow=arrow(length=unit(3,"pt"),type="closed")) +
  geom_text(data=stats,aes(x=pmax(m0,m1)+.06,y=arrow_y,label=sep_txt),hjust=0,size=2.5,colour="#0b0b0b",fontface="bold") +
  annotate("text",x=median(ref0)+.03,y=max(ypos)+1.06,hjust=0,vjust=1,size=2.3,colour=INK2,label="PRONIA nontransition median") +
  annotate("segment",x=.5,xend=.64,y=1.62,yend=1.62,colour=OUT_COL["Nontransition"],linewidth=1.6,alpha=.7) +
  annotate("segment",x=.5,xend=.64,y=1.46,yend=1.46,colour=OUT_COL["Transition"],linewidth=1.6,alpha=.7) +
  annotate("text",x=.7,y=c(1.62,1.46),label=c("Nontransition","Transition"),hjust=0,size=2.4,colour="#0b0b0b") +
  scale_y_continuous(breaks=stats$y,labels=stats$lab,limits=c(.85,max(ypos)+1.1),expand=c(0,0)) +
  scale_x_continuous(limits=xr,breaks=seq(-1,1,.5),expand=c(0,0)) +
  labs(x="MRI risk score (SVM decision value)",y=NULL) +
  theme_figure + theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(hjust=1,lineheight=.95,vjust=.1))
# ---- Figure 2B: observed vs case-mix-expected AUC ---------------------------------------
c1 <- read.csv(file.path(out_dir,"tables","C1_casemix_standardized_reference.csv")) %>% filter(covariates=="Age, sex")
obs <- bind_rows(lapply(cohort_levels[-1],function(co){x<-filter(d,cohort==co);a<-delong(x$Mean_Score,x$EXP_LABEL)
  tibble(cohort=co,N=nrow(x),ev=sum(x$EXP_LABEL),AUC=a[2],lo=a[1],hi=a[3])}))
pr <- filter(d,cohort=="PRONIA",Studygroup=="CHR"); pr_auc <- delong(pr$Mean_Score,pr$EXP_LABEL)
napls <- filter(d,cohort=="NAPLS-3"); youngest <- min(pr$age)
split3 <- bind_rows(lapply(list(inside=filter(napls,age>=youngest),below=filter(napls,age<youngest)),function(x){a<-delong(x$Mean_Score,x$EXP_LABEL)
  tibble(N=nrow(x),ev=sum(x$EXP_LABEL),AUC=a[2],lo=a[1],hi=a[3])}),.id="part")
pos3 <- c(PRONIA=7.6,ZInEP=6.4,FePsy=5.4,`NAPLS-3`=4.4,inside=3.65,below=3.05,`MUC-FRUE`=1.9); OFF <- .14
std <- inner_join(transmute(c1,cohort=external_cohort,std=reference_AUC_casemix_standardized,std_lo=lo,std_hi=hi),obs,by="cohort") %>%
  mutate(y=pos3[cohort],grp=status[cohort],lab=sprintf("%s (%d/%d)",cohort,N,ev),
         txt=sprintf("%.2f vs %.2f (%s)",AUC,std,sub("-","−",sprintf("%.2f",AUC-std))))
split3 <- mutate(split3,y=pos3[part],lab=ifelse(part=="inside",sprintf("NAPLS-3, within PRONIA age range (%d/%d)",N,ev),
                                                sprintf("NAPLS-3, younger than PRONIA (%d/%d)",N,ev)),txt=sprintf("%.2f",AUC))
ref3 <- tibble(y=pos3[["PRONIA"]],AUC=pr_auc[2],lo=pr_auc[1],hi=pr_auc[3],lab=sprintf("PRONIA CHR-P, unweighted (%d/%d)",nrow(pr),sum(pr$EXP_LABEL)),txt=sprintf("%.2f",pr_auc[2]))
labs3 <- bind_rows(select(ref3,y,lab),select(std,y,lab),select(split3,y,lab)); XR <- c(.25,1.06)
p3 <- ggplot() +
  annotate("rect",xmin=-Inf,xmax=Inf,ymin=pos3[["below"]]-.4,ymax=pos3[["NAPLS-3"]]+.45,fill="#f6f5f2") +
  annotate("segment",x=.5,xend=.5,y=.75,yend=8,colour=INK2,linewidth=.3,linetype=2) +
  geom_errorbar(data=ref3,aes(y=y,xmin=lo,xmax=hi),width=0,linewidth=.6,colour=COL2[["Development"]],orientation="y") +
  geom_point(data=ref3,aes(y=y,x=AUC),shape=15,size=2.2,colour=COL2[["Development"]]) +
  geom_segment(data=std,aes(y=y,yend=y,x=std,xend=AUC),colour="#b5b4ae",linewidth=1.6,alpha=.6) +
  geom_errorbar(data=std,aes(y=y+OFF,xmin=std_lo,xmax=std_hi),width=0,linewidth=.45,colour=INK2,orientation="y") +
  geom_errorbar(data=std,aes(y=y-OFF,xmin=lo,xmax=hi,colour=grp),width=0,linewidth=.6,orientation="y") +
  geom_point(data=std,aes(y=y,x=std),shape=22,size=2.4,fill="white",colour="#0b0b0b",stroke=.6) +
  geom_point(data=std,aes(y=y,x=AUC,colour=grp,shape=grp),size=2.3) +
  geom_errorbar(data=split3,aes(y=y,xmin=lo,xmax=hi),width=0,linewidth=.45,colour=EXT_COL,alpha=.7,orientation="y") +
  geom_point(data=split3,aes(y=y,x=AUC),shape=21,size=1.8,colour=EXT_COL,fill="white",stroke=.6) +
  geom_text(data=bind_rows(select(ref3,y,txt),select(std,y,txt),select(split3,y,txt)),aes(y=y,x=XR[2],label=txt),hjust=1,size=2.3,colour=INK2) +
  annotate("text",x=XR[2],y=8.35,label="Observed vs expected\nAUC (difference)",hjust=1,vjust=0,size=2.3,colour="#0b0b0b",lineheight=.9) +
  annotate("point",x=.283,y=8.55,shape=16,size=2.2,colour=EXT_COL) +
  annotate("text",x=.31,y=8.55,label="Observed AUC in the external cohort (colored line, DeLong 95% CI)",hjust=0,size=2.3) +
  annotate("point",x=.283,y=8.15,shape=22,size=2.4,fill="white",colour="#0b0b0b",stroke=.6) +
  annotate("text",x=.31,y=8.15,label="Expected AUC, PRONIA CHR-P reweighted to the cohort's age and sex (gray line, bootstrap 95% CI)",hjust=0,size=2.3) +
  scale_colour_manual(values=COL_EXT,guide="none") + scale_shape_manual(values=SHAPE_EXT,guide="none") +
  scale_y_continuous(breaks=labs3$y,labels=labs3$lab,limits=c(.75,8.9),expand=c(0,0)) +
  scale_x_continuous(limits=XR,breaks=seq(.3,.9,.1)) + labs(x="AUC (95% CI)",y=NULL) +
  theme_figure + theme(axis.text.y=element_text(hjust=1))
fig2 <- p2a / p3 + plot_layout(heights=c(3.9,3.6)) + plot_annotation(tag_levels="A")
ggsave(file.path(fig_dir,"Figure2.pdf"),fig2,width=7.2,height=7.6,device=cairo_pdf)
ggsave(file.path(fig_dir,"Figure2.png"),fig2,width=7.2,height=7.6,dpi=600)

write.csv(forest,file.path(fig_dir,"Figure1A_source.csv"),row.names=FALSE)
write.csv(cal_cohort,file.path(fig_dir,"Figure1B_cohort_source.csv"),row.names=FALSE)
write.csv(cal_dec,file.path(fig_dir,"Figure1B_decile_source.csv"),row.names=FALSE)
write.csv(select(stats,-lab,-arrow_y,-sep_txt),file.path(fig_dir,"Figure2_separation_source.csv"),row.names=FALSE)
write.csv(bind_rows(select(std,-lab,-txt),select(split3,-lab,-txt)),file.path(fig_dir,"Figure2B_source.csv"),row.names=FALSE)

# ---- eFigure 1: decision curves, external MRI model and LOCO clinical models ---------
dca_all <- read.csv(file.path(out_dir,"tables","A3_decision_curves.csv")) %>%
  mutate(panel=ifelse(grepl("^External",set),"External cohorts, MRI model (685/106)","Symptom cohorts, leave-one-cohort-out (658/84)"),
         strategy=recode(strategy,`MRI model, original probabilities`="Original probabilities",
                         `MRI model, held-out recalibration`="Held-out recalibration",`Clinical + MRI score`="Clinical + MRI score"))
EDCA_COL <- c(DCA_COL,`Clinical only`="#2a78d6",`Clinical + MRI score`="#eb6834")
EDCA_LT <- c(DCA_LT,`Clinical only`="solid",`Clinical + MRI score`="longdash")
edca_panel <- function(p) {
  x <- filter(dca_all,panel==p)
  ggplot(x,aes(threshold,net_benefit,colour=strategy,linetype=strategy)) +
    geom_hline(yintercept=0,colour=GRID,linewidth=.25) + geom_line(linewidth=.6) +
    scale_colour_manual(values=EDCA_COL) + scale_linetype_manual(values=EDCA_LT) +
    coord_cartesian(xlim=c(.05,.5),ylim=c(-.16,.12)) + scale_x_continuous(breaks=seq(.1,.5,.1),labels=function(v)paste0(v*100,"%")) +
    labs(title=p,x="Threshold probability",y="Net benefit") + guides(colour=guide_legend(nrow=2),linetype=guide_legend(nrow=2)) +
    theme_figure + theme(plot.title=element_text(size=8,face="bold"))
}
# calibration by cohort after intercept updating (external_auc_synthesis_and_calibration.R)
cal_pts <- read.csv(file.path(out_dir,"tables","D2_calibration_points.csv"))
cal_curve <- read.csv(file.path(out_dir,"tables","D2_calibration_curves.csv"))
cal_fit <- read.csv(file.path(out_dir,"tables","D2_calibration_by_cohort.csv"))
cal_panel <- function(co) {
  f <- filter(cal_fit,cohort==co)
  ttl <- sprintf("%s (%d/%d)\nslope %.2f (%.2f to %.2f)",co,f$N,f$events,f$slope,f$slope_lo,f$slope_hi)
  ggplot(filter(cal_pts,cohort==co),aes(predicted,observed)) +
    geom_abline(slope=1,intercept=0,colour=INK2,linewidth=.3,linetype=2) +
    geom_line(data=filter(cal_curve,cohort==co,p<=.8),aes(p,smooth),colour="#2a78d6",linewidth=.6) +
    geom_errorbar(aes(ymin=lo,ymax=hi),width=0,colour="#0b0b0b",linewidth=.4) +
    geom_point(size=1.6,colour="#0b0b0b") +
    coord_cartesian(xlim=c(0,.8),ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,.8,.4),labels=function(v)paste0(v*100,"%")) +
    scale_y_continuous(breaks=seq(0,1,.25),labels=function(v)paste0(v*100,"%")) +
    labs(title=ttl,x="Predicted risk",y=if(co=="NAPLS-3")"Observed proportion" else NULL) +
    theme_figure + theme(plot.title=element_text(size=7,face="bold"))
}
top <- edca_panel("External cohorts, MRI model (685/106)") + edca_panel("Symptom cohorts, leave-one-cohort-out (658/84)")
bottom <- cal_panel("NAPLS-3") + cal_panel("MUC-FRUE") + cal_panel("ZInEP") + cal_panel("FePsy") + plot_layout(nrow=1)
efig1 <- (top / bottom) + plot_layout(heights=c(1.25,1)) + plot_annotation(tag_levels="A")
ggsave(file.path(fig_dir,"eFigure1.pdf"),efig1,width=7.4,height=6.4,device=cairo_pdf)
ggsave(file.path(fig_dir,"eFigure1.png"),efig1,width=7.4,height=6.4,dpi=600)
message("Figures written to ",fig_dir)

# ---- eFigure (supplement): external-cohort AUCs with random-effects summary ---------------
pi_re <- predict(re)  # Hartung-Knapp CI and 95% prediction interval on the logit scale
ext_rows <- filter(forest,label%in%c("ZInEP","FePsy","NAPLS-3","MUC-FRUE")) %>% mutate(pos=4:1)
summ <- tibble(label="Random-effects summary",AUC=plogis(pi_re$pred),lo=plogis(pi_re$ci.lb),hi=plogis(pi_re$ci.ub),
               pi_lo=plogis(pi_re$pi.lb),pi_hi=plogis(pi_re$pi.ub),pos=-.2)
pe <- ggplot() +
  geom_vline(xintercept=.5,colour=INK2,linewidth=.3,linetype=2) +
  geom_errorbar(data=ext_rows,aes(y=pos,xmin=lo,xmax=hi,colour=group),width=0,linewidth=.6,orientation="y") +
  geom_point(data=ext_rows,aes(y=pos,x=AUC,colour=group,shape=group),size=2.2) +
  geom_errorbar(data=summ,aes(y=pos,xmin=pi_lo,xmax=pi_hi),width=0,linewidth=.4,linetype=3,colour=INK2,orientation="y") +
  geom_errorbar(data=summ,aes(y=pos,xmin=lo,xmax=hi),width=0,linewidth=.8,colour="#3d3c39",orientation="y") +
  geom_point(data=summ,aes(y=pos,x=AUC),shape=18,size=3.4,colour="#3d3c39") +
  geom_text(data=bind_rows(ext_rows,summ),aes(y=pos,x=.97,label=sprintf("%.2f (%.2f-%.2f)",AUC,lo,hi)),hjust=1,size=2.3,colour=INK2) +
  annotate("text",x=.97,y=4.75,label="AUC (95% CI)",hjust=1,size=2.3,colour="#0b0b0b") +
  scale_colour_manual(values=COL,guide="none") + scale_shape_manual(values=SHAPE,guide="none") +
  scale_y_continuous(breaks=c(ext_rows$pos,summ$pos),labels=c(ext_rows$text,"Random-effects summary (4 cohorts)"),limits=c(-.7,5)) +
  scale_x_continuous(limits=c(.25,.97),breaks=seq(.3,.9,.1)) + labs(x="AUC",y=NULL) + theme_figure
ggsave(file.path(fig_dir,"eFigure_random_effects.png"),pe,width=7.2,height=7.2*3995928/6080760,dpi=600)
write.csv(transmute(summ,AUC,lo,hi,pi_lo,pi_hi,I2=re$I2),file.path(fig_dir,"eFigure_random_effects_source.csv"),row.names=FALSE)

# ---- eFigure (supplement): MRI score x SIPS-N moderation, NAPLS-3 reference ---------------
sv <- c("Mean_Score","COGDIS_score","SIPS_Positiv_Gesamt","SIPS_Negativ_Gesamt")
dm <- filter(d,cohort%in%c("PRONIA","NAPLS-3","ZInEP")) %>% filter(if_all(all_of(sv),~!is.na(.x))) %>%
  mutate(across(all_of(sv),~as.numeric(scale(.x)),.names="{.col}_z"),y=EXP_LABEL,cohort_f=relevel(factor(cohort),"NAPLS-3"))
m1 <- glm(y~Mean_Score_z*(COGDIS_score_z+SIPS_Positiv_Gesamt_z+SIPS_Negativ_Gesamt_z)+cohort_f,data=dm,family=binomial())
nd <- tidyr::expand_grid(Mean_Score_z=seq(quantile(dm$Mean_Score_z,.025),quantile(dm$Mean_Score_z,.975),length.out=160),
                         SIPS_Negativ_Gesamt_z=c(-1,0,1),COGDIS_score_z=0,SIPS_Positiv_Gesamt_z=0,
                         cohort_f=factor("NAPLS-3",levels=levels(dm$cohort_f)))
pp <- predict(m1,nd,type="link",se.fit=TRUE)
nd <- mutate(nd,p=plogis(pp$fit),lo=plogis(pp$fit-qnorm(.975)*pp$se.fit),hi=plogis(pp$fit+qnorm(.975)*pp$se.fit),
             level=factor(SIPS_Negativ_Gesamt_z,levels=c(-1,0,1),labels=c("SIPS-N −1 SD","SIPS-N mean","SIPS-N +1 SD")))
LEV_COL <- c(`SIPS-N −1 SD`="#f2b49b",`SIPS-N mean`="#eb6834",`SIPS-N +1 SD`="#8a3412")  # one hue, light to dark
pm <- ggplot(nd,aes(Mean_Score_z,p,colour=level,fill=level)) +
  geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.12,colour=NA) + geom_line(linewidth=.7) +
  scale_colour_manual(values=LEV_COL) + scale_fill_manual(values=LEV_COL) +
  scale_y_continuous(limits=c(0,NA),labels=function(v)paste0(v*100,"%")) +
  labs(x="MRI risk score (standardized)",y="Fitted transition probability") + theme_figure
ggsave(file.path(fig_dir,"eFigure_moderation.png"),pm,width=5.2,height=4,dpi=600)
co <- summary(m1)$coefficients["Mean_Score_z:SIPS_Negativ_Gesamt_z",]
message(sprintf("moderation check: n %d, transitions %d, beta %.3f (%.3f to %.3f), P %.3f",nrow(dm),sum(dm$y),co[1],co[1]-1.96*co[2],co[1]+1.96*co[2],co[4]))
message(sprintf("random effects: AUC %.3f, HK CI %.3f-%.3f, PI %.3f-%.3f, I2 %.1f",plogis(pi_re$pred),plogis(pi_re$ci.lb),plogis(pi_re$ci.ub),plogis(pi_re$pi.lb),plogis(pi_re$pi.ub),re$I2))
