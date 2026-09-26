#!/usr/bin/env Rscript
# Render statistical figures from a completed run's aggregate exports.
suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(trailingOnly=TRUE)
value <- function(flag,default=NULL) {i<-match(flag,args);if(!is.na(i)&&i<length(args))args[[i+1L]]else default}
if('--help'%in%args){cat('Rscript R/render_figures.R --input-dir RUN --output-dir FIGURES [--figure NAME]\n');quit(status=0)}
input<-value('--input-dir'); out<-value('--output-dir'); which<-value('--figure','all')
if(is.null(input)||!dir.exists(input)) stop('--input-dir must be the supplied current run output directory')
if(is.null(out)) stop('--output-dir is required')
if(normalizePath(out,mustWork=FALSE)==normalizePath(input)) stop('Output must differ from run input')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,'plotted_values'),showWarnings=FALSE)
read_result<-function(rel){p<-file.path(input,rel);if(!file.exists(p))stop('Missing run input: ',rel);read.csv(p,check.names=FALSE)}
allowed<-c('primary_1d_landscape_six_panel','sips_n_brainage_local_auc_surface','canonical_linear_sips_n_moderation','external_auc_random_effects_forest','external_calibration_original_vs_apparent','gam_sips_n_uncertainty_supported_regions','component_selected_statistic_heatmap')
if(!which%in%c('all',allowed))stop('Unknown --figure. Choose: ',paste(allowed,collapse=', '))
emit<-function(p,data,stem,width,height){
  mode<-jsonlite::read_json(file.path(input,'run_manifest.json'))$mode
  if(mode=='synthetic')p<-p+labs(caption='SYNTHETIC TEST DATA — does not reproduce study results')
  for(ext in c('png','svg','pdf')){
    device<-switch(ext,svg=grDevices::svg,pdf=grDevices::cairo_pdf,png='png')
    ggsave(file.path(out,paste0(stem,'.',ext)),p,width=width,height=height,dpi=320,device=device,bg='white')
  }
  write.csv(data,file.path(out,'plotted_values',paste0(stem,'.csv')),row.names=FALSE)
  message('PASS: render ',stem,' (',nrow(data),' rows)')
}
if(which%in%c('all',allowed[1])){
  one_d<-read_result('figure_data/primary_1d_window_details_figure_ready.csv')
  labels<-c(age='Age',brainage_sbc='BrainAGE',cogdis='COGDIS',sips_p='SIPS-P',sips_n='SIPS-N',psychosoz='Functioning')
  one_d$display_moderator<-factor(one_d$display_moderator,levels=unname(labels))
  stopifnot(all(one_d$bootstrap_valid[one_d$supported]>0),all(one_d$bootstrap_invalid==0L))
  B<-unique(one_d$bootstrap_valid[one_d$supported]);stopifnot(length(B)==1)
  ann<-one_d%>%group_by(display_moderator)%>%summarise(center=min(center),panel_note=first(panel_note),.groups='drop')
  p<-ggplot(one_d,aes(center,AUC,group=segment_id))+
    geom_ribbon(aes(ymin=AUC_bootstrap_lo,ymax=AUC_bootstrap_hi),fill='#5B8DB8',alpha=.22)+geom_line(color='#24557A',linewidth=.7)+
    geom_point(aes(size=N,fill=events),shape=21,color='#17384F',stroke=.25)+geom_hline(yintercept=.5,linetype=2,color='grey55')+
    geom_text(data=ann,aes(x=center,y=.98,label=panel_note),inherit.aes=FALSE,hjust=0,vjust=1,size=3.1)+
    facet_wrap(~display_moderator,ncol=2,scales='free_x')+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25))+
    scale_fill_viridis_c(option='C',end=.85,name='Transitions in window')+scale_size_continuous(range=c(1.8,5.5),name='Participants in window')+
    labs(x='Moderator window centre (midpoint of fixed observed-scale window)',y=paste0('Local AUC (',B,'-bootstrap 95% interval)'),
      subtitle='Supported windows require at least 7 transitions and 7 non-transitions; nominal within-cohort p and six-test BH q are shown')+
    theme_minimal(base_size=11)+theme(legend.position='bottom',panel.grid.minor=element_blank(),strip.text=element_text(face='bold'))
  emit(p,one_d,allowed[1],10.5,12)
}
if(which%in%c('all',allowed[2])){
  d<-read_result('figure_data/primary_2d_cell_details.csv')%>%filter(pair=='brainage_sbc_x_sips_n')
  if(!nrow(d))stop('Expected canonical BrainAGE x SIPS-N surface is absent; refusing automatic axis swap')
  stopifnot(all(d$AUC[d$supported]>=0&d$AUC[d$supported]<=1))
  p<-ggplot(d,aes(center_1,center_2,fill=AUC))+geom_tile(color='white',linewidth=.15)+geom_text(data=filter(d,supported),aes(label=paste0('n=',N,'\ne=',events)),size=2)+scale_fill_gradient2(low='#3b0f70',mid='white',high='#d7191c',midpoint=.5,limits=c(0,1),na.value='grey88')+labs(x='BrainAGE window center (years)',y='SIPS-N window center',fill='Local AUC')+theme_minimal(base_size=11)+theme(panel.grid=element_blank())
  emit(p,d,allowed[2],7,5.4)
}
if(which%in%c('all',allowed[3])){
  nd<-read_result('figure_data/canonical_linear_moderation_curve.csv')
  nd$SIPS_N<-factor(nd$SIPS_N,levels=c('Low (-1 SD)','Mean','High (+1 SD)'))
  p<-ggplot(nd,aes(Mean_Score_z,probability,color=SIPS_N,fill=SIPS_N))+geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.14,color=NA)+geom_line(linewidth=.9)+labs(x='MRI risk score (standardized within analysis set)',y='Fitted transition probability',subtitle='NAPLS-3 reference; COGDIS and SIPS-P fixed at analysis-set means',color='SIPS-N',fill='SIPS-N')+theme_minimal(base_size=11)+theme(legend.position='bottom')
  emit(p,nd,allowed[3],7,5)
}
if(which%in%c('all',allowed[4])){
  d<-read_result('figure_data/external_auc_random_effects_forest_data.csv')
  cohort_levels<-rev(d$cohort)
  p<-ggplot(d,aes(y=factor(cohort,levels=cohort_levels),x=AUC,xmin=AUC_lo,xmax=AUC_hi))+
    geom_vline(xintercept=.5,linetype=2,color='grey50')+geom_errorbar(orientation='y',width=.18,color='#1f4e79')+
    geom_point(aes(shape=cohort=='Random-effects summary',size=cohort=='Random-effects summary'),color='#1f4e79')+
    scale_shape_manual(values=c(`FALSE`=16,`TRUE`=18),guide='none')+scale_size_manual(values=c(`FALSE`=2.7,`TRUE`=4),guide='none')+
    scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.2))+labs(x='AUC (95% CI)',y=NULL,title='External-cohort discrimination')+theme_minimal(base_size=11)+theme(panel.grid.major.y=element_blank())
  emit(p,d,allowed[4],7,4.6)
}
if(which%in%c('all',allowed[5])){
  d<-read_result('figure_data/external_original_vs_apparent_calibration_curve.csv')
  p<-ggplot(d,aes(mean_predicted,observed_fraction,color=method))+geom_abline(slope=1,intercept=0,linetype=2,color='grey50')+geom_line()+geom_point(aes(size=events))+coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(x='Mean predicted risk',y='Observed transition fraction',color='Probability',size='Transitions',subtitle='Apparent updates were fitted and assessed in the same external data')+guides(color=guide_legend(nrow=2,byrow=TRUE),size=guide_legend(nrow=1))+theme_minimal(base_size=11)+theme(legend.position='bottom',legend.box='vertical')
  emit(p,d,allowed[5],8.4,6.5)
}
if(which%in%c('all',allowed[6])){
  seg<-read_result('figure_data/gam_sips_n_uncertainty_segmented.csv')
  seg$SIPS_N<-factor(seg$SIPS_N,levels=c('Low (-1 SD)','Mean','High (+1 SD)'))
  supported_multi<-filter(seg,supported_region,segment_size>=2); supported_single<-filter(seg,supported_region,singleton_segment)
  p<-ggplot(seg,aes(Mean_Score_z,probability,color=SIPS_N,fill=SIPS_N))+
    geom_ribbon(data=supported_multi,aes(ymin=lo,ymax=hi,group=segment_id),alpha=.14,color=NA)+
    geom_line(data=filter(seg,segment_size>=2),aes(group=segment_id,linetype=supported_region),linewidth=.9)+
    geom_errorbar(data=supported_single,aes(ymin=lo,ymax=hi),width=0,linewidth=.6)+geom_point(data=supported_single,size=1.8)+
    scale_linetype_manual(values=c(`TRUE`='solid',`FALSE`='dotted'),labels=c(`TRUE`='Pooled local support ≥10',`FALSE`='Pooled local support <10'),name='Local support')+
    labs(x='MRI risk score (standardized within analysis set)',y='Fitted transition probability',color='SIPS-N',fill='SIPS-N',subtitle='Pointwise 95% intervals; support is pooled N within ±0.25 SD of both displayed score and SIPS-N values')+
    guides(fill='none',linetype=guide_legend(nrow=1),color=guide_legend(nrow=1))+theme_minimal(base_size=11)+theme(legend.position='bottom',legend.box='vertical')
  emit(p,seg,allowed[6],8.6,6.2)
}
if(which%in%c('all',allowed[7])){
  d<-read_result('tables/component_selected_statistic_key.csv')
  d<-d%>%mutate(display_context=factor(display_context,levels=unique(display_context)),display_component=factor(display_component,levels=paste('Comp',1:4)))
  p<-ggplot(d,aes(display_component,display_context,fill=rho))+geom_tile(color='white',linewidth=.3)+geom_text(aes(label=paste0(display_metric,'\n',sprintf('ρ %.2f',rho))),size=2.1)+scale_fill_gradient2(low='#2166ac',mid='white',high='#b2182b',midpoint=0,limits=c(-1,1))+facet_grid(dimension~.,scales='free_y',space='free_y')+labs(x='NeuroMiner component (mapping in accompanying key)',y=NULL,fill='Selected ρ',subtitle='Each cell selects the largest |ρ| among six summaries; q values use all 24 tests per context')+theme_minimal(base_size=10)+theme(panel.grid=element_blank(),strip.text.y=element_text(angle=0))
  emit(p,d,allowed[7],9,8.5)
}
