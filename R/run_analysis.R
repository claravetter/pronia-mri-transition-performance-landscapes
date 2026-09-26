#!/usr/bin/env Rscript
pipeline_main <- function(args=commandArgs(trailingOnly=TRUE)) {
  value<-function(flag,default=NULL){i<-match(flag,args);if(is.na(i))return(default);if(i==length(args)||startsWith(args[i+1],"--"))stop("Missing value for ",flag);args[i+1]}
  if("--help"%in%args||!length(args)){cat("Usage: Rscript R/run_analysis.R --mode synthetic|full --output NEW_DIR --config config/study.yml [--input AUTHORIZED.csv --correction RECORD.json --authorize-private-analysis] [--cores 1] [--component-map-dir MAP_DIR] [--check-input]\n--check-input validates full-mode input without inference or an output directory.\nUse python/run_pipeline.py for the complete statistical and rendering workflow.\n");return(invisible(NULL))}
  allowed<-c("--mode","--output","--config","--input","--correction","--authorize-private-analysis","--check-input","--cores","--component-map-dir")
  if(length(setdiff(args[startsWith(args,"--")],allowed)))stop("Unknown argument")
  script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]);root<-normalizePath(file.path(dirname(script),".."))
  mode<-value("--mode","full");if(!mode%in%c("full","synthetic"))stop("Unknown mode")
  if(mode=="full"&&!"--check-input"%in%args&&!"--authorize-private-analysis"%in%args)stop("Protected inference requires explicit --authorize-private-analysis")
  env<-environment()
  for(f in c("common","contract","calibration"))source(file.path(root,"R/lib",paste0(f,".R")),local=env)
  cfg_path<-normalizePath(value("--config",file.path(root,"config/study.yml")));cfg<-read_study_config(cfg_path)
  options(pronia.seed=cfg$permutations$seed_base)
  counts<-cfg[[mode]];n_perm<-counts$permutations;n_boot_local<-counts$local_auc_bootstrap;n_boot_cal<-counts$pooled_calibration_bootstrap
  n_boot_cohort<-counts$source_cohort_heldout_bootstrap;n_boot_loco<-counts$loco_bootstrap;n_fisher<-counts$fisher_simulations;n_kcheck<-counts$gam_kcheck;n_gam_sim<-counts$gam_simulations
  n_cores<-as.integer(value("--cores","1"));if(is.na(n_cores)||n_cores<1)stop("Invalid cores")
  map_root<-value("--component-map-dir","")
  correction_path<-value("--correction");correction<-if(is.null(correction_path))NULL else jsonlite::read_json(correction_path,simplifyVector=TRUE)
  input_path<-value("--input")
  if(mode=="synthetic"&&!is.null(input_path))stop("Synthetic mode refuses external input")
  if(mode=="full"&&(is.null(input_path)||!file.exists(input_path)))stop("Authorized input missing")
  if("--check-input"%in%args){d<-canonicalize_data(input_path,correction=correction,expected=cfg$cohorts);cat("Input validation complete; no inference. N=",nrow(d),"\n");return(invisible(NULL))}
  run_dir<-value("--output");if(is.null(run_dir)||file.exists(run_dir))stop("--output must be a new directory")
  run_dir<-external_output_path(run_dir,root)
  dir.create(run_dir,recursive=TRUE);run_dir<-normalizePath(run_dir);run_id<-basename(run_dir)
  tables_dir<-file.path(run_dir,"tables");figdata_dir<-file.path(run_dir,"figure_data");figures_dir<-file.path(run_dir,"figures");private_dir<-file.path(run_dir,"private")
  for(path in c(tables_dir,figdata_dir,figures_dir,private_dir))dir.create(path,recursive=TRUE)
  finished<-FALSE;on.exit(if(!finished)writeLines("status: FAILED\nPartial output; do not treat as a completed analysis.",file.path(run_dir,"RUN_STATUS.txt")),add=TRUE)
  writeLines(paste("status: RUNNING\nmode:",mode),file.path(run_dir,"RUN_STATUS.txt"))
  if(mode=="synthetic"){input_path<-file.path(private_dir,"synthetic_input.csv");write.csv(synthetic_data(cfg),input_path,row.names=FALSE)}
  d<-canonicalize_data(input_path,private_dir,correction,if(mode=="full")cfg$cohorts else NULL)
  file.copy(cfg_path,file.path(run_dir,"consumed_study.yml"))
  manifest<-list(run_id=run_id,mode=mode,empirical=mode=="full",input_sha256=digest::digest(file=input_path,algo="sha256"),config_sha256=digest::digest(file=cfg_path,algo="sha256"),correction_sha256=if(is.null(correction_path))NULL else digest::digest(file=correction_path,algo="sha256"),command=commandArgs(),resamples=counts,R_version=R.version.string,upstream_model_training=FALSE)
  jsonlite::write_json(manifest,file.path(run_dir,"run_manifest.json"),auto_unbox=TRUE,pretty=TRUE)
  writeLines(capture.output(sessionInfo()),file.path(run_dir,"sessionInfo.txt"))
  resampling_registry<-list()
add_population <- function(id,label,x,represented=NULL,required="") {
  a <- bind_rows(x %>% group_by(cohort) %>% summarise(N=n(),events=sum(EXP_LABEL),non_events=N-events,.groups="drop"),
                x %>% summarise(cohort="ALL",N=n(),events=sum(EXP_LABEL),non_events=N-events))
  if (is.null(represented)) a <- a %>% mutate(represented_N=NA_integer_,represented_events=NA_integer_)
  else {
    b <- bind_rows(represented %>% group_by(cohort) %>% summarise(represented_N=n(),represented_events=sum(EXP_LABEL),.groups="drop"),
                  represented %>% summarise(cohort="ALL",represented_N=n(),represented_events=sum(EXP_LABEL)))
    a <- a %>% left_join(b,by="cohort")
  }
  a %>% mutate(analysis_id=id,population=label,required_variables=required,.before=1)
}
mods <- c(age="age",brainage_sbc="BrainAGE_SBC_corr",cogdis="COGDIS_score",
          sips_p="SIPS_Positiv_Gesamt",sips_n="SIPS_Negativ_Gesamt",psychosoz="Psychosoz_aequiv")
surface_vars <- c(age="age",brainage_sbc="BrainAGE_SBC_corr",cogdis="COGDIS_score",sips_n="SIPS_Negativ_Gesamt")
pairs <- combn(names(surface_vars),2,simplify=FALSE)
population_rows <- list(
  add_population("POP_ALL","All canonical participants",d,required="outcome; MRI risk score; inherited class; original probability"),
  add_population("POP_EXTERNAL","Four external cohorts",filter(d,cohort!="PRONIA"),required="outcome; MRI risk score; inherited class; original probability"))
external <- filter(d,cohort!="PRONIA")
modules <- c("descriptive","validation","landscapes","moderation","gam","calibration",
             "cohort_calibration","loco","correlations","components","outputs")
for (module in modules) {
  stage_start<-Sys.time()
  message(format(stage_start,"%Y-%m-%d %H:%M:%S")," Running statistical module: ",module)
  sys.source(file.path(root,"R","modules",paste0(module,".R")),envir=env)
  message("Completed ",module," in ",round(as.numeric(difftime(Sys.time(),stage_start,units="secs")),1)," seconds")
}
writeLines(c("status: STATISTICS_COMPLETE","python_summaries_and_rendering: pending pipeline runner"),
           file.path(run_dir,"RUN_STATUS.txt"))
message("Statistical outputs complete: ",run_dir)

finished<-TRUE
}
if(sys.nframe()==0L)pipeline_main()
