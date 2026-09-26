suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(mgcv)
  library(metafor)
})

RNG_KIND <- "L'Ecuyer-CMRG"
RNG_NORMAL_KIND <- "Inversion"
RNG_SAMPLE_KIND <- "Rejection"

# Indexed random streams preserve reproducibility across worker schedules.
stream_cache <- new.env(parent=emptyenv())
stream_states <- function(analysis_id, n, base_seed=getOption("pronia.seed",20260913L)) {
  stopifnot(length(analysis_id)==1L, nzchar(analysis_id), n >= 1L)
  key <- paste(base_seed,analysis_id,sep=":")
  if (!exists(key,stream_cache,inherits=FALSE)) {
    old_kind <- RNGkind(); old <- if(exists(".Random.seed",.GlobalEnv)) .Random.seed else NULL
    on.exit({do.call(RNGkind,as.list(old_kind));if(is.null(old)) {if(exists(".Random.seed",.GlobalEnv)) rm(".Random.seed",envir=.GlobalEnv)} else assign(".Random.seed",old,.GlobalEnv)})
    seed <- strtoi(substr(digest::digest(key,algo="sha256",serialize=FALSE),1,7),16L)
    set.seed(seed,kind=RNG_KIND,normal.kind=RNG_NORMAL_KIND,sample.kind=RNG_SAMPLE_KIND)
    assign(key,list(.Random.seed),stream_cache)
  }
  states <- get(key,stream_cache)
  while(length(states)<n) states[[length(states)+1L]] <- parallel::nextRNGStream(states[[length(states)]])
  assign(key,states,stream_cache); states[seq_len(n)]
}
set_replicate_seed <- function(analysis_id, replicate, base_seed=getOption("pronia.seed",20260913L)) {
  state <- stream_states(analysis_id,replicate,base_seed)[[replicate]]
  assign(".Random.seed",state,.GlobalEnv)
  invisible(paste(state,collapse=":"))
}
deterministic_rows <- function(n, analysis_id, fun, cores=1L, base_seed=getOption("pronia.seed",20260913L), order=seq_len(n)) {
  stopifnot(n>=1L, cores>=1L, setequal(order,seq_len(n)), !anyDuplicated(order))
  states <- stream_states(analysis_id,n,base_seed)
  old_kind <- RNGkind(); old <- if(exists(".Random.seed",.GlobalEnv)) .Random.seed else NULL
  on.exit({do.call(RNGkind,as.list(old_kind));if(is.null(old)) {if(exists(".Random.seed",.GlobalEnv)) rm(".Random.seed",envir=.GlobalEnv)} else assign(".Random.seed",old,.GlobalEnv)})
  worker <- function(i) {assign(".Random.seed",states[[i]],.GlobalEnv); fun(i)}
  ans <- if(cores<=1L) lapply(order,worker) else parallel::mclapply(order,worker,mc.cores=cores,mc.set.seed=FALSE)
  if(any(vapply(ans,inherits,logical(1),"try-error"))) stop("Resampling worker failed: ",analysis_id)
  do.call(rbind,ans[match(seq_len(n),order)])
}

strict_binary <- function(x) {
  text <- as.character(x)
  if(any(!is.na(text) & !text %in% c("0","1"))) stop("Outcome/class values must be explicitly coded 0/1")
  as.integer(text)
}
reconcile_rows <- function(d, correction=NULL) {
  if(anyNA(d$Cases)||anyNA(d$cohort)||any(!nzchar(as.character(d$Cases)))) stop("Missing participant/cohort key")
  duplicates <- duplicated(d[c("cohort","Cases")]) | duplicated(d[c("cohort","Cases")],fromLast=TRUE)
  if(!any(duplicates)) return(d)

  # One external record authorizes one two-row clinical-field conflict. All
  # other fields, including the outcome and fixed model outputs, must agree.
  required <- c("key_sha256", "cohort", "field", "conflicting_values", "retained_value")
  scalar_text <- function(x) is.character(x) && length(x)==1L && !is.na(x) && nzchar(x)
  if(!is.list(correction) || anyDuplicated(names(correction)) ||
     !setequal(names(correction), required))
    stop("Duplicate correction requires key_sha256, cohort, field, conflicting_values and retained_value")
  if(!scalar_text(correction$key_sha256) || !grepl("^[a-f0-9]{64}$", correction$key_sha256) ||
     !scalar_text(correction$cohort) || !scalar_text(correction$field))
    stop("Invalid duplicate-correction identity or field")
  clinical_fields <- c("age", "sex", "COGDIS_score", "SIPS_Positiv_Gesamt",
                       "SIPS_Negativ_Gesamt", "Psychosoz_aequiv")
  field <- correction$field
  values <- correction$conflicting_values
  retained <- correction$retained_value
  if(!field %in% clinical_fields || !field %in% names(d))
    stop("Duplicate correction must name an available clinical field; identifiers, outcomes and model outputs are fixed")
  if(!is.numeric(values) || length(values)!=2L || any(!is.finite(values)) || anyDuplicated(values) ||
     !is.numeric(retained) || length(retained)!=1L || !is.finite(retained) || !retained %in% values)
    stop("Duplicate correction requires two distinct finite values and a retained value from that pair")
  block <- d[duplicates,]
  if(nrow(block)!=2L || length(unique(block$Cases))!=1L || any(block$cohort!=correction$cohort) ||
     digest::digest(as.character(block$Cases[[1]]),algo="sha256",serialize=FALSE)!=correction$key_sha256 ||
     !is.numeric(block[[field]]) || !setequal(block[[field]], values))
    stop("Duplicate does not match the correction's identity, cohort or conflicting values")
  fixed <- setdiff(names(block), field)
  if(!all(vapply(block[fixed],function(z) length(unique(z))==1L,logical(1))))
    stop("Duplicate rows disagree outside the authorized clinical field")
  keep <- !duplicates | (duplicates & d[[field]]==retained)
  d[which(keep),,drop=FALSE]
}
validate_pronia_groups <- function(d) {
  groups <- d$Studygroup[which(d$cohort == "PRONIA")]
  if(anyNA(groups) || any(!groups %in% c("CHR", "ROP")))
    stop("PRONIA Studygroup must be explicitly coded CHR or ROP without missing values")
  if(length(groups) && !any(groups == "CHR"))
    stop("PRONIA contains no CHR records for the configured sensitivity analysis")
  invisible(d)
}

chr_sensitivity_population <- function(d) {
  validate_pronia_groups(d)
  d[which(d$cohort != "PRONIA" | d$Studygroup == "CHR"), , drop=FALSE]
}

canonicalize_data <- function(input_path, private_dir=NULL, correction=NULL, expected=NULL) {
  required <- c("Cases","cohort","EXP_LABEL","PRED_LABEL","Mean_Score","Probs_platt","age","sex","Studygroup",
    "COGDIS_score","SIPS_Positiv_Gesamt","SIPS_Negativ_Gesamt","Psychosoz_aequiv","BrainAGE_corr",
    "BrainAGE_SBC_corr","PredictedAge_SBC_corr",paste0("Comp",1:4))
  d <- read.csv(input_path,check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(Cases="character"))
  if(anyDuplicated(names(d))||length(setdiff(required,names(d)))) stop("Missing or duplicated input columns")
  d$cohort[d$cohort=="NAPLES"] <- "NAPLS-3"
  # Standardized covariates are computed after analysis-population selection.
  stale_scales <- c("cogdis_z","sips_n_z","sips_p_z","psychosoz_z","decision_z",
                    "age_z","brainage_z","brainage_SBC_z",paste0("comp",1:4,"_z"))
  d <- d[setdiff(names(d),stale_scales)]
  d$EXP_LABEL <- strict_binary(d$EXP_LABEL); d$PRED_LABEL <- strict_binary(d$PRED_LABEL)
  numeric_cols <- setdiff(required,c("Cases","cohort","Studygroup"))
  if(!all(vapply(d[numeric_cols],is.numeric,logical(1)))) stop("Numeric schema mismatch; harmonize using the input dictionary")
  if(anyNA(d$EXP_LABEL)||anyNA(d$PRED_LABEL)||any(!is.finite(d$Mean_Score))||
     any(!is.finite(d$Probs_platt))||any(d$Probs_platt<0|d$Probs_platt>1)||
     any(!is.finite(as.matrix(d[paste0("Comp",1:4)])))) stop("Invalid fixed input fields")
  if(any(!is.na(d$sex)&!d$sex%in%c(1,2))) stop("Sex must follow documented 1/2 coding")
  d <- as_tibble(reconcile_rows(d,correction))
  validate_pronia_groups(d)
  d <- d[order(d$cohort,d$Cases),] # stable execution independent of input row order
  if(!is.null(expected)) {
    actual <- d %>% group_by(cohort) %>% summarise(N=n(),events=sum(EXP_LABEL),.groups="drop")
    target <- bind_rows(lapply(names(expected),function(co) tibble(cohort=co,N=expected[[co]]$N,events=expected[[co]]$events)))
    comparison <- full_join(actual,target,by="cohort",suffix=c("_actual","_expected"))
    if(!is.null(private_dir)) {dir.create(private_dir,recursive=TRUE,showWarnings=FALSE);write.csv(comparison,file.path(private_dir,"population_input_reconciliation.csv"),row.names=FALSE)}
    if(anyNA(comparison)||any(comparison$N_actual!=comparison$N_expected | comparison$events_actual!=comparison$events_expected))
      stop("Cohort counts do not match the configured population; inspect population_input_reconciliation.csv when an output directory is supplied")
  }
  if(!is.null(private_dir)) saveRDS(d,file.path(private_dir,"canonical_participants.rds"))
  d
}

eligible <- function(x, variables=character()) {
  cols <- unique(c("EXP_LABEL","Mean_Score",variables))
  ok <- rep(TRUE,nrow(x))
  for(v in cols) ok <- ok & !is.na(x[[v]]) & if(is.numeric(x[[v]])) is.finite(x[[v]]) else TRUE
  if(any(!is.na(x$EXP_LABEL)&!x$EXP_LABEL%in%c(0,1))) stop("Nonbinary outcome")
  ok
}
sample_indices <- function(ii,replace=FALSE) if(!length(ii)) integer() else ii[sample.int(length(ii),length(ii),replace=replace)]
safe_cor <- function(x,y,method="spearman") {
  ok <- is.finite(x)&is.finite(y)
  if(sum(ok)<3||length(unique(x[ok]))<2||length(unique(y[ok]))<2) return(NA_real_)
  cor(x[ok],y[ok],method=method)
}
scale_fit <- function(x,vars,id) {
  if(any(!is.finite(as.matrix(x[vars])))) stop("Select complete analysis population before scaling")
  params <- map_dfr(vars,function(v)tibble(scaling_id=id,variable=v,N=nrow(x),mean=mean(x[[v]]),SD=sd(x[[v]])))
  if(any(!is.finite(params$SD)|params$SD<=0)) stop("Zero/nonfinite training scale")
  data <- x
  for(v in vars) {p<-filter(params,variable==v);data[[paste0(v,"_z")]]<-(x[[v]]-p$mean)/p$SD}
  list(data=data,params=params)
}
segment_support <- function(supported) cumsum(c(TRUE,diff(as.integer(supported))!=0))

auc_fast <- function(score, y) {
  y <- strict_binary(y); ok <- is.finite(score) & !is.na(y); score <- score[ok]; y <- y[ok]
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (n1 < 1L || n0 < 1L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

bacc_fast <- function(pred, y) {
  ok <- !is.na(pred) & !is.na(y); pred <- pred[ok]; y <- y[ok]
  if (!any(y == 1L) || !any(y == 0L)) return(NA_real_)
  (mean(pred[y == 1L] == 1L) + mean(pred[y == 0L] == 0L)) / 2
}

source_files <- vapply(sys.frames(),function(e)if(is.null(e$ofile))""else e$ofile,character(1))
source(file.path(dirname(tail(source_files[nzchar(source_files)],1)),"windows.R"),local=TRUE)

cohort_composition <- function(x) {
  lev <- c("PRONIA", "NAPLS-3", "ZInEP", "FePsy", "MUC-FRUE")
  tab <- table(factor(x, levels = lev))
  paste(paste(names(tab), as.integer(tab), sep = "="), collapse = ";")
}

bootstrap_auc <- function(score, y, analysis_id, B = 2000L) {
  ok <- is.finite(score)&!is.na(y); score<-score[ok];y<-strict_binary(y[ok])
  pos <- which(y == 1L); neg <- which(y == 0L)
  if(!length(pos)||!length(neg)) return(c(lo=NA,hi=NA,valid=0,invalid=B))
  values <- deterministic_rows(B, analysis_id, function(i) {
    ii <- c(sample_indices(pos,TRUE), sample_indices(neg,TRUE))
    auc_fast(score[ii], y[ii])
  }, cores = 1L)[, 1]
  c(lo = unname(quantile(values, .025, na.rm = TRUE)), hi = unname(quantile(values, .975, na.rm = TRUE)),
    valid = sum(is.finite(values)), invalid = sum(!is.finite(values)))
}

window_stats_1d <- function(x, variable, width_fraction = .20, step_fraction = .05,
                            min_events = 7L, min_non_events = 7L, bootstrap_B = 0L,
                            analysis_id = "window") {
  x <- x[eligible(x,variable),,drop=FALSE]
  grid <- window_grid(x[[variable]], width_fraction, step_fraction)
  rows <- map_dfr(seq_len(nrow(grid)), function(i) {
    ii <- which(in_window(x[[variable]], grid, i)); z <- x[ii, ]; events <- sum(z$EXP_LABEL == 1L); non_events <- nrow(z) - events
    supported <- events >= max(1,min_events) && non_events >= max(1,min_non_events)
    auc <- if (supported) auc_fast(z$Mean_Score, z$EXP_LABEL) else NA_real_
    boot <- if (supported && bootstrap_B > 0L) bootstrap_auc(z$Mean_Score, z$EXP_LABEL, paste0(analysis_id, ".", i), bootstrap_B) else c(lo = NA, hi = NA, valid = 0, invalid = 0)
    tibble(window_id = i, start = grid$start[[i]], end = grid$end[[i]], center = grid$center[[i]],
           right_closed = grid$right_closed[[i]], N = nrow(z), events = events, non_events = non_events,
           supported = supported, cohort_composition = cohort_composition(z$cohort), moderator_mean = if(nrow(z)) mean(z[[variable]]) else NA_real_,
           AUC = auc, AUC_bootstrap_lo = boot[["lo"]], AUC_bootstrap_hi = boot[["hi"]],
           bootstrap_valid = as.integer(boot[["valid"]]), bootstrap_invalid = as.integer(boot[["invalid"]]))
  })
  if(!nrow(rows)) return(list(rows=rows,Dmax=NA_real_,Dmean=NA_real_,rho=NA_real_))
  good <- rows %>% filter(supported, is.finite(AUC))
  list(rows = rows, Dmax = if(nrow(good) >= 3) diff(range(good$AUC)) else NA_real_,
       Dmean = if(nrow(good) >= 3) mean(abs(good$AUC - mean(good$AUC))) else NA_real_,
       rho = if(nrow(good) >= 3) safe_cor(good$AUC,good$moderator_mean) else NA_real_)
}

window_dmax <- function(x, variable, moderator_values, grid, min_events = 7L, min_non_events = 7L) {
  aucs <- numeric(0)
  for (i in seq_len(nrow(grid))) {
    ii <- which(in_window(moderator_values, grid, i)); y <- x$EXP_LABEL[ii]
    if (sum(y == 1L) >= min_events && sum(y == 0L) >= min_non_events) aucs <- c(aucs, auc_fast(x$Mean_Score[ii], y))
  }
  if (length(aucs) >= 3L) diff(range(aucs)) else NA_real_
}

surface_stats <- function(x, variable_1, variable_2, width_fraction = .20, step_fraction = .05,
                          min_events = 7L, min_non_events = 7L, bootstrap_B = 0L,
                          analysis_id = "surface") {
  x <- x[eligible(x,c(variable_1,variable_2)),,drop=FALSE]
  g1 <- window_grid(x[[variable_1]], width_fraction, step_fraction)
  g2 <- window_grid(x[[variable_2]], width_fraction, step_fraction)
  rows <- map_dfr(seq_len(nrow(g1)), function(i) map_dfr(seq_len(nrow(g2)), function(j) {
    ii <- which(in_window(x[[variable_1]], g1, i) & in_window(x[[variable_2]], g2, j)); z <- x[ii, ]
    events <- sum(z$EXP_LABEL == 1L); non_events <- nrow(z) - events; supported <- events >= max(1,min_events) && non_events >= max(1,min_non_events)
    auc <- if (supported) auc_fast(z$Mean_Score, z$EXP_LABEL) else NA_real_
    boot <- if (supported && bootstrap_B > 0L) bootstrap_auc(z$Mean_Score, z$EXP_LABEL, paste0(analysis_id, ".", i, ".", j), bootstrap_B) else c(lo = NA, hi = NA, valid = 0, invalid = 0)
    tibble(cell_i = i, cell_j = j, start_1 = g1$start[[i]], end_1 = g1$end[[i]], center_1 = g1$center[[i]], right_closed_1 = g1$right_closed[[i]],
           start_2 = g2$start[[j]], end_2 = g2$end[[j]], center_2 = g2$center[[j]], right_closed_2 = g2$right_closed[[j]],
           N = nrow(z), events = events, non_events = non_events, supported = supported,
           cohort_composition = cohort_composition(z$cohort), AUC = auc,
           AUC_bootstrap_lo = boot[["lo"]], AUC_bootstrap_hi = boot[["hi"]],
           bootstrap_valid = as.integer(boot[["valid"]]), bootstrap_invalid = as.integer(boot[["invalid"]]))
  }))
  if(!nrow(rows)) return(list(rows=rows,Dmax=NA_real_,Dmean=NA_real_,rho=NA_real_))
  good <- rows %>% filter(supported, is.finite(AUC))
  list(rows = rows, Dmax = if(nrow(good) >= 3) diff(range(good$AUC)) else NA_real_,
       Dmean = if(nrow(good) >= 3) mean(abs(good$AUC - mean(good$AUC))) else NA_real_, g1 = g1, g2 = g2)
}

surface_dmax <- function(x, values_1, values_2, g1, g2, min_events = 7L, min_non_events = 7L) {
  aucs <- numeric(0)
  for (i in seq_len(nrow(g1))) for (j in seq_len(nrow(g2))) {
    ii <- which(in_window(values_1, g1, i) & in_window(values_2, g2, j)); y <- x$EXP_LABEL[ii]
    if (sum(y == 1L) >= min_events && sum(y == 0L) >= min_non_events) aucs <- c(aucs, auc_fast(x$Mean_Score[ii], y))
  }
  if (length(aucs) >= 3L) diff(range(aucs)) else NA_real_
}

permute_within <- function(values, strata) {
  out <- values
  for (s in unique(strata)) { ii <- which(strata == s); out[ii] <- values[sample_indices(ii)] }
  out
}

null_summary <- function(values, observed, analysis_id, statistic, tail, reference = NA_real_) {
  values <- as.numeric(values)
  finite_mask <- is.finite(values)
  finite_values <- values[finite_mask]
  attempted_count <- length(values)
  valid_count <- length(finite_values)
  invalid_count <- attempted_count - valid_count
  observed_defined <- length(observed) == 1L && is.finite(observed)
  reference_defined <- tail != "absolute_reference" || (length(reference) == 1L && is.finite(reference))

  if (valid_count == 0L) {
    exceedance_count <- NA_integer_
    p_value <- NA_real_
    null_mean_value <- null_sd_value <- null_q025_value <- null_q50_value <- null_q975_value <- NA_real_
    estimability_status <- "NOT ESTIMABLE: no finite null draws"
  } else {
    null_mean_value <- mean(finite_values)
    null_sd_value <- if (valid_count >= 2L) stats::sd(finite_values) else NA_real_
    null_quantiles <- stats::quantile(finite_values, c(.025, .5, .975), names = FALSE)
    null_q025_value <- null_quantiles[[1L]]
    null_q50_value <- null_quantiles[[2L]]
    null_q975_value <- null_quantiles[[3L]]
    if (!observed_defined || !reference_defined) {
      exceedance_count <- NA_integer_
      p_value <- NA_real_
      estimability_status <- "NOT ESTIMABLE: observed statistic or required reference undefined"
    } else {
      exceedance_count <- if (tail == "upper") {
        sum(finite_values >= observed)
      } else if (tail == "absolute_zero") {
        sum(abs(finite_values) >= abs(observed))
      } else if (tail == "absolute_reference") {
        sum(abs(finite_values - reference) >= abs(observed - reference))
      } else stop("unknown tail")
      p_value <- (1 + exceedance_count) / (1 + valid_count)
      estimability_status <- "ESTIMABLE"
    }
  }

  tibble(analysis_id = analysis_id, statistic = statistic, tail = tail, reference = reference,
         observed = observed, attempted = attempted_count, valid = valid_count, invalid = invalid_count,
         exceedances = exceedance_count, p_add_one = p_value,
         validity_policy="conditional_on_finite_draws", method_review_required=invalid_count>0,
         null_mean = null_mean_value, null_sd = null_sd_value,
         null_q025 = null_q025_value, null_q50 = null_q50_value, null_q975 = null_q975_value,
         estimability_status = estimability_status)
}

component_stats <- function(v, y) {
  pos <- v[y == 1L]; neg <- v[y == 0L]
  c(mad_all = mad(v, constant = 1), mad_neg = mad(neg, constant = 1), mad_pos = mad(pos, constant = 1),
    mean_pos = mean(pos), mean_neg = mean(neg), delta_mean = mean(pos) - mean(neg))
}

component_1d <- function(x, variable, moderator_values, grid, comp_cols = paste0("Comp", 1:4)) {
  keep <- eligible(x,comp_cols)&is.finite(moderator_values); x<-x[keep,];moderator_values<-moderator_values[keep]
  rows <- map_dfr(seq_len(nrow(grid)), function(i) {
    ii <- which(in_window(moderator_values, grid, i)); y <- x$EXP_LABEL[ii]; events <- sum(y == 1L); non_events <- sum(y == 0L)
    base <- tibble(window_id = i, start = grid$start[[i]], end = grid$end[[i]], center = grid$center[[i]], right_closed = grid$right_closed[[i]],
                   N = length(ii), events = events, non_events = non_events, moderator_mean = if(length(ii)) mean(moderator_values[ii]) else NA_real_,
                   AUC = if(events >= 7L && non_events >= 7L) auc_fast(x$Mean_Score[ii], y) else NA_real_)
    metric_names <- unlist(lapply(tolower(comp_cols),function(cc)paste0(cc,"_",c("mad_all","mad_neg","mad_pos","mean_pos","mean_neg","delta_mean"))))
    if (events < 7L || non_events < 7L) return(bind_cols(base,as_tibble_row(setNames(rep(NA_real_,24),metric_names))))
    vals <- unlist(map(comp_cols, function(cc) setNames(component_stats(x[[cc]][ii], y), paste0(tolower(cc), "_", names(component_stats(x[[cc]][ii], y))))))
    bind_cols(base, as_tibble_row(vals))
  })
  good <- rows %>% filter(is.finite(AUC)); metrics <- grep("^comp[1-4]_", names(good), value = TRUE)
  rho <- vapply(metrics, function(m) suppressWarnings(safe_cor(good$AUC,good[[m]])), numeric(1))
  list(rows = rows, rho = rho)
}

component_2d <- function(x, values_1, values_2, g1, g2, comp_cols = paste0("Comp", 1:4)) {
  keep <- eligible(x,comp_cols)&is.finite(values_1)&is.finite(values_2);x<-x[keep,];values_1<-values_1[keep];values_2<-values_2[keep]
  metric_names <- unlist(lapply(tolower(comp_cols), function(cc) paste0(cc, "_", c("mad_all","mad_neg","mad_pos","mean_pos","mean_neg","delta_mean"))))
  records <- vector("list", nrow(g1) * nrow(g2)); k <- 0L
  for (i in seq_len(nrow(g1))) for (j in seq_len(nrow(g2))) {
    k <- k + 1L; ii <- which(in_window(values_1, g1, i) & in_window(values_2, g2, j)); y <- x$EXP_LABEL[ii]
    events <- sum(y == 1L); non_events <- sum(y == 0L)
    vals <- rep(NA_real_, length(metric_names)); names(vals) <- metric_names
    auc <- NA_real_
    if (events >= 7L && non_events >= 7L) {
      auc <- auc_fast(x$Mean_Score[ii], y)
      vals <- unlist(lapply(comp_cols, function(cc) component_stats(x[[cc]][ii], y)), use.names = FALSE); names(vals) <- metric_names
    }
    records[[k]] <- bind_cols(tibble(cell_i=i,cell_j=j,center_1=g1$center[[i]],center_2=g2$center[[j]],N=length(ii),events=events,non_events=non_events,AUC=auc), as_tibble_row(vals))
  }
  rows <- bind_rows(records); good <- rows %>% filter(is.finite(AUC))
  rho <- vapply(metric_names, function(m) suppressWarnings(safe_cor(good$AUC,good[[m]])), numeric(1))
  list(rows = rows, rho = rho)
}

# Allocation-light variants used inside permutation loops. They return the
# same 24 correlations without constructing per-window/cell tibbles.
component_rho_1d_fast <- function(x, moderator_values, grid, comp_cols = paste0("Comp", 1:4)) {
  keep<-eligible(x,comp_cols)&is.finite(moderator_values);x<-x[keep,];moderator_values<-moderator_values[keep]
  values <- matrix(NA_real_, nrow(grid), 25L)
  for (i in seq_len(nrow(grid))) {
    ii <- which(in_window(moderator_values, grid, i)); y <- x$EXP_LABEL[ii]
    if (sum(y == 1L) < 7L || sum(y == 0L) < 7L) next
    values[i, 1L] <- auc_fast(x$Mean_Score[ii], y)
    values[i, -1L] <- unlist(lapply(comp_cols, function(cc) component_stats(x[[cc]][ii], y)), use.names = FALSE)
  }
  good <- is.finite(values[, 1L])
  out <- vapply(2:ncol(values), function(j) suppressWarnings(safe_cor(values[good,1L],values[good,j])), numeric(1))
  names(out) <- unlist(lapply(tolower(comp_cols), function(cc) paste0(cc, "_", c("mad_all","mad_neg","mad_pos","mean_pos","mean_neg","delta_mean"))))
  out
}

component_rho_2d_fast <- function(x, values_1, values_2, g1, g2, comp_cols = paste0("Comp", 1:4)) {
  keep<-eligible(x,comp_cols)&is.finite(values_1)&is.finite(values_2);x<-x[keep,];values_1<-values_1[keep];values_2<-values_2[keep]
  values <- matrix(NA_real_, nrow(g1) * nrow(g2), 25L); k <- 0L
  for (i in seq_len(nrow(g1))) for (j in seq_len(nrow(g2))) {
    k <- k + 1L; ii <- which(in_window(values_1, g1, i) & in_window(values_2, g2, j)); y <- x$EXP_LABEL[ii]
    if (sum(y == 1L) < 7L || sum(y == 0L) < 7L) next
    values[k, 1L] <- auc_fast(x$Mean_Score[ii], y)
    values[k, -1L] <- unlist(lapply(comp_cols, function(cc) component_stats(x[[cc]][ii], y)), use.names = FALSE)
  }
  good <- is.finite(values[, 1L])
  out <- vapply(2:ncol(values), function(j) suppressWarnings(safe_cor(values[good,1L],values[good,j])), numeric(1))
  names(out) <- unlist(lapply(tolower(comp_cols), function(cc) paste0(cc, "_", c("mad_all","mad_neg","mad_pos","mean_pos","mean_neg","delta_mean"))))
  out
}

# Resolve a new output path through existing ancestors, including symlinks.
external_output_path <- function(path, root) {
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(getwd(), path)
  suffix <- character()
  ancestor <- path
  while (!file.exists(ancestor)) {
    suffix <- c(basename(ancestor), suffix)
    parent <- dirname(ancestor)
    if (identical(parent, ancestor)) stop("Cannot resolve output parent")
    ancestor <- parent
  }
  resolved <- normalizePath(ancestor, winslash="/", mustWork=TRUE)
  for (part in suffix) {
    if (part == "..") resolved <- dirname(resolved)
    else if (part != ".") resolved <- file.path(resolved, part)
    if (file.exists(resolved)) resolved <- normalizePath(resolved, winslash="/", mustWork=TRUE)
  }
  root <- normalizePath(root, winslash="/", mustWork=TRUE)
  if (identical(resolved, root) || startsWith(resolved, paste0(root,"/")))
    stop("Output must be outside the code checkout")
  resolved
}
