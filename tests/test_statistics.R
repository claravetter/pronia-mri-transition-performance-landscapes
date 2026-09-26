#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)[1]
root <- normalizePath(file.path(dirname(sub("^--file=", "", arg)), ".."))
source(file.path(root,"R","lib","common.R"))
direct_summary <- function(values, observed, tail, reference = NA_real_) {
  finite_values <- as.numeric(values)[is.finite(values)]
  n <- length(finite_values)
  if (!n) return(c(attempted=length(values), valid=0, invalid=length(values), exceedances=NA,
                   p_add_one=NA, null_mean=NA, null_sd=NA, null_q025=NA, null_q50=NA, null_q975=NA))
  exceed <- if (tail == "upper") sum(finite_values >= observed) else
    if (tail == "absolute_zero") sum(abs(finite_values) >= abs(observed)) else
      sum(abs(finite_values-reference) >= abs(observed-reference))
  c(attempted=length(values), valid=n, invalid=length(values)-n, exceedances=exceed,
    p_add_one=(1+exceed)/(1+n), null_mean=base::mean(finite_values),
    null_sd=if(n >= 2) stats::sd(finite_values) else NA_real_,
    null_q025=stats::quantile(finite_values,.025,names=FALSE),
    null_q50=stats::quantile(finite_values,.5,names=FALSE),
    null_q975=stats::quantile(finite_values,.975,names=FALSE))
}

cases <- list(
  ordinary=list(values=c(-1,0,1,2,3), observed=1.5, tail="upper", reference=NA_real_),
  mixed_nonfinite=list(values=c(NA,-2,Inf,0,4,-Inf), observed=1, tail="absolute_zero", reference=0),
  constant_finite=list(values=rep(.25,8), observed=.25, tail="upper", reference=NA_real_)
)

for (case_name in names(cases)) {
  z <- cases[[case_name]]
  got <- null_summary(z$values,z$observed,case_name,"test",z$tail,z$reference)
  expected <- direct_summary(z$values,z$observed,z$tail,z$reference)
  for (field in names(expected)) {
    actual <- as.numeric(got[[field]])
    target <- unname(expected[[field]])
    stopifnot((is.na(actual) && is.na(target)) || isTRUE(all.equal(actual,target,tolerance=1e-14)))
  }
  stopifnot(got$attempted == got$valid + got$invalid)
}

constant <- null_summary(rep(7,5),7,"constant","test","upper")
stopifnot(identical(constant$null_sd,0), constant$estimability_status == "ESTIMABLE")
zero <- null_summary(c(NA,Inf,-Inf),1,"zero","test","upper")
stopifnot(is.na(zero$p_add_one), zero$estimability_status == "NOT ESTIMABLE: no finite null draws")
undefined <- null_summary(1:3,NA_real_,"undefined","test","upper")
stopifnot(is.na(undefined$p_add_one), is.finite(undefined$null_mean), grepl("NOT ESTIMABLE",undefined$estimability_status))
cat("PASS: null_summary unit tests\n")

for (w in c(.15,.20,.25,.30)) {
  g <- window_grid(c(-2.3, 11.7),w,w/4)
  stopifnot(abs(tail(g$end,1)-11.7)<1e-12, tail(g$right_closed,1),
            all(!head(g$right_closed,-1)), in_window(11.7,g,nrow(g)))
  stopifnot(!in_window(g$end[1],g,1))
}
stopifnot(auc_fast(c(0,1,2,3),c(0,0,1,1))==1,
          auc_fast(rep(1,4),c(0,0,1,1))==.5,
          is.na(auc_fast(1:4,rep(0,4))))
f <- function(i) runif(3)
a <- deterministic_rows(8,"unit_test_rng",f,cores=1)
b <- deterministic_rows(8,"unit_test_rng",f,cores=2)
stopifnot(identical(a,b),identical(a,deterministic_rows(8,"unit_test_rng",f)))
v <- c(2,3,7,12,19,21); y <- c(0,0,0,1,1,1)
st <- component_stats(v,y)
stopifnot(st['delta_mean']==mean(v[y==1])-mean(v[y==0]),st['mad_all']==mad(v,constant=1))
cat("PASS: window boundaries, rank AUC, component summaries and serial/multicore keyed RNG\n")
