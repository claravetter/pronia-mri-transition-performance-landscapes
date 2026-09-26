window_grid <- function(x, width_fraction = .20, step_fraction = .05) {
  if(!is.finite(width_fraction)||width_fraction<=0||width_fraction>1||!is.finite(step_fraction)||step_fraction<=0||step_fraction>width_fraction) stop("Invalid window width/step")
  x <- x[is.finite(x)]
  if(!length(x)) return(tibble(window_id=integer(),start=double(),end=double(),center=double(),right_closed=logical(),width_fraction=double(),step_fraction=double()))
  lo <- min(x); hi <- max(x); rg <- hi-lo
  if(rg==0) return(tibble(window_id=1L,start=lo,end=hi,center=lo,right_closed=TRUE,width_fraction=width_fraction,step_fraction=step_fraction))
  width <- rg * width_fraction; step <- rg * step_fraction; target <- hi - width
  starts <- seq(lo, target, by = step)
  tolerance <- max(1, abs(hi), abs(lo)) * 1e-12
  if (abs(tail(starts, 1) - target) > tolerance) starts <- c(starts, target) else starts[length(starts)] <- target
  starts <- starts[c(TRUE, diff(starts) > tolerance)]
  out <- tibble(window_id = seq_along(starts), start = starts, end = pmin(starts + width, hi),
                center = (start + end) / 2, right_closed = seq_along(starts) == length(starts),
                width_fraction = width_fraction, step_fraction = step_fraction)
  stopifnot(abs(tail(out$end, 1) - hi) <= tolerance, isTRUE(tail(out$right_closed, 1)))
  out
}

in_window <- function(x, grid, i) {
  if (grid$right_closed[[i]]) is.finite(x) & x >= grid$start[[i]] & x <= grid$end[[i]]
  else is.finite(x) & x >= grid$start[[i]] & x < grid$end[[i]]
}

