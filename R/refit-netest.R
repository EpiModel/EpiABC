
#' @title Refit TERGM with Biased Edges Dissolution Approximation
#'
#' @description Refits a TERGM estimated with \code{netest} using the edges
#'              dissolution approximation that had residual biased estimates.
#'
#' @param est An object of class \code{netest}, from the
#'        \code{\link{netest}} function.
#' @param nsims Number of simulations per sequential ABC wave (more = fewer
#'              waves required for convergence).
#' @param ncores Number of cores for each wave of simulations.
#' @param nsteps Number of time steps to simulate the TERGM for each refitting.
#' @param coefs.to.fit Integer vector of coefficient position to refit. Implicit
#'                     default is the full coefficient vector, but use this to
#'                     subset fitting to selected coefficients.
#' @param prior.min Absolute lower bound of adjustment to be made across all
#'                  model coefficients.
#' @param prior.max Absolute upper bound of adjustment to be made across all
#'                  model coefficients.
#' @param p_acc_min Convergence threshold for ABC waves, as a proportion of
#'                  simulations accepted (range from 0 to 1, with lower =
#'                  greater precision).
#'
#' @details
#' Fitting an TERGM in \code{\link{netest}} can be accelerated with the
#' edges-dissolution approximation method (\code{edapprox = TRUE}), but this
#' sometimes results in bias (poor fit to the target statistics) when the mean
#' network degree is high, mean duration is short, or other combinations. All of
#' this poor fit is revealed with model diagnostics in \code{\link{netdx}}. If
#' the approximation method does not work, then the immediate next step is to
#' try fitting a full TERGM without the approximation (\code{edapprox = FALSE}).
#' This is often successful must has challenges with larger and more complex
#' network models, particularly those with longer durations. Therefore, this
#' function provides an alternative approach using approximate Bayesian
#' computation (ABC) methods. This function simply wraps a commonly used
#' sequential ABC approach and does the additional coefficient adjustment
#' automatically. It may not work for all models, and there is no guarantee of
#' model convergence or even improvement from the original fit; one should
#' evaluate that fit with an additional round of model diagnostics after the
#' ABC algorithm completes.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Original fit with edges dissolution approximation: High mean degree and
#' # short duration often leads to poor performance of the approximation method
#' nw <- network_initialize(1000)
#' formation <- ~edges + concurrent
#' target.stats <- c(1*(1000/2), 250)
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
#' est <- netest(nw, formation, target.stats, coef.diss)
#'
#' # Diagnostics will show bias here
#' dx <- netdx(est, nsims = 10, nsteps = 300, ncores = 4)
#' print(dx)
#'
#' # Refit the model with Bayesian ABC, with uniform prior adjustments applying
#' # to all model coefficients
#' est_new <- netest_refit_abc(est, nsims = 10, ncores = 5, nsteps = 300,
#'                             prior.min = -0.5, prior.max = 0)
#'
#' # Raw output from EasyABC
#' est_new$refit
#'
#' # Compare old and new coefficient values
#' est$coef.form
#' est_new$coef.form
#'
#' # Rerun diagnostics to evaluate fit (no guarantee that ABC works)
#' dx_new <- netdx(est_new, nsims = 10, nsteps = 300, ncores = 4)
#' dx_new
#'}
#'
netest_refit_abc <- function(est, nsims, ncores, nsteps,
                             coefs.to.fit,
                             prior.min = 0, prior.max = 0,
                             p_acc_min = 0.1) {

  est_orig <- est
  if (missing(nsims)) {
    stop("Specify nsims for sequential ABC wave")
  }
  if (missing(nsteps)) {
    stop("Specify nsteps for netdx resimulation", call. = FALSE)
  }
  if (missing(ncores)) {
    ncores <- 1
  }
  if (ncores > parallel::detectCores()) {
    ncores <- parallel::detectCores()
  }
  if (missing(coefs.to.fit)) {
    coefs.to.fit <- seq_along(est_orig$coef.form)
  }

  save(est_orig, nsteps, coefs.to.fit, file = "temp-refit-abc.rda")

  myfunc <- function(x) {
    set.seed(x[1])
    # require(EpiModel)
    load("temp-refit-abc.rda")
    est_temp <- est_orig
    est_temp$coef.form[coefs.to.fit] <- est_temp$coef.form[coefs.to.fit] +
      x[2:(length(coefs.to.fit) + 1)]
    dx <- EpiModel::netdx(est_temp, nsims = 1, nsteps = nsteps, verbose = FALSE)
    out <- EpiModel::get_nwstats(dx)
    out <- out[, which(!names(out) %in% c("time", "sim")), drop = FALSE]
    out <- colMeans(out)[coefs.to.fit]
    return(out)
  }

  targets <- est_orig$target.stats[coefs.to.fit]
  n_targets <- length(targets)
  priors <- list()
  for (ii in seq_len(n_targets)) {
    priors[[ii]] <- c("unif", prior.min, prior.max)
  }

  refit <- EasyABC::ABC_sequential(
    method = "Lenormand",
    model = myfunc,
    prior = priors,
    nb_simul = nsims,
    summary_stat_target = targets,
    p_acc_min = p_acc_min,
    progress_bar = TRUE,
    verbose = FALSE,
    n_cluster = ncores,
    use_seed = TRUE
  )

  if (n_targets == 1) {
    coef.adj <- sum(refit$param * refit$weights)
  } else {
    coef.adj <- rep(NA, n_targets)
    for (jj in seq_len(n_targets)) {
      coef.adj[jj] <- sum(refit$param[, jj] * refit$weights)
    }
  }

  est_new <- est_orig
  est_new$coef.form[coefs.to.fit] <- est_new$coef.form[coefs.to.fit] + coef.adj
  est_new$refit <- refit

  unlink("temp-refit-abc.rda")

  return(est_new)
}
