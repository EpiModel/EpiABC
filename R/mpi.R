
#' Approximate Bayesian Computation with Sequential Monte Carlo Sampling on
#' High-Performance Computing Clusters
#'
#' This function implements the Lenormand algorithm for ABC-SMC designed from the
#' \code{EasyABC} R package on arbitrarily user-defined clusters.
#'
#' @param model a \code{R} function implementing the model to be simulated. It
#'        must take as arguments a vector of model parameter values and it
#'        must return a vector of summary statistics. When using the option
#'       \code{use_seed=TRUE}, \code{model} must take as arguments a vector
#'       containing a seed value and the model parameter values. A tutorial
#'       is provided in the package's vignette to dynamically link a binary
#'       code to a \code{R} function. Users may alternatively wish to wrap
#'       their binary executables using the provided functions
#'       \code{\link{binary_model}} and \code{\link{binary_model_cluster}}.
#'       The use of these functions is associated with slightly different
#'       constraints on the design of the binary code (see
#'       \code{\link{binary_model}} and \code{\link{binary_model_cluster}}).
#' @param prior a list of prior information. Each element of the list
#'        corresponds to a model parameter. The list element must be a vector
#'        whose first argument determines the type of prior distribution:
#'        possible values are \code{"unif"} for a uniform distribution on
#'        a segment, \code{"normal"} for a normal distribution,
#'        \code{"lognormal"} for a lognormal distribution or \code{"exponential"}
#'        for an exponential distribution.
#' @param nb_simul a positive integer equal to the desired number of simulations
#'        of the model below the tolerance threshold when \code{method} is
#'        \code{"Beaumont"}, \code{"Drovandi"} and \code{"Delmoral"}.
#'        When \code{method} is \code{"Lenormand"}, the number of simulations
#'        below the tolerance threshold is equal to \code{nb_simul * alpha}.
#'        See the package's vignette and Lenormand et al. (2012) for details.
#' @param summary_stat_target a vector containing the targeted (observed)
#'        summary statistics.
#' @param prior_test a string expressing the constraints between model
#'        parameters. This expression will be evaluated as a logical expression,
#'        you can use all the logical operators including \code{"<"}, \code{">"},
#'        \ldots{}. Each parameter should be designated with \code{"X1"}, \code{"X2"},
#'        \ldots{} in the same order as in the prior definition. If not provided,
#'        no constraint will be applied.
#' @param n_cluster a positive integer. If larger than 1 (the default value),
#'        \code{ABC_sequential} will launch \code{model} simulations in parallel on
#'        \code{n_cluster} cores of the computer.
#' @param verbose If \code{TRUE}, \code{ABC_sequential} writes in the current
#'        directory intermediary results at the end of each step of the algorithm
#'        various files.
#' @param dist_weights a vector containing the weights to apply to the distance
#'        between the computed and the targeted statistics. These weights can
#'        be used to give more importance to a summary statistisc for example.
#'        The weights will be normalized before applying them. If not provided,
#'        no weights will be applied.
#' @param cl cluster object, such as that produced by \code{snow::makeCluster}.
#' @param \dots Additional arguments can be passed depending on the choosen
#'        method (see below).
#'
#' @return The returned value is a list containing the following components:
#'    \item{param}{The model parameters used in the \code{model} simulations.}
#'    \item{stats}{The summary statistics obtained at the end of the \code{model}
#'          simulations.}
#'    \item{weights}{The weights of the different \code{model} simulations.}
#'    \item{stats_normalization}{The standard deviation of the summary statistics
#'          across the \code{model} simulations of the initial step. These values
#'          are used to normalize the summary statistics before the computation
#'          of the Euclidean distance between simulations and data.}
#'     \item{epsilon}{The final maximal distance between simulations and data in
#'           the retained sample of particles.}
#'     \item{nsim}{The number of \code{model} simulations performed.}
#'     \item{computime}{The computing time to perform the simulations.}
#'     \item{intermediary}{Intermediary results stored when the option
#'           \code{"verbose=TRUE"} is chosen. Each element of this list
#'           corresponds to a different step. See the argument \code{verbose}
#'           above for more details on the information stored.}
#'
#' @export
#'
abc_smc_cluster <- function(model,
                            prior,
                            nb_simul,
                            summary_stat_target,
                            prior_test = NULL,
                            n_cluster = 1,
                            verbose = FALSE,
                            dist_weights = NULL,
                            cl,
                            ...) {

  use_seed <- TRUE

  ## checking errors in the inputs

  if (missing(model))
    stop("model is missing")
  if (missing(prior))
    stop("prior is missing")
  data <- .wrap_constants_in_model(prior, model, use_seed)
  prior <- data$new_prior
  model <- data$new_model
  prior <- .process_prior(prior)
  if (!is.null(prior_test))
    .check_prior_test(length(prior), prior_test)
  if (missing(nb_simul))
    stop("nb_simul is missing")
  if (missing(summary_stat_target))
    stop("summary_stat_target is missing")
  if (!is.vector(nb_simul) || length(nb_simul) > 1 || nb_simul < 1)
    stop("nb_simul must be a number larger than 1.")
  nb_simul <- floor(nb_simul)
  if (!is.vector(summary_stat_target))
    stop("'summary_stat_target' has to be a vector.")
  if (!is.vector(n_cluster) || length(n_cluster) > 1 || n_cluster < 1)
    stop("'n_cluster' has to be a positive number.")
  n_cluster <- floor(n_cluster)
  if (!is.logical(verbose))
    stop("verbose has to be boolean")
  if (!is.null(dist_weights) && length(dist_weights) != length(summary_stat_target)) {
    stop("'dist_weights' has to be the same length than 'summary_stat_target'")
  }
  sequential <- NULL
  if (n_cluster == 1) {
    stop("This version of ABC-SMC designed for multi-core, set n_cluster > 1")
  } else {
    sequential <- abc_lenormand_cluster(model = model, prior = prior, prior_test = prior_test,
                                        nb_simul = nb_simul, summary_stat_target = summary_stat_target,
                                        n_cluster = n_cluster, verbose = verbose,
                                        dist_weights = dist_weights, cl = cl, ...)
  }
  return(sequential)
}

abc_lenormand_cluster <- function(model,
                                  prior,
                                  prior_test,
                                  nb_simul,
                                  summary_stat_target,
                                  n_cluster,
                                  verbose,
                                  alpha = 0.5,
                                  p_acc_min = 0.05,
                                  dist_weights = NULL,
                                  cl,
                                  seed_count = 0,
                                  inside_prior = TRUE,
                                  progress_bar = TRUE,
                                  max_pick = 10000) {

  if (!is.vector(alpha) || length(alpha) > 1)
    stop("alpha has to be a number")
  if (alpha <= 0 | alpha >= 1)
    stop("alpha has to be between 0 and 1")
  if (!is.vector(p_acc_min) || length(p_acc_min) > 1)
    stop("p_acc_min has to be a number")
  if (p_acc_min <= 0 | p_acc_min >= 1)
    stop("p_acc_min has to be between 0 and 1")
  if (!is.vector(seed_count) || length(seed_count) > 1 || seed_count < 0)
    stop("seed_count has to be a positive number")
  seed_count <- floor(seed_count)
  if (!is.logical(inside_prior))
    stop("inside_prior has to be boolean")
  start <- Sys.time()
  if (progress_bar == TRUE) {
    print("## Lenormand et al. (2012)'s algorithm ##")
  }

  seed_count_ini <- seed_count
  nparam <- length(prior)
  nstat <- length(summary_stat_target)
  if (!.all_unif(prior)) {
    stop("Prior distributions must be uniform to use the Lenormand et al. (2012)'s algorithm.")
  }
  n_alpha <- ceiling(nb_simul * alpha)

  ## step 1 ABC rejection step with LHS
  tab_ini <- abc_rejection_lhs_cluster(model,
                                       prior,
                                       prior_test,
                                       nb_simul,
                                       seed_count,
                                       n_cluster,
                                       cl)
  # initially, weights are equal
  tab_weight <- array(1, n_alpha)
  seed_count <- seed_count + nb_simul
  # determination of the normalization constants in each dimension associated to
  # each summary statistic, this normalization will not change during all the
  # algorithm
  sd_simul <- sapply(as.data.frame(tab_ini[, (nparam + 1):(nparam + nstat)]), sd,
                    na.rm = TRUE)
  # selection of the alpha quantile closest simulations
  simul_below_tol <- NULL
  simul_below_tol <- rbind(simul_below_tol,
                          .selec_simul_alpha(summary_stat_target,
                                             tab_ini[, 1:nparam],
                                             tab_ini[, (nparam + 1):(nparam + nstat)], sd_simul,
                                             alpha, dist_weights = dist_weights))
  # to be sure that there are not two or more simulations at a distance equal
  #   to the tolerance determined by the quantile
  simul_below_tol <- simul_below_tol[1:n_alpha, ]
  tab_dist <- .compute_dist(summary_stat_target,
                            as.matrix(as.matrix(simul_below_tol)[, (nparam + 1):(nparam + nstat)]),
                            sd_simul, dist_weights = dist_weights)
  if (!is.null(dist_weights)) {
    tab_dist <- tab_dist * (dist_weights/sum(dist_weights))
  }
  tol_next <- max(tab_dist)
  if (progress_bar) {
    print("Step 1 Completed ...")
  }

  ## following steps
  p_acc <- p_acc_min + 1
  nb_simul_step <- nb_simul - n_alpha
  it <- 1

  while (p_acc > p_acc_min) {
    it <- it + 1
    simul_below_tol2 <- NULL
    param_previous_step <- as.matrix(as.matrix(simul_below_tol)[, 1:nparam])
    tab_inic <- abc_launcher_not_uniformc_cluster(model = model,
                                                  prior = prior,
                                                  param_previous_step = param_previous_step,
                                                  tab_weight = tab_weight/sum(tab_weight),
                                                  nb_simul = nb_simul_step,
                                                  seed_count = seed_count,
                                                  inside_prior = inside_prior,
                                                  n_cluster = n_cluster,
                                                  cl = cl,
                                                  max_pick = max_pick)
    tab_ini <- as.matrix(tab_inic[[1]])
    tab_ini <- as.numeric(tab_ini)
    dim(tab_ini) <- c(nb_simul_step, (nparam + nstat))
    seed_count <- seed_count + nb_simul_step
    if (!inside_prior) {
      tab_weight2 <- .compute_weightb(as.matrix(as.matrix(as.matrix(tab_ini)[, 1:nparam])),
                                      as.matrix(as.matrix(as.matrix(simul_below_tol)[, 1:nparam])),
                                      tab_weight/sum(tab_weight), prior)
    } else {
      tab_weight2 <- tab_inic[[2]] * (.compute_weightb(as.matrix(as.matrix(as.matrix(tab_ini)[,1:nparam])),
                                                       as.matrix(as.matrix(as.matrix(simul_below_tol)[, 1:nparam])),
                                                       tab_weight/sum(tab_weight), prior))
    }
    simul_below_tol2 <- rbind(as.matrix(simul_below_tol), as.matrix(tab_ini))
    tab_weight <- c(tab_weight, tab_weight2)
    tab_dist2 <- .compute_dist(summary_stat_target,
                               as.matrix(as.matrix(tab_ini)[, (nparam + 1):(nparam + nstat)]),
                               sd_simul, dist_weights = dist_weights)
    if (!is.null(dist_weights)) {
      tab_dist2 <- tab_dist2 * (dist_weights/sum(dist_weights))
    }
    p_acc <- length(tab_dist2[!is.na(tab_dist2) & tab_dist2 <= tol_next])/nb_simul_step
    tab_dist <- c(tab_dist, tab_dist2)
    tol_next <- sort(tab_dist)[n_alpha]
    simul_below_tol2 <- simul_below_tol2[!is.na(tab_dist) & tab_dist <= tol_next, ]
    tab_weight <- tab_weight[!is.na(tab_dist) & tab_dist <= tol_next]
    tab_weight <- tab_weight[1:n_alpha]
    tab_dist <- tab_dist[!is.na(tab_dist) & tab_dist <= tol_next]
    odist <- order(tab_dist, decreasing = FALSE)[1:n_alpha]
    tab_dist_new <- tab_dist
    simul_below_tol <- matrix(0, n_alpha, (nparam + nstat))
    for (i1 in 1:n_alpha) {
      tab_dist_new[i1] <- tab_dist[odist[i1]]
      for (i2 in 1:(nparam + nstat)) {
        simul_below_tol[i1, i2] <- as.numeric(simul_below_tol2[odist[i1], i2])
      }
    }
    tab_dist <- tab_dist_new[1:n_alpha]
    if (progress_bar) {
      print(paste("Step ", it, " Completed: p_acc = ", p_acc, sep = ""))
    }
  }

  final_res <- NULL
  final_res <- list(param = as.matrix(as.matrix(simul_below_tol)[, 1:nparam]),
                   stats = as.matrix(as.matrix(simul_below_tol)[, (nparam + 1):(nparam + nstat)]),
                   weights = tab_weight/sum(tab_weight), stats_normalization = as.numeric(sd_simul),
                   epsilon = max(.compute_dist(summary_stat_target,
                                               as.matrix(as.matrix(simul_below_tol)[, (nparam + 1):(nparam + nstat)]),
                                               sd_simul, dist_weights = dist_weights)),
                   nsim = (seed_count - seed_count_ini), computime = as.numeric(difftime(Sys.time(),
                                                                                         start, units = "secs")))

  return(final_res)
}

abc_rejection_lhs_cluster <- function(model,
                                      prior,
                                      prior_test,
                                      nb_simul,
                                      seed_count,
                                      n_cluster,
                                      cl) {

  tab_simul_summarystat <- NULL
  tab_param <- NULL
  list_param <- list(NULL)
  # n_end <- nb_simul
  nparam <- length(prior)
  l <- length(prior)
  random_tab <- NULL
  all_unif_prior <- .all_unif(prior)
  if (all_unif_prior) {
    random_tab <- randomLHS(nb_simul, nparam)
  }

  list_param <- list(NULL)
  for (i in 1:nb_simul) {
    param <- array(0, l)
    if (!all_unif_prior) {
      param <- .sample_prior(prior, prior_test)
    } else {
      for (j in 1:l) {
        param[j] <- as.numeric(prior[[j]]$sampleArgs[2]) +
          (as.numeric(prior[[j]]$sampleArgs[3]) - as.numeric(prior[[j]]$sampleArgs[2])) *
          random_tab[i, j]
      }
    }
    param <- c((seed_count + i), param)
    list_param[[i]] <- param
    tab_param <- rbind(tab_param, param[2:(l + 1)])
  }
  seed_count <- seed_count + nb_simul
  list_simul_summarystat <- parLapplyLB(cl, list_param, model)
  for (i in 1:nb_simul) {
    tab_simul_summarystat <- rbind(tab_simul_summarystat, as.numeric(list_simul_summarystat[[i]]))
  }

  options(scipen = 0)
  out <- cbind(tab_param, tab_simul_summarystat)
  return(out)
}

abc_launcher_not_uniformc_cluster <- function(model,
                                              prior,
                                              param_previous_step,
                                              tab_weight,
                                              nb_simul,
                                              seed_count,
                                              inside_prior,
                                              n_cluster,
                                              cl,
                                              max_pick = 10000) {
  tab_simul_summarystat <- NULL
  tab_param <- NULL
  k_acc <- 0
  list_param <- list(NULL)

  list_param <- list(NULL)
  for (i in 1:nb_simul) {
    l <- dim(param_previous_step)[2]
    counter <- 0
    repeat {
      k_acc <- k_acc + 1
      counter <- counter + 1
      # pick a particle
      param_picked <- .particle_pick(param_previous_step, tab_weight)
      # move it
      # only variable parameters are moved, computation of a WEIGHTED variance
      param_moved <- .move_particle(as.numeric(param_picked),
                                   2 * cov.wt(as.matrix(as.matrix(param_previous_step)),
                                              as.vector(tab_weight))$cov)
      if ((!inside_prior) || (.is_included(param_moved, prior)) || (counter >= max_pick)) {
        break
      }
    }
    if (counter == max_pick) {
      stop("The proposal jumps outside of the prior distribution too often -
                   consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
    }
    param <- param_previous_step[1, ]
    param <- param_moved
    param <- c((seed_count + i), param)
    list_param[[i]] <- param
    tab_param <- rbind(tab_param, param[2:(l + 1)])
  }
  seed_count <- seed_count + nb_simul
  list_simul_summarystat <- parLapplyLB(cl, list_param, model)
  for (i in 1:nb_simul) {
    tab_simul_summarystat <- rbind(tab_simul_summarystat, as.numeric(list_simul_summarystat[[i]]))
  }

  out <- list(cbind(tab_param, tab_simul_summarystat), nb_simul/k_acc)
  return(out)
}
