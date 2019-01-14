
#' @export
abc_smc_prep <- function(model,
                         prior,
                         nb_simul,
                         summary_stat_target,
                         prior_test = NULL,
                         n_cluster = 2,
                         dist_weights = NULL,
                         alpha = 0.5,
                         p_acc_min = 0.1,
                         ...) {

  ## checking errors in the inputs
  if (missing(model))
    stop("'model' is missing")
  if (missing(prior))
    stop("'prior' is missing")
  data = .wrap_constants_in_model(prior, model, use_seed = TRUE)
  prior = data$new_prior
  model = data$new_model
  prior = .process_prior(prior)
  if (!is.null(prior_test))
    .check_prior_test(length(prior), prior_test)
  if (missing(nb_simul))
    stop("'nb_simul' is missing")
  if (missing(summary_stat_target))
    stop("'summary_stat_target' is missing")
  if (!is.vector(nb_simul))
    stop("'nb_simul' has to be a number.")
  if (length(nb_simul) > 1)
    stop("'nb_simul' has to be a number.")
  if (nb_simul < 1)
    stop("'nb_simul' must be a number larger than 1.")
  nb_simul = floor(nb_simul)
  if (!is.vector(summary_stat_target))
    stop("'summary_stat_target' has to be a vector.")
  if (!is.vector(n_cluster))
    stop("'n_cluster' has to be a number.")
  if (length(n_cluster) > 1)
    stop("'n_cluster' has to be a number.")
  if (n_cluster < 1)
    stop("'n_cluster' has to be a positive number.")
  n_cluster = floor(n_cluster)
  if (!is.null(dist_weights) && length(dist_weights) != length(summary_stat_target)) {
    stop("'dist_weights' has to be the same length than 'summary_stat_target'")
  }
  if (!is.vector(alpha))
    stop("'alpha' has to be a number.")
  if (length(alpha) > 1)
    stop("'alpha' has to be a number.")
  if (alpha <= 0)
    stop("'alpha' has to be between 0 and 1.")
  if (alpha >= 1)
    stop("'alpha' has to be between 0 and 1.")
  if (!is.vector(p_acc_min))
    stop("'p_acc_min' has to be a number.")
  if (length(p_acc_min) > 1)
    stop("'p_acc_min' has to be a number.")
  if (p_acc_min <= 0)
    stop("'p_acc_min' has to be between 0 and 1.")
  if (p_acc_min >= 1)
    stop("'p_acc_min' has to be between 0 and 1.")

  out <- list(model = model, prior = prior, prior_test = prior_test,
              nb_simul = nb_simul, summary_stat_target = summary_stat_target,
              n_cluster = n_cluster, dist_weights = dist_weights,
              alpha = alpha, p_acc_min = p_acc_min)

  return(out)
}



#' @export
abc_smc_wave <- function(input, wave, batch) {

  if (wave == 0) {

    # fixed params/settings
    model <- input$model
    prior <- input$prior
    prior_test <- input$prior_test
    nb_simul <- input$nb_simul
    summary_stat_target <- input$summary_stat_target
    n_cluster <- input$n_cluster
    # dist_weights <- input$dist_weights
    alpha <- input$alpha
    p_acc_min <- input$p_acc_min

    seed_count <- 0
    inside_prior <- input$inside_prior <- TRUE
    max_pick <- input$max_pick <- 10000

    nparam <- input$nparam <- length(prior)
    nstat <- input$nstat <- length(summary_stat_target)
    if (!.all_unif(prior)) {
      stop("Prior distributions must be uniform")
    }
    n_alpha <- input$n_alpha <- ceiling(nb_simul * alpha)

    tab_ini <- abc_wave0(model,
                         prior,
                         prior_test,
                         nb_simul,
                         seed_count,
                         n_cluster,
                         batch = batch)

    out <- list(init = input, seed_count = seed_count, tab_ini = tab_ini)
    return(out)
  }


  if (wave > 0) {

    # fixed
    p_acc_min <- input$init$p_acc_min
    n_alpha <- input$init$n_alpha
    nb_simul <- input$init$nb_simul
    model <- input$init$model
    prior <- input$init$prior
    nparam <- input$init$nparam
    inside_prior <- input$init$inside_prior
    max_pick <- input$init$max_pick
    nstat <- input$init$nstat
    n_cluster <- input$init$n_cluster
    nb_simul_step <- input$init$nb_simul_step <- nb_simul - n_alpha

    #prior wave
    simul_below_tol <- input$pwave$simul_below_tol
    tab_weight <- input$pwave$tab_weight
    seed_count <- input$pwave$seed_count

    tab_inic <- abc_waveN(model = model,
                          prior = prior,
                          param_previous_step = as.matrix(as.matrix(simul_below_tol)[, 1:nparam]),
                          tab_weight = tab_weight/sum(tab_weight),
                          nb_simul = nb_simul_step,
                          seed_count = seed_count,
                          inside_prior = inside_prior,
                          n_cluster = n_cluster,
                          max_pick = max_pick,
                          batch = batch)

    out <- list(init = input$init, pwave = input$pwave, tab_inic = tab_inic)
    return(out)
  }

}

#' @export
abc_smc_process <- function(input, wave) {

  if (wave == 0) {

    # fixed
    n_alpha <- input$init$n_alpha
    nparam <- input$init$nparam
    nstat <- input$init$nstat
    nb_simul <- input$init$nb_simul
    summary_stat_target <- input$init$summary_stat_target
    alpha <- input$init$alpha
    dist_weights <- input$init$dist_weights

    # current wave simulation
    seed_count <- input$seed_count
    tab_ini <- input$tab_ini

    # initially, weights are equal
    tab_weight <- array(1, n_alpha)
    seed_count <- seed_count + nb_simul
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

    out <- list(init = input$init,
                pwave = list(tab_weight = tab_weight, seed_count = seed_count,
                             simul_below_tol = simul_below_tol, tab_dist = tab_dist,
                             tol_next = tol_next, sd_simul = as.numeric(sd_simul)))
    return(out)
  }

  if (wave > 0) {

    # fixed
    prior <- input$init$prior
    p_acc_min <- input$init$p_acc_min
    n_alpha <- input$init$n_alpha
    nparam <- input$init$nparam
    nstat <- input$init$nstat
    nb_simul <- input$init$nb_simul
    summary_stat_target <- input$init$summary_stat_target
    dist_weights <- input$init$dist_weights
    inside_prior <- input$init$inside_prior
    nb_simul_step <- input$init$nb_simul_step

    # prior wave
    seed_count <- input$pwave$seed_count
    simul_below_tol <- input$pwave$simul_below_tol
    tab_weight <- input$pwave$tab_weight
    sd_simul <- input$pwave$sd_simul
    tol_next <- input$pwave$tol_next
    tab_dist <- input$pwave$tab_dist

    # current wave simulation
    tab_inic <- input$tab_inic

    if (wave == 1) {
      p_acc <- p_acc_min + 1
    }

    simul_below_tol2 <- NULL

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

    out <- list(init = input$init,
                pwave = list(tab_weight = tab_weight, seed_count = seed_count,
                             simul_below_tol = simul_below_tol, tab_dist = tab_dist,
                             tol_next = tol_next, sd_simul = sd_simul, p_acc = p_acc))
    return(out)
  }

  if (wave == Inf) {
    final_res <- list(param = as.matrix(as.matrix(simul_below_tol)[, 1:nparam]),
                     stats = as.matrix(as.matrix(simul_below_tol)[, (nparam + 1):(nparam + nstat)]),
                     weights = tab_weight/sum(tab_weight), stats_normalization = as.numeric(sd_simul),
                     epsilon = max(.compute_dist(summary_stat_target,
                                                 as.matrix(as.matrix(simul_below_tol)[, (nparam + 1):(nparam + nstat)]),
                                                 sd_simul, dist_weights = dist_weights)))

    return(final_res)
  }

}


#' @export
abc_wave0 <- function(model,
                      prior,
                      prior_test,
                      nb_simul,
                      seed_count,
                      n_cluster,
                      batch) {

  tab_param <- NULL
  list_param <- list(NULL)
  nparam <- length(prior)
  l <- length(prior)
  random_tab <- NULL
  all_unif_prior <- .all_unif(prior)
  if (all_unif_prior) {
    random_tab <- randomLHS(nb_simul, nparam)
  }

  list_param = list(NULL)
  for (i in 1:nb_simul) {
    param = array(0, l)
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

  if (!is.null(batch)) {
    simset <- batch_to_sims(n_cluster, batch)
    list_param <- list_param[simset]
    tab_param <- tab_param[simset, , drop = FALSE]
  }

  cl <- makeCluster(n_cluster)
  list_simul_summarystat <- parLapplyLB(cl, list_param, model)
  stopCluster(cl)

  tab_simul_summarystat <- do.call("rbind", list_simul_summarystat)

  options(scipen = 0)

  out <- cbind(tab_param, tab_simul_summarystat)
  return(out)
}

#' @export
abc_waveN <- function(model,
                      prior,
                      param_previous_step,
                      tab_weight,
                      nb_simul,
                      seed_count,
                      inside_prior,
                      n_cluster,
                      max_pick = 10000,
                      batch) {

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

  if (!is.null(batch)) {
    simset <- batch_to_sims(n_cluster, batch)
    list_param <- list_param[simset]
  }

  cl <- makeCluster(n_cluster)
  list_simul_summarystat = parLapplyLB(cl, list_param, model)
  stopCluster(cl)

  tab_simul_summarystat <- do.call("rbind", list_simul_summarystat)

  out <- list(cbind(tab_param, tab_simul_summarystat), nb_simul/k_acc)
  return(out)
}

batch_to_sims <- function(batchSize, batchNum) {
  ((batchSize*batchNum) - (batchSize - 1)):(batchSize*batchNum)
}
