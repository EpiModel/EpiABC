
#' @export
abc_smc_prep <- function(model,
                         prior,
                         nsims,
                         summary_stat_target,
                         prior_test = NULL,
                         ncores = 2,
                         dist_weights = NULL,
                         alpha = 0.5,
                         ...) {

  p_acc_min <- 0.1

  ## checking errors in the inputs
  if (missing(model))
    stop("model is missing")
  if (missing(prior))
    stop("prior is missing")
  data <- .wrap_constants_in_model(prior, model, use_seed = TRUE)
  prior <- data$new_prior
  model <- data$new_model
  prior <- .process_prior(prior)
  if (!is.null(prior_test))
    .check_prior_test(length(prior), prior_test)
  if (missing(nsims))
    stop("nsims is missing")
  if (missing(summary_stat_target))
    stop("summary_stat_target is missing")
  if (!is.vector(nsims) || length(nsims) > 1 || nsims < 1)
    stop("nsims must be a number larger than 1.")
  nsims <- floor(nsims)
  if (!is.vector(summary_stat_target))
    stop("'summary_stat_target' has to be a vector.")
  if (!is.vector(ncores) || length(ncores) > 1 || ncores < 1)
    stop("'ncores' has to be a positive number.")
  ncores <- floor(ncores)
  if (!is.null(dist_weights) && length(dist_weights) != length(summary_stat_target)) {
    stop("'dist_weights' has to be the same length than 'summary_stat_target'")
  }

  # batch sizes for wave 0 and waves 1+
  batchSize <- c(ceiling(nsims/ncores),
                 ceiling((nsims - ceiling(nsims * alpha))/ncores))

  out <- list(model = model, prior = prior, prior_test = prior_test,
              nsims = nsims, batchSize = batchSize,
              summary_stat_target = summary_stat_target,
              ncores = ncores, dist_weights = dist_weights,
              alpha = alpha, p_acc_min = p_acc_min)

  return(out)
}


#' @export
abc_smc_wave <- function(input = "data/", wave, batch, save = TRUE, outdir = "data/") {

  if (class(input) == "character") {
    file <- list.files(input, pattern = paste0("wave", wave - 1, ".rda"), full.names = TRUE)
    input <- readRDS(file)
  }

  if (wave == 0) {

    # fixed params/settings
    model <- input$model
    prior <- input$prior
    prior_test <- input$prior_test
    nsims <- input$nsims
    summary_stat_target <- input$summary_stat_target
    ncores <- input$ncores
    alpha <- input$alpha

    seed_count <- 0
    input$inside_prior <- TRUE
    input$max_pick <- 10000

    input$nparam <- length(prior)
    input$nstat <- length(summary_stat_target)
    if (!.all_unif(prior)) {
      stop("Prior distributions must be uniform")
    }
    n_alpha <- input$n_alpha <- ceiling(nsims * alpha)
    input$nsims_step <- nsims - n_alpha

    tab_ini <- abc_wave0(model,
                         prior,
                         prior_test,
                         nsims,
                         seed_count,
                         ncores,
                         batch = batch)

    out <- list(init = input, seed_count = seed_count, tab_ini = tab_ini)

    if (save == TRUE) {
      saveRDS(out, file = paste0(outdir, "abc.wave0.batch",
                                 stringr::str_pad(batch, 4, pad = "0"), ".rda"))
    } else {
      return(out)
    }
  }

  if (wave > 0) {
    tab_inic <- abc_waveN(input = input, batch = batch)
    out <- list(init = input$init, pwave = input$pwave, tab_inic = tab_inic)
    if (save == TRUE) {
      saveRDS(out, file = paste0(outdir, "abc.wave", wave, ".batch",
                                 stringr::str_pad(batch, 4, pad = "0"), ".rda"))
    } else {
      return(out)
    }
  }
}

#' @export
abc_smc_process <- function(input = "data/", wave, save = TRUE, outdir = "data/") {

  if (class(input) == "character") {
    file <- list.files(input, pattern = paste0("wave", wave, ".rda"), full.names = TRUE)
    if (length(file) == 0) return(NULL)
    input <- readRDS(file)
  }

  if (wave == 0) {

    # fixed
    prior <- input$init$prior
    n_alpha <- input$init$n_alpha
    nparam <- input$init$nparam
    nstat <- input$init$nstat
    nsims <- input$init$nsims
    nsims_step <- input$init$nsims_step
    summary_stat_target <- input$init$summary_stat_target
    alpha <- input$init$alpha
    dist_weights <- input$init$dist_weights
    inside_prior <- input$init$inside_prior
    max_pick <- input$init$max_pick


    # current wave simulation
    seed_count <- input$seed_count
    tab_ini <- input$tab_ini

    # initially, weights are equal
    tab_weight <- array(1, n_alpha)
    seed_count <- seed_count + nsims
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


    # particle selection for wave 1 -------------------------------------------

    param_previous_step <- as.matrix(as.matrix(simul_below_tol)[, 1:nparam])
    tab_weight <- tab_weight/sum(tab_weight)

    tab_param <- NULL
    k_acc <- 0
    list_param <- list(NULL)

    list_param <- list(NULL)
    for (i in 1:nsims_step) {
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
    seed_count <- seed_count + nsims_step

    out <- list(init = input$init,
                pwave = list(tab_weight = tab_weight, seed_count = seed_count,
                             simul_below_tol = simul_below_tol, tab_dist = tab_dist,
                             tol_next = tol_next, sd_simul = as.numeric(sd_simul)),
                cwave = list(list_param = list_param, tab_param = tab_param, k_acc = k_acc))

    if (save == TRUE) {
      saveRDS(out, file = paste0(outdir, "abc.wave", wave, ".rda"))
    } else {
      return(out)
    }
  }

  if (wave > 0) {

    # fixed
    prior <- input$init$prior
    p_acc_min <- input$init$p_acc_min
    n_alpha <- input$init$n_alpha
    nparam <- input$init$nparam
    nstat <- input$init$nstat
    nsims <- input$init$nsims
    summary_stat_target <- input$init$summary_stat_target
    dist_weights <- input$init$dist_weights
    inside_prior <- input$init$inside_prior
    nsims_step <- input$init$nsims_step
    max_pick <- input$init$max_pick

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
    dim(tab_ini) <- c(nsims_step, (nparam + nstat))
    seed_count <- seed_count + nsims_step
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
    p_acc <- length(tab_dist2[!is.na(tab_dist2) & tab_dist2 <= tol_next])/nsims_step
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


    # particle selection for wave N+1 -------------------------------------------

    param_previous_step = as.matrix(as.matrix(simul_below_tol)[, 1:nparam])
    tab_weight = tab_weight/sum(tab_weight)

    tab_param <- NULL
    k_acc <- 0
    list_param <- list(NULL)

    list_param <- list(NULL)
    for (i in 1:nsims_step) {
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
    seed_count <- seed_count + nsims_step

    out <- list(init = input$init,
                pwave = list(tab_weight = tab_weight, seed_count = seed_count,
                             simul_below_tol = simul_below_tol, tab_dist = tab_dist,
                             tol_next = tol_next, sd_simul = sd_simul, p_acc = p_acc),
                cwave = list(list_param = list_param, tab_param = tab_param, k_acc = k_acc))

    if (save == TRUE) {
      saveRDS(out, file = paste0(outdir, "abc.wave", wave, ".rda"))
    } else {
      return(out)
    }
  }

}


#' @export
out_abc <- function(wave, input = "data/") {

  if (class(input) == "character") {
    file <- list.files(input, pattern = paste0("wave", wave, ".rda"), full.names = TRUE)
    if (length(file) == 0) return(NULL)
    input <- readRDS(file)
  }

  # fixed
  nparam <- input$init$nparam
  nstat <- input$init$nstat
  summary_stat_target <- input$init$summary_stat_target
  dist_weights <- input$init$dist_weights

  # prior wave
  simul_below_tol <- input$pwave$simul_below_tol
  tab_weight <- input$pwave$tab_weight
  sd_simul <- input$pwave$sd_simul
  p_acc <- input$pwave$p_acc

  final_res <- list(param = as.matrix(as.matrix(simul_below_tol)[, 1:nparam]),
                    stats = as.matrix(as.matrix(simul_below_tol)[, (nparam + 1):(nparam + nstat)]),
                    target = summary_stat_target,
                    weights = tab_weight/sum(tab_weight), stats_normalization = as.numeric(sd_simul),
                    epsilon = max(.compute_dist(summary_stat_target,
                                                as.matrix(as.matrix(simul_below_tol)[, (nparam + 1):(nparam + nstat)]),
                                                sd_simul, dist_weights = dist_weights)),
                    wave = wave, p_acc = p_acc)

  return(final_res)

}

#' @export
summary_abc <- function(input) {

  cat("Wave:", input$wave)
  cat("\np_acc:", input$p_acc)

  paramSumm <- apply(input$param, 2, summary)
  colnames(paramSumm) <- rep("", ncol(paramSumm))
  statsSumm <- apply(input$stats, 2, summary)
  colnames(statsSumm) <- rep("", ncol(statsSumm))

  cat("\n\nTarget Stats: \n", input$target)
  cat("\nSimulated Targets")
  cat("\n-----------------")
  print(statsSumm)

  cat("\nEsimated Parameters")
  cat("\n-------------------")
  print(paramSumm)

}


#' @export
abc_wave0 <- function(model,
                      prior,
                      prior_test,
                      nsims,
                      seed_count,
                      ncores,
                      batch) {

  tab_param <- NULL
  list_param <- list(NULL)
  nparam <- length(prior)
  l <- length(prior)
  random_tab <- NULL
  all_unif_prior <- .all_unif(prior)
  if (all_unif_prior) {
    random_tab <- randomLHS(nsims, nparam)
  }

  list_param = list(NULL)
  for (i in 1:nsims) {
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
  seed_count <- seed_count + nsims

  if (!is.null(batch)) {
    simset <- batch_to_sims(ncores, batch, nsims)
    list_param <- list_param[simset]
    tab_param <- tab_param[simset, , drop = FALSE]
  }

  cl <- makeCluster(ncores)
  list_simul_summarystat <- parLapplyLB(cl, list_param, model)
  stopCluster(cl)

  tab_simul_summarystat <- do.call("rbind", list_simul_summarystat)

  options(scipen = 0)

  out <- cbind(tab_param, tab_simul_summarystat)
  return(out)
}

#' @export
abc_waveN <- function(input, batch) {

  # Fixed
  model <- input$init$model
  nsims <- input$init$nsims_step
  ncores <- input$init$ncores

  # Current wave
  list_param <- input$cwave$list_param
  tab_param <- input$cwave$tab_param
  k_acc <- input$cwave$k_acc

  if (!is.null(batch)) {
    simset <- batch_to_sims(ncores, batch, nsims)
    list_param <- list_param[simset]
    tab_param <- tab_param[simset, , drop = FALSE]
  }

  cl <- makeCluster(ncores)
  list_simul_summarystat = parLapplyLB(cl, list_param, model)
  stopCluster(cl)

  tab_simul_summarystat <- do.call("rbind", list_simul_summarystat)

  out <- list(cbind(tab_param, tab_simul_summarystat), nsims/k_acc)

  return(out)
}

batch_to_sims <- function(batchSize, batchNum, totSims) {
  ((batchSize*batchNum) - (batchSize - 1)):(min(batchSize*batchNum, totSims))
}

#' @export
merge_abc <- function(wave, indir = "data/", outdir = "data/") {

  files <- list.files(indir, pattern = paste0("wave", wave, ".batch"), full.names = TRUE)

  if (wave == 0) {
    for (i in 1:length(files)) {
      if (i == 1) {
        dat <- readRDS(files[i])
        nBatches <- ceiling(dat$init$nsims/dat$init$ncores)
        if (length(files) < nBatches) return(NULL)
      } else {
        temp.dat <- readRDS(files[i])
        dat$tab_ini <- rbind(dat$tab_ini, temp.dat$tab_ini)
      }
    }
  }

  if (wave > 0) {
    for (i in 1:length(files)) {
      if (i == 1) {
        dat <- readRDS(files[i])
        nBatches <- ceiling((dat$init$nsims - ceiling(dat$init$nsims * dat$init$alpha))/dat$init$ncores)
        if (length(files) < nBatches) return(NULL)
      } else {
        temp.dat <- readRDS(files[i])
        dat$tab_inic[[1]] <- rbind(dat$tab_inic[[1]], temp.dat$tab_inic[[1]])
      }
    }
  }

  saveRDS(dat, file = paste0(outdir, "abc.wave", wave, ".rda"))
  unlink(files)

}
