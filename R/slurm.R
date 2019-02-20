
#' ABC-SMC Slurm Model Preparation
#'
#' This is the first step in the ABC-SMC Slurm workflow.
#'
#' @inheritParams abc_smc_cluster
#'
#' @param ncores Number of cores per node (defines a batch size).
#' @param alpha Number of simulations to retain in waves 1+.
#'
#' @export
#'
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
    stop("summary_stat_target must be a vector.")
  if (!is.vector(ncores) || length(ncores) > 1 || ncores < 1)
    stop("ncores must be a positive number.")
  ncores <- floor(ncores)
  if (!is.null(dist_weights) && length(dist_weights) != length(summary_stat_target)) {
    stop("dist_weights must be the same length than 'summary_stat_target'")
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


#' Simulates the ABC-SMC Slurm Model
#'
#' This is the second step in the ABC-SMC Slurm workflow.
#'
#' @param input A character string containing the directory of the output file
#'        from \code{\link{abc_smc_prep}} which should be saved as an RDS file
#'        with the name \code{abc.wave0.rda} (the default), or the object itself.
#' @param wave SMC wave number, where the initial wave = 0. In the standard Slurm
#'        workflow, this would get passed in as an environmental variable
#'        \code{wave} from the master bash script.
#' @param batch Batch number for the simulation set, which corresponds to the
#'        array number passed in as an environmental Slurm variable
#'        \code{SLURM_ARRAY_TASK_ID}.
#' @param save If \code{TRUE}, writes output to an RDS file with the name
#'        \code{abc.waveX.batchY.rda} in the directory specified by \code{outdir},
#'        where \code{X} is the value of \code{wave} and \code{Y} is the value
#'        of \code{batch}.
#' @param outdir Path to save the output RDS file if \code{save=TRUE}.
#'
#' @export
#'
abc_smc_wave <- function(input = "data/",
                         wave,
                         batch,
                         save = TRUE,
                         outdir = "data/") {

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


#' Processes ABC-SMC Slurm Simulation Output and Selects Next Wave Particles
#'
#' This is the fourth step in the ABC-SMC Slurm workflow.
#'
#' @inheritParams abc_smc_wave
#'
#' @export
#'
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


    ## particle selection for wave 1

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


    ## particle selection for wave N+1

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


#' Extracts Posterior Distribution for Parameters and Summary Statistics
#'
#' Once a simulation wave is complete, this function processes the output and
#' stores it in a format that is useful for numerical analysis and plotting.
#'
#' @param wave If \code{input} is a character string, the wave file that should
#'        be read from that directory.
#' @param input Either a character string with the directory to read the wave
#'        files created with \code{\link{abc_smc_wave}} from, or the direct object
#'        itself.
#'
#' @export
#'
get_posterior <- function(wave, input = "data/") {

  if (class(input) == "character") {
    file <- list.files(input, pattern = paste0("wave", wave, ".rda"), full.names = TRUE)
    if (length(file) == 0) {
      stop("No files in ", input, " named ", paste0("abc.wave", wave, ".rda"), call. = FALSE)
    }
    input <- readRDS(file)
  }

  # fixed
  nparam <- input$init$nparam
  nstat <- input$init$nstat
  summary_stat_target <- input$init$summary_stat_target
  dist_weights <- input$init$dist_weights

  priors <- list()
  for (i in 1:length(input$init$prior)) {
    priors[[i]] <- input$init$prior[[i]]$sampleArgs[2:3]
  }

  # prior wave
  simul_below_tol <- input$pwave$simul_below_tol
  tab_weight <- input$pwave$tab_weight
  sd_simul <- input$pwave$sd_simul
  p_acc <- input$pwave$p_acc

  out <- list(param = as.matrix(as.matrix(simul_below_tol)[, 1:nparam]),
              priors = priors,
              stats = as.matrix(as.matrix(simul_below_tol)[, (nparam + 1):(nparam + nstat)]),
              target = summary_stat_target,
              weights = tab_weight/sum(tab_weight), stats_normalization = as.numeric(sd_simul),
              epsilon = max(.compute_dist(summary_stat_target,
                                          as.matrix(as.matrix(simul_below_tol)[, (nparam + 1):(nparam + nstat)]),
                                          sd_simul, dist_weights = dist_weights)),
              wave = wave, p_acc = p_acc)

  class(out) <- "abcsmc"
  return(out)
}


#' Numerical Summary of Posterior Distribution of Parameters and Summary Statistics
#'
#' @param object Output from \code{\link{get_posterior}}.
#' @param digits Significant digits to print in output.
#' @param ... Additional arguments based to generic \code{summary}.
#'
#' @method summary abcsmc
#' @export
#'
summary.abcsmc <- function(object, digits = 3, ...) {

  cat("ABC-SMC Model Summary")
  cat("\n==========================")
  cat("\nWave:", object$wave)
  cat("\np_acc:", object$p_acc)

  statsSumm <- apply(object$stats, 2, summary)
  statsSumm <- rbind(object$target, statsSumm)
  colnames(statsSumm) <- rep("", ncol(statsSumm))
  rownames(statsSumm)[1] <- "Target"
  statsSumm <- round(statsSumm, digits)

  paramSumm <- apply(object$param, 2, summary)

  allPriors <- matrix(nrow = 2, ncol = length(object$priors))
  for (i in 1:length(object$priors)) {
    allPriors[, i] <- as.numeric(object$priors[[i]])
  }
  paramSumm <- rbind(allPriors, paramSumm)
  colnames(paramSumm) <- rep("", ncol(paramSumm))
  rownames(paramSumm)[1:2] <- c("Prior L", "Prior H")
  paramSumm <- round(paramSumm, digits)

  cat("\n\nModel Statistics")
  cat("\n-------------------")
  print(statsSumm)

  cat("\nModel Parameters")
  cat("\n-------------------")
  print(paramSumm)

}


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


#' Combines Individual Batch Simulation Files for Processing in Next Wave
#'
#' This is the third step in the ABC-SMC Slurm workflow.
#'
#' @inheritParams abc_smc_process
#' @param indir Directory containing the batch files output from
#'        \code{\link{abc_smc_wave}} for merging into a single file.
#'
#' @export
#'
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


#' @title Create sbatch Bash Shell Script ABC-SMC Workflow
#'
#' @description Creates a master-level SLURM::sbatch script given a ABC-SMC
#'              model prepared and output from \code{abc_smc_prep}.
#'
#' @param input Output object from \code{\link{abc_smc_prep}}.
#' @param nwaves Number of waves of SMC to run, not including wave = 0 that is the
#'        starting wave based on ABC rejection.
#' @param master.file Name of the output bash shell script file to write. If
#'        \code{""}, then will print to console.
#' @param runsim.file Name of the bash shell script file that runs the R simulation
#'        job.
#' @param ckpt If \code{TRUE}, use the checkpoint queue to submit jobs. If
#'        numeric, will specify the first X jobs on the grid as non-backfill.
#' @param append If \code{TRUE}, will append lines to a previously created shell
#'        script. New simno will either start with value of \code{simno.start}
#'        or the previous value if missing.
#' @param mem Amount of memory needed per node within each Slurm job.
#' @param walltime Amount of clock time needed per Slurm job.
#' @param partition.main Name of primary HPC partition (passed to -p).
#' @param partition.ckpt Name of checkpoint HPC partition (passed to -p).
#' @param account.main Name of primary account (passed to -A).
#' @param account.ckpt Name of checkpoint account (passed to -A).
#'
#' @export
#'
sbatch_master_abc <- function(input,
                              nwaves,
                              master.file = "master.sh",
                              runsim.file = "runsim.sh",
                              ckpt = FALSE,
                              append = FALSE,
                              mem = "55G",
                              walltime = "1:00:00",
                              partition.main = "csde",
                              partition.ckpt = "ckpt",
                              account.main = "csde",
                              account.ckpt = "csde-ckpt") {



  if (append == TRUE) {
    tab <- readLines(master.file, skipNul = TRUE)
    tab.last <- tail(tab, 1)
    tab.split <- strsplit(tab.last, " ")[[1]]
    job <- tab.split[grep("--job-name", tab.split)]
    last.job <- as.numeric(strsplit(job, "wave")[[1]][2])
    start.job <- last.job + 1
  } else {
    start.job <- 0
  }

  ncores <- input$ncores
  batchSize <- input$batchSize

  if (ckpt == TRUE) {
    pA <- paste("-p", partition.ckpt, "-A", account.ckpt)
  } else {
    pA <- paste("-p", partition.main, "-A", account.main)
  }

  if (append == FALSE) {
    cat("#!/bin/bash\n", file = master.file)
  }

  for (i in start.job:nwaves) {

    if (i == 0) {
      array <- paste0("--array=1-", batchSize[1])
    } else {
      array <- paste0("--array=1-", batchSize[2])
    }
    job <- paste0("--job-name=wave", i)
    expt <- paste0("--export=ALL,wave=", i)
    if (i == start.job) {
      depend <- ""
    } else {
      depend <- paste0("--depend=afterany:$(squeue --noheader --format %i --name wave", i - 1, ")")
    }
    ntasks <- paste0("--ntasks-per-node=", ncores)
    memout <- paste0("--mem=", mem)
    time <- paste0("--time=", walltime)

    cat("\nsbatch", pA, array,
        job, expt, depend, ntasks,
        memout, time, runsim.file,
        file = master.file, append = TRUE, sep = " ")
  }
  cat("\n", file = master.file, append = TRUE)

}


#' Kernel Density Plot for Posterior Distribution of Parameters and Summary Statistics
#'
#' @param x Output from \code{\link{get_posterior}}.
#' @param type Character string of \code{type="stats"} for model statistics
#'        or \code{type="param"} for model parameters.
#' @param stats Numeric vector of statistics (or parameters) by position in
#'        \code{x} to plot.
#' @param mean.line If \code{TRUE}, plot the mean of the stat/parameter distribution
#'        as a blue line.
#' @param stats.positive If \code{TRUE}, constrain the lower bound of the kernel
#'        density at 0 for summary statistics plots.
#' @param ... Additional arguments based to generic \code{plot}.
#'
#' @method plot abcsmc
#' @export
#'
plot.abcsmc <- function(x, type, stats, mean.line = TRUE,
                        stats.positive = FALSE, ...) {

  summstats <- as.data.frame(x$stats)
  param <- as.data.frame(x$param)
  if (type == "stats") {
    nstats <- ncol(summstats)
  } else if (type == "param") {
    nstats <- ncol(param)
  }

  if (missing(stats)) {
    stats <- 1:nstats
  } else {
    if (max(stats) > nstats) {
      stop("Maximum for stats is ", nstats, call. = FALSE)
    }
    nstats <- length(stats)
  }

  if (nstats == 1) dimens <- c(1, 1)
  if (nstats == 2) dimens <- c(1, 2)
  if (nstats == 3) dimens <- c(1, 3)
  if (nstats == 4) dimens <- c(2, 2)
  if (nstats == 5) dimens <- c(2, 3)
  if (nstats == 6) dimens <- c(2, 3)
  if (nstats %in% 7:9) dimens <- c(3, 3)
  if (nstats %in% 10:12) dimens <- c(4, 3)
  if (nstats %in% 13:16) dimens <- c(4, 4)
  if (nstats > 16) dimens <- rep(ceiling(sqrt(nstats)), 2)

  # Pull graphical parameters
  ops <- list(mar = par()$mar, mfrow = par()$mfrow, mgp = par()$mgp)
  par(mar = c(2.2,3,2,1), mgp = c(2, 1, 0), mfrow = dimens)

  if (type == "stats") {
    if (stats.positive == FALSE) {
      for (i in stats) {
        plot(density(summstats[, i]), main = paste0("Statistic ", i), lwd = 2)
        grid()
        abline(v = x$target[i], col = "red", lty = 2, lwd = 2)
        if (mean.line == TRUE) {
          abline(v = colMeans(summstats[, i, drop = FALSE]), lty = 1, lwd = 1.5, col = "blue")
        }
      }
    } else {
      for (i in stats) {
        plot(density(summstats[, i], from = 0), main = paste0("Statistic ", i), lwd = 2)
        grid()
        abline(v = x$target[i], col = "red", lty = 2, lwd = 2)
        if (mean.line == TRUE) {
          abline(v = colMeans(summstats[, i, drop = FALSE]), lty = 1, lwd = 1.5, col = "blue")
        }
      }
    }
  }

  if (type == "param") {
    for (i in stats) {
      plot(density(param[, i], from = as.numeric(x$priors[[i]][1]),
                   to = as.numeric(x$priors[[i]][2])),
           main = paste0("Parameter ", i), lwd = 2,
           xlim = as.numeric(x$priors[[i]]))
      grid()
      abline(v = as.numeric(x$priors[[i]]), col = "red", lty = 2, lwd = 2)
      if (mean.line == TRUE) {
        abline(v = colMeans(param[, i, drop = FALSE]), lty = 1, lwd = 1.5, col = "blue")
      }
    }
  }

  on.exit(par(ops))
}

#' Boxplot for Posterior Distribution of Parameters and Summary Statistics
#'
#' @param x Output from \code{\link{get_posterior}}.
#' @param type Character string of \code{type="stats"} for model statistics
#'        or \code{type="param"} for model parameters.
#' @param stats Numeric vector of statistics (or parameters) by position in
#'        \code{x} to plot.
#' @param ... Additional arguments based to generic \code{plot}.
#'
#' @method boxplot abcsmc
#' @export
#'
boxplot.abcsmc <- function(x, type, stats, ...) {

  summstats <- as.data.frame(x$stats)
  param <- as.data.frame(x$param)
  if (type == "stats") {
    nstats <- ncol(summstats)
  } else if (type == "param") {
    nstats <- ncol(param)
  }

  if (missing(stats)) {
    stats <- 1:nstats
  } else {
    if (max(stats) > nstats) {
      stop("Maximum for stats is ", nstats, call. = FALSE)
    }
    nstats <- length(stats)
  }

  if (nstats == 1) dimens <- c(1, 1)
  if (nstats == 2) dimens <- c(1, 2)
  if (nstats == 3) dimens <- c(1, 3)
  if (nstats == 4) dimens <- c(2, 2)
  if (nstats == 5) dimens <- c(2, 3)
  if (nstats == 6) dimens <- c(2, 3)
  if (nstats %in% 7:9) dimens <- c(3, 3)
  if (nstats %in% 10:12) dimens <- c(4, 3)
  if (nstats %in% 13:16) dimens <- c(4, 4)
  if (nstats > 16) dimens <- rep(ceiling(sqrt(nstats)), 2)

  # Pull graphical parameters
  ops <- list(mar = par()$mar, mfrow = par()$mfrow, mgp = par()$mgp)
  par(mar = c(2,2.5,2,1), mgp = c(2, 1, 0), mfrow = dimens)

  if (type == "stats") {
    for (i in stats) {
      boxplot(summstats[, i], col = adjustcolor("steelblue", alpha.f = 0.6),
              main = paste0("Statistic ", i))
      points(x$target[i], pch = 20, cex = 2, col = "red")
    }
  }

  if (type == "param") {
    for (i in stats) {
      boxplot(param[, i], col = adjustcolor("steelblue", alpha.f = 0.6),
              main = paste0("Parameter ", i))
    }
  }

  on.exit(par(ops))
}
