
.ABC_sequential_cluster <- function(method, model, prior, prior_test, nb_simul, summary_stat_target,
                                    n_cluster, use_seed, verbose, dist_weights=NULL, cl, ...) {
  if (use_seed == FALSE) {
    stop("For parallel implementations, you must specify the option 'use_seed=TRUE' and modify your model accordingly -
             see the package's vignette for more details.")
  }
  options(scipen = 50)
  return(switch(EXPR = method,
                Beaumont = .ABC_PMC_cluster(model, prior, prior_test,
                                            nb_simul, summary_stat_target, n_cluster, verbose,
                                            dist_weights = dist_weights, cl, ...),
                Drovandi = .ABC_Drovandi_cluster(model,
                                                 prior, nb_simul, summary_stat_target, n_cluster, verbose,
                                                 dist_weights = dist_weights, cl,  ...),
                Delmoral = .ABC_Delmoral_cluster(model,
                                                 prior, prior_test, nb_simul, summary_stat_target, n_cluster,
                                                 verbose, dist_weights = dist_weights, cl, ...),
                Lenormand = .ABC_Lenormand_cluster(model = model, prior = prior, prior_test = prior_test,
                                                   nb_simul = nb_simul, summary_stat_target = summary_stat_target,
                                                   n_cluster = n_cluster, verbose = verbose,
                                                   dist_weights = dist_weights, cl = cl, ...)))
  options(scipen = 0)
}


## sequential algorithm of Lenormand et al. 2012
.ABC_Lenormand_cluster <- function(model, prior, prior_test, nb_simul, summary_stat_target,
                                   n_cluster, verbose, alpha = 0.5, p_acc_min = 0.05,
                                   dist_weights=NULL, cl, seed_count = 0,
                                   inside_prior = TRUE, progress_bar = TRUE, max_pick=10000) {

  ## checking errors in the inputs
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
  if (!is.vector(seed_count))
    stop("'seed_count' has to be a number.")
  if (length(seed_count) > 1)
    stop("'seed_count' has to be a number.")
  if (seed_count < 0)
    stop("'seed_count' has to be a positive number.")
  seed_count = floor(seed_count)
  if (!is.logical(inside_prior))
    stop("'inside_prior' has to be boolean.")
  start = Sys.time()
  if (progress_bar) {
    print("## Lenormand et al. (2012)'s algorithm ##")
  }

  seed_count_ini = seed_count
  nparam = length(prior)
  nstat = length(summary_stat_target)
  if (!.all_unif(prior)) {
    stop("Prior distributions must be uniform to use the Lenormand et al. (2012)'s algorithm.")
  }
  n_alpha = ceiling(nb_simul * alpha)

  ## step 1 ABC rejection step with LHS
  tab_ini = .ABC_rejection_lhs_cluster(model, prior, prior_test, nb_simul, seed_count, n_cluster, cl)
  # initially, weights are equal
  tab_weight = array(1, n_alpha)
  seed_count = seed_count + nb_simul
  # determination of the normalization constants in each dimension associated to
  # each summary statistic, this normalization will not change during all the
  # algorithm
  sd_simul = sapply(as.data.frame(tab_ini[, (nparam + 1):(nparam + nstat)]), sd,
                    na.rm = TRUE)
  # selection of the alpha quantile closest simulations
  simul_below_tol = NULL
  simul_below_tol = rbind(simul_below_tol,
                          .selec_simul_alpha(summary_stat_target,
                                             tab_ini[, 1:nparam],
                                             tab_ini[, (nparam + 1):(nparam + nstat)], sd_simul,
                                             alpha, dist_weights = dist_weights))
  # to be sure that there are not two or more simulations at a distance equal
  #   to the tolerance determined by the quantile
  simul_below_tol = simul_below_tol[1:n_alpha, ]
  tab_dist = .compute_dist(summary_stat_target,
                           as.matrix(as.matrix(simul_below_tol)[, (nparam + 1):(nparam + nstat)]),
                           sd_simul, dist_weights = dist_weights)
  if (!is.null(dist_weights)) {
    tab_dist = tab_dist * (dist_weights/sum(dist_weights))
  }
  tol_next = max(tab_dist)
  intermediary_steps = list(NULL)
  if (progress_bar) {
    print("Step 1 Completed ...")
  }

  ## following steps
  p_acc = p_acc_min + 1
  nb_simul_step = nb_simul - n_alpha
  it = 1

  while (p_acc > p_acc_min) {
    it = it + 1
    simul_below_tol2 = NULL
    tab_inic = .ABC_launcher_not_uniformc_cluster(model = model,
                                                  prior = prior,
                                                  param_previous_step = as.matrix(as.matrix(simul_below_tol)[, 1:nparam]),
                                                  tab_weight = tab_weight/sum(tab_weight),
                                                  nb_simul = nb_simul_step,
                                                  seed_count = seed_count,
                                                  inside_prior = inside_prior,
                                                  n_cluster = n_cluster,
                                                  cl = cl,
                                                  max_pick = max_pick)
    tab_ini = as.matrix(tab_inic[[1]])
    tab_ini = as.numeric(tab_ini)
    dim(tab_ini) = c(nb_simul_step, (nparam + nstat))
    seed_count = seed_count + nb_simul_step
    if (!inside_prior) {
      tab_weight2 = .compute_weightb(as.matrix(as.matrix(as.matrix(tab_ini)[, 1:nparam])),
                                     as.matrix(as.matrix(as.matrix(simul_below_tol)[, 1:nparam])),
                                     tab_weight/sum(tab_weight), prior)
    } else {
      tab_weight2 = tab_inic[[2]] * (.compute_weightb(as.matrix(as.matrix(as.matrix(tab_ini)[,1:nparam])),
                                                      as.matrix(as.matrix(as.matrix(simul_below_tol)[, 1:nparam])),
                                                      tab_weight/sum(tab_weight), prior))
    }
    simul_below_tol2 = rbind(as.matrix(simul_below_tol), as.matrix(tab_ini))
    tab_weight = c(tab_weight, tab_weight2)
    tab_dist2 = .compute_dist(summary_stat_target,
                              as.matrix(as.matrix(tab_ini)[, (nparam + 1):(nparam + nstat)]),
                              sd_simul, dist_weights = dist_weights)
    if (!is.null(dist_weights)) {
      tab_dist2 = tab_dist2 * (dist_weights/sum(dist_weights))
    }
    p_acc = length(tab_dist2[!is.na(tab_dist2) & tab_dist2 <= tol_next])/nb_simul_step
    tab_dist = c(tab_dist, tab_dist2)
    tol_next = sort(tab_dist)[n_alpha]
    simul_below_tol2 = simul_below_tol2[!is.na(tab_dist) & tab_dist <= tol_next, ]
    tab_weight = tab_weight[!is.na(tab_dist) & tab_dist <= tol_next]
    tab_weight = tab_weight[1:n_alpha]
    tab_dist = tab_dist[!is.na(tab_dist) & tab_dist <= tol_next]
    odist = order(tab_dist, decreasing = FALSE)[1:n_alpha]
    tab_dist_new = tab_dist
    simul_below_tol = matrix(0, n_alpha, (nparam + nstat))
    for (i1 in 1:n_alpha) {
      tab_dist_new[i1] = tab_dist[odist[i1]]
      for (i2 in 1:(nparam + nstat)) {
        simul_below_tol[i1, i2] = as.numeric(simul_below_tol2[odist[i1], i2])
      }
    }
    tab_dist = tab_dist_new[1:n_alpha]
    if (progress_bar) {
      print(paste("Step ", it, " Completed: p_acc = ", p_acc, sep = ""))
    }
  }

  final_res = NULL
  final_res = list(param = as.matrix(as.matrix(simul_below_tol)[, 1:nparam]),
                   stats = as.matrix(as.matrix(simul_below_tol)[, (nparam + 1):(nparam + nstat)]),
                   weights = tab_weight/sum(tab_weight), stats_normalization = as.numeric(sd_simul),
                   epsilon = max(.compute_dist(summary_stat_target, as.matrix(as.matrix(simul_below_tol)[, (nparam + 1):(nparam + nstat)]),
                                               sd_simul, dist_weights = dist_weights)),
                   nsim = (seed_count - seed_count_ini), computime = as.numeric(difftime(Sys.time(), start, units = "secs")))

  final_res
}


## function to sample in the prior distributions using a Latin Hypercube sample
.ABC_rejection_lhs_cluster <- function(model, prior, prior_test, nb_simul, seed_count,
                                       n_cluster, cl) {

  tab_simul_summarystat = NULL
  tab_param = NULL
  list_param = list(NULL)
  # n_end = nb_simul
  nparam = length(prior)
  l = length(prior)
  random_tab = NULL
  all_unif_prior = .all_unif(prior)
  if (all_unif_prior) {
    random_tab = randomLHS(nb_simul, nparam)
  }

  list_param = list(NULL)
  for (i in 1:nb_simul) {
    param = array(0, l)
    if (!all_unif_prior) {
      param = .sample_prior(prior, prior_test)
    } else {
      for (j in 1:l) {
        param[j] = as.numeric(prior[[j]]$sampleArgs[2]) +
          (as.numeric(prior[[j]]$sampleArgs[3]) - as.numeric(prior[[j]]$sampleArgs[2])) *
          random_tab[i, j]
      }
    }
    param = c((seed_count + i), param)
    list_param[[i]] = param
    tab_param = rbind(tab_param, param[2:(l + 1)])
  }
  seed_count = seed_count + nb_simul
  list_simul_summarystat = parLapplyLB(cl, list_param, model)
  for (i in 1:nb_simul) {
    tab_simul_summarystat = rbind(tab_simul_summarystat, as.numeric(list_simul_summarystat[[i]]))
  }

  options(scipen = 0)
  cbind(tab_param, tab_simul_summarystat)
}


## function to perform ABC simulations from a non-uniform prior (derived from a
## set of particles)
.ABC_launcher_not_uniformc_cluster <- function(model, prior, param_previous_step,
                                               tab_weight, nb_simul, seed_count,
                                               inside_prior, n_cluster, cl, max_pick = 10000) {
  tab_simul_summarystat = NULL
  tab_param = NULL
  k_acc = 0
  list_param = list(NULL)

  list_param = list(NULL)
  for (i in 1:nb_simul) {
    l = dim(param_previous_step)[2]
    counter = 0
    repeat {
      k_acc = k_acc + 1
      counter = counter + 1
      # pick a particle
      param_picked = .particle_pick(param_previous_step, tab_weight)
      # move it
      # only variable parameters are moved, computation of a WEIGHTED variance
      param_moved = .move_particle(as.numeric(param_picked),
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
    param = param_previous_step[1, ]
    param = param_moved
    param = c((seed_count + i), param)
    list_param[[i]] = param
    tab_param = rbind(tab_param, param[2:(l + 1)])
  }
  seed_count = seed_count + nb_simul
  list_simul_summarystat = parLapplyLB(cl, list_param, model)
  for (i in 1:nb_simul) {
    tab_simul_summarystat = rbind(tab_simul_summarystat, as.numeric(list_simul_summarystat[[i]]))
  }

  list(cbind(tab_param, tab_simul_summarystat), nb_simul/k_acc)
}
