
library("methods")
library("EasyABCMPI")

toy_model_parallel <- function(x){
  set.seed(x[1])
  a <- rnorm(1e7)
  # a <- 0
  2 * x[2] + 5 + rnorm(1,mean(a),0.1)
}
sum_stat_obs <- 6.5
toy_prior <- list(c("unif",0,1))

prep <- abc_smc_prep(model = toy_model_parallel,
                     prior = toy_prior,
                     nb_simul = 500,
                     summary_stat_target = sum_stat_obs,
                     n_cluster = 4,
                     alpha = 0.5,
                     p_acc_min = 0.45)
prep

wavedist <- prep
wave <- 0
for (i in 1:20) {
  wavedat <- abc_smc_wave(input = wavedist, wave = wave, batch = 1)
  wavedist <- abc_smc_process(input = wavedat, wave = wave)
  cat("\n", wavedist$pwave$p_acc)
  wave <- wave + 1
}
