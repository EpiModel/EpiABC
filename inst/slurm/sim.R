
library("methods")
library("EasyABCMPI")

toy_model_parallel <- function(x){
  set.seed(x[1])
  # a <- rnorm(1e7)
  a <- 0
  2 * x[2] + 5 + rnorm(1,mean(a),0.1)
}
sum_stat_obs <- 6.5
toy_prior <- list(c("unif",0,1))

prep <- abc_smc_prep(model = toy_model_parallel,
                     prior = toy_prior,
                     nb_simul = 1000,
                     summary_stat_target = sum_stat_obs,
                     n_cluster = 4,
                     alpha = 0.2,
                     p_acc_min = 0.1)
prep

wavedist <- prep
wave <- 0
for (i in 1:20) {
  wavedat <- abc_smc_wave(input = wavedist, wave = wave, batch = 1)
  wavedist <- abc_smc_process(input = wavedat, wave = wave)
  cat("\n", wavedist$pwave$p_acc)
  wave <- wave + 1
}

## TODO:
## batch size = n_cluster
## need fx to query sim range from batch number
## split work only in abc_wave0 and abc_waveN to sim range
## abc_smc_merge fx to combine output from abc_smc_wave, following same order
## abc_smc_process runs on merged data
## master sim script does smart file names based on batch
## helper function to write sbatch scripts?
