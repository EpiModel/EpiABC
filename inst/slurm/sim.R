
library("methods")
library("EpiABC")

rm(list = ls())
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
                     nb_simul = 500,
                     summary_stat_target = sum_stat_obs,
                     n_cluster = 4,
                     alpha = 0.2,
                     p_acc_min = 0.1)
prep

wavedist <- prep
wave <- 0
for (i in 1:20) {
  wavedat <- abc_smc_wave(input = wavedist, wave = wave, batch = NULL)
  wavedist <- abc_smc_process(input = wavedat, wave = wave)
  cat("\n", wavedist$pwave$p_acc)
  wave <- wave + 1
}




# batch workflow ----------------------------------------------------------

# interactive
prep <- abc_smc_prep(model = toy_model_parallel,
                     prior = toy_prior,
                     nb_simul = 100,
                     summary_stat_target = sum_stat_obs,
                     n_cluster = 8,
                     alpha = 0.5,
                     p_acc_min = 0.1)
prep

# write bash script with this
nBatches <- ceiling(prep$nb_simul/prep$n_cluster)
nBatches

# run in batch model
for (batch in 1:nBatches) {
  wavedat <- abc_smc_wave(input = prep, wave = 0, batch = batch, save = TRUE, outdir = "dat/")
  cat("\n Batch:", batch)
}

# interactive, or could run in batch mode by looking up file size done
merge_abc(wave = 0, indir = "dat/", outdir = "dat/")
abc_smc_process(input = "dat/", wave = 0, save = TRUE, outdir = "dat/")

# write a new batch script
nBatches <- ceiling(wavedat$init$nb_simul_step/wavedat$init$n_cluster)

# run in batch mode
for (batch in 1:nBatches) {
  wavedat <- abc_smc_wave(input = wavedist, wave = 1, batch = batch, save = TRUE, outdir = "dat/")
  cat("\n Batch:", batch)
}

# interactive
merge_abc(wave = 1, indir = "dat/", outdir = "dat/")
wavedist <- abc_smc_process(input = "dat/", wave = 1)
wavedist$pwave$p_acc

# run in batch mode
for (batch in 1:nBatches) {
  wavedat <- abc_smc_wave(input = wavedist, wave = 2, batch = batch)
  saveRDS(wavedat, file = paste0("dat/abc.wave2.batch", stringr::str_pad(batch, 5, pad = "0"), ".rda"))
  cat("\n Batch:", batch)
}

# interactive
merge_abc(wave = 2, indir = "dat/", outdir = "dat/")
wavedat <- readRDS("dat/abc.wave2.rda")
wavedist <- abc_smc_process(input = wavedat, wave = 1)
wavedist$pwave$p_acc


## TODO:
## batch size = n_cluster # done
## abc_smc_merge fx to combine output from abc_smc_wave, following same order
## abc_smc_process runs on merged data
## master sim script does smart file names based on batch
## helper function to write sbatch scripts?

