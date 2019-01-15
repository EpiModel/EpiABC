
library("methods")
suppressMessages(library("EpiABC"))
suppressMessages(library("EpiModel"))

prep <- readRDS("data/abc.prep.rda")

batch <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
wave <- as.numeric(Sys.getenv("wave"))

# Wave 0
if (wave == 0) {
  abc_smc_wave(input = prep, wave = wave, batch = batch)
} else {
  abc_smc_wave(wave = wave, batch = batch)
}

merge_abc(wave = wave)
abc_smc_process(wave = wave)
