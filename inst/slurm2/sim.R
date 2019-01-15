
library("methods")
suppressMessages(library("EpiABC"))
suppressMessages(library("EpiModel"))

prep <- readRDS("data/abc.prep.rda")

batch <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(batch)
wave <- as.numeric(Sys.getenv("wave"))
print(wave)

# Wave 0
if (wave == 0) {
  abc_smc_wave(input = prep, wave = 0, batch = batch)
} else {
  abc_smc_wave(wave = 1, batch = batch)
}

merge_abc(wave = wave)
abc_smc_process(wave = wave)
