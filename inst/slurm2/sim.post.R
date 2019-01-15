
library("methods")
suppressMessages(library("EpiABC"))
suppressMessages(library("EpiModel"))

prep <- readRDS(file = "data/abc.prep.rda")

out <- out_abc(wave = 5)
summary_abc(out)
