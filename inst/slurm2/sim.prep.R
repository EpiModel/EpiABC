
library("methods")
suppressMessages(library("EpiABC"))
suppressMessages(library("EpiModel"))


# Model Setup -------------------------------------------------------------

n <- 1000
nw <- network.initialize(n = n, directed = FALSE)
formation <- ~edges
target.stats <- 0.75*(n/2)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
save(est, file = "inst/slurm2/est.rda")



# Main Model Fx -----------------------------------------------------------

myfunc <- function(x) {
  set.seed(x[1])
  require(EpiModel)
  load("inst/slurm2/est.rda")
  param <- param.net(inf.prob = x[2], rec.rate = x[3])
  init <- init.net(i.num = 50, status.rand = FALSE)
  control <- control.net(type = "SIS", nsteps = 300, nsims = 1, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  df <- as.data.frame(mod)
  out <- mean(tail(df$i.num/df$num, control$nsteps/10))
  return(out)
}


# ABC Priors and Target Stats ---------------------------------------------

priors <- list(c("unif", 0.2, 0.4),
               c("unif", 0.05, 0.25))
prev.targ <- 0.25



# Run ABC Prep ------------------------------------------------------------

prep <- abc_smc_prep(model = myfunc,
                     prior = priors,
                     nsims = 500,
                     summary_stat_target = prev.targ,
                     ncores = 16,
                     alpha = 0.2)
prep
saveRDS(prep, file = "inst/slurm2/data/abc.prep.rda")

debug(sbatch_master_abc)
sbatch_master_abc(prep,
                  master.file = "inst/slurm2/master.sh",
                  nwaves = 10, ckpt = TRUE)

sbatch_master_abc(prep,
                  master.file = "inst/slurm2/master.sh",
                  nwaves = 20, ckpt = TRUE, append = TRUE)

