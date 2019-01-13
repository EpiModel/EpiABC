

library("methods")
suppressMessages(library("EasyABCMPI"))
suppressMessages(library("EpiModel"))

n <- 1000
nw <- network.initialize(n = n, directed = FALSE)
formation <- ~edges
target.stats <- 0.75*(n/2)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
save(est, file = "est.rda")

myfunc <- function(x) {
  set.seed(x[1])
  require(EpiModel)
  load("est.rda")
  param <- param.net(inf.prob = x[2], rec.rate = x[3])
  init <- init.net(i.num = 50, status.rand = FALSE)
  control <- control.net(type = "SIS", nsteps = 200, nsims = 1, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  df <- as.data.frame(mod)
  out <- mean(tail(df$i.num/df$num, control$nsteps/10))
  return(out)
}

priors <- list(c("unif", 0.2, 0.4),
               c("unif", 0.05, 0.25))
prev.targ <- 0.25

prep <- abc_smc_prep(model = myfunc,
                     prior = priors,
                     nb_simul = 20,
                     summary_stat_target = prev.targ,
                     n_cluster = 4,
                     alpha = 0.5,
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
