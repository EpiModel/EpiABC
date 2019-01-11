
library("methods")
library("Rmpi")
library("snow")
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
  control <- control.net(type = "SIS", nsteps = 250, nsims = 1, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  df <- as.data.frame(mod)
  out <- mean(tail(df$i.num/df$num, control$nsteps/10))
  return(out)
}

priors <- list(c("unif", 0.2, 0.4),
               c("unif", 0.05, 0.25))
prev.targ <- 0.25

uni <- mpi.universe.size() - 1
cluster <- makeCluster(uni)

a <- ABC_sequential(method = "Lenormand",
                    model = myfunc,
                    prior = priors,
                    nb_simul = 200,
                    summary_stat_target = prev.targ,
                    p_acc_min = 0.05,
                    n_cluster = uni,
                    use_seed = TRUE,
                    cl = cluster)
save(a, file = "sisfit.rda")

stopCluster(cluster)
mpi.exit()



library("snow")
suppressMessages(library("EpiModel"))
suppressMessages(library("EasyABCMPI"))

uni = 2
uni = parallel::detectCores()
uni
# cluster <- makeCluster(uni)
a <- EasyABC::ABC_sequential(method = "Lenormand",
                    model = myfunc,
                    prior = priors,
                    nb_simul = 200,
                    summary_stat_target = prev.targ,
                    p_acc_min = 0.05,
                    n_cluster = uni,
                    progress_bar = TRUE,
                    use_seed = TRUE)
# stopCluster(cluster)


nb_simul <- 200
n_cluster <- 95
( npar = floor(nb_simul/(100 * n_cluster)) )
( n_end = nb_simul - (npar * 100 * n_cluster) )
