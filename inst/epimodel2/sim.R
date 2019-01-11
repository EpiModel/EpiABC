
n <- 1000
nw <- network.initialize(n = n, directed = FALSE)
formation <- ~edges
target.stats <- 1*(n/2)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
save("est", file = "temp.est.rda")

dx <- netdx(est, nsims = 10, nsteps = 500, ncores = 4, verbose = FALSE)
dx

plot(dx)

myfunc <- function(x) {
  set.seed(x[1])
  require(EpiModel)
  load("temp.est.rda")
  est$coef.form <- est$coef.form + x[2]
  dx <- netdx(est, nsims = 1, nsteps = 200, verbose = FALSE)
  out <- mean(tail(dx$stats[[1]][, "edges"], 50))
  return(out)
}

priors <- list(c("unif", -0.2, 0))
prev.targ <- 500

a <- ABC_sequential(method = "Lenormand",
                    model = myfunc,
                    prior = priors,
                    nb_simul = 200,
                    summary_stat_target = prev.targ,
                    p_acc_min = 0.02,
                    progress_bar = TRUE,
                    verbose = FALSE,
                    n_cluster = 16,
                    use_seed = TRUE)

a
save(a, file = "carnegie.rda")

load("carnegie.rda")
par(mfrow = c(1, 2))
plot(density(a$param))
plot(density(a$stats))

coef.adj <- sum(a$param * a$weights)
coef.adj

est$coef.form <- est$coef.form + coef.adj
dx <- netdx(est, nsims = 10, nsteps = 500, ncores = 4, verbose = FALSE)
dx

par(mfrow = c(1,1))
plot(dx)
