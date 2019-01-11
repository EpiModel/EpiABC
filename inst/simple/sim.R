
library("methods")
library("Rmpi")
library("snow")
library("EasyABCMPI")

uni <- mpi.universe.size()-1
print(uni)
cluster <- makeCluster(uni)

toy_model_parallel <- function(x){
  set.seed(x[1])
  2 * x[2] + 5 + rnorm(1,0,0.1) }
sum_stat_obs <- 6.5
toy_prior <- list(c("unif",0,1))

fit <- ABC_sequential(method = "Lenormand",
                      model = toy_model_parallel,
                      prior = toy_prior,
                      nb_simul = 1000,
                      summary_stat_target = sum_stat_obs,
                      p_acc_min = 0.1,
                      use_seed = TRUE,
                      n_cluster = uni,
                      cl = cluster)
save(fit, file = "fit.rda")

stopCluster(cluster)
mpi.exit()
