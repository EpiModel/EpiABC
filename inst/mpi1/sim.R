
library("methods")
library("snow")
library("EpiABC")

if (!require("Rmpi", character.only = TRUE)) {
  uni <- parallel::detectCores() - 1
} else {
  uni <- mpi.universe.size()-1
}

cluster <- makeCluster(uni)
cluster

toy_model_parallel <- function(x){
  set.seed(x[1])
  a <- rnorm(1e6)
  2 * x[2] + 5 + rnorm(1,mean(a), 0.1)
}
sum_stat_obs <- 6.5
toy_prior <- list(c("unif", 0, 1))

fit <- abc_smc_cluster(model = toy_model_parallel,
                       prior = toy_prior,
                       summary_stat_target = sum_stat_obs,
                       nsims = 200,
                       alpha = 0.5,
                       p_acc_min = 0.45,
                       cl = cluster)

stopCluster(cluster)

