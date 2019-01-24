
library("methods")
library("Rmpi")
library("snow")
library("EpiABC")

uni <- mpi.universe.size()-1
print(uni)
uni = 4
cluster <- makeCluster(uni)
cluster

toy_model_parallel <- function(x){
  set.seed(x[1])
  a <- rnorm(1e7)
  2 * x[2] + 5 + rnorm(1,mean(a),0.1)
}
sum_stat_obs <- 6.5
toy_prior <- list(c("unif",0,1))

fit <- ABC_sequential(method = "Lenormand",
                      model = toy_model_parallel,
                      prior = toy_prior,
                      summary_stat_target = sum_stat_obs,
                      nb_simul = 200,
                      alpha = 0.5,
                      p_acc_min = 0.45,
                      use_seed = TRUE,
                      n_cluster = uni,
                      cl = cluster)

stopCluster(cluster)

