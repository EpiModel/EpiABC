
library("methods")
library("Rmpi")
library("snow")
library("EasyABCMPI")

uni <- mpi.universe.size()-1
print(uni)
uni = 8
cluster <- makeCluster(uni)
cluster

toy_model_parallel <- function(x){
  set.seed(x[1])
  a <- rnorm(1e7)
  2 * x[2] + 5 + rnorm(1,mean(a),0.1)
}
sum_stat_obs <- 6.5
toy_prior <- list(c("unif",0,1))

system.time({
  fit <- ABC_sequential(method = "Lenormand",
                        model = toy_model_parallel,
                        prior = toy_prior,
                        nb_simul = 100,
                        summary_stat_target = sum_stat_obs,
                        p_acc_min = 0.2,
                        use_seed = TRUE,
                        n_cluster = uni,
                        cl = cluster)
})

# save(fit, file = "fit.rda")

stopCluster(cluster)
# mpi.exit()

# uni 2 st 110, 142
# uni 4 st 86, 68
# uni 8 st 59, 52
