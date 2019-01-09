
sum_stat_obs <- c(100, 2.5, 20, 30000)
toy_model_parallel <- function(x){ 
  set.seed(x[1])
  2 * x[2] + 5 + rnorm(1,0,0.1) }
sum_stat_obs <- 6.5
toy_prior <- list(c("unif",0,1))

fit <- ABC_sequential(method = "Lenormand", 
                      model = toy_model_parallel, 
                      prior = toy_prior,
                      nb_simul = 100, 
                      summary_stat_target = sum_stat_obs, 
                      p_acc_min = 0.05, 
                      use_seed = TRUE, 
                      n_cluster = 95)
fit
hist(fit$param, breaks = 20)
hist(fit$stats, breaks = 20)

# Stack
ABC_sequential
.ABC_sequential_cluster
.ABC_Lenormand_cluster

.ABC_rejection_lhs_cluster # first step
cl <- makeCluster(getOption("cl.cores", n_cluster))
parLapplyLB(cl, list_param, model)
stopCluster(cl)

.ABC_launcher_not_uniformc_cluster
makeCluster(getOption("cl.cores", n_cluster))
parLapplyLB(cl, list_param, model)
stopCluster(cl)

undebug(ABC_sequential)
debug(EasyABC:::.ABC_Lenormand_cluster)

library(snow)

model <- function(x) {
  a <- runif(1e8)*x
  return(mean(a))
}

ncores <- parallel::detectCores()/2
system.time({
  cl <- makeCluster(ncores)
  print(cl)
  d <- parLapplyLB(cl, 1:11, model)
  stopCluster(cl)
})
d
