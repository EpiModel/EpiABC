
system("scp hyak:~/sisfit.rda inst/epimodel1")

load("inst/epimodel1/sisfit.rda")
plot(density(a$stats))

par(mfrow = c(1,2))
plot(density(a$param[, 1]))
plot(density(a$param[, 2]))

library(MASS)
kde1 <- kde2d(a$param[,1], a$param[,2], n = 250)
par(mar = c(3,3,2,1), mgp = c(2,1,0), mfrow = c(1,1))
image(kde1, col = heat.colors(100), ylab = "Recovery Rate",
      xlab = "Infection Probability")

kde1.max <- kde1$z == max(kde1$z)
kde1.row <- which.max(rowSums(kde1.max))
kde1.col <- which.max(colSums(kde1.max))
bestfit1 <- c(kde1$x[kde1.row], kde1$y[kde1.col])
bestfit1
