# SDE HMM example 
library(smoothSDE)

formulas <- list(mu = ~ s(x, bs = "cs", k = 5), 
                 sigma = ~ s(x, bs = "cs", k = 5))

# format data 
data <- data.frame(ID = 1, 
                   y = do.call(c, sim_Bm_object$y), 
                   coarse_time = rep(1:50, each = 200),
                   x = rep(seq(0, 1, length = 200), 50),
                   time = 1:(200 * 50))

type <- "BM_HMM"
my_sde <- SDE_HMM$new(formulas = formulas, data = data, type = type, response = "y", nstates = 2)
my_sde$fit(silent = FALSE)
my_sde$par(t = 10, state = 1)
my_sde$par(t = 10, state = 2)

par(mfrow=c(2,2))
plot(data$time[1:200], my_sde$par(t = 1:200, state = 1)[,1], type = "l", lwd = 1.5, xlab = "t", ylab = "mu state 1")
plot(data$time[1:200], my_sde$par(t = 1:200, state = 2)[,1], type = "l", lwd = 1.5, xlab = "t", ylab = "mu state 2")
plot(data$time[1:200], my_sde$par(t = 1:200, state = 1)[,2], type = "l", lwd = 1.5, xlab = "t", ylab = "sd state 1")
plot(data$time[1:200], my_sde$par(t = 1:200, state = 2)[,2], type = "l", lwd = 1.5, xlab = "t", ylab = "sd state 2")


