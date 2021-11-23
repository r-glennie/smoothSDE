# SDE HMM example 
library(smoothSDE)

formulas <- list(mu = ~ 1, 
                 sigma = ~ 1)

# format data 
data <- data.frame(ID = 1, 
                   y = do.call(c, sim_Bm_object$y), 
                   coarse_time = rep(1:50, each = 200), 
                   time = 1:(200 * 50))

type <- "BM_HMM"
my_sde <- SDE_HMM$new(formulas = formulas, data = data, type = type, response = "y", nstates = 2)
my_sde$fit(silent = FALSE)
my_sde$plot_par("x1", n_post = 100)


