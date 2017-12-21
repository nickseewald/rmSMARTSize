###   Design II Simulations   ###
### All assumptions satisifed ###

library(MASS)
library(doParallel)
library(doRNG)
library(RPushbullet)
# if (!is.loaded("mpi_initialize")) {
#   library("Rmpi")
# }

source("functions.R", echo = FALSE)

# Set up cluster for parallel computation
clus <- makeCluster(8)
# clusterEvalQ(clus, source("longitudinal-simulations-auxfunctions.R"))
clusterEvalQ(clus, {
  source("functions.R")
  library(MASS)
  library(RPushbullet)
  })
registerDoParallel(clus)

##### Simulation Setup #####

times <- c(0, 1, 2)
spltime <- 1

niter <- 5000
maxiter.solver <- 1000
tol <- 1e-8
seed <- 1001

pbPost("note", "Starting Simulations",
       paste("All systems go so far!\nRunning on", length(clus), "cores."),
       recipients = c('pixel', 'spectre'))

##### All assumptions satisfied #####
sigma <- 6
sigma.r1 <- 6.1
sigma.r0 <- 6.1

### Effect size: 0.3
gammas <- c(33.5, -0.8, 0.9, -0.8, 0.4, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

for (corr in c(0, .3, .6, .8)) {
  for (resp in c(0.4, 0.6)) {
    for (sharp in c(FALSE, TRUE)) {
      r1 <- r0 <- resp
      set.seed(seed)
      assign(paste0("d2small.r", resp * 10, ".exch", corr * 10, ifelse(sharp, ".sharp", "")),
             sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                 alpha = .05, power = .8, delta = 0.3, design = 2, conservative = !sharp,
                 sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                 constant.var.time = FALSE, L.eos = c(0, 0, 2, 0, 2, 2, 0),
                 corstr = "exch", rho = corr, niter = niter, notify = TRUE, pbDevice = c("pixel", "spectre"),
                 pbIdentifier = paste0("all assumptions ok\n",
                                       ifelse(sharp, "sharp n", "conservative n"))), envir = .GlobalEnv)
      save.image(file = "simsDesign2.RData")
    }
  }
}

pbPost()

### Effect size: 0.5
gammas <- c(33.5, -0.8, 1.1, -0.8, 0.8, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

for (corr in c(0, .3, .6, .8)) {
  for (resp in c(0.4, 0.6)) {
    for (sharp in c(FALSE, TRUE)) {
      r1 <- r0 <- resp
      set.seed(seed)
      assign(paste0("d2med.r", resp * 10, ".exch", corr * 10, ifelse(sharp, ".sharp", "")), 
             sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                 alpha = .05, power = .8, delta = 0.5, design = 2, conservative = !sharp,
                 sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                 constant.var.time = FALSE, L.eos = c(0, 0, 2, 0, 2, 2, 0),
                 corstr = "exch", rho = corr, niter = niter, notify = TRUE, pbDevice = c("pixel", "spectre"),
                 pbIdentifier = paste0("all assumptions ok\n",
                                       ifelse(sharp, "sharp n", "conservative n"))), envir = .GlobalEnv)
      save.image(file = "simsDesign2_med.RData")
    }
  }
}

pbPost("note", "All done!",
       paste("All simulations are complete for effect size 0.5."),
       recipients = c('pixel', 'spectre'))

stopCluster(clus)
