###   Design II Simulations   ###

library(MASS)
library(doParallel)
library(doRNG)
library(RPushbullet)
# if (!is.loaded("mpi_initialize")) {
#   library("Rmpi")
# }

source("functions.R", echo = FALSE)

# Set up cluster for parallel computation
clus <- makeCluster(10)
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

## Effect size: 0.3
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
      save(file = "simsDesign2-delta3-noViol.RData",
           list = c(grep("d2small", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                    "gammas", "lambdas", "seed", "times", "spltime"))
    }
  }
}

pbPost("note", "All done!",
       paste("All simulations are complete for effect size 0.3, no assumption violations"),
       recipients = c('pixel', 'spectre'))

## Effect size: 0.5
gammas <- c(33.5, -0.8, 1.1, -0.8, 0.8, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

load('simsDesign2-delta5-noViol-test.RData')
for (corr in c(.6, .8)) {
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
      save(file = "simsDesign2-delta5-noViol.RData",
           list = c(grep("d2med", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                    "gammas", "lambdas", "seed", "times", "spltime"))
    }
  }
}

pbPost("note", "All done!",
       paste("All simulations are complete for effect size 0.5, no assumption violations"),
       recipients = c('pixel', 'spectre'))


# ##### Violation of S1(a) #####
# ### Effect size: 0.3
gammas <- c(33.5, -0.8, 0.9, -0.8, 0.4, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

sigma <- 6
sigma.r1 <- 5.2
sigma.r0 <- 5.2

lapply(list(c(0, 0, 0), c(0.3, 0.31, .32), c(.6, .62, .63), c(.8, .82, .83)), function(corr) {
  for (resp in c(0.4, 0.6)) {
    for (sharp in c(FALSE, TRUE)) {
      r1 <- r0 <- resp
      set.seed(seed)
      assign(paste0("d2small.r", resp * 10, ".exch", corr * 10, ifelse(sharp, ".sharp", ""), ".viol1"),
             sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                 alpha = .05, power = .8, delta = 0.3, design = 2, conservative = !sharp,
                 sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                 constant.var.time = FALSE, L.eos = c(0, 0, 2, 0, 2, 2, 0),
                 corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3],
                 niter = niter, notify = TRUE, pbDevice = c("pixel", "spectre"),
                 pbIdentifier = paste0("Violate S1(a)\n",
                                       ifelse(sharp, "sharp n", "conservative n"))), envir = .GlobalEnv)
      save(file = "simsDesign2-delta3-viol1.RData",
           list = c(grep("d2small", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                    "gammas", "lambdas", "seed", "times", "spltime"))
    }
  }
})

pbPost("note", "All done!",
       paste("All simulations are complete for effect size 0.3 with the violation of S1(a)."),
       recipients = c('pixel', 'spectre'))

sigma <- 6
sigma.r1 <- 4.9
sigma.r0 <- 4.9

lapply(list(c(0, 0, 0), c(0.3, 0.37, .37), c(.6, .74, .74), c(.8, .98, .98)), function(corr) {
  for (resp in c(0.4, 0.6)) {
    for (sharp in c(FALSE, TRUE)) {
      r1 <- r0 <- resp
      set.seed(seed)
      assign(paste0("d2small.r", resp * 10, ".exch", corr * 10, ifelse(sharp, ".sharp", ""), ".viol1"),
             sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                 alpha = .05, power = .8, delta = 0.3, design = 2, conservative = !sharp,
                 sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                 constant.var.time = FALSE, L.eos = c(0, 0, 2, 0, 2, 2, 0),
                 corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3],
                 niter = niter, notify = TRUE, pbDevice = c("pixel", "spectre"),
                 pbIdentifier = paste0("Violate S1(a)\n",
                                       ifelse(sharp, "sharp n", "conservative n"))), envir = .GlobalEnv)
      save(file = "simsDesign2-delta3-viol2.RData",
           list = c(grep("d2small", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                    "gammas", "lambdas", "seed", "times", "spltime"))
    }
  }
})

pbPost("note", "All done!",
       paste("All simulations are complete for effect size 0.3 with the violation of S1(a)."),
       recipients = c('pixel', 'spectre'))

stopCluster(clus)
