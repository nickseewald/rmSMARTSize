###   Design II Simulations   ###
#####   Violation of S1(a)  #####
library(MASS)
library(doParallel)
library(doRNG)
library(slackr)
# if (!is.loaded("mpi_initialize")) {
#   library("Rmpi")
# }

slackrSetup(config_file = "d3slack.dcf")
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

startString <- paste("All systems go so far!\nRunning on", length(clus), "cores.")

slackr_bot(startString)

### Effect size: 0.3
gammas <- c(33.5, -0.8, 0.9, -0.8, 0.4, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

sigma <- 6
sigma.r1 <- 5.9
sigma.r0 <- 5.75

save(file = "simsDesign2-delta3-violS1a-1.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

lapply(list(c(0, 0, 0), c(0.3, 0.31, .32), c(.6, .62, .63), c(.8, .82, .83)), function(corr) {
  for (resp in c(0.4, 0.6)) {
    for (sharp in c(FALSE, TRUE)) {
      r1 <- r0 <- resp
      set.seed(seed)
      assign(paste0("d2small.r", resp * 10, ".exch", corr * 10, ".viol1", ifelse(sharp, ".sharp", "")),
             sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                 alpha = .05, power = .8, delta = 0.3, design = 2, conservative = !sharp,
                 sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                 constant.var.time = FALSE, L.eos = c(0, 0, 2, 0, 2, 2, 0),
                 corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3],
                 niter = niter, notify = TRUE,
                 postIdentifier = paste0("Violate S1(a)\n",
                                       ifelse(sharp, "sharp n", "conservative n"))), envir = .GlobalEnv)
      save(file = "simsDesign2-delta3-violS1a-1.RData",
           list = c(grep("^d2small.*viol1", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                    "gammas", "lambdas", "seed", "times", "spltime"),
           safe = TRUE, precheck = TRUE)
    }
  }
})

slackr_bot("All simulations are complete for effect size 0.3 with the violation of S1(a).")

sigma <- 6
sigma.r1 <- 5.2
sigma.r0 <- 5.2

save(file = "simsDesign2-delta3-violS1a-2.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

lapply(list(c(0, 0, 0), c(0.3, 0.35, .35), c(.6, .7, .7), c(.8, .94, .94)), function(corr) {
  for (resp in c(0.4, 0.6)) {
    for (sharp in c(FALSE, TRUE)) {
      r1 <- r0 <- resp
      set.seed(seed)
      assign(paste0("d2small.r", resp * 10, ".exch", corr * 10, ".viol2", ifelse(sharp, ".sharp", "")),
             sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                 alpha = .05, power = .8, delta = 0.3, design = 2, conservative = !sharp,
                 sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                 constant.var.time = FALSE, L.eos = c(0, 0, 2, 0, 2, 2, 0),
                 corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3],
                 niter = niter, notify = TRUE, pbDevice = c("pixel", "spectre"),
                 pbIdentifier = paste0("Violate S1(a)\n",
                                       ifelse(sharp, "sharp n", "conservative n"))), envir = .GlobalEnv)
      save(file = "simsDesign2-delta3-violS1a-2.RData",
           list = c(grep("^d2small.*viol2", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                    "gammas", "lambdas", "seed", "times", "spltime"),
           safe = TRUE, precheck = TRUE)
    }
  }
})

slackr_bot("All simulations are complete for effect size 0.3 with the stronger violation of S1(a).")

