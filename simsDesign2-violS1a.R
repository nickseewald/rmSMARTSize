###   Design II Simulations   ###
#####   Violation of S1(a)  #####
library(MASS)
library(doParallel)
library(doRNG)
library(slackr)
# if (!is.loaded("mpi_initialize")) {
#   library("Rmpi")
# }

notify <- TRUE

if(notify) {
  slackrSetup(config_file = "d3slack.dcf")
  source("functions.R", echo = FALSE)
}

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

sigma <- 6
sigma.r0 <- 5.92
sigma.r1 <- 5.8

simGrid <- expand.grid(list(
  sharp = c(FALSE, TRUE),
  resp = c(.2, .4, .6, .8),
  corr = list(rep(0, 3),
              c(0.3, 0.3 * sigma / sigma.r1, 0.3 * sigma / sigma.r0),
              c(0.6, 0.6 * sigma / sigma.r1, 0.6 * sigma / sigma.r0),
              c(0.8, 0.8 * sigma / sigma.r1, 0.8 * sigma / sigma.r0))
))

startString <- paste("All systems go so far!\nRunning on", length(clus), "cores.")

if(notify)
  slackr_bot(startString)

### Effect size: 0.3
gammas <- c(33.5, -0.8, 0.9, -0.8, 0.4, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

save(file = "simsDesign2-delta3-violS1a-low.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

for (i in 1:nrow(simGrid)) {
  # Extract simulation conditions from simGrid
  resp <- r1 <- r0 <- simGrid$resp[i]
  sharp <- simGrid$sharp[i]
  corr <- simGrid$corr[[i]]
  set.seed(seed)
  assign(paste0("d2small.r", resp * 10, ".exch", corr[1] * 10, ifelse(sharp, ".sharp", "")),
         try(sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                 alpha = .05, power = .8, delta = 0.3, design = 2, conservative = !sharp,
                 sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                 constant.var.time = FALSE, L = c(0, 0, 2, 0, 2, 2, 0),
                 corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3], 
                 niter = niter, notify = notify, pbDevice = c("pixel", "spectre"),
                 postIdentifier = paste0("Violation of S1(a)\n",
                                         ifelse(sharp, "sharp n", "conservative n")))), envir = .GlobalEnv)
  save(file = "simsDesign2-delta3-violS1a-low.RData",
       list = c(grep("d2small", ls(), value = T), "sigma", "sigma.r1", "sigma.r0", "corr",
                "gammas", "lambdas", "seed", "times", "spltime"), precheck = TRUE)
}

### Effect size 0.5
gammas <- c(33.5, -0.8, 1.1, -0.8, 0.8, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

save(file = "simsDesign2-delta5-violS1a-low.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

for (i in 1:nrow(simGrid)) {
  # Extract simulation conditions from simGrid
  resp <- r1 <- r0 <- simGrid$resp[i]
  sharp <- simGrid$sharp[i]
  corr <- simGrid$corr[[i]]
  set.seed(seed)
  assign(paste0("d2med.r", resp * 10, ".exch", corr * 10, ifelse(sharp, ".sharp", "")),
         try(sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                 alpha = .05, power = .8, delta = 0.5, design = 2, conservative = !sharp,
                 sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                 constant.var.time = FALSE, L = c(0, 0, 2, 0, 2, 2, 0),
                 corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3], 
                 niter = niter, notify = notify, pbDevice = c("pixel", "spectre"),
                 postIdentifier = paste0("Violation of S1(a)\n",
                                         ifelse(sharp, "sharp n", "conservative n")))), envir = .GlobalEnv)
  save(file = "simsDesign2-delta5-violS1a-low.RData",
       list = c(grep("d2med", ls(), value = T), "sigma", "sigma.r1", "sigma.r0", "corr",
                "gammas", "lambdas", "seed", "times", "spltime"), precheck = TRUE)
}

slackr_bot("All simulations are complete for effect size 0.5 with small violation of S1(a).")

### Effect size 0.8
gammas <- c(33.5, -0.8, 1.5, -0.8, 1.3, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

save(file = "simsDesign2-delta8-violS1a-low.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

for (i in 1:nrow(simGrid)) {
  # Extract simulation conditions from simGrid
  resp <- r1 <- r0 <- simGrid$resp[i]
  sharp <- simGrid$sharp[i]
  corr <- simGrid$corr[[i]]
  set.seed(seed)
  assign(paste0("d2large.r", resp * 10, ".exch", corr * 10, ifelse(sharp, ".sharp", "")),
         try(sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                 alpha = .05, power = .8, delta = 0.8, design = 2, conservative = !sharp,
                 sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                 constant.var.time = FALSE, L = c(0, 0, 2, 0, 2, 2, 0),
                 corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3], 
                 niter = niter, notify = notify, pbDevice = c("pixel", "spectre"),
                 postIdentifier = paste0("Violation of S1(a)\n",
                                         ifelse(sharp, "sharp n", "conservative n")))), envir = .GlobalEnv)
  save(file = "simsDesign2-delta8-violS1a-low.RData",
       list = c(grep("d2large", ls(), value = T), "sigma", "sigma.r1", "sigma.r0", "corr",
                "gammas", "lambdas", "seed", "times", "spltime"), precheck = TRUE)
}

slackr_bot("All simulations are complete for effect size 0.8 with small violation of S1(a).")


##### Large violation #####
sigma <- 6
sigma.r1 <- sigma.r0 <- 5

simGrid <- expand.grid(list(
  sharp = c(FALSE, TRUE),
  resp = c(.2, .4, .6, .8),
  corr = list(rep(0, 3),
              c(0.3, 0.3 * sigma / sigma.r1, 0.3 * sigma / sigma.r0),
              c(0.6, 0.6 * sigma / sigma.r1, 0.6 * sigma / sigma.r0),
              c(0.8, 0.8 * sigma / sigma.r1, 0.8 * sigma / sigma.r0))
))

### Effect size: 0.3
gammas <- c(33.5, -0.8, 0.9, -0.8, 0.4, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

save(file = "simsDesign2-delta3-violS1a-high.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

for (i in 1:nrow(simGrid)) {
  # Extract simulation conditions from simGrid
  resp <- r1 <- r0 <- simGrid$resp[i]
  sharp <- simGrid$sharp[i]
  corr <- simGrid$corr[[i]]
  set.seed(seed)
  assign(paste0("d2small.r", resp * 10, ".exch", corr * 10, ifelse(sharp, ".sharp", "")),
         try(sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                 alpha = .05, power = .8, delta = 0.3, design = 2, conservative = !sharp,
                 sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                 constant.var.time = FALSE, L = c(0, 0, 2, 0, 2, 2, 0),
                 corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3], 
                 niter = niter, notify = notify, pbDevice = c("pixel", "spectre"),
                 postIdentifier = paste0("Violation of S1(a)\n",
                                         ifelse(sharp, "sharp n", "conservative n")))), envir = .GlobalEnv)
  save(file = "simsDesign2-delta3-violS1a-high.RData",
       list = c(grep("d2small", ls(), value = T), "sigma", "sigma.r1", "sigma.r0", "corr",
                "gammas", "lambdas", "seed", "times", "spltime"), precheck = TRUE)
}

### Effect size 0.5
gammas <- c(33.5, -0.8, 1.1, -0.8, 0.8, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

save(file = "simsDesign2-delta5-violS1a-high.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

for (i in 1:nrow(simGrid)) {
  # Extract simulation conditions from simGrid
  resp <- r1 <- r0 <- simGrid$resp[i]
  sharp <- simGrid$sharp[i]
  corr <- simGrid$corr[[i]]
  set.seed(seed)
  assign(paste0("d2med.r", resp * 10, ".exch", corr * 10, ifelse(sharp, ".sharp", "")),
         try(sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                 alpha = .05, power = .8, delta = 0.5, design = 2, conservative = !sharp,
                 sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                 constant.var.time = FALSE, L = c(0, 0, 2, 0, 2, 2, 0),
                 corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3], 
                 niter = niter, notify = notify, pbDevice = c("pixel", "spectre"),
                 postIdentifier = paste0("Violation of S1(a)\n",
                                         ifelse(sharp, "sharp n", "conservative n")))), envir = .GlobalEnv)
  save(file = "simsDesign2-delta5-violS1a-high.RData",
       list = c(grep("d2med", ls(), value = T), "sigma", "sigma.r1", "sigma.r0", "corr",
                "gammas", "lambdas", "seed", "times", "spltime"), precheck = TRUE)
}

slackr_bot("All simulations are complete for effect size 0.5 with large violation of S1(a).")

### Effect size 0.8
gammas <- c(33.5, -0.8, 1.5, -0.8, 1.3, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

save(file = "simsDesign2-delta8-violS1a-high.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

for (i in 1:nrow(simGrid)) {
  # Extract simulation conditions from simGrid
  resp <- r1 <- r0 <- simGrid$resp[i]
  sharp <- simGrid$sharp[i]
  corr <- simGrid$corr[[i]]
  set.seed(seed)
  assign(paste0("d2large.r", resp * 10, ".exch", corr * 10, ifelse(sharp, ".sharp", "")),
         try(sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                 alpha = .05, power = .8, delta = 0.8, design = 2, conservative = !sharp,
                 sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                 constant.var.time = FALSE, L = c(0, 0, 2, 0, 2, 2, 0),
                 corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3], 
                 niter = niter, notify = notify, pbDevice = c("pixel", "spectre"),
                 postIdentifier = paste0("Violation of S1(a)\n",
                                         ifelse(sharp, "sharp n", "conservative n")))), envir = .GlobalEnv)
  save(file = "simsDesign2-delta8-violS1a-high.RData",
       list = c(grep("d2large", ls(), value = T), "sigma", "sigma.r1", "sigma.r0", "corr",
                "gammas", "lambdas", "seed", "times", "spltime"), precheck = TRUE)
}

slackr_bot("All simulations are complete for effect size 0.8 with large violation of S1(a).")


stopCluster(clus)