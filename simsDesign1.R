###    Design I Simulations   ###
### All assumptions satisifed ###

library(MASS)
library(doParallel)
library(doRNG)
library(slackr)

### SET FLAG TO NOTIFY USER OF SIMULATION
### PROGRESS VIA SLACK (requires config.)
### Do not use if not properly configured
notify <- TRUE

if (notify)
  slackrSetup(config_file = "d3slack.dcf")

source("functions.R", echo = FALSE)

# Set up cluster for parallel computation
clus <- makeCluster(15)
clusterEvalQ(clus, {
  source("functions.R")
  library(MASS)
})
registerDoParallel(clus)

##### Simulation Setup #####

times <- c(0, 1, 2)
spltime <- 1

niter <- 5000
maxiter.solver <- 1000
tol <- 1e-8
seed <- 1001

if (notify) {
  startString <- paste("All systems go so far!\nRunning on", length(clus), "cores.")
  slackr_bot(startString)
}

sigma <- 6
sigma.r1 <- 6
sigma.r0 <- 6

simGrid <- expand.grid(list(
  sharp = c(FALSE, TRUE),
  resp = c(.2, .4, .6, .8),
  corr = list(rep(0, 3),
              c(0.3, 0.3 * sigma / sigma.r1, 0.3 * sigma / sigma.r0),
              c(0.6, 0.6 * sigma / sigma.r1, 0.6 * sigma / sigma.r0),
              c(0.8, 0.8 * sigma / sigma.r1, 0.8 * sigma / sigma.r0))
))

#### Effect size: 0.3 #####
gammas <- c(35, -4, 2.2, -1.6, -1.3, 1, -1, 1, 1)
lambdas <- c(1, -2)

save(file = "simsDesign1-delta3-noViol.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

for (i in 1:nrow(simGrid)) {
  # Extract simulation conditions from simGrid
  resp <- r1 <- r0 <- simGrid$resp[i]
  sharp <- simGrid$sharp[i]
  corr <- simGrid$corr[[i]]
  
  # When corr is zero, sharp and conservative sample sizes are identical
  # Skip the one that we label "sharp" to avoid doing the same sim twice
  if (sum(corr == 0) == length(times) & sharp)
    next
  
  # Set the seed for every unique simulation
  set.seed(seed)
  
  # Simulate
  assign(paste0("d1small.r", resp * 10, ".exch", corr[1] * 10, ifelse(sharp, ".sharp", "")),
         sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
             alpha = .05, power = .8, delta = 0.3, design = 1, conservative = !sharp,
             sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
             constant.var.time = FALSE, L = c(0, 0, 2, 0, 2, 2, 0),
             corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3],
             niter = niter, notify = notify, pbDevice = c("pixel", "spectre"),
             postIdentifier = paste0("all assumptions ok\n",
                                     ifelse(sharp, "sharp n", "conservative n"))), envir = .GlobalEnv)
  
  # Save the result
  save(file = "simsDesign1-delta3-noViol.RData",
       list = c(grep("d1small", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                "gammas", "lambdas", "seed", "times", "spltime"), precheck = TRUE)
}

if (notify)
  slackr_bot("All simulations are complete for effect size 0.3, no assumption violations")

## Effect size: 0.5
gammas <- c(33.5, -0.8, 1.1, -0.8, 0.8, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

save(file = "simsDesign1-delta5-noViol.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

for (i in 1:nrow(simGrid)) {
  resp <- r1 <- r0 <- simGrid$resp[i]
  sharp <- simGrid$sharp[i]
  corr <- simGrid$corr[[i]]
  set.seed(seed)
  assign(paste0("d1med.r", resp * 10, ".exch", corr * 10, ifelse(sharp, ".sharp", "")),
         sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
             alpha = .05, power = .8, delta = 0.5, design = 1, conservative = !sharp,
             sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
             constant.var.time = FALSE, L = c(0, 0, 2, 0, 2, 2, 0),
             corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3],
             niter = niter, notify = notify, pbDevice = c("pixel", "spectre"),
             postIdentifier = paste0("all assumptions ok\n",
                                     ifelse(sharp, "sharp n", "conservative n"))), envir = .GlobalEnv)
  save(file = "simsDesign1-delta5-noViol.RData",
       list = c(grep("d1med", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                "gammas", "lambdas", "seed", "times", "spltime"), precheck = TRUE)
}

if (notify)
  slackr_bot("All simulations are complete for effect size 0.5, no assumption violations")

## Effect size: 0.8
gammas <- c(33.5, -0.8, 1.5, -0.8, 1.3, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

save(file = "simsDesign1-delta8-noViol.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

for (i in 1:nrow(simGrid)) {
  # Extract simulation conditions from simGrid
  resp <- r1 <- r0 <- simGrid$resp[i]
  sharp <- simGrid$sharp[i]
  corr <- simGrid$corr[[i]]
  set.seed(seed)
  assign(paste0("d1large.r", resp * 10, ".exch", corr * 10, ifelse(sharp, ".sharp", "")),
         sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
             alpha = .05, power = .8, delta = 0.8, design = 1, conservative = !sharp,
             sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
             constant.var.time = FALSE, L = c(0, 0, 2, 0, 2, 2, 0),
             corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3], 
             niter = niter, notify = notify, pbDevice = c("pixel", "spectre"),
             postIdentifier = paste0("all assumptions ok\n",
                                     ifelse(sharp, "sharp n", "conservative n"))), envir = .GlobalEnv)
  save(file = "simsDesign1-delta8-noViol.RData",
       list = c(grep("d1large", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                "gammas", "lambdas", "seed", "times", "spltime"), precheck = TRUE)
}

if (notify)
  slackr_bot("All simulations are complete for effect size 0.8, no assumption violations")

stopCluster(clus)
