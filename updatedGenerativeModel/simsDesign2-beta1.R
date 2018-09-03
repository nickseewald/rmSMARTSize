###    Design II Simulations   ###
### Thresholded Response Model ###

library(here)

notify <- TRUE
seed <- 1001
source(here("init.R"), echo = FALSE)

##### Simulation Setup #####

times <- c(0, 1, 2)
spltime <- 1

niter <- 5000
maxiter.solver <- 1000
tol <- 1e-8

sigma <- 6
sigma.r1 <- 6
sigma.r0 <- 6

simGrid <- expand.grid(list(
  sharp = c(FALSE, TRUE),
  resp = c(.4, .6),
  corr = c(0, .3, .6, .8)
))
simGrid$corr.r1 <- simGrid$corr * sigma / sigma.r1
simGrid$corr.r0 <- simGrid$corr * sigma / sigma.r0

startString <- paste("All systems go so far!\nRunning on", nWorkers, "cores.")

if (notify)
  slackr_bot(startString)

#### Effect size: 0.3 #####
gammas <- c(33.5, -0.8, 0.9, -0.8, 0.4, -0.4, 0.1)
lambdas <- c(0, 0)

save(file = "Results/simsDesign2-delta3-noViol-beta.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

for (i in 1:nrow(simGrid)) {
  # Extract simulation conditions from simGrid
  resp <- r1 <- r0 <- simGrid$resp[i]
  sharp <- simGrid$sharp[i]
  corr <- c(simGrid$corr[i], simGrid$corr.r1[i], simGrid$corr.r0[i])
  
  # When corr is zero, sharp and conservative sample sizes are identical
  # Skip the one that we label "sharp" to avoid doing the same sim twice
  if (sum(corr == 0) == length(times) & sharp)
    next
  
  # Set the seed for every unique simulation
  set.seed(seed)
  
  # Simulate
  assign(paste0("d2small.r", resp * 10, ".exch", corr[1] * 10, ifelse(sharp, ".sharp", "")),
         try(
           simulateSMART(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                         alpha = .05, power = .8, delta = 0.3, design = 2, conservative = !sharp,
                         sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                         pool.time = TRUE, L = c(0, 0, 2, 0, 2, 2, 0),
                         corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3],
                         niter = niter, notify = notify, respModel = "beta",
                         postIdentifier = paste0("Design 2\nbeta response\n",
                                                 ifelse(sharp, "sharp n", "conservative n")))),
         envir = .GlobalEnv)
  
  # Save the result
  save(file = "Results/simsDesign2-delta3-noViol-twoT.RData",
       list = c(grep("d2small", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                "gammas", "lambdas", "seed", "times", "spltime", "simGrid"), precheck = TRUE)
}

if (notify)
  slackr_bot("All simulations are complete for effect size 0.3, no assumption violations, beta")

## Effect size: 0.5
gammas <- c(33.5, -0.8, 1.1, -0.8, 0.8, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

save(file = "Results/simsDesign2-delta5-noViol-beta.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

for (i in 1:nrow(simGrid)) {
  resp <- r1 <- r0 <- simGrid$resp[i]
  sharp <- simGrid$sharp[i]
  corr <- c(simGrid$corr[i], simGrid$corr.r1[i], simGrid$corr.r0[i])
  set.seed(seed)
  
  assign(paste0("d2med.r", resp * 10, ".exch", corr[1] * 10, ifelse(sharp, ".sharp", "")),
         try(
           simulateSMART(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                         alpha = .05, power = .8, delta = 0.5, design = 2, conservative = !sharp,
                         sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                         pool.time = TRUE, L = c(0, 0, 2, 0, 2, 2, 0),
                         corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3],
                         niter = niter, notify = notify, respModel = "beta",
                         postIdentifier = paste0("Design 2\nAll assumptions ok\n",
                                                 ifelse(sharp, "sharp n", "conservative n")))
         ), envir = .GlobalEnv)
  save(file = "Results/simsDesign2-delta5-noViol-beta.RData",
       list = c(grep("d2med", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                "gammas", "lambdas", "seed", "times", "spltime", "simGrid"), precheck = TRUE)
}

if (notify)
  slackr_bot("All simulations are complete for effect size 0.5, no assumption violations, beta")

## Effect size: 0.8
gammas <- c(33.5, -0.8, 1.5, -0.8, 1.3, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

save(file = "Results/simsDesign2-delta8-noViol-beta.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

for (i in 1:nrow(simGrid)) {
  # Extract simulation conditions from simGrid
  resp <- r1 <- r0 <- simGrid$resp[i]
  sharp <- simGrid$sharp[i]
  corr <- c(simGrid$corr[i], simGrid$corr.r1[i], simGrid$corr.r0[i])
  set.seed(seed)
  assign(paste0("d2large.r", resp * 10, ".exch", corr[1] * 10, ifelse(sharp, ".sharp", "")),
         try(
           simulateSMART(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                         alpha = .05, power = .8, delta = 0.8, design = 2, conservative = !sharp,
                         sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                         pool.time = TRUE, L = c(0, 0, 2, 0, 2, 2, 0),
                         corstr = "exch", rho = corr[1], rho.r1 = corr[2], rho.r0 = corr[3], 
                         niter = niter, notify = notify, respModel = "beta",
                         postIdentifier = paste0("Design 2\nAll assumptions ok\n",
                                                 ifelse(sharp, "sharp n", "conservative n")))
         ), envir = .GlobalEnv)
  save(file = "Results/simsDesign2-delta8-noViol-beta.RData",
       list = c(grep("d2large", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                "gammas", "lambdas", "seed", "times", "spltime", "simGrid"), precheck = TRUE)
}

if (notify)
  slackr_bot("All simulations are complete for effect size 0.8, no assumption violations, beta")

if (Sys.info()["sysname"] != "Windows") {
  closeCluster(clus)
  mpi.quit()
} else {
  stopCluster(clus)
}
