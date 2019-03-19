# simsDesign1-nrCorrDown50.R
# Copyright 2018 Nicholas J. Seewald
#
# This file is part of rmSMARTsize.
# 
# rmSMARTsize is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# rmSMARTsize is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with rmSMARTsize.  If not, see <https://www.gnu.org/licenses/>.

### ---------------------------------------------- ###
###             Simulations for Design 1           ###
### Violation of Conditional Covariance Assumption ###
###   50% Relative Reduction in NRs Correlation    ###
### ---------------------------------------------- ###

library(here)

notify <- TRUE
seed <- 1001
source(here("init.R"), echo = FALSE)

##### Simulation Setup #####

times <- c(0, 1, 2)
spltime <- 1

niter <- 3000
maxiter.solver <- 1000
tol <- 1e-8

sigma <- 6
sigma.r1 <- 6
sigma.r0 <- 6

# Generate a grid of simulation scenarios
simGrid <- expand.grid(
  list(
    sharp = c(FALSE, TRUE),
    r0 = c(.4, .6),
    r1 = c(.4, .6),
    corr = c(.3, .6, .8),
    respFunction = "response.indep",
    oldModel = TRUE
  ),
  stringsAsFactors = F
)

simGrid$corr.r1[simGrid$corr == .3 & simGrid$r1 == .4] <- .525
simGrid$corr.r0[simGrid$corr == .3 & simGrid$r0 == .4] <- .525
simGrid$corr.r1[simGrid$corr == .3 & simGrid$r1 == .6] <- .400
simGrid$corr.r0[simGrid$corr == .3 & simGrid$r0 == .6] <- .400

simGrid$corr.r1[simGrid$corr == .6 & simGrid$r1 == .4] <- .747
simGrid$corr.r0[simGrid$corr == .6 & simGrid$r0 == .4] <- .747
simGrid$corr.r1[simGrid$corr == .6 & simGrid$r1 == .6] <- .747
simGrid$corr.r0[simGrid$corr == .6 & simGrid$r0 == .6] <- .747

simGrid$corr.r1[simGrid$corr == .8 & simGrid$r1 == .4] <- .874
simGrid$corr.r0[simGrid$corr == .8 & simGrid$r0 == .4] <- .874
simGrid$corr.r1[simGrid$corr == .8 & simGrid$r1 == .6] <- .874
simGrid$corr.r0[simGrid$corr == .8 & simGrid$r0 == .6] <- .874


# Send initial notification that simulations are about to start
if(notify) {
  startString <- paste("All systems go so far!\nRunning on", nWorkers, "cores.")
  slackr_bot(startString)
  rm(startString)
}


##### Effect size: 0.3 #####

gammas <- c(35, -4, 2.2, -1.6, -1.3, 0.4, -0.4, 0.4, 0.4)
lambdas <- c(0.3, 0.4)

save(file = here("Results", "simsDesign1-delta3-nrCorrDown50.RData"),
     list = c("sigma", "simGrid", "gammas", "lambdas",
              "seed", "times", "spltime"))

for (scenario in 1:nrow(simGrid)) {
  # Extract simulation conditions from simGrid
  r0 <- simGrid$r0[scenario]
  r1 <- simGrid$r1[scenario]
  sharp <- simGrid$sharp[scenario]
  corr <- c(simGrid$corr[scenario], simGrid$corr.r1[scenario],
            simGrid$corr.r0[scenario])
  respFunc.name <- simGrid$respFunction[scenario]
  old <- simGrid$oldModel[scenario]
  
  postID <- paste0(
    "Scenario ", scenario, " of ", nrow(simGrid), "\n",
    "nrCorrDown50 simulation setup\n",
    "Effect size: 0.3\n",
    ifelse(sharp, "sharp n",
           "conservative n")
  )
  
  x <- simGrid[scenario, ]
  simName <- paste0("d1_delta3.",
                    ifelse(x$r0 == x$r1, paste0("r", x$r0* 10),
                           paste0("r0_", x$r0*10, ".r1_", x$r1*10)),
                    ".exch", x$corr * 10,
                    ifelse(x$sharp, ".sharp", ""),
                    ".nrCorrDown50")
  rm(x)

  # Set the seed for every unique simulation
  set.seed(seed)
  
  # Simulate
  if (notify) slackr_bot(simName)
  assign(simName,
         try(simulateSMART(
           gammas = gammas,
           lambdas = lambdas,
           r1 = r1,
           r0 = r0,
           times = times,
           spltime = spltime,
           alpha = .05,
           power = .8,
           delta = 0.3,
           design = 1,
           conservative = !sharp,
           sigma = sigma,
           sigma.r1 = sigma.r1,
           sigma.r0 = sigma.r0,
           L = c(0, 0, 2, 0, 2, 2, 2, 0, 0),
           corstr = "exch",
           rho = corr[1],
           rho.r1 = corr[2],
           rho.r0 = corr[3],
           respFunction = get(unlist(respFunc.name)),
           niter = niter,
           notify = notify,
           old = old,
           postIdentifier = postID
         )),
         envir = .GlobalEnv)
  
  # Save the result
  save(file = here("Results", "simsDesign1-delta3-nrCorrDown50.RData"),
       list = c(grep("d1_delta3", ls(), value = T), "sigma", "sigma.r1",
                "sigma.r0", "simGrid", "gammas", "lambdas", "seed", "times",
                "spltime"),
       precheck = TRUE)
}

if (notify) {
  x <- paste("All simulations are complete for Design 1,", 
             "effect size 0.3\n for nrCorrDown50 scenarios.")
  slackr_bot(x)
  rm(x)
}

rm(list = grep("d1_delta3", ls(), value = T))


##### Effect size: 0.5 #####

gammas <- c(35, -4, 2.2, -1.6, -0.7, 0.4, -0.4, 0.4, 0.4)
lambdas <- c(0.3, 0.4)

save(file = here("Results", "simsDesign1-delta5-nrCorrDown50.RData"),
     list = c("sigma", "simGrid", "gammas", "lambdas",
              "seed", "times", "spltime"))

for (scenario in 1:nrow(simGrid)) {
  # Extract simulation conditions from simGrid
  r0 <- simGrid$r0[scenario]
  r1 <- simGrid$r1[scenario]
  sharp <- simGrid$sharp[scenario]
  corr <- c(simGrid$corr[scenario], simGrid$corr.r1[scenario],
            simGrid$corr.r0[scenario])
  respFunc.name <- simGrid$respFunction[scenario]
  old <- simGrid$oldModel[scenario]
  
  postID <- paste0(
    "Scenario ", scenario, " of ", nrow(simGrid), "\n",
    "nrCorrDown50 simulation setup\n",
    "Effect size: 0.5\n",
    ifelse(sharp, "sharp n",
           "conservative n")
  )
  
  x <- simGrid[scenario, ]
  simName <- paste0("d1_delta5.",
                    ifelse(x$r0 == x$r1, paste0("r", x$r0* 10),
                           paste0("r0_", x$r0*10, ".r1_", x$r1*10)),
                    ".exch", x$corr * 10,
                    ifelse(x$sharp, ".sharp", ""),
                    ".nrCorrDown50")
  rm(x)
  
  # Set the seed for every unique simulation
  set.seed(seed)
  
  assign(simName,
         try(simulateSMART(
           gammas = gammas,
           lambdas = lambdas,
           r1 = r1,
           r0 = r0,
           times = times,
           spltime = spltime,
           alpha = .05,
           power = .8,
           delta = 0.5,
           design = 1,
           conservative = !sharp,
           sigma = sigma,
           sigma.r1 = sigma.r1,
           sigma.r0 = sigma.r0,
           L = c(0, 0, 2, 0, 2, 2, 2, 0, 0),
           corstr = "exch",
           rho = corr[1],
           rho.r1 = corr[2],
           rho.r0 = corr[3],
           respFunction = get(unlist(respFunc.name)),
           niter = niter,
           notify = notify,
           old = old,
           postIdentifier = postID
         )),
         envir = .GlobalEnv)
  
  save(file = here("Results", "simsDesign1-delta5-nrCorrDown50.RData"),
       list = c(grep("d1_delta5", ls(), value = T), "sigma", "sigma.r1",
                "sigma.r0", "simGrid", "gammas", "lambdas", "seed", "times",
                "spltime"),
       precheck = TRUE)
}

if (notify) {
  x <- paste("All simulations are complete for Design 1,", 
             "effect size 0.5\n for nrCorrDown50 scenarios.")
  slackr_bot(x)
  rm(x)
}

rm(list = grep("d1_delta5", ls(), value = T))

if(check.dompi) {
  closeCluster(clus)
  mpi.quit()
} else {
  stopCluster(clus)
}
