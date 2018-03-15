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
clus <- makeCluster(15)
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

notify <- T
startString <- paste("All systems go so far!\nRunning on", length(clus), "cores.")

slackr_bot(startString)

### Effect size: 0.3
gammas <- c(33.5, -0.8, 0.9, -0.8, 0.4, -0.4, 0.1)
lambdas <- c(0.1, -0.5)

sigma <- 6
sigma.r1 <- 5.6
sigma.r0 <- 5.6

save(file = "simsDesign2-delta3-violS1b-test.RData",
     list = c("sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"))

## 10% violation of S1(b)
simGrid <- expand.grid(list(sharp = c(FALSE, TRUE),
                            resp = c(.4),
                            corr = list(c(0.3, .23, .23))
))

for (i in 1:dim(simGrid)[1]) {
  r1 <- r0 <- simGrid$resp[i]
  set.seed(seed)
  with(simGrid, 
       assign(paste0("d2small.r", resp[i] * 10, ".exch",
                     corr[[i]][1] * 10, ".violS1b.1", ifelse(sharp[i], ".sharp", "")),
              sim(gammas = gammas, lambdas = lambdas, r1 = r1, r0 = r0, times = times, spltime = spltime,
                  alpha = .05, power = .8, delta = 0.3, design = 2, conservative = !sharp[i],
                  sigma = sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                  constant.var.time = FALSE, L.eos = c(0, 0, 2, 0, 2, 2, 0),
                  corstr = "exch", rho = corr[[i]][1], rho.r1 = corr[[i]][2], rho.r0 = corr[[i]][3],
                  niter = niter, notify = notify, pbDevice = c("pixel", "spectre"),
                  postIdentifier = paste0("Violate S1(b)\n", ifelse(sharp[i], "sharp n", "conservative n"))),
              envir = .GlobalEnv))
  save(file = "simsDesign2-delta3-violS1b-test.RData",
       list = c(grep("^d2small.*violS1b.1", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
                "gammas", "lambdas", "seed", "times", "spltime"), precheck = TRUE, envir = environment())
}

save(file = "simsDesign2-delta3-violS1b-test.RData",
     list = c(grep("^d2small.*violS1b.1", ls(), value = T), "sigma", "sigma.r1", "sigma.r0",
              "gammas", "lambdas", "seed", "times", "spltime"), precheck = TRUE)

slackr_bot("All simulations are complete for effect size 0.3.")

stopCluster(clus)
