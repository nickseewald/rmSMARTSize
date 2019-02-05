# generativeFunctions.R
# Data Generation Helper Functions
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

generateStage1 <- function(n,
                           times,
                           spltime,
                           r1,
                           r0,
                           gammas,
                           sigma,
                           corstr = c("identity", "exchangeable", "ar1"),
                           rho = 0,
                           respFunction,
                           respDirection = NULL,
                           balanceRand = FALSE,
                           empirical = FALSE,
                           old = FALSE) {
  
  if (old) 
    warning("Option 'old' currently does nothing in generateStage1()")
  if (is.null(n))
    stop("There's a problem with n being generated")
  
  corstr <- match.arg(corstr)
  
  ## Initialize data frame
  d <- data.frame("id" = 1:n)
  
  ## Generate baseline outcome
  d$Y0 <- gammas[1] + rnorm(n, 0, sigma)
  
  ## Generate observed treatment assignment
  if (balanceRand) {
    s <- sample(1:n, size = floor(n/2), replace = FALSE)
    d$A1[s]  <- 1
    d$A1[-s] <- -1
  } else {
    d$A1 <- 2 * rbinom(n, 1, 0.5) - 1
  }
  
  ## Generate stage 1 potential outcomes
  if (corstr %in% c("exchangeable", "identity")) {
    # FIXME: This only works for time 1 right now! Any other time will not have 
    # proper correlation structure
    d[, paste0("Y", times[times > 0 & times <= spltime], ".0")] <-
      (1-rho)*gammas[1] + rho*d$Y0 + gammas[2] - gammas[3] + 
      rnorm(n, 0, sqrt((1-rho^2))*sigma)
    d[, paste0("Y", times[times > 0 & times <= spltime], ".1")] <-
      (1-rho)*gammas[1] + rho*d$Y0 + gammas[2] + gammas[3] + 
      rnorm(n, 0, sqrt((1-rho^2))*sigma)
  } else if (corstr == "ar1") {
    # FIXME: This only works for time 1 right now! Any other time will not have 
    # proper correlation structure
    d[[paste0("Y", time, ".0")]] <- (1-rho)*gammas[1] + rho*d$Y0 + gammas[2] -
      gammas[3] + rnorm(n, 0, sqrt((1-rho^2))*sigma)
    d[[paste0("Y", time, ".1")]] <- (1-rho)*gammas[1] + rho*d$Y0 + gammas[2] +
      gammas[3] + rnorm(n, 0, sqrt((1-rho^2))*sigma)
  } else {
    warning("Y1 not created")
  }
  
  ## Generate response status
  resp <- respFunction(d, gammas, r1, r0, respDirection, sigma, causal = T)
  d <- resp$data
  r1 <- resp$r1
  r0 <- resp$r0
  
  ## Select potential Y1 value to observe based on randomization
  d$Y1 <- NA
  d$Y1[d$A1 == 1]  <- d$Y1.1[d$A1 == 1]
  d$Y1[d$A1 == -1] <- d$Y1.0[d$A1 == -1]
  
  ## Select potential R value to observe based on randomization
  d$R <- NA
  d$R[d$A1 == 1]  <- d$R.1[d$A1 == 1]
  d$R[d$A1 == -1] <- d$R.0[d$A1 == -1]
  
  output <- list("data" = d, "r0" = r0, "r1" = r1, "rho" = rho)
  class(output) <- c("generateStage1", class(output))
  output
}

generateStage2.means <-
  function(stage1,
           times,
           spltime,
           gammas,
           lambdas,
           design,
           corstr = c("identity", "exchangeable", "ar1")) {
    
  if (!("generateStage1") %in% class(stage1))
    stop("'stage1' must be generated from 'generateStage1'")
  
  corstr <- match.arg(corstr)
  
  # Extract stage-1 data and parameters
  d   <- stage1$data
  r1  <- stage1$r1
  r0  <- stage1$r0
  rho <- stage1$rho
  
  n <- nrow(d)
  nDTR <- switch(design, 8, 4, 3)
  DTRs <- dtrIndex(design)
  dtrTriples <- dtrNames(design)
  
  dtrTxts <- do.call(rbind, DTRs)
  
  # Replicate the "potential" response status variables for each DTR,
  # so that each DTR gets the appropriate potential response
  respMatrix <- do.call(cbind, lapply((1+dtrTxts[1, ])/2, 
                                      function(a1) d[[paste0("R.", a1)]]))
  
  # Create a vector of response probabilities, each entry of which corresponds 
  # to the appropriate value (r0 or r1) for each DTR
  respProbVec <- sapply(paste0("r", (1+dtrTxts[1, ])/2), get,
                        envir = sys.frame(sys.parent(0)))
  
  # Create matrices of the "adjustments" added to the mean model that are 
  # specific to responders and non-responders
  if (design == 1) {
    designMeanAdj.R <- matrix(rep(
      dtrTxts[2,] * (gammas[6] + gammas[8] * dtrTxts[1,]) / respProbVec +
        (1 - respProbVec) * (lambdas[1] + lambdas[2] * dtrTxts[1,]),
      n
    ),
    nrow = n,
    byrow = T)
    designMeanAdj.NR <- matrix(rep(
      dtrTxts[3,] * (gammas[7] + gammas[9] * dtrTxts[1,]) / (1 - respProbVec) -
        respProbVec * (lambdas[1] + lambdas[2] * dtrTxts[1,]),
      n
    ),
    nrow = n,
    byrow = T)
  } else if (design == 2) {
    designMeanAdj.R <- 
      matrix(rep((1 - respProbVec) * (lambdas[1] + lambdas[2] * dtrTxts[1,]),
                 n), nrow = n, byrow = T)
    designMeanAdj.NR <- matrix(rep(
      dtrTxts[3,] * (gammas[6] + gammas[7] * dtrTxts[1,]) / (1 - respProbVec) -
        respProbVec * (lambdas[1] + lambdas[2] * dtrTxts[1,]),
      n
    ),
    nrow = n,
    byrow = T)
  } else if (design == 3) {
    
  } else
    stop("'design' must be one of 1, 2, or 3.")
  
  if (corstr %in% c("identity", "exchangeable")) {
    # FIXME: This only works for time 2 right now! Any other time will not have 
    # proper correlation structure
    d[, paste0("Y", times[times > spltime], ".", dtrTriples)] <- 
      # Start with a matrix of just the appropriate linear combinations of 
      # gammas (only up to the stage 1 parts of the time 2 model)
      matrix(rep(
        ((1 - rho) / (1 + rho)) * gammas[1] + 
          (gammas[2] + gammas[3] * dtrTxts[1,]) / (1 + rho) +
          gammas[4] + gammas[5] * dtrTxts[1,],
        n
      ),
      nrow = n,
      byrow = T) + 
      # Add Y0 and the appropriate causal Y1
      (rho/(1+rho)) * (d$Y0 + 
                         matrix(c(rep(d$Y1.1, sum(dtrTxts[1, ] == 1)),
                                      rep(d$Y1.0, sum(dtrTxts[1, ] == -1))),
                                nrow = n)
                       ) +
      # Add response adjustment
      respMatrix * designMeanAdj.R + (1 - respMatrix) * designMeanAdj.NR
  } else if (corstr == "ar1") {
    stop("generateStage2 not yet implemented for corstr = 'ar1'")
  }
  d
}

testVarianceInput <- function(x, a, d, r0, r1, rho) {
  #TODO: is this design-agnostic?
  #TODO: What happens when this is always positive?
  y2string <- paste0("Y2.", paste0(a, collapse = ""))
  rstring <- paste0("R.", a[1])
  rprobstring <- paste0("r", a[1])
  sigma.nr <- 
    # sqrt(max(0.1,
             (sigma^2 - get(rprobstring) * x^2 - 
                get(paste0("r", a[1]))*(1-get(paste0("r", a[1]))) * 
                with(d, mean(get(y2string)[get(rstring) == 1]) -
                       mean(get(y2string)[get(rstring) == 0]))^2) /
               (1 - get(rprobstring))
    # ))
  sigma.nr - (rho / (1 + rho))^2 * with(subset(d, get(rstring) == 0), 
                                          var(Y0 + get(paste0("Y1.", a[1]))))
  
  # ((sigma^2 - r0 * sigma.r0^2 - r0*(1-r0) * 
  #     with(d, mean(Y2.000[R.0 == 1]) - mean(Y2.000[R.0 == 0]))^2) /
  #     (1 - r0))
}

testGenerativeVariances <- function(simGrid, times, spltime, gammas, sigma,
                                    corstr, design, balanceRand = F, 
                                    empirical = F, seed = 6781) {
  
  if (!is.null(seed))
    set.seed(seed)
  
  sg <- foreach(i = 1:nrow(simGrid), .combine = rbind, .inorder = TRUE, 
          .export = ls(envir = .GlobalEnv)) %dorng% {

    # Extract parameters from simGrid
    rho <- simGrid$corr[i]
    respFunction <- get(unlist(simGrid$respFunction[i]))
    
    r0 <- simGrid$r0[i]
    r1 <- simGrid$r1[i]
    
    # Generate a big data frame to estimate some of these things
    s1 <- generateStage1(1e5, times, spltime, r1, r0, gammas, sigma, corstr,
                         rho = rho, respFunction = respFunction,
                         respDirection = "high", balanceRand = F,
                         empirical = F)
    
    s <- generateStage2.means(s1, times, spltime, gammas, lambdas, 
                              design = design, corstr = corstr)
    
    # Extract (potentially-modified) response probabilities from s1
    r0  <- s1$r0
    r1  <- s1$r1
    rho <- s1$rho
    
    dtrTriples <- dtrNames(design)
    
    sigma.r0.LB <- 
      sqrt(with(s,
                max(
                  sapply(dtrTriples[substr(dtrTriples, 1, 1) == 0],
                         function(dtr) {
                           max(0.1, sigma^2 + (1-r0)/r0 *
                                 (mean(get(paste0("Y2.", dtr))[R.0 == 0]) -
                                    mean(get(paste0("Y2.", dtr))))^2 -
                                 (mean(get(paste0("Y2.", dtr))[R.0 == 1]) -
                                    mean(get(paste0("Y2.", dtr))[R.0 == 0]))^2)
                         })
                )
      ))
    
    # simGrid$sigma.r0.LB[i] <- 
    #   with(s, max(sqrt(sigma^2 + (1-r0)/(r0) * 
    #                      (mean(Y2.000[R.0 == 0]) - mean(Y2.000))^2 - 
    #                      (mean(Y2.000[R.0 == 1]) - mean(Y2.000[R.0 == 0]))^2),
    #               sqrt(sigma^2 + (1-r0)/(r0) * 
    #                      (mean(Y2.001[R.0 == 0]) - mean(Y2.001))^2 -
    #                      (mean(Y2.001[R.0 == 1]) - mean(Y2.001[R.0 == 0]))^2))
    #   )
    
    sigma.r0.UB<- 
      min(uniroot(testVarianceInput, interval = c(-1, 2*sigma),
                  a = c(0, 0, 1), d = s, r0 = r0, r1 = r1, rho = rho,
                  extendInt = "yes")$root,
          uniroot(testVarianceInput, interval = c(-1, 2*sigma),
                  a = c(0, 0, 0), d = s, r0 = r0, r1 = r1, rho = rho,
                  extendInt = "yes")$root)
    
    sigma.r0 <- mean(c(sigma.r0.LB, sigma.r0.UB))
    
    sigma.r1.LB <- 
      sqrt(with(s,
                max(
                  sapply(dtrTriples[substr(dtrTriples, 1, 1) == 1],
                         function(dtr) {
                           max(0.1, sigma^2 + (1-r1)/r1 *
                                 (mean(get(paste0("Y2.", dtr))[R.1 == 0]) -
                                    mean(get(paste0("Y2.", dtr))))^2 -
                                 (mean(get(paste0("Y2.", dtr))[R.1 == 1]) -
                                    mean(get(paste0("Y2.", dtr))[R.1 == 0]))^2)
                         })
                )
      ))
    
    # simGrid$sigma.r1.LB[i] <- 
    #   with(s, max(sqrt(sigma^2 + (1-r1)/(r1) * 
    #                      (mean(Y2.100[R.1 == 0]) - mean(Y2.100))^2 - 
    #                      (mean(Y2.100[R.1 == 1]) - mean(Y2.100[R.1 == 0]))^2),
    #               sqrt(sigma^2 + (1-r1)/(r1) * 
    #                      (mean(Y2.101[R.1 == 0]) - mean(Y2.101))^2 -
    #                      (mean(Y2.101[R.1 == 1]) - mean(Y2.101[R.1 == 0]))^2))
    #   )
    # 
    sigma.r1.UB <- 
      min(uniroot(testVarianceInput, interval = c(-1, 2*sigma), 
                  a = c(1, 0, 1), d = s, r0 = r0, r1 = r1, rho = rho,
                  extendInt = "yes")$root,
          uniroot(testVarianceInput, interval = c(-1, 2*sigma),
                  a = c(1, 0, 0), d = s, r0 = r0, r1 = r1, rho = rho,
                  extendInt = "yes")$root)
    
    sigma.r1 <- mean(c(sigma.r1.LB, sigma.r1.UB))
    
    cbind(simGrid[i, ], "sigma.r0.LB" = sigma.r0.LB, 
          "sigma.r0.UB" = sigma.r0.UB, "sigma.r0" = sigma.r0,
          "sigma.r1.LB" = sigma.r1.LB, "sigma.r1.UB" = sigma.r1.UB,
          "sigma.r1" = sigma.r1)
          }
  sg
}
