# generateSMART.R
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


#' Generate data from sequential, multiple-assignment, randomized trial with a longitudinal outcome
#'
#' @param n number of participants in the trial to generate
#' @param times vector of observation times at which the repeated-measures outcome is collected
#' @param spltime the observation time immediately before assessment of response/non-response and subsequent re-randomization
#' @param r1 probability of response to first-stage treatment A1 = 1
#' @param r0 probability of response to first-stage treatment A1 = -1
#' @param gammas vector of parameters indexing the marginal mean model
#' @param lambdas 
#' @param design type of SMART design to simulate, one of (1, 2, 3). See Details.
#' @param balanceRand logical. 
#' @param sigma 
#' @param sigma.r1 
#' @param sigma.r0 
#' @param corstr true correlation structure. one of "identity", "exchangeable", or "ar1".
#' @param uneqsdDTR 
#' @param uneqsd
#' @param varmats
#' @param respModel  
#' @param rho 
#' @param rho.r1 
#' @param rho.r0 
#' @param empirical 
#'
#' @return A data.frame
#' @export
#'
#' @examples


### TODO: assign class "generateSMART" to output list??
### - idea would be to include some design parameters as list elements to avoid having to take so many arguments to functions
 
generateSMART <- function(n, times, spltime, r1, r0, gammas, lambdas, design, balanceRand = FALSE,
                          sigma, sigma.r1, sigma.r0, corstr = c("identity", "exchangeable", "ar1"),
                          rho = NULL, rho.r1 = rho, rho.r0 = rho,
                          uneqsdDTR = NULL, uneqsd = NULL, varmats = NULL,
                          respDirection = NULL, respFunction,
                          empirical = FALSE, old = FALSE) {
  
  call <- match.call()
  
  ## Input Checks
  if (sum(times > spltime) == 0) stop("Currently spltime must be less than max(times)")
  
  ## Handle correlation structure
  corstr <- match.arg(corstr)
  if (corstr == "identity") {
    rho <- rho.r1 <- rho.r0 <- 0
  }
  
  ## Make sure design input is valid
  # if (!(design %in% 1:3)) stop("Invalid design choice. Must be one of 1-3.")
  if (design != 2) stop("New generative model structure only implemented for design 2")
  
  #### LEGACY GENERATIVE MODEL ####
  if (old) {
    
    d <- data.frame("id"   = 1:n,
                    "A1"   = 0)
    
    ## Generate first-stage treatment allocations
    if (balanceRand) {
      s <- sample(1:n, size = floor(n/2), replace = FALSE)
      d$A1[s]  <- 1
      d$A1[-s] <- -1
    } else {
      d$A1 <- 2 * rbinom(n, 1, 0.5) - 1
    }
    
    warning(paste0("Specifying old == TRUE forces response to be independent of previous outcomes.\n",
                   "respFunction and respDirection will be ignored."))
    
    ## Generate response status independent of previous outcomes
    
    d$R <- NA
    d$R[d$A1 ==  1] <- rbinom(sum(d$A1 ==  1), 1, r1)
    d$R[d$A1 == -1] <- rbinom(sum(d$A1 == -1), 1, r0)
    
    ## Generate second-stage treatment allocations.
    ## NOTE: this is design-specific
    
    d$A2NR <- d$A2R <- rep(0, n)
    
    if (design == 1) {
      if (balanceRand) {
        A1Rcombos <- expand.grid("A1" = c(1, -1), "R" = c(0, 1))
        for (i in 1:nrow(A1Rcombos)) {
          s <- sample(which(d$A1 == A1Rcombos$A1[i] & d$R == A1Rcombos$R[i]),
                      floor(sum(d$A1 == A1Rcombos$A1[i] & d$R == A1Rcombos$R[i]) / 2),
                      replace = FALSE)
          if (A1Rcombos$R[i] == 1) {
            d$A2R[s] <- 1
            d$A2R[setdiff(which(d$A1 == A1Rcombos$A1[i] & d$R == A1Rcombos$R[i]), s)] <- -1
          } else {
            d$A2NR[s] <- 1
            d$A2NR[setdiff(which(d$A1 == A1Rcombos$A1[i] & d$R == A1Rcombos$R[i]), s)] <- -1
          }
        }
      } else {
        d$A2R[d$A1 == 1  & d$R == 1] <- 2 * rbinom(sum(d$A1 == 1  & d$R == 1), 1, .5) - 1
        d$A2R[d$A1 == -1 & d$R == 1] <- 2 * rbinom(sum(d$A1 == -1 & d$R == 1), 1, .5) - 1
        d$A2NR[d$A1 == 1  & d$R == 0] <- 2 * rbinom(sum(d$A1 == 1  & d$R == 0), 1, .5) - 1
        d$A2NR[d$A1 == -1 & d$R == 0] <- 2 * rbinom(sum(d$A1 == -1 & d$R == 0), 1, .5) - 1
      }
    } else if (design == 2) {
      d$A2NR <- d$A2R <- rep(0, n)
      if (balanceRand) {
        for (i in c(-1, 1)) {
          s <- sample(which(d$A1 == i & d$R == 0), floor(sum(d$A1 == i & d$R == 0) / 2))
          d$A2NR[s] <- 1
          d$A2NR[setdiff(which(d$A1 == i & d$R == 0), s)] <- -1
        }
      } else {
        d$A2NR[d$A1 == 1  & d$R == 0]  <- 2 * rbinom(sum(d$A1 == 1  & d$R == 0), 1, .5) - 1
        d$A2NR[d$A1 == -1 & d$R == 0]  <- 2 * rbinom(sum(d$A1 == -1 & d$R == 0), 1, .5) - 1
      }
    } else if (design == 3) {
      d$A2NR[d$A1 == 1 & d$R == 0] <- 2 * rbinom(sum(d$A1 == 1 & d$R == 0), 1, .5) - 1
    }
    
    ## Generate Y's at their conditional means
    
    d <- do.call("rbind", lapply(split.SMART(d), function(z) {
      x <- do.call("data.frame", lapply(times, function(time) {
        conditional.model(z$A1, z$R, z$A2R, z$A2NR,
                          time, spltime, design,
                          get(paste0("r", (1 + unique(z$A1)) / 2)),
                          gammas, lambdas)
      }))
      names(x) <- sapply(times, function(t) paste0("Y", t))
      z <- cbind(z, x)
    }))
    d <- d[order(d$id), ]
    rownames(d) <- d$id
    
    ## Compute variance matrices if not already provided
    
    if (is.null(varmats)) {
      varmats <- conditionalVarmat(times, spltime, design, r1, r0,
                                   corstr, sigma, sigma.r1, sigma.r0,
                                   uneqsd = uneqsd, uneqsdDTR = uneqsdDTR,
                                   rho, rho.r1, rho.r0,
                                   gammas, lambdas)
    }
    
    ## Add noise with variance according to varmats
    
    d <- do.call("rbind",
                 lapply(split.SMART(d), function(x) {
                   varIndex <- as.character(paste(unique(x$A1), unique(x$R),
                                                  unique(x$A2R), unique(x$A2NR)))
                   x[, grepl("Y", names(x))] <- x[, grepl("Y", names(x))] +
                     mvrnorm(
                       n = nrow(x),
                       mu = rep(0, length(times)),
                       Sigma = varmats[[varIndex]],
                       empirical = empirical
                     )
                   x
                 })
    )
  } else {
    #### UPDATED GENERATIVE MODEL ####
    
    ## Initialize data frame
    d <- generateStage1(n = n, times = times, spltime = spltime, 
                        r1 = r1, r0 = r0,
                        gammas = gammas,
                        sigma = sigma, corstr = corstr,
                        rho = rho, 
                        respFunction = respFunction,
                        respDirection, 
                        balanceRand,
                        empirical)
    
    d <- generateStage2.means(d, times, spltime, gammas, lambdas, design, corstr)
    
    if (design == 1) {
      ## Time 2 variances
      
      # Check variance assumption:
      sigma.r1.LB <- sqrt(sigma^2 + ((1 - r1) / r1) * (mean(d$Y2.111[d$R.1 == 0]) - mean(d$Y2.111))^2 -
                            (mean(d$Y2.111[d$R.1 == 1]) - mean(d$Y2.111[d$R.1 == 0]))^2)
      sigma.r0.LB <- sqrt(sigma^2 + ((1 - r0) / r0) * (mean(d$Y2.000[d$R.0 == 0]) - mean(d$Y2.000))^2 -
                            (mean(d$Y2.000[d$R.0 == 1]) - mean(d$Y2.000[d$R.0 == 0]))^2)
      
      # Notation: 
      ## v2.(stage1).(R/NR).(stage2 given response) refers to *residual* variance
      ##  i.e., v2.1.R.1 is the additional variance needed for responders to A1=1 who are randomized to A2R=1
      
      v2.1.R.1   <- sigma.r1^2 - (rho / (1 + rho))^2 * with(subset(d, R.1 == 1), var(Y0 + Y1.1))
      sigma.nr11 <- sqrt((sigma^2 - r1 * sigma.r1^2 - r1*(1-r1) * with(d, mean(Y2.111[R.1 == 1]) - mean(Y2.111[R.1 == 0]))^2) / (1 - r1))
      v2.1.NR.1  <- sigma.nr11^2 - (rho / (1 + rho))^2 * with(subset(d, R.1 == 0), var(Y0 + Y1.1))
      sigma.r10  <- sqrt((sigma^2 - (1 - r1) * sigma.nr11^2 - r1*(1-r1) * with(d, mean(Y2.101[R.1 == 1]) - mean(Y2.101[R.1 == 0]))^2) / r1)
      v2.1.R.0   <- sigma.r10^2 - (rho / (1 + rho))^2 * with(subset(d, R.1 == 1), var(Y0 + Y1.1))
      sigma.nr10 <- sqrt((sigma^2 - r1 * sigma.r10^2 - r1*(1-r1) * with(d, mean(Y2.100[R.1 == 1]) - mean(Y2.100[R.1 == 0]))^2) / (1 - r1))
      v2.1.NR.0  <- sigma.nr10^2 - (rho / (1 + rho))^2 * with(subset(d, R.1 == 0), var(Y0 + Y1.1))
      
      v2.0.R.1   <- sigma.r0^2 - (rho / (1 + rho))^2 * with(subset(d, R.0 == 1), var(Y0 + Y1.0))
      sigma.nr01 <- sqrt((sigma^2 - r0 * sigma.r0^2 - r0*(1-r0) * with(d, mean(Y2.011[R.0 == 1]) - mean(Y2.011[R.0 == 0]))^2) / (1 - r0))
      v2.0.NR.1  <- sigma.nr01^2 - (rho / (1 + rho))^2 * with(subset(d, R.0 == 0), var(Y0 + Y1.0))
      sigma.r00  <- sqrt((sigma^2 - (1 - r0) * sigma.nr01^2 - r0*(1-r0) * with(d, mean(Y2.001[R.0 == 1]) - mean(Y2.001[R.0 == 0]))^2) / r0)
      v2.0.R.0   <- sigma.r00^2 - (rho / (1 + rho))^2 * with(subset(d, R.0 == 1), var(Y0 + Y1.0))
      sigma.nr00 <- sqrt((sigma^2 - r0 * sigma.r00^2 - r0*(1-r0) * with(d, mean(Y2.000[R.1 == 1]) - mean(Y2.000[R.1 == 0]))^2) / (1 - r0))
      v2.0.NR.0  <- sigma.nr00^2 - (rho / (1 + rho))^2 * with(subset(d, R.0 == 0), var(Y0 + Y1.0))
      
      ## Add residuals at time 2
      e2.1.R.1 <- rnorm(sum(d$R.1 == 1), 0, sqrt(v2.1.R.1))
      e2.1.R.0 <- rnorm(sum(d$R.1 == 1), 0, sqrt(v2.1.R.0))
      e2.0.R.1 <- rnorm(sum(d$R.0 == 1), 0, sqrt(v2.0.R.1))
      e2.0.R.0 <- rnorm(sum(d$R.0 == 1), 0, sqrt(v2.0.R.0))
      
      e2.1.NR.1 <- rnorm(sum(d$R.1 == 0), 0, sqrt(v2.1.NR.1))
      e2.1.NR.0 <- rnorm(sum(d$R.1 == 0), 0, sqrt(v2.1.NR.0))
      e2.0.NR.1 <- rnorm(sum(d$R.0 == 0), 0, sqrt(v2.0.NR.1))
      e2.0.NR.0 <- rnorm(sum(d$R.0 == 0), 0, sqrt(v2.0.NR.0))
      
      d$Y2.000[d$R.0 == 1] <- d$Y2.000[d$R.0 == 1] + e2.0.R.0
      d$Y2.000[d$R.0 == 0] <- d$Y2.000[d$R.0 == 0] + e2.0.NR.0
      
      d$Y2.001[d$R.0 == 1] <- d$Y2.001[d$R.0 == 1] + e2.0.R.0
      d$Y2.001[d$R.0 == 0] <- d$Y2.001[d$R.0 == 0] + e2.0.NR.1
      
      d$Y2.010[d$R.0 == 1] <- d$Y2.010[d$R.0 == 1] + e2.0.R.1
      d$Y2.010[d$R.0 == 0] <- d$Y2.010[d$R.0 == 0] + e2.0.NR.0
      
      d$Y2.011[d$R.0 == 1] <- d$Y2.011[d$R.0 == 1] + e2.0.R.1
      d$Y2.011[d$R.0 == 0] <- d$Y2.011[d$R.0 == 0] + e2.0.NR.1
      
      d$Y2.100[d$R.1 == 1] <- d$Y2.100[d$R.1 == 1] + e2.1.R.0
      d$Y2.100[d$R.1 == 0] <- d$Y2.100[d$R.1 == 0] + e2.1.NR.0
      
      d$Y2.101[d$R.1 == 1] <- d$Y2.101[d$R.1 == 1] + e2.1.R.0
      d$Y2.101[d$R.1 == 0] <- d$Y2.101[d$R.1 == 0] + e2.1.NR.1
      
      d$Y2.110[d$R.1 == 1] <- d$Y2.110[d$R.1 == 1] + e2.1.R.1
      d$Y2.110[d$R.1 == 0] <- d$Y2.110[d$R.1 == 0] + e2.1.NR.0
      
      d$Y2.111[d$R.1 == 1] <- d$Y2.111[d$R.1 == 1] + e2.1.R.1
      d$Y2.111[d$R.1 == 0] <- d$Y2.111[d$R.1 == 0] + e2.1.NR.1
      
      # Second-stage randomization
      d$A2NR <- d$A2R  <- 0
      d$A2R[d$R == 1]  <- 2 * rbinom(sum(d$R == 1), 1, .5) - 1
      d$A2NR[d$R == 0] <- 2 * rbinom(sum(d$R == 0), 1, .5) - 1
      
      d$weight <- 4
      
      d$dtr1 <- as.numeric(with(d, A1 ==  1 & (A2R ==  1 | A2NR ==  1)))
      d$dtr2 <- as.numeric(with(d, A1 ==  1 & (A2R ==  1 | A2NR == -1)))
      d$dtr3 <- as.numeric(with(d, A1 ==  1 & (A2R == -1 | A2NR ==  1)))
      d$dtr4 <- as.numeric(with(d, A1 ==  1 & (A2R == -1 | A2NR == -1)))
      d$dtr5 <- as.numeric(with(d, A1 == -1 & (A2R ==  1 | A2NR ==  1)))
      d$dtr6 <- as.numeric(with(d, A1 == -1 & (A2R ==  1 | A2NR == -1)))
      d$dtr7 <- as.numeric(with(d, A1 == -1 & (A2R == -1 | A2NR ==  1)))
      d$dtr8 <- as.numeric(with(d, A1 == -1 & (A2R == -1 | A2NR == -1)))
    } else if (design == 2) {
      ## DESIGN 2
      ## Time 2 variances
      
      # Check variance assumption:
      sigma.r1.LB <- sqrt(max(0, sigma^2 + ((1 - r1) / r1) * (mean(d$Y2.101[d$R.1 == 0]) - mean(d$Y2.101))^2 -
                            (mean(d$Y2.101[d$R.1 == 1]) - mean(d$Y2.101[d$R.1 == 0]))^2))
      sigma.r0.LB <- sqrt(max(0, sigma^2 + ((1 - r0) / r0) * (mean(d$Y2.000[d$R.0 == 0]) - mean(d$Y2.000))^2 -
                            (mean(d$Y2.000[d$R.0 == 1]) - mean(d$Y2.000[d$R.0 == 0]))^2))
      
      if (sigma.r1.LB > sigma.r1 | sigma.r0.LB > sigma.r0) {
        message("Conditional variation working assumption violated.")
        condVarAssump <- 1
      } else {
        condVarAssump <- 0
      }
      
      sigma.nr00 <- sqrt((sigma^2 - r0 * sigma.r0^2 - r0*(1-r0) * with(d, mean(Y2.000[R.0 == 1]) - mean(Y2.000[R.0 == 0]))^2) / (1 - r0))
      sigma.nr01 <- sqrt((sigma^2 - r0 * sigma.r0^2 - r0*(1-r0) * with(d, mean(Y2.001[R.0 == 1]) - mean(Y2.001[R.0 == 0]))^2) / (1 - r0))
      sigma.nr10 <- sqrt((sigma^2 - r1 * sigma.r1^2 - r1*(1-r1) * with(d, mean(Y2.100[R.1 == 1]) - mean(Y2.100[R.1 == 0]))^2) / (1 - r1))
      sigma.nr11 <- sqrt((sigma^2 - r1 * sigma.r1^2 - r1*(1-r1) * with(d, mean(Y2.101[R.1 == 1]) - mean(Y2.101[R.1 == 0]))^2) / (1 - r1))
      
      v2.0.R <- sigma.r0^2 - (rho / (1 + rho))^2 * with(subset(d, R.0 == 1), var(Y0 + Y1.0))
      v2.0.NR.0 <- sigma.nr00^2 - (rho / (1 + rho))^2 * with(subset(d, R.0 == 0), var(Y0 + Y1.0))
      v2.0.NR.1 <- sigma.nr01^2 - (rho / (1 + rho))^2 * with(subset(d, R.0 == 0), var(Y0 + Y1.0))
      
      v2.1.R <- sigma.r1^2 - (rho / (1 + rho))^2 * with(subset(d, R.1 == 1), var(Y0 + Y1.1))
      v2.1.NR.0 <- sigma.nr10^2 - (rho / (1 + rho))^2 * with(subset(d, R.1 == 0), var(Y0 + Y1.1))
      v2.1.NR.1 <- sigma.nr11^2 - (rho / (1 + rho))^2 * with(subset(d, R.1 == 0), var(Y0 + Y1.1))
      
      ## Add variance
      
      e2.0.R <- rnorm(sum(d$R.0), 0, sqrt(v2.0.R))
      e2.1.R <- rnorm(sum(d$R.1), 0, sqrt(v2.1.R))
      
      d$Y2.000[d$R.0 == 1] <- d$Y2.000[d$R.0 == 1] + e2.0.R
      d$Y2.000[d$R.0 == 0] <- d$Y2.000[d$R.0 == 0] + rnorm(sum(1 - d$R.0), 0, sqrt(v2.0.NR.0))
      
      d$Y2.001[d$R.0 == 1] <- d$Y2.001[d$R.0 == 1] + e2.0.R
      d$Y2.001[d$R.0 == 0] <- d$Y2.001[d$R.0 == 0] + rnorm(sum(1 - d$R.0), 0, sqrt(v2.0.NR.1))
      
      d$Y2.100[d$R.1 == 1] <- d$Y2.100[d$R.1 == 1] + e2.1.R
      d$Y2.100[d$R.1 == 0] <- d$Y2.100[d$R.1 == 0] + rnorm(sum(1 - d$R.1), 0, sqrt(v2.1.NR.0))
      
      d$Y2.101[d$R.1 == 1] <- d$Y2.101[d$R.1 == 1] + e2.1.R
      d$Y2.101[d$R.1 == 0] <- d$Y2.101[d$R.1 == 0] + rnorm(sum(1 - d$R.1), 0, sqrt(v2.1.NR.1))
      
     
      
      # Randomize stage 2
      d$A2R <- 0
      d$A2NR <- 0
      d$A2NR[d$R == 0] <- 2 * rbinom(sum(d$R == 0), 1, .5) - 1
      
      # Construct weights and DTR indicators
      d$weight <- 2*(d$R + 2 * (1 - d$R))
      d$dtr1 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR ==  1))))
      d$dtr2 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR == -1))))
      d$dtr3 <- as.numeric(with(d, (A1 == -1) * (R + (1 - R) * (A2NR ==  1))))
      d$dtr4 <- as.numeric(with(d, (A1 == -1) * (R + (1 - R) * (A2NR == -1))))
      
      # Select potential Y2 value to observe based on randomization
      d$Y2 <- NA
      d$Y2[d$dtr1 == 1] <- d$Y2.101[d$dtr1 == 1]
      d$Y2[d$dtr2 == 1 & d$R == 0] <- d$Y2.100[d$dtr2 == 1 & d$R == 0]
      d$Y2[d$dtr3 == 1] <- d$Y2.001[d$dtr3 == 1]
      d$Y2[d$dtr4 == 1 & d$R == 0] <- d$Y2.000[d$dtr4 == 1 & d$R == 0]
      
      # Compute weights
      d$weight <- 2 * (d$R + 2 * (1 - d$R))
      
      d$dtr1 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR ==  1))))
      d$dtr2 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR == -1))))
      d$dtr3 <- as.numeric(with(d, (A1 == -1) * (R + (1 - R) * (A2NR ==  1))))
      d$dtr4 <- as.numeric(with(d, (A1 == -1) * (R + (1 - R) * (A2NR == -1))))
      # }
    } else {
      
      d$weight <- 2 + 2 * (1 - d$R) * (d$A1 == 1)
      
      d$dtr1 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR ==  1))))
      d$dtr2 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR == -1))))
      d$dtr3 <- as.numeric(with(d, (A1 == -1)))
    }
  }
  
  d <- d[order(d$id), ]
  rownames(d) <- 1:n
  d.full <- d
  d <- subset(d, select = c("id", "Y0", "A1", "Y1", "R", "A2R",
                            "A2NR", "Y2", "weight",
                            grep("dtr", names(d), value = T)))
  
  class(d) <- c("SMART", class(d))
  
  # Check validity of treatment allocations
  if (!validTrial(d, design)) {
    return(list("data" = d,  "valid" = F, 
                "params" = list("times" = times, "spltime" = spltime, "r1" = r1, "r0" = r0,
                                "gammas" = gammas, "lambdas" = lambdas, "design" = design),
                "assumptions" = list("conditionalVariation" = NULL)))
  }
  
  output <- list("data" = d, "potentialData" = d.full, "valid" = T, 
                 "params" = list("times" = times, "spltime" = spltime, "r1" = r1, "r0" = r0,
                                 "gammas" = gammas, "lambdas" = lambdas, "design" = design),
                 "assumptions" = list("conditionalVariation" = condVarAssump))
  class(output) <- c('generateSMART', class(output))
  return(output)
}
