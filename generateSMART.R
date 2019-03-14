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


#' Generate data from sequential, multiple-assignment, randomized trial with a
#'  longitudinal outcome
#'
#' @param n number of participants in the trial to generate
#' @param times vector of observation times at which the repeated-measures 
#' outcome is collected
#' @param spltime the observation time immediately before assessment of 
#' response/non-response and subsequent re-randomization
#' @param r1 probability of response to first-stage treatment A1 = 1
#' @param r0 probability of response to first-stage treatment A1 = -1
#' @param gammas vector of parameters indexing the marginal mean model
#' @param lambdas 
#' @param design type of SMART design to simulate, one of (1, 2, 3). 
#' See Details.
#' @param balanceRand logical. 
#' @param sigma 
#' @param sigma.r1 
#' @param sigma.r0 
#' @param corstr true correlation structure. one of "identity", "exchangeable",
#'  or "ar1".
#' @param uneqsdDTR 
#' @param uneqsd
#' @param varmats
#' @param variances a named vector of variance components to be used in the
#'  updated generative model. Because the new generative model needs to have 
#'  variance parameters defined numerically, these MUST be provided externally
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

generateSMART <- function(n, times, spltime, r1, r0, gammas, lambdas, design,
                          balanceRand = FALSE,
                          sigma, sigma.r1, sigma.r0, 
                          corstr = c("identity", "exchangeable", "ar1"),
                          rho = NULL, rho.r1 = rho, rho.r0 = rho,
                          uneqsdDTR = NULL, uneqsd = NULL,
                          varmats = NULL, variances = NULL,
                          respDirection = c("high", "low"), respFunction,
                          empirical = FALSE, old = FALSE) {
  
  call <- match.call()
  respDirection <- match.arg(respDirection)
  
  ## Input Checks
  if (sum(times > spltime) == 0) 
    stop("Currently spltime must be less than max(times)")
  if (!old & is.null(variances))
    stop(paste("If using the updated generative model, you currently need to",
               "provide the variances argument."))
  
  ## Handle correlation structure
  corstr <- match.arg(corstr)
  if (corstr == "identity") {
    rho <- rho.r1 <- rho.r0 <- 0
  }
  
  ## Make sure design input is valid
  if (!(design %in% 1:3)) stop("Invalid design choice. Must be one of 1-3.")
  # if (design != 2) 
  #   stop("New generative model structure only implemented for design 2")
  
  dtrTriples <- dtrNames(design)

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
    
    warning(paste("Specifying old == TRUE forces response to be independent of", 
                  "previous outcomes.\n", "respFunction and respDirection will", 
                  "be ignored."))
    
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
                      floor(sum(d$A1 == A1Rcombos$A1[i] & 
                                  d$R == A1Rcombos$R[i]) / 2),
                      replace = FALSE)
          if (A1Rcombos$R[i] == 1) {
            d$A2R[s] <- 1
            d$A2R[setdiff(which(d$A1 == A1Rcombos$A1[i] & 
                                  d$R == A1Rcombos$R[i]), s)] <- -1
          } else {
            d$A2NR[s] <- 1
            d$A2NR[setdiff(which(d$A1 == A1Rcombos$A1[i] & 
                                   d$R == A1Rcombos$R[i]), s)] <- -1
          }
        }
      } else {
        d$A2R[d$A1 == 1  & d$R == 1] <-
          2 * rbinom(sum(d$A1 == 1  & d$R == 1), 1, .5) - 1
        d$A2R[d$A1 == -1 & d$R == 1] <-
          2 * rbinom(sum(d$A1 == -1 & d$R == 1), 1, .5) - 1
        d$A2NR[d$A1 == 1  & d$R == 0] <-
          2 * rbinom(sum(d$A1 == 1  & d$R == 0), 1, .5) - 1
        d$A2NR[d$A1 == -1 & d$R == 0] <-
          2 * rbinom(sum(d$A1 == -1 & d$R == 0), 1, .5) - 1
      }
    } else if (design == 2) {
      d$A2NR <- d$A2R <- rep(0, n)
      if (balanceRand) {
        for (i in c(-1, 1)) {
          s <- sample(which(d$A1 == i & d$R == 0), 
                      floor(sum(d$A1 == i & d$R == 0) / 2))
          d$A2NR[s] <- 1
          d$A2NR[setdiff(which(d$A1 == i & d$R == 0), s)] <- -1
        }
      } else {
        d$A2NR[d$A1 == 1  & d$R == 0] <- 
          2 * rbinom(sum(d$A1 == 1  & d$R == 0), 1, .5) - 1
        d$A2NR[d$A1 == -1 & d$R == 0] <-
          2 * rbinom(sum(d$A1 == -1 & d$R == 0), 1, .5) - 1
      }
    } else if (design == 3) {
      d$A2NR[d$A1 == 1 & d$R == 0] <-
        2 * rbinom(sum(d$A1 == 1 & d$R == 0), 1, .5) - 1
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
                   varIndex <- as.character(paste(unique(x$A1),
                                                  unique(x$R),
                                                  unique(x$A2R),
                                                  unique(x$A2NR)))
                   x[, grepl("Y", names(x))] <-
                     x[, grepl("Y", names(x))] +
                     mvrnorm(
                       n = nrow(x),
                       mu = rep(0, length(times)),
                       Sigma = varmats[[varIndex]],
                       empirical = empirical
                     )
                   x
                 })
    )
    
    # Add DTR indicators
    d <- createDTRIndicators(d, design)
    
  } else {
    #### UPDATED GENERATIVE MODEL ####
    ## Initialize data frame
    d <- generateStage1(n = n, times = times, spltime = spltime, 
                        r1 = r1, r0 = r0,
                        gammas = gammas,
                        sigma = sigma, corstr = corstr,
                        rho = rho, 
                        design = design,
                        respFunction = respFunction,
                        respDirection, 
                        balanceRand,
                        empirical)
    
    # Extract parameters from stage 1
    r0 <- d$r0
    r1 <- d$r1
    rho <- d$rho
    
    d <- generateStage2.means(d, times, spltime, gammas, lambdas,
                              design, corstr)
    
    if (design == 1) {
      ## Time 2 variances
      
      # Notation: 
      ## v2.(stage1).(R/NR).(stage2 given response) refers to *residual* 
      ## variance; i.e., v2.1.R.1 is the additional variance needed for 
      ## responders to A1=1 who are randomized to A2R=1
      
      
      ## Add residuals at time 2
      e2.1.R.1 <- rnorm(sum(d$R.1 == 1), 0, sqrt(variances$v2.1.R.1))
      e2.1.R.0 <- rnorm(sum(d$R.1 == 1), 0, sqrt(variances$v2.1.R.0))
      e2.0.R.1 <- rnorm(sum(d$R.0 == 1), 0, sqrt(variances$v2.0.R.1))
      e2.0.R.0 <- rnorm(sum(d$R.0 == 1), 0, sqrt(variances$v2.0.R.0))
      
      e2.1.NR.1 <- rnorm(sum(d$R.1 == 0), 0, sqrt(variances$v2.1.NR.1))
      e2.1.NR.0 <- rnorm(sum(d$R.1 == 0), 0, sqrt(variances$v2.1.NR.0))
      e2.0.NR.1 <- rnorm(sum(d$R.0 == 0), 0, sqrt(variances$v2.0.NR.1))
      e2.0.NR.0 <- rnorm(sum(d$R.0 == 0), 0, sqrt(variances$v2.0.NR.0))
      
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
      
      # Add DTR indicators
      d <- createDTRIndicators(d, design)
      
      # Select potential Y2 value to observe based on randomization
      d$Y2 <- NA
      
      # Responders, dtr1 & dtr2
      d$Y2[d$dtr1 == 1 & d$R == 1] <- d$Y2.111[d$dtr1 == 1 & d$R == 1]
      # Responders, dtr3 & dtr4
      d$Y2[d$dtr4 == 1 & d$R == 1] <- d$Y2.100[d$dtr4 == 1 & d$R == 1]
      # Non-responders, dtr1 & dtr3
      d$Y2[d$dtr3 == 1 & d$R == 0] <- d$Y2.101[d$dtr3 == 1 & d$R == 0]
      # Non-responders, dtr2 & dtr4
      d$Y2[d$dtr2 == 1 & d$R == 0] <- d$Y2.110[d$dtr2 == 1 & d$R == 0]
      
      # Responders, dtr5 & dtr6
      d$Y2[d$dtr5 == 1 & d$R == 1] <- d$Y2.011[d$dtr5 == 1 & d$R == 1]
      # Responders, dtr7 & dtr8
      d$Y2[d$dtr8 == 1 & d$R == 1] <- d$Y2.000[d$dtr8 == 1 & d$R == 1]
      # Non-responders, dtr5 & dtr7
      d$Y2[d$dtr7 == 1 & d$R == 0] <- d$Y2.001[d$dtr7 == 1 & d$R == 0]
      # Non-responders, dtr6 & dtr8
      d$Y2[d$dtr6 == 1 & d$R == 0] <- d$Y2.010[d$dtr6 == 1 & d$R == 0]
      
    } else if (design == 2) {
      ## DESIGN 2
      ## Time 2 variances
      
      e2.0.R <- rnorm(sum(d$R.0), 0, sqrt(variances$v2.0.R.0))
      e2.1.R <- rnorm(sum(d$R.1), 0, sqrt(variances$v2.1.R.1))
      
      d$Y2.000[d$R.0 == 1] <- d$Y2.000[d$R.0 == 1] + e2.0.R
      d$Y2.000[d$R.0 == 0] <- d$Y2.000[d$R.0 == 0] + 
        rnorm(sum(1 - d$R.0), 0, sqrt(variances$v2.0.NR.0))
      
      d$Y2.001[d$R.0 == 1] <- d$Y2.001[d$R.0 == 1] + e2.0.R
      d$Y2.001[d$R.0 == 0] <- d$Y2.001[d$R.0 == 0] + 
        rnorm(sum(1 - d$R.0), 0, sqrt(variances$v2.0.NR.1))
      
      d$Y2.100[d$R.1 == 1] <- d$Y2.100[d$R.1 == 1] + e2.1.R
      d$Y2.100[d$R.1 == 0] <- d$Y2.100[d$R.1 == 0] + 
        rnorm(sum(1 - d$R.1), 0, sqrt(variances$v2.1.NR.0))
      
      d$Y2.101[d$R.1 == 1] <- d$Y2.101[d$R.1 == 1] + e2.1.R
      d$Y2.101[d$R.1 == 0] <- d$Y2.101[d$R.1 == 0] + 
        rnorm(sum(1 - d$R.1), 0, sqrt(variances$v2.1.NR.1))
      
      DTRs <- dtrIndex(design)
      dtrTriples <- as.data.frame((do.call(cbind, DTRs) + 1) / 2)
      dtrTriples[dtrTriples == 0.5] <- 0
      dtrTriples <- apply(dtrTriples, 1, function(x) paste(x, collapse = ""))
      
      # Randomize stage 2
      d$A2R <- 0
      d$A2NR <- 0
      d$A2NR[d$R == 0] <- 2 * rbinom(sum(d$R == 0), 1, .5) - 1
      
      # Add DTR indicators
      d <- createDTRIndicators(d, design)
      
      # Select potential Y2 value to observe based on randomization
      d$Y2 <- NA
      d$Y2[d$dtr1 == 1] <- d$Y2.101[d$dtr1 == 1]
      d$Y2[d$dtr2 == 1 & d$R == 0] <- d$Y2.100[d$dtr2 == 1 & d$R == 0]
      d$Y2[d$dtr3 == 1] <- d$Y2.001[d$dtr3 == 1]
      d$Y2[d$dtr4 == 1 & d$R == 0] <- d$Y2.000[d$dtr4 == 1 & d$R == 0]
    } else {
      ## DESIGN 3
      ## Time 2 variances
      
      # Generate noise for responders to a1 = 1
      e2.1.R <- rnorm(sum(d$R.1 == 1), 0, sqrt(variances$v2.1.R.1))
      
      d$Y2.101[d$R.1 == 1] <- d$Y2.101[d$R.1 == 1] + e2.1.R
      d$Y2.101[d$R.1 == 0] <- d$Y2.101[d$R.1 == 0] + 
        rnorm(sum(1 - d$R.1), 0, sqrt(variances$v2.1.NR.1))
      
      d$Y2.100[d$R.1 == 1] <- d$Y2.100[d$R.1 == 1] + e2.1.R
      d$Y2.100[d$R.1 == 0] <- d$Y2.100[d$R.1 == 0] + 
        rnorm(sum(1 - d$R.1), 0, sqrt(variances$v2.1.NR.0))
      
      d$Y2.000[d$R.0 == 1] <- d$Y2.000[d$R.0 == 1] +
        rnorm(sum(d$R.0), 0, sqrt(variances$v2.0.R.0))
      d$Y2.000[d$R.0 == 0] <- d$Y2.000[d$R.0 == 0] + 
        rnorm(sum(1 - d$R.0), 0, sqrt(variances$v2.0.NR.0))
      
      # Check conditional variation assumption
      DTRs <- dtrIndex(design)
      
      dtrTriples <- as.data.frame((do.call(cbind, DTRs) + 1) / 2)
      dtrTriples[dtrTriples == 0.5] <- 0
      dtrTriples <- apply(dtrTriples, 1, function(x) paste(x, collapse = ""))
      
      # Randomize stage 2
      d$A2R <- 0
      d$A2NR <- 0
      d$A2NR[d$R == 0 & d$A1 == 1] <-
        2 * rbinom(sum(d$R == 0 & d$A1 == 1), 1, .5) - 1
      
      # Add DTR indicators
      d <- createDTRIndicators(d, design)
      
      # Select potential Y2 value to observe based on randomization
      d$Y2 <- NA
      d$Y2[d$dtr1 == 1] <- d$Y2.101[d$dtr1 == 1]
      d$Y2[d$dtr2 == 1 & d$R == 0] <- d$Y2.100[d$dtr2 == 1 & d$R == 0]
      d$Y2[d$dtr3 == 1] <- d$Y2.000[d$dtr3 == 1]
    }
  }
  # 
  # varCheck <- 
  #   do.call(c, 
  #           lapply(dtrTriples, function(dtr) {
  #             index <- with(d, get(paste0("R.", substr(dtr, 1, 1))) == 0)
  #             x <- with(d, mean((get(paste0("Y2.", dtr))[index] -
  #                                  mean(get(paste0("Y2.", dtr))))^2))
  #           }))
  # 
  # if (any(varCheck > sigma^2)) {
  #   condVarAssump <- 1
  # } else {
  #   condVarAssump <- 0
  # }
  
  # Assign weights
  if (design == 1) {
    d$weight <- 4
  } else if (design == 2) {
    d$weight <- 2*(d$R + 2 * (1 - d$R))
  } else {
    d$weight <- 2 + 2 * (1 - d$R) * (d$A1 == 1)  
  }
  
  
  d <- d[order(d$id), ]
  rownames(d) <- 1:n
  d.full <- d
  if (old) d.full <- NULL
  d <- subset(d, select = c("id", "Y0", "A1", "Y1", "R", "A2R",
                            "A2NR", "Y2", "weight",
                            grep("dtr", names(d), value = T)))
  
  class(d) <- c("SMART", class(d))
  
  # Check validity of treatment allocations
  if (!validTrial(d, design)) {
    return(list(
      "data" = d,
      "valid" = F,
      "params" = list(
        "times" = times,
        "spltime" = spltime,
        "r1" = r1,
        "r0" = r0,
        "gammas" = gammas,
        "lambdas" = lambdas,
        "design" = design,
        "sigma.r0" = sigma.r0,
        "sigma.r1" = sigma.r1
      )
      # "assumptions" = list("conditionalVariation" = NULL)
    ))
  }
  
  output <- list(
    "data" = d,
    "potentialData" = d.full,
    "valid" = T,
    "params" = list(
      "times" = times,
      "spltime" = spltime,
      "r1" = r1,
      "r0" = r0,
      "gammas" = gammas,
      "lambdas" = lambdas,
      "design" = design,
      "sigma.r0" = sigma.r0,
      "sigma.r1" = sigma.r1
    )
    # "assumptions" = list("conditionalVariation" = condVarAssump)
  )
  class(output) <- c('generateSMART', class(output))
  return(output)
}
