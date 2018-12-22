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


#' Generate data from sequential, multiple-assignment, randomized trial
#'
#' @param n number of participants in the trial
#' @param times vector of observation times at which the repeated-measures outcome is collected
#' @param spltime the observation time immediately before re-randomization
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
                          respDirection = NULL,
                          respFunction,
                          empirical = FALSE, causal = FALSE) {
  
  # times:    Vector of times at which measurements are collected (currently limited to length three)
  # spltime:  Time (contained in times) at which re-randomization occurs
  # sigma:    Marginal variance of Y (assumed constant over time and DTR)
  # sigma.r1: Vector of conditional variances of Y for responders to treatment A1 = 1 in the second stage 
  # sigma.r0: Vector of conditional variances of Y for responders to treatment A1 = -1 in the second stage
  # corstr:   True correlation structure
  
  # uneqvar:  Vector of same length as uneqsdDTR. If not all DTRs have the same variance (i.e., sigma^2),
  #           this is the sigma
  # rho:      If corstr is either "exch" or "exchangeable", the within-person correlation used to compute V 
  #           (assumed constant across DTRs)
  # rho.r1:   Vector of conditional correlations between Ys for responders to treatment A1 = 1 in the second stage
  # rho.r0:   Exchangeable correlation for responders to A1 == -1
  
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
  
  ## Generate treatment allocations and response
  d <- data.frame("id"   = 1:n,
                  "Y0"   = NA,
                  "A1"   = 0)
  
  d$Y0 <- gammas[1] + rnorm(n, 0, sigma)
  
  if (balanceRand) {
    s <- sample(1:n, size = floor(n/2), replace = FALSE)
    d$A1[s]  <- 1
    d$A1[-s] <- -1
  } else {
    d$A1 <- 2 * rbinom(n, 1, 0.5) - 1
  }
  
  d$Y1 <- (1 - rho) * gammas[1] + rho * d$Y0 + gammas[2] + 
    gammas[3] * d$A1 + rnorm(n, 0, sqrt((1 - rho^2) * sigma^2))
  
  d$R <- rep(NA, n)
  
  d <- respFunction(d, gammas, r1, r0, respDirection)
  
  # respModel <- match.arg(respModel)
  # respDirection <- match.arg(respDirection)
  # 
  # if (respModel == "independent") {
  #   d$R[d$A1 ==  1] <- rbinom(sum(d$A1 == 1), 1, r1)
  #   d$R[d$A1 == -1] <- rbinom(sum(d$A1 == -1), 1, r0)
  # } else if (respModel == "oneThreshold") {
  #   if (design != 2) stop("oneThreshold is not yet implemented for any design other than 2")
  #   tail <- switch(respDirection, "high" = F, "low" = T)
  #   upsilon <- qnorm(r1, as.numeric(c(rep(1, 3), rep(0, 4)) %*% gammas), sigma, lower.tail = tail)
  #   d$R <- as.numeric(d$Y1 >= upsilon)
  #   r0temp <- pnorm(upsilon, sum(gammas[1:2] - gammas[3]), sigma, lower.tail = FALSE)
  #   if (r0temp != r0) {
  #     warning(paste("Overwriting the provided value of r0 to accomodate the oneThreshold respModel.",
  #                   "The provided value is", r0, "and the new value is", round(r0temp, 3)))
  #     r0 <- r0temp
  #   }
  # } else if (respModel == "twoThreshold") {
  #   d <- response.twoT(d, gammas, r1, r0, respDirection)
  # } else if (respModel == "beta") {
  #   shape1 <- r1 / (1 - r1)
  #   shape0 <- r0 / (1 - r0)
  #   x <- pnorm(d$Y1, mean = gammas[1] + gammas[2] + gammas[3]*d$A1, sd = sigma)
  #   respProb <- vector("numeric", n)
  #   respProb[d$A1 ==  1] <- qbeta(x[d$A1 ==  1], shape1 = shape1, shape2 = 1)
  #   respProb[d$A1 == -1] <- qbeta(x[d$A1 == -1], shape1 = shape0, shape2 = 1)
  #   d$respProb <- respProb
  #   d$R <- sapply(1:n, function(i) rbinom(1, 1, respProb[i]))
  # }
  
  # ## Generate Y's at their (conditional) means
  # d <- do.call("rbind", lapply(split.SMART(d), function(z) {
  #   x <- do.call("data.frame", lapply(times, function(time) {
  #     conditional.model(z$A1, z$R, z$A2R, z$A2NR,
  #                       time, spltime, design, 
  #                       get(paste0("r", (1 + unique(z$A1)) / 2)), 
  #                       gammas, lambdas)
  #   }))
  #   names(x) <- sapply(times, function(t) paste0("Y", t))
  #   z <- cbind(z, x)
  # }))
  # d <- d[order(d$id), ]
  # rownames(d) <- d$id
  # 
  # ## Compute variance matrices if not already provided
  # if (is.null(varmats)) {
  #   varmats <- conditionalVarmat(times, spltime, design, r1, r0,
  #                                corstr, sigma, sigma.r1, sigma.r0,
  #                                uneqsd = uneqsd, uneqsdDTR = uneqsdDTR,
  #                                rho, rho.r1, rho.r0,
  #                                gammas, lambdas)
  # }
  # 
  # ## Add noise with variance according to varmats
  # d <- do.call("rbind",
  #              lapply(split.SMART(d), function(x) {
  #                varIndex <- as.character(paste(unique(x$A1), unique(x$R), 
  #                                  unique(x$A2R), unique(x$A2NR)))
  #                x[, grepl("Y", names(x))] <- x[, grepl("Y", names(x))] +
  #                  mvrnorm(
  #                    n = nrow(x),
  #                    mu = rep(0, length(times)),
  #                    Sigma = varmats[[varIndex]],
  #                    empirical = empirical
  #                  )
  #                x
  #              })
  # )
  
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
    
    # generate
    sigmaCorrection <- 2 * sigma^2 * (rho^2 / (1 + rho))
    
    ### FIXME: THIS DOESN'T INCORPORATE SPLTIME!!
    d$Y2 <- (1 - rho)/(1 + rho)*gammas[1] + spltime/(1 + rho)*(gammas[2] + gammas[3]*d$A1) +
      rho/(1 + rho)*(d$Y0 + d$Y1) + 
      (gammas[4] + gammas[5]*d$A1 + 
                                (gammas[6]*d$A2NR + gammas[7]*d$A1*d$A2NR)/(1 - ifelse(d$A1 == 1, r1, r0)) + 
                                (d$R - ifelse(d$A1 == 1, r1, r0))*(lambdas[1] + lambdas[2]*d$A1))
    
    d$Y2[d$A1 == 1 & d$R == 1] <- d$Y2[d$A1 == 1 & d$R == 1] +
      rnorm(sum(d$A1 == 1 & d$R == 1), 0, sd = sqrt(sigma.r1^2 - sigmaCorrection))
    d$Y2[d$A1 == -1 & d$R == 1] <- d$Y2[d$A1 == -1 & d$R == 1] +
      rnorm(sum(d$A1 == -1 & d$R == 1), 0, sd = sqrt(sigma.r0^2 - sigmaCorrection))
    
    sigma.nr11 <- sqrt((1/(1 - r1)) * (sigma^2 - r1*sigma.r1^2 - (1 - r1) * sigmaCorrection -
                                         r1/(1 - r1)*(sum(gammas[6:7])/(1 - r1) - sum(lambdas))^2))
    sigma.nr10 <- sqrt((1/(1 - r1)) * (sigma^2 - r1*sigma.r1^2 - (1 - r1) * sigmaCorrection - 
                                         r1/(1 - r1)*((-sum(gammas[6:7]))/(1 - r1) - sum(lambdas))^2))
    
    d$Y2[d$A1 == 1 & d$R == 0 & d$A2NR == 1] <- d$Y2[d$A1 == 1 & d$R == 0 & d$A2NR == 1] +
      rnorm(sum(d$A1 == 1 & d$R == 0 & d$A2NR == 1), 0, sigma.nr11)
    d$Y2[d$A1 == 1 & d$R == 0 & d$A2NR == -1] <- d$Y2[d$A1 == 1 & d$R == 0 & d$A2NR == -1] +
      rnorm(sum(d$A1 == 1 & d$R == 0 & d$A2NR == -1), 0, sigma.nr10)
    
    sigma.nr01 <- sqrt((1/(1 - r0)) * (sigma^2 - r0*sigma.r0^2 - (1 - r0)*sigmaCorrection -
                                         r0/(1 - r0)*((gammas[6] - gammas[7])/(1 - r0) - lambdas[1] + lambdas[2])^2))
    sigma.nr00 <- sqrt((1/(1 - r0)) * (sigma^2 - r0*sigma.r0^2 - (1 - r0)*sigmaCorrection -
                                         r0/(1 - r0)*((-gammas[6] + gammas[7])/(1 - r0) - lambdas[1] + lambdas[2])^2))
    
    d$Y2[d$A1 == -1 & d$R == 0 & d$A2NR == 1] <- d$Y2[d$A1 == -1 & d$R == 0 & d$A2NR == 1] +
      rnorm(sum(d$A1 == -1 & d$R == 0 & d$A2NR == 1), 0, sigma.nr01)
    d$Y2[d$A1 == -1 & d$R == 0 & d$A2NR == -1] <- d$Y2[d$A1 == -1 & d$R == 0 & d$A2NR == -1] +
      rnorm(sum(d$A1 == -1 & d$R == 0 & d$A2NR == -1), 0, sigma.nr00)
    
    d$weight <- 2 * (d$R + 2 * (1 - d$R))
    
    d$dtr1 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR ==  1))))
    d$dtr2 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR == -1))))
    d$dtr3 <- as.numeric(with(d, (A1 == -1) * (R + (1 - R) * (A2NR ==  1))))
    d$dtr4 <- as.numeric(with(d, (A1 == -1) * (R + (1 - R) * (A2NR == -1))))
  } else {
    d$A2NR[d$A1 == 1 & d$R == 0] <- 2 * rbinom(sum(d$A1 == 1 & d$R == 0), 1, .5) - 1
    d$weight <- 2 + 2 * (1 - d$R) * (d$A1 == 1)
    
    d$dtr1 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR ==  1))))
    d$dtr2 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR == -1))))
    d$dtr3 <- as.numeric(with(d, (A1 == -1)))
  }
  
  class(d) <- c("SMART", class(d))
  
  # Check validity of treatment allocations
  if (!validTrial(d, design)) {
    return(list("data" = d, "valid" = F, "times" = times, "spltime" = spltime))
  }
  d <- d[order(d$id), ]
  rownames(d) <- 1:n
  return(list("data" = d, "valid" = T, "times" = times, "spltime" = spltime))
}
