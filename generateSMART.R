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
generateSMART <- function(n, times, spltime, r1, r0 = NULL, gammas, lambdas, design, balanceRand = FALSE,
                          sigma, sigma.r1, sigma.r0, corstr = c("identity", "exchangeable", "ar1"),
                          uneqsdDTR = NULL, uneqsd = NULL, varmats = NULL,
                          respModel = c("independent", "oneThreshold", "twoThreshold"),
                          rho = NULL, rho.r1 = rho, rho.r0 = rho, empirical = FALSE) {
  
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
                  "Y0"   = gammas[1] + rnorm(n, 0, sigma),
                  "A1"   = 0,
                  "R"    = NA,
                  "A2R"  = 0,
                  "A2NR" = 0)
  
  if (balanceRand) {
    s <- sample(1:n, size = floor(n/2), replace = FALSE)
    d$A1[s]  <- 1
    d$A1[-s] <- -1
  } else {
    d$A1 <- 2 * rbinom(n, 1, 0.5) - 1
  }
  
  d$Y1 <- (1 - rho) * gammas[1] + rho * d$Y0 + gammas[2] + 
    gammas[3] * d$A1 + rnorm(n, 0, sqrt((1 - rho^2) * sigma^2))
  
  respModel <- match.arg(respModel)
  if (respModel == "independent") {
    d$R[d$A1 ==  1] <- rbinom(sum(d$A1 == 1), 1, r1)
    d$R[d$A1 == -1] <- rbinom(sum(d$A1 == -1), 1, r0)
  } else if (respModel == "oneThreshold") {
    if (design != 2) stop("oneThreshold is not yet implemented for any design other than 2")
    upsilon <- qnorm(r1, as.numeric(c(rep(1, 3), rep(0, 4)) %*% gammas), sigma, lower.tail = F)
    d$R <- as.numeric(d$Y1 >= upsilon)
    warning("Overwriting the provided value of r0 to accomodate the oneThreshold respModel.")
    r0 <- pnorm(upsilon, sum(gammas[1:2] - gammas[3]), sigma, lower.tail = FALSE)
  } else if (respModel == "twoThreshold") {
    if (design != 2) stop("twoThreshold is not yet implemented for any design other than 2")
    upsilon1 <- qnorm(r1, as.numeric(c(rep(1, 3), rep(0, 4)) %*% gammas), sigma, lower.tail = F)
    upsilon0 <- qnorm(r0, as.numeric(c(c(1, 1, -1), rep(0, 4)) %*% gammas), sigma, lower.tail = F)
    d$R[d$A1 ==  1] <- as.numeric(d$Y1[d$A1 ==  1] >= upsilon1)
    d$R[d$A1 == -1] <- as.numeric(d$Y1[d$A1 == -1] >= upsilon0)
  }
  
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
     
    d$Y2 <- (1 - rho)/(1 + rho) * gammas[1] + 1/(1 + rho)*(gammas[2] + gammas[3]*d$A1) +
      rho/(1 + rho)*d$Y0 + rho/(1 + rho) * d$Y1 + gammas[4] + gammas[5] * d$A1 +
      (gammas[6] * d$A2NR + gammas[7] * d$A1 * d$A2NR) /
      (1 - ifelse(d$A1 == 1, r1, r0))
    
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
  
  ## Check validity of treatment allocations
  if (!validTrial(d, design)) {
    return(list("data" = d, "valid" = F))
  }
  
  sigmaCorrection <- 2 * (rho^2 / (1 + rho)) * sigma^2
  
  d$Y2[d$A1 == 1 & d$R == 1] <- d$Y2[d$A1 == 1 & d$R == 1] +
    rnorm(sum(d$A1 == 1 & d$R == 1), 0, sd = sqrt(sigma.r1^2 - sigmaCorrection))
  d$Y2[d$A1 == -1 & d$R == 1] <- d$Y2[d$A1 == -1 & d$R == 1] +
    rnorm(sum(d$A1 == -1 & d$R == 1), 0, sd = sqrt(sigma.r0^2 - sigmaCorrection))
  
  sigma.nr11 <- sqrt((1/(1 - r1)) * (sigma^2 - r1*sigma.r1^2 - (1 - r1) * sigmaCorrection - r1/(1 - r1)*(sum(gammas[6:7])^2)))
  sigma.nr10 <- sqrt((1/(1 - r1)) * (sigma^2 - r1*sigma.r1^2 - (1 - r1) * sigmaCorrection - r1/(1 - r1)*(gammas[6] - gammas[7])^2))
  d$Y2[d$A1 == 1 & d$R == 0 & d$A2NR == 1] <- d$Y2[d$A1 == 1 & d$R == 0 & d$A2NR == 1] +
    rnorm(sum(d$A1 == 1 & d$R == 0 & d$A2NR == 1), 0, sigma.nr11)
  d$Y2[d$A1 == 1 & d$R == 0 & d$A2NR == -1] <- d$Y2[d$A1 == 1 & d$R == 0 & d$A2NR == -1] +
    rnorm(sum(d$A1 == 1 & d$R == 0 & d$A2NR == -1), 0, sigma.nr10)
  
  sigma.nr01 <- sqrt((1/(1 - r0)) * (sigma^2 - r0*sigma.r0^2 - (1 - r0)*sigmaCorrection - r0/(1 - r0)*(-sum(gammas[6:7]))^2))
  sigma.nr00 <- sqrt((1/(1 - r0)) * (sigma^2 - r0*sigma.r0^2 - (1 - r0)*sigmaCorrection - r0/(1 - r0)*(-gammas[6] + gammas[7])^2))
  d$Y2[d$A1 == -1 & d$R == 0 & d$A2NR == 1] <- d$Y2[d$A1 == -1 & d$R == 0 & d$A2NR == 1] +
    rnorm(sum(d$A1 == -1 & d$R == 0 & d$A2NR == 1), 0, sigma.nr01)
  d$Y2[d$A1 == -1 & d$R == 0 & d$A2NR == -1] <- d$Y2[d$A1 == -1 & d$R == 0 & d$A2NR == -1] +
    rnorm(sum(d$A1 == -1 & d$R == 0 & d$A2NR == -1), 0, sigma.nr00)
  
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
  
  d <- d[order(d$id), ]
  rownames(d) <- 1:n
  return(list("data" = d, "valid" = T))
}
