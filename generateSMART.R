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
#' @param rho 
#' @param rho.r1 
#' @param rho.r0 
#' @param empirical 
#'
#' @return A data.frame
#' @export
#'
#' @examples
generateSMART <- function(n, times, spltime, r1, r0, gammas, lambdas, design, balanceRand = FALSE,
                          sigma, sigma.r1, sigma.r0, corstr = c("identity", "exchangeable", "ar1"),
                          uneqsdDTR = NULL, uneqsd = NULL, varmats = NULL,
                          respModel = NULL,
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
  if (!(design %in% 1:3)) stop("Invalid design choice. Must be one of 1-3.")
  
  ## Generate treatment allocations and response
  d <- data.frame("id"   = 1:n,
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
  
  # if (!is.null(respModel)) {
  #   if (length(respModel) != 3)
  #     stop("If providing a respModel model for response, it must be a vector of length 3.")
  #   respEst <- respModel[1] + respModel[2] * d$Y0
  # } else {
    d$R[d$A1 ==  1] <- rbinom(sum(d$A1 == 1), 1, r1)
    d$R[d$A1 == -1] <- rbinom(sum(d$A1 == -1), 1, r0)
  # }
  
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
    d$A2NR[d$A1 == 1  & d$R == 0]  <- 2 * rbinom(sum(d$A1 == 1  & d$R == 0), 1, .5) - 1
    d$A2NR[d$A1 == -1 & d$R == 0]  <- 2 * rbinom(sum(d$A1 == -1 & d$R == 0), 1, .5) - 1
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
  
  ## Check validity of treatment allocations
  if (!validTrial(d, design)) {
    return(list("data" = d, "valid" = F))
  }
  
  class(d) <- c("SMART", class(d))
  
  ## Generate Y's at their (conditional) means
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
  
  d <- d[order(d$id), ]
  rownames(d) <- 1:n
  return(list("data" = d, "valid" = T))
}
