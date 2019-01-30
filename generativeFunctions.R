## Data Generation Helper Functions

generateStage1 <- function(n, times, spltime, r1, r0, gammas,
                            sigma, corstr = c("identity", "exchangeable", "ar1"),
                            rho = 0, 
                            respFunction, respDirection = NULL, 
                            balanceRand = FALSE,
                            empirical = FALSE, old = FALSE) {
  if (old) 
    warning("Option 'old' currently does nothing in generateStage1()")
  
  corstr <- match.arg(corstr)
  
  ## Initialize data frame
  try(assign("d", data.frame("id" = 1:n)))
  
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
    # FIXME: This only works for time 1 right now! Any other time will not have proper correlation structure
    d[, paste0("Y", times[times > 0 & times <= spltime], ".0")] <-
      (1-rho)*gammas[1] + rho*d$Y0 + gammas[2] - gammas[3] + rnorm(n, 0, sqrt((1-rho^2))*sigma)
    d[, paste0("Y", times[times > 0 & times <= spltime], ".1")] <-
      (1-rho)*gammas[1] + rho*d$Y0 + gammas[2] + gammas[3] + rnorm(n, 0, sqrt((1-rho^2))*sigma)
  } else if (corstr == "ar1") {
    # FIXME: This only works for time 1 right now! Any other time will not have proper correlation structure
    d[[paste0("Y", time, ".0")]] <- (1-rho)*gammas[1] + rho*d$Y0 + gammas[2] - gammas[3] + rnorm(n, 0, sqrt((1-rho^2))*sigma)
    d[[paste0("Y", time, ".1")]] <- (1-rho)*gammas[1] + rho*d$Y0 + gammas[2] + gammas[3] + rnorm(n, 0, sqrt((1-rho^2))*sigma)
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
  
  class(d) <- c("stage1SMART", class(d))
  d
}

generateStage2.means <- function(d, times, spltime, gammas, lambdas,
                                 design, corstr = c("identity", "exchangeable", "ar1")) {
  if (!("stage1SMART") %in% class(d))
    stop("'d' must be of class 'stage1SMART'")
  
  n <- nrow(d)
  nDTR <- switch(design, 8, 4, 3)
  DTRs <- dtrIndex(design)
  
  dtrTriples <- as.data.frame((do.call(cbind, DTRs) + 1) / 2)
  dtrTriples[dtrTriples == 0.5] <- 0
  dtrTriples <- apply(dtrTriples, 1, function(x) paste(x, collapse = ""))
  
  dtrTxts <- do.call(rbind, DTRs)
  
  # Replicate the "potential" response status variables for each DTR,
  # so that each DTR gets the appropriate potential response
  respMatrix <- do.call(cbind, lapply((1+dtrTxts[1, ])/2, function(a1) d[[paste0("R.", a1)]]))
  
  # Create a vector of response probabilities, each entry of which corresponds to 
  # the appropriate value (r0 or r1) for each DTR
  respProbVec <- sapply(paste0("r", (1+dtrTxts[1, ])/2), get)
  
  # Create matrices of the "adjustments" added to the mean model that are 
  # specific to responders and non-responders
  if (design == 1) {
    designMeanAdj.R <- matrix(rep(
      dtrTxts[2, ] * (gammas[6] + gammas[8] * dtrTxts[1, ]) / respProbVec +
        (1 - respProbVec) * (lambdas[1] + lambdas[2] * dtrTxts[1, ]),
      n), nrow = n, byrow = T)
    designMeanAdj.NR <- matrix(rep(
      dtrTxts[3, ] * (gammas[7] + gammas[9] * dtrTxts[1, ]) / (1 - respProbVec) -
        respProbVec * (lambdas[1] + lambdas[2] * dtrTxts[1, ]),
      n), nrow = n, byrow = T)
  } else if (design == 2) {
    designMeanAdj.R <- matrix(rep(
      (1 - respProbVec) * (lambdas[1] + lambdas[2] * dtrTxts[1, ]),
      n), nrow = n, byrow = T)
    designMeanAdj.NR <- matrix(rep(
      dtrTxts[3, ] * (gammas[6] + gammas[7] * dtrTxts[1, ]) / (1 - respProbVec) -
        respProbVec * (lambdas[1] + lambdas[2] * dtrTxts[1, ]),
      n), nrow = n, byrow = T)
  } else if (design == 3) {
    
  } else stop("'design' must be one of 1, 2, or 3.")
  
  if (corstr %in% c("identity", "exchangeable")) {
    # FIXME: This only works for time 2 right now! Any other time will not have proper correlation structure
    d[, paste0("Y", times[times > spltime], ".", dtrTriples)] <- 
      # Start with a matrix of just the appropriate linear combinations of gammas
      # (only up to the stage 1 parts of the time 2 model)
      matrix(rep(
        ((1-rho)/(1+rho)) * gammas[1] + (gammas[2] + gammas[3] * dtrTxts[1, ])/(1+rho) +
          gammas[4] + gammas[5] * dtrTxts[1, ],
        n), nrow = n, byrow = T) + 
      # Add Y0 and the appropriate causal Y1
      (rho/(1+rho)) * (d$Y0 + 
                         matrix(c(rep(d$Y1.1, sum(dtrTxts[1, ] == 1)),
                                      rep(d$Y1.0, sum(dtrTxts[1, ] == -1))), nrow = n)
                       ) +
      # Add response adjustment
      respMatrix * designMeanAdj.R + (1 - respMatrix) * designMeanAdj.NR
  } else if (corstr == "ar1") {
    stop("generateStage2 not yet implemented for corstr = 'ar1'")
  }
  d
}