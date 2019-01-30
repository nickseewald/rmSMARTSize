## Data Generation Helper Functions

generateStage1 <- function(n, times, spltime, r1, r0, gammas,
                            sigma, corstr = c("identity", "exchangeable", "ar1"),
                            rho = 0, 
                            respFunction, respDirection = NULL, 
                            balanceRand = FALSE,
                            empirical = FALSE, old = FALSE) {
  
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
  
  d
}

generateStage2.means <- function(d, times, spltime, gammas, lambdas,
                                 design, corstr = c("identity", "exchangeable", "ar1")) {
  
  nDTR <- switch(design, 8, 4, 3)
  DTRs <- dtrIndex(design)
  
  dtrTriples <- as.data.frame((do.call(cbind, DTRs) + 1) / 2)
  dtrTriples[dtrTriples == 0.5] <- 0
  dtrTriples <- apply(dtrTriples, 1, function(x) paste(x, collapse = ""))
  
  dtrTxts <- do.call(rbind, DTRs)
  
  respMatrix <- do.call(cbind, lapply((1+dtrTxts[1, ])/2, function(a1) d[[paste0("R.", a1)]]))
  
  # Create a vector of response probabilities, each entry of which corresponds to 
  # the appropriate value (r0 or r1) for each DTR
  respProbVec <- sapply(paste0("r", (1+dtrTxts[1, ])/2), get)
  
  if (design == 1) {
    designMeanAdj.R <- matrix(rep(
      dtrTxts[2, ] * (gammas[6] + gammas[8] * dtrTxts[1, ]) / respProbVec +
        (1 - respProbVec) * (lambdas[1] + lambdas[2] * dtrTxts[1, ]),
      n), nrow = n, byrow = T)
    designMeanAdj.NR <- matrix(rep(
      dtrTxts[3, ] * (gammas[7] + gammas[9] * dtrTxts[1, ]) / (1 - respProbVec) -
        respProbVec * (lambdas[1] + lambdas[2] * dtrTxts[1, ]),
      n), nrow = n, byrow = T)
  }
  
  if (corstr %in% c("identity", "exchangeable")) {
    # FIXME: This only works for time 2 right now! Any other time will not have proper correlation structure
    d[, paste0("Y", times[times > spltime], ".", dtrTriples)] <- 
      # Start with a matrix of just the appropriate linear combinations of gammas
      # (only up to the stage 1 parts of the time 2 model)
      matrix(rep(
        ((1-rho) / (1+rho)) * gammas[1] + (gammas[2] + gammas[3] * dtrTxts[1, ])/(1+rho) +
          gammas[4] + gammas[5] * dtrTxts[1, ],
        n), nrow = n, byrow = T) + 
      # Add Y0 and the appropriate causal Y1
      (rho/(1+rho)) * (d$Y0 + matrix(c(rep(d$Y1.1, 4), rep(d$Y1.0, 4)), nrow = n)) +
      # Add 
      respMatrix * designMeanAdj.R + (1 - respMatrix) * designMeanAdj.NR
    
    

    
  } else if (corstr == "ar1") {
    stop("generateStage2 not yet implemented for corstr = 'ar1'")
  }
  
  
  # for (i in 1:nDTR){
  #   d[[paste0("Y2.", dtrTriples[1])]] <- 
  # }
  
  
  ## Time 2 means
  d$Y2.11 <- (1-rho)/(1+rho) * gammas[1] + (gammas[2] + gammas[3]) / (1+rho) +
    rho/(1+rho) * (d$Y0 + d$Y1.1) + gammas[4] + gammas[5] + (1-d$R.1)*(gammas[6] + gammas[7])/(1-r1) +
    (d$R.1 - r1) * (lambdas[1] + lambdas[2]) 
  
  d$Y2.10 <- (1-rho)/(1+rho) * gammas[1] + (gammas[2] + gammas[3]) / (1+rho) +
    rho/(1+rho) * (d$Y0 + d$Y1.1) + gammas[4] + gammas[5] + (1-d$R.1)*(-gammas[6] - gammas[7])/(1-r1) +
    (d$R.1 - r1) * (lambdas[1] + lambdas[2])
  
  d$Y2.01 <- (1-rho)/(1+rho) * gammas[1] + (gammas[2] - gammas[3]) / (1+rho) +
    rho/(1+rho) * (d$Y0 + d$Y1.0) + gammas[4] - gammas[5] + (1-d$R.0)*(gammas[6] - gammas[7])/(1-r0) +
    (d$R.0 - r0) * (lambdas[1] - lambdas[2])
  
  d$Y2.00 <- (1-rho)/(1+rho) * gammas[1] + spltime * (gammas[2] - gammas[3]) / (1+rho) +
    rho/(1+rho) * (d$Y0 + d$Y1.0) + gammas[4] - gammas[5] + (1-d$R.0)*(-gammas[6] + gammas[7])/(1-r0) +
    (d$R.0 - r0) * (lambdas[1] - lambdas[2])
  
  
}