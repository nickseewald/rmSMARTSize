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