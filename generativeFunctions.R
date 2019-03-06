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
      
      if (length(gammas) != 9) 
        stop("For design 1, gammas must be of length 9.")
      
      designMeanAdj.R <- matrix(rep(
        dtrTxts[2,] * (gammas[6] + gammas[8] * dtrTxts[1, ]) / respProbVec +
          (1 - respProbVec) * (lambdas[1] + lambdas[2] * dtrTxts[1, ]),
        n
      ),
      nrow = n,
      byrow = T)
      designMeanAdj.NR <- matrix(rep(
        dtrTxts[3,] * (gammas[7] + gammas[9] * dtrTxts[1, ]) / (1 - respProbVec) -
          respProbVec * (lambdas[1] + lambdas[2] * dtrTxts[1, ]),
        n
      ),
      nrow = n,
      byrow = T)
    } else if (design == 2) {
      
      if (length(gammas) != 7)
        stop("For design 2, gammas must be of length 7.")
      
      designMeanAdj.R <-  matrix(rep(
        (1 - respProbVec) * (lambdas[1] + lambdas[2] * dtrTxts[1, ]),
        n
        ), nrow = n, byrow = T)
      
      designMeanAdj.NR <- matrix(rep(
        dtrTxts[3,] * (gammas[6] + gammas[7] * dtrTxts[1, ]) / (1 - respProbVec) -
          respProbVec * (lambdas[1] + lambdas[2] * dtrTxts[1, ]),
        n
      ),
      nrow = n,
      byrow = T)
    } else if (design == 3) {
      
      if (length(gammas) != 6)
        stop("For design 3, gammas must be of length 6.")
      
      designMeanAdj.R <- matrix(rep(
        (1 - respProbVec) * (lambdas[1] + lambdas[2] * dtrTxts[1, ]),
        n
        ), nrow = n, byrow = T)
      
      designMeanAdj.NR <- matrix(rep(
        gammas[6] * as.numeric(dtrTxts[1, ] == 1) * dtrTxts[3, ] / 
          (1 - respProbVec) - 
          respProbVec * (lambdas[1] + lambdas[2] * dtrTxts[1, ]),
        n
      ), nrow = n, byrow = T)
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


generateStage2.var <- function() {
  sigma.nr00 <- 
    sqrt((sigma^2 - r0 * sigma.r0^2 - r0*(1-r0) * 
            with(d, mean(Y2.000[R.0 == 1]) -
                   mean(Y2.000[R.0 == 0]))^2) / (1 - r0))
  sigma.nr01 <- 
    sqrt((sigma^2 - r0 * sigma.r0^2 - r0*(1-r0) *
            with(d, mean(Y2.001[R.0 == 1]) - 
                   mean(Y2.001[R.0 == 0]))^2) / (1 - r0))
  sigma.nr10 <- 
    sqrt((sigma^2 - r1 * sigma.r1^2 - 
            r1*(1-r1) * with(d, mean(Y2.100[R.1 == 1]) -
                               mean(Y2.100[R.1 == 0]))^2) / (1 - r1))
  sigma.nr11 <- 
    sqrt((sigma^2 - r1 * sigma.r1^2 - 
            r1*(1-r1) * with(d, mean(Y2.101[R.1 == 1]) -
                               mean(Y2.101[R.1 == 0]))^2) / (1 - r1))
}


testVarianceInput <- function(x, a, d, r0, r1, rho, resp = c("nr", "r")) {
  #TODO: is this design-agnostic?
  #TODO: What happens when this is always positive?
  
  resp <- tolower(resp)
  resp <- match.arg(resp)
  
  y1string <- paste0("Y1.", a[1])
  y2string <- paste0("Y2.", paste0(a, collapse = ""))
  rstring <- paste0("R.", a[1])
  rprobstring <- paste0("r", a[1])
  
  if (resp == "nr") {
    sigma.nr <- 
      (sigma^2 - get(rprobstring) * x^2 - 
         get(rprobstring)*(1-get(rprobstring)) * 
         with(d, mean(get(y2string)[get(rstring) == 1]) -
                mean(get(y2string)[get(rstring) == 0]))^2) /
      (1 - get(rprobstring))
    # ))
    sigma.nr - (rho / (1 + rho))^2 * with(subset(d, get(rstring) == 0), 
                                          var(Y0 + get(paste0("Y1.", a[1]))))
  } else {
    x^2 - 
      (rho / (1 + rho))^2 * 
      with(subset(d, get(rstring) == 1), var(Y0 + get(y1string)))
  }
}


#' Compute variance bounds for simulation scenario
#'
#' @param a1 
#' @param d 
#' @param design 
#' @param sigma 
#' @param r 
#' @param rho 
#' @param times 
#' @param corstr 
#' @param bound 
#' @param boundType A `character` string, either `feasibility` or `assumption`.
#' Ignored for `bound == "upper"`. See Details.
#'
#' @return
#' @export
#'
#' @examples
computeVarBound <- function(a1, d, design, sigma, r, rho = 0, times,
                            corstr = c("identity", "exchangeable", "ar1"),
                            bound = c("lower", "upper"),
                            boundType = c("feasibility", "assumption")) {
 
  bound     <- match.arg(bound)
  boundType <- match.arg(boundType)
  corstr    <- match.arg(corstr)
  
  if (corstr == "identity")
    corstr <- "exchangeable"
  
  sigma <- reshapeSigma(sigma, times, design)
  
  dtrs <- dtrNames(design)
  
  responderData <- subset(d, get(paste0("R.", a1)) == 1)
  nonResponderData <- subset(d, get(paste0("R.", a1)) == 0)
  
  # For each design and first-stage treatment, specify "reference" DTRs from 
  # which we can get mean outcomes for responders and non-responders which are
  # consistent with the means needed for the variance bounds
  if (design == 1) {
    if (a1 == 1) {
      dtr1 <- "111"
      dtr2 <- "101"
    } else if (a1 %in% c(0, -1)) {
      # FIXME: I think....
      dtr1 <- "000"
      dtr2 <- "001"
    } else {
      stop("Invalid choice of a1: must be one of 0, 1, or -1")
    }
  } else if (design == 2) {
    if (a1 == 1) {
      dtr1 <- "101"
      dtr2 <- "100"
    } else if (a1 %in% c(0, -1)) {
      dtr1 <- "000"
      dtr2 <- "001"
    } else {
      stop("Invalid choice of a1: must be one of 0, 1, or -1")
    }
  } else if (design == 3) {
    if (a1 == 1) {
      dtr1 <- "101"
      dtr2 <- "100"
    } else if (a1 %in% c(0, -1)) {
      dtr1 <- "000"
      dtr2 <- "000"
    }
  }
  
  dtrCol1 <- match(dtr1, dtrs)
  dtrCol2 <- match(dtr2, dtrs)
  
  c0 <- rho * sigma[2, dtrCol1] * sigma[3] / 
    (sigma[1] * (rho * sigma[1] + sigma[2]))
  c1 <- rho * sigma[3] / (rho * sigma[1] + sigma[2])
  
  
  stage1var <- with(responderData, var(c0 * Y0 + c1 * get(paste0("Y1.", a1))))
  
  Y2.dtr1    <- with(d,                get(paste0("Y2.", dtr1)))
  Y2.R.dtr1  <- with(responderData,    get(paste0("Y2.", dtr1)))
  Y2.NR.dtr1 <- with(nonResponderData, get(paste0("Y2.", dtr1)))
  
  Y2.dtr2    <- with(d,                get(paste0("Y2.", dtr2)))
  Y2.R.dtr2  <- with(responderData,    get(paste0("Y2.", dtr2)))
  Y2.NR.dtr2 <- with(nonResponderData, get(paste0("Y2.", dtr2)))
  
  if (bound == "lower") {
    if (boundType == "feasibility") {
      sqrt(max(c(
        0.1,
        stage1var,
        (1 - r) * ((mean(Y2.R.dtr2) - mean(Y2.NR.dtr2))^2 -
                     (mean(Y2.R.dtr1) - mean(Y2.NR.dtr1))^2 + stage1var)
      )))
    } else if (boundType == "assumption") {
      sqrt(max(c(
        0.1,
        sigma^2 + ((1 - r) / r) * (mean(Y2.NR.dtr1) - mean(Y2.dtr1))^2 -
          (1 - r) * (mean(Y2.R.dtr1) - mean(Y2.NR.dtr1))^2,
        sigma^2 + ((1 - r) / r) * (mean(Y2.NR.dtr2) - mean(Y2.dtr2))^2 -
          (1 - r) * (mean(Y2.R.dtr2) - mean(Y2.NR.dtr2))^2
      )))
    }
  } else {
    sqrt(max(
      c(0.1,
        sigma^2 / r - (1 - r) * 
          max(c(
            (mean(Y2.R.dtr1) - mean(Y2.NR.dtr1))^2, 
            (mean(Y2.R.dtr2) - mean(Y2.NR.dtr2))^2)) - stage1var)
    ))
  }
}

#' Title
#'
#' @param simGrid 
#' @param times 
#' @param spltime 
#' @param gammas 
#' @param sigma 
#' @param corstr 
#' @param design 
#' @param balanceRand 
#' @param varCombine A function of one vector argument which describes how to
#' choose sigma.rX from the lower and upper bounds. It's assumed that the
#' function will be able to handle vectors of length two in which the elements
#' are c(`<lower bound>`, `<upper bound>`). Defaults to `mean`. 
#' @param empirical 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
computeVarGrid <- function(simGrid, times, spltime, gammas, sigma,
                           corstr, design, balanceRand = F,
                           varCombine = mean,
                           empirical = F, seed = 6781) {
  
  # if (design == 3) stop("computeVarGrid doesn't yet work for design 3")
  
  if (is.list(varCombine) & length(varCombine) == 2) {
    varCombine.11 <- varCombine[[1]]
    varCombine.00 <- varCombine[[2]]
  } else if (is.list(varCombine) & length(varCombine) == 1) {
    varCombine.11 <- varCombine.00 <- varCombine[[1]]
  } else {
    varCombine.11 <- varCombine.00 <- varCombine
  }
  
  # Get names for variances to create
  dtrTriples <- dtrNames(design)
  nrNames <- unique(sapply(dtrTriples, function(x) 
    paste0(substr(x, 1, 1), substr(x, 3, 3))))
  dtrs <- do.call(rbind, dtrIndex(design))
  
  if (design == 1) {
    rNames <- unique(sapply(dtrNames(1), function(x) paste0(substr(x, 1, 2))))
    varOrder <- data.frame("genPair" = c("11", "10", "10", "00", "01", "01"),
                           "refPair" = c(rep("11", 3), rep("00", 3)),
                           "genR"    = rep(c("nr", "r", "nr"), 2))
  } else if (design == 2) {
    rNames <- c("11", "00")
    varOrder <- data.frame("genPair" = c("11", "10", "01", "00"),
                           "refPair" = c("11", "11", "00", "00"),
                           "genR"    = rep("nr", 4))
  } else if (design == 3) {
    rNames <- c("11", "00")
    varOrder <- data.frame("genPair" = c("11", "10", "00"),
                           "refPair" = c("11", "11", "00"),
                           "genR"    = rep("nr", 3))
  }
  
  # Create an environment to store variances created in loops
  # (This is to simplify scoping)
  varEnv <- new.env()
  
  # Reconstruct simGrid, adding appropriate values
  sg <- foreach(i = 1:nrow(simGrid), .combine = rbind, .inorder = TRUE, 
                .export = ls(envir = .GlobalEnv)
  ) %dorng% {
    
    if (!is.null(seed)) set.seed(seed)
    
    # cat(paste0("Starting iteration ", i, ".\n"))
    
    # Extract parameters from simGrid
    rho <- simGrid$corr[i]
    respFunction <- get(unlist(simGrid$respFunction[i]))
    if ("respDirection" %in% names(simGrid)){
      respDir <- simGrid$respDirection[i]
    } else {
      respDir <- "high"
    }
    
    r0 <- simGrid$r0[i]
    r1 <- simGrid$r1[i]
    
    # Generate a big data frame to estimate some of these things
    s1 <- generateStage1(2e5, times, spltime, r1, r0, gammas, sigma, corstr,
                         rho = rho, respFunction = respFunction,
                         respDirection = respDir, balanceRand = F,
                         empirical = F)
    
    s <- generateStage2.means(s1, times, spltime, gammas, lambdas, 
                              design = design, corstr = corstr)
    
    # Extract (potentially-modified) response probabilities from s1
    r0  <- s1$r0
    r1  <- s1$r1
    rho <- s1$rho
    
    # Compute bounds on manipulable variances.
    # Note that the lower bound is for the violation of the conditional
    # variation assumption -- it's generally possible to go below it.
    # On the other hand, the upper bound is a hard limit -- going over this
    # value will produce errors
    sigma.r0.LBf <- computeVarBound(0, s, design, sigma, r0, rho, times, corstr,
                                   bound = "lower", boundType = "feasibility")
    sigma.r0.LBa <- computeVarBound(0, s, design, sigma, r0, rho, times, corstr,
                                    bound = "lower", boundType = "assumption")
    sigma.r0.UB <- computeVarBound(0, s, design, sigma, r0, rho, times, corstr, 
                                   bound = "upper")
    assign("sigma.r00", varCombine.00(c(sigma.r0.LBa, sigma.r0.LBf, sigma.r0.UB)),
           envir = varEnv)
    
    sigma.r1.LBf <- computeVarBound(1, s, design, sigma, r0, rho, times, corstr,
                                    bound = "lower", boundType = "feasibility")
    sigma.r1.LBa <- computeVarBound(1, s, design, sigma, r0, rho, times, corstr,
                                    bound = "lower", boundType = "assumption")
    sigma.r1.UB <- computeVarBound(1, s, design, sigma, r1, rho, times, corstr, 
                                   bound = "upper")
    assign("sigma.r11", varCombine.11(c(sigma.r1.LBa, sigma.r1.LBf, sigma.r1.UB)),
           envir = varEnv)
    
    # If the upper bound is below the lower bound, choose sigma.rX to be 
    # below the upper bound to avoid errors (note that this violates the
    # conditional variation assumption)
    if (varEnv$sigma.r00 > sigma.r0.UB) 
      assign("sigma.r00", sigma.r0.UB * .95, envir = varEnv)
    if (varEnv$sigma.r11 > sigma.r1.UB) 
      assign("sigma.r11", sigma.r1.UB * .95, envir = varEnv)
    
    
    for (index in 1:nrow(varOrder)) {
      r <- get(paste0("r", substr(varOrder$genPair[index], 1, 1)))
      R <- s[[paste0("R.", substr(varOrder$genPair[index], 1, 1))]]
      
      if (varOrder$genR[index] == "nr") {
        dtr <- with(varOrder, paste0(substr(genPair[index], 1, 1),
                                     0,
                                     substr(genPair[index], 2, 2)))
        Y2.R  <- s[[paste0("Y2.", dtr)]][R == 1]
        Y2.NR <- s[[paste0("Y2.", dtr)]][R == 0]
        sigma.r <- get(paste0("sigma.r", varOrder$refPair[index]), 
                       envir = varEnv)
        
        assign(paste0("sigma.nr", varOrder$genPair[index]),
               sqrt((sigma^2 - r * sigma.r^2 - 
                       r * (1 - r) * (mean(Y2.R) - mean(Y2.NR))^2) / (1 - r)),
               envir = varEnv)
      } else {
        dtr <- with(varOrder, paste0(substr(genPair[index], 1, 1),
                                     0,
                                     substr(refPair[index], 2, 2)))
        Y2.R  <- s[[paste0("Y2.", dtr)]][R == 1]
        Y2.NR <- s[[paste0("Y2.", dtr)]][R == 0]
        
        sigma.nr <- get(paste0("sigma.nr", varOrder$refPair[index]), 
                        envir = varEnv)
        assign(paste0("sigma.r", varOrder$genPair[index]), 
               sqrt((sigma^2 - (1 - r) * sigma.nr^2 - 
                       r * (1 - r) * (mean(Y2.R) - mean(Y2.NR))^2) / r),
               envir = varEnv)
      }
    }
    
    # Compute residual variances for use in generateSMART()
    for (index in nrNames) {
      sigma.nr <- get(paste0("sigma.nr", index), envir = varEnv)
      r <- get(paste0("r", substr(index, 1, 1)))
      Y1 <- with(s, get(paste0("Y1.", substr(index, 1, 1))))
      R <- with(s, get(paste0("R.", substr(index, 1, 1))))
      Y2 <- with(s, get(paste0("Y2.",
                               substr(index, 1, 1), 0, substr(index, 2, 2))))
      
      assign(paste0("v2.", substr(index, 1, 1), ".NR.", substr(index, 2, 2)),
             sigma.nr^2 -
               (rho / (1 + rho))^2 * var(s$Y0[R == 0] + Y1[R == 0]),
             envir = varEnv)
    }
    for (index in rNames) {
      Y1 <- with(s, get(paste0("Y1.", substr(index, 1, 1))))
      R <- with(s, get(paste0("R.", substr(index, 1, 1))))
      sigma.r <- get(paste0("sigma.r", index), env = varEnv)
      
      assign(paste0("v2.", substr(index, 1, 1), ".R.", substr(index, 2, 2)),
             sigma.r^2 - (rho / (1 + rho))^2 * var(s$Y0[R == 1] + Y1[R == 1]),
             envir = varEnv)
    }
    
    v <- do.call(cbind, lapply(names(varEnv), function(v) get(v, envir = varEnv)))
    colnames(v) <- names(varEnv)
    
    cbind(simGrid[i, ], "sigma.r0.LBf" = sigma.r0.LBf,
          "sigma.r0.LBa" = sigma.r0.LBa,
          "sigma.r0.UB" = sigma.r0.UB,
          "sigma.r1.LBf" = sigma.r1.LBf,
          "sigma.r1.LBa" = sigma.r1.LBa,
          "sigma.r1.UB" = sigma.r1.UB, v)
  }
  
  sg
}

createDTRIndicators <- function(d, design) {
  if (sum(c("A1", "A2R", "A2NR") %in% names(d)) != 3)
    stop("Must provide data from a SMART.")
  if (design == 1) {
    d$dtr1 <- as.numeric(with(d, A1 ==  1 & (A2R ==  1 | A2NR ==  1)))
    d$dtr2 <- as.numeric(with(d, A1 ==  1 & (A2R ==  1 | A2NR == -1)))
    d$dtr3 <- as.numeric(with(d, A1 ==  1 & (A2R == -1 | A2NR ==  1)))
    d$dtr4 <- as.numeric(with(d, A1 ==  1 & (A2R == -1 | A2NR == -1)))
    d$dtr5 <- as.numeric(with(d, A1 == -1 & (A2R ==  1 | A2NR ==  1)))
    d$dtr6 <- as.numeric(with(d, A1 == -1 & (A2R ==  1 | A2NR == -1)))
    d$dtr7 <- as.numeric(with(d, A1 == -1 & (A2R == -1 | A2NR ==  1)))
    d$dtr8 <- as.numeric(with(d, A1 == -1 & (A2R == -1 | A2NR == -1)))
  } else if (design == 2) {
    d$dtr1 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR ==  1))))
    d$dtr2 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR == -1))))
    d$dtr3 <- as.numeric(with(d, (A1 == -1) * (R + (1 - R) * (A2NR ==  1))))
    d$dtr4 <- as.numeric(with(d, (A1 == -1) * (R + (1 - R) * (A2NR == -1))))
  } else if (design == 3) {
    d$dtr1 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR ==  1))))
    d$dtr2 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR == -1))))
    d$dtr3 <- as.numeric(with(d, (A1 == -1)))
  } else stop("design must be one of 1, 2, 3.")
  
  return(d)
}

checkVarGridValidity <- function(varGrid) {
  invalidSims <- subset(varGrid, ((sigma.r0.LBa <= sigma.r0.LBf) |
                                           (sigma.r1.LBa <= sigma.r1.LBf) | 
                                           (sigma.r0.UB <= sigma.r0.LBf) |
                                           (sigma.r1.UB <= sigma.r1.LBf)),
                        select = "simName")
  invalidSims$noViol1 <- invalidSims$noViol0 <- F
  invalidSims$noPositiveVar1 <- invalidSims$noPositiveVar0 <- F
  
  for (i in 1:nrow(invalidSims)) {
    sgrow <- varGrid[varGrid$simName == invalidSims$simName[i], ]
    if (sgrow$sigma.r0.LBa <= sgrow$sigma.r0.LBf) 
      invalidSims$noViol0[i] <- T
    if (sgrow$sigma.r1.LBa <= sgrow$sigma.r1.LBf) 
      invalidSims$noViol1[i] <- T
    if (sgrow$sigma.r0.UB <= sgrow$sigma.r0.LBf)
      invalidSims$noPositiveVar0[i] <- T
    if (sgrow$sigma.r1.UB <= sgrow$sigma.r1.LBf)
      invalidSims$noPositiveVar1[i] <- T
  }
  
  invalidSims
}
