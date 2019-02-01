# functions.R
# Functions for simulations of SMARTs with continuous repeated-measures outcomes
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

### Check assumptions in generative model
checkSMART <- function(smart) {
  # Argument 'smart' must be of class 'generateSMART'
  if (!("generateSMART" %in% class(smart))) stop("Argument 'smart' must be of class 'generateSMART'")
  
  # 1. Check response correlation
  respCor <- estimate.respCor(smart$data, smart$params$design, smart$params$time,
                              smart$params$spltime, smart$params$gammas)
  
  # 2. Check 
}

### Custom combine function for the foreach %dopar% loop in sim()
combine.results <- function(list1, list2) {
  pval       <- c(list1$pval, list2$pval)
  param.hat  <- rbind(list1$param.hat, list2$param.hat)
  sigma2.hat <- Reduce(function(x, y) {x[is.na(x)] <- 0; y[is.na(y)] <- 0; x + y},
                       list(list1$sigma2.hat, list2$sigma2.hat))
  param.var  <- Reduce(function(x, y) {x[is.na(x)] <- 0; y[is.na(y)] <- 0; x + y},
                       list(list1$param.var, list2$param.var))
  rho.hat    <- c(list1$rho.hat, list2$rho.hat)
  coverage   <- list1$coverage + list2$coverage
  iter       <- c(list1$iter, list2$iter)
  valid      <- list1$valid + list2$valid
  condVars   <- Map(function(x, y) {x[is.na(x)] <- 0; y[is.na(y)] <- 0; x + y},
                    list1$condVars, list2$condVars)
  respCor    <- Reduce(function(x, y) {x[is.na(x)] <- 0; y[is.na(y)] <- 0; x + y},
                       list(list1$respCor, list2$respCor))
  condCov    <- Reduce(function(x, y) {x[is.na(x)] <- 0; y[is.na(y)] <- 0; x + y},
                       list(list1$condCov, list2$condCov))
  assumptionViolations <- Map(function(x, y) {x[is.na(x)] <- 0; y[is.na(y)] <- 0; x + y},
                              list1$assumptionViolations, list2$assumptionViolations)
  
  result <- list("pval" = pval, "param.hat" = param.hat, "sigma2.hat" = sigma2.hat, "iter" = iter,
                 "param.var" = param.var, "rho.hat" = rho.hat, "valid" = valid, "coverage" = coverage,
                 "condVars" = condVars, "respCor" = respCor, "condCov" = condCov, 
                 "assumptionViolations" = assumptionViolations)
  
  if ("data" %in% names(list1)) {
    dataList <- c(list1$data, list2$data)
    result[["data"]] <- dataList
  }
  
  return(result)
}

compute.effectSize <- function(gammas, sigma, design, L) {
  
}

conditionalIndex <- function(design) {
  dtrs <- dtrIndex(design)
  grid <- expand.grid("a1" = c(1, -1), "r" = c(0, 1),
                      "a2R" = c(1, 0, -1), "a2NR" = c(1, 0, -1))
  if (design == 1) {
    index <- subset(grid, (r == 1 & (a2R != 0 & a2NR == 0)) |
                      (r == 0 & (a2NR != 0 & a2R == 0)) )
  } else if (design == 2) {
    index <- subset(grid, (r == 1 & (a2R == 0 & a2NR == 0)) |
                      (r == 0 & (a2NR != 0 & a2R == 0)))
  } else if (design == 3) {
    index <- subset(grid, (a1 == -1 & a2R == 0 & a2NR == 0) |
                      (r == 1 & a1 == 1 & a2R == 0 & a2NR == 0) |
                      (r == 0 & a1 == 1 & a2NR != 0 & a2R == 0))
  }
  index <- index[order(index$a1, index$r, decreasing = T), ]
  rownames(index) <- 1:nrow(index)
  index
}

### Compute conditional mean Y at a particular time for a particular cell,
### per the specified model
conditional.model <- function(a1, r, a2R, a2NR, obsTime, spltime, design, rprob, gammas, lambdas) {
  if (obsTime <= spltime) {
    mu <- gammas[1] + obsTime * (gammas[2] + gammas[3] * a1)
  } else {
    if (design == 1) {
      if (length(gammas) != 9) stop("for design 1 gammas must be length 9.")
      mu <- (obsTime - spltime) * (gammas[4] + gammas[5] * a1 + r * (gammas[6] * a2R + gammas[8] * a1 * a2R) / rprob + 
                                     (1 - r) * (gammas[7] * a2NR + gammas[9] * a1 * a2NR) / (1 - rprob) +
                                     (r - rprob) * (lambdas[1] + lambdas[2] * a1))
    } else if (design == 2) {
      if (length(gammas) != 7) stop("for design 2, gammas must be length 7.")
      mu <- (obsTime - spltime) * 
        (gammas[4] + gammas[5] * a1 + 
           (gammas[6] * (1 - r) * a2NR + gammas[7] * a1 * (1 - r) * a2NR)/(1 - rprob)) +
        (obsTime - spltime) * (r - rprob) * (lambdas[1] + lambdas[2] * a1)
    } else if (design == 3) {
      if (length(gammas) != 6) stop("for design 3, gammas must be length 6.")
      mu <- (obsTime - spltime) * (gammas[4] + gammas[5] * a1 + gammas[6] * (a1 == 1) * a2NR * (1 - r) / (1 - rprob)) +
        (obsTime - spltime) * (r - rprob) * (lambdas[1] + lambdas[2] * a1)
    }
    mu <- gammas[1] + spltime * (gammas[2] + gammas[3] * a1) + mu
  }
  mu
}

conditionalVarmat <- function(times, spltime, design, r1, r0,
                              corstr = c("identity", "exchangeable", "ar1"),
                              sigma, sigma.r1, sigma.r0,
                              uneqsd = NULL, uneqsdDTR = NULL,
                              rho, rho.r1, rho.r0,
                              gammas, lambdas) {
  
  nDTR <- switch(design, 8, 4, 3)
  varEnv <- new.env()
  corstr <- match.arg(corstr)
  A <- conditionalIndex(design)
  dtrs <- dtrIndex(design)
  dtrMatrix <- do.call(cbind, dtrs)
  
  if (length(sigma) == 1) {
    sigma <- matrix(rep(sigma, nDTR * length(times)), ncol = length(times))
  } else if (length(sigma) == length(times)) {
    sigma <- matrix(rep(sigma, nDTR), nrow = nDTR, byrow = T)
  } else if (length(sigma) != length(times) & length(sigma) != nDTR * length(times)) {
    stop("sigma must either be length-1, the same length as times, or must be an nDTR-by-length(times) matrix.")
  }
  
  ## Handle sigma.r1, sigma.r0 inputs in case sum(times > spltime) > 1 
  ## (i.e., more than one measurement after re-randomization)
  if (length(sigma.r1) == 1 & sum(times > spltime) > 1) {
    sigma.r1 <- rep(sigma.r1, sum(times > spltime))
  }
  if (length(sigma.r0) == 1 & sum(times > spltime) > 1) {
    sigma.r0 <- rep(sigma.r0, sum(times > spltime))
  }
  if (length(sigma.r1) != length(times)) {
    sigma.r1 <- c(sigma[min(which(dtrs$a1 == 1)), times <= spltime], sigma.r1)
  }
  if (length(sigma.r0) != length(times)) {
    sigma.r0 <- c(sigma[min(which(dtrs$a1 == -1)), times <= spltime], sigma.r0)
  }
  
  # Construct a matrix of (X,Y) coordinates which need to be adjusted in the conditional cor matrices
  # i.e., identify which cells of the conditional cor matrix correspond to correlations between the 
  # second stage and times before the second randomization.
  # The rows of m are the coordinates (row, column)
  m <- expand.grid(which(times > spltime), which(times <= spltime))
  m <- rbind(as.matrix(m), as.matrix(m[, 2:1]))
  m <- rbind(m, expand.grid(which(times > spltime), which(times > spltime)))
  m <- subset(m, Var1 != Var2)
  
  # Construct "base" correlation matrices for responders 
  # (these are the marginal correlation matrices and will be modified below)
  cormat.r1 <- cormat.r0 <- cormat(rho, length(times), corstr)
  
  # Create responders' correlation matrices by looping over m and computing appropriate entries
  for (i in 1:nrow(m)) {
    dif <- abs(m[i, 1] - m[i, 2])
    cormat.r1[m[i, 1], m[i, 2]] <- ifelse(corstr == "ar1", rho.r1^dif, rho.r1)
    cormat.r0[m[i, 1], m[i, 2]] <- ifelse(corstr == "ar1", rho.r0^dif, rho.r0)
  }
  
  if (design == 1) {
    dtrs <- do.call("rbind", dtrIndex(design))
    if (!is.null(uneqsdDTR)) {
      # Determine which DTRs are the unequal variance ones
      uneqsdDTRindex <- sapply(1:2, function(i) 
        which((dtrs[1, ] == uneqsdDTR[[i]][1]) + 
                (dtrs[2, ] == uneqsdDTR[[i]][2]) +
                (dtrs[3, ] == uneqsdDTR[[i]][3]) == 3)
      )
      # Compute the standard deviation for times after spltime in the DTR with unequal variance
      uneqsd <- lapply(uneqsdDTR, function(dtr) {
        sapply(times, function(obsTime) {
          generate.sd(sigma[which(times == obsTime)], a1 = dtr[1], a2R = dtr[2], a2NR = dtr[3], obsTime = obsTime, spltime = spltime, 
                      design, rprob = get(paste0("r", (dtr[1] + 1) / 2)), gammas, lambdas)
        })
      })
    }
    assign("sigma.r11", sigma.r1, envir = varEnv)
    assign("sigma.r01", sigma.r0, envir = varEnv)
    
    #Loop over indices for first-stage treatment and compute conditional standard deviations
    invisible(sapply(c(0, 1), function(a1ind) {
      # Use DTR (a1ind, 1, 1) to 
      assign(paste0("sigma.nr", a1ind, "1"),
             sapply(1:length(times), function(i) {
               generate.cond.sd(a1 = 2 * a1ind - 1, r = 0, a2R = 1, a2NR = 1, obsTime = times[i], spltime = spltime,
                                design = 1, rprob = get(paste0("r", a1ind)),
                                sigma = ifelse(!is.null(uneqsdDTR), 
                                               ifelse(sum(sapply(uneqsdDTR, function(dtr) sum(dtr == c(2 * a1ind - 1, 1, 1)) == 3)) > 0,
                                                      uneqsdDTR.sd[[which(sapply(uneqsdDTR, function(dtr) 
                                                        sum(dtr == c(2 * a1ind -1, 1, 1)) == 3))]][i],
                                                      sigma[i]),
                                               sigma[i]),
                                sigma.cond = get(paste0("sigma.r", a1ind, "1"), envir = varEnv)[i], gammas, lambdas)
             }), envir = varEnv)
      assign(paste0("sigma.r", a1ind, "0"),
             sapply(1:length(times), function(i) {
               generate.cond.sd(a1 = 2 * a1ind - 1, r = 1, a2R = -1, a2NR = 1, obsTime = times[i], spltime = spltime,
                                design = 1, rprob = get(paste0("r", a1ind)),
                                sigma = ifelse(!is.null(uneqsdDTR),
                                               ifelse(sum(sapply(uneqsdDTR, function(dtr) sum(dtr == c(2 * a1ind - 1, -1, 1)) == 3)) > 0,
                                                      uneqsdDTR.sd[[which(sapply(uneqsdDTR, function(dtr) 
                                                        sum(dtr == c(2 * a1ind -1, -1, 1)) == 3))]][i],
                                                      sigma[i]),
                                               sigma[i]),
                                sigma.cond = get(paste0("sigma.nr", a1ind, "1"), envir = varEnv)[i], gammas, lambdas)
             }), envir = varEnv)
      assign(paste0("sigma.nr", a1ind, "0"),
             sapply(1:length(times), function(i) {
               generate.cond.sd(a1 = 2 * a1ind - 1, r = 0, a2R = -1, a2NR = -1, obsTime = times[i], spltime = spltime,
                                design = 1, rprob = get(paste0("r", a1ind)),
                                sigma = ifelse(!is.null(uneqsdDTR),
                                               ifelse(sum(sapply(uneqsdDTR, function(dtr) sum(dtr == c(2 * a1ind - 1, -1, -1)) == 3)) > 0,
                                                      uneqsdDTR.sd[[which(sapply(uneqsdDTR, function(dtr) 
                                                        sum(dtr == c(2 * a1ind -1, -1, -1)) == 3))]][i],
                                                      sigma[i]),
                                               sigma[i]),
                                sigma.cond = get(paste0("sigma.r", a1ind, "0"), envir = varEnv)[i], gammas, lambdas)
             }), envir = varEnv)
    }))
    
    assign("cormat.r11", cormat.r1, env = varEnv)
    assign("cormat.r01", cormat.r0, env = varEnv)
    
    for (a1ind in c(0, 1)) {
      x1 <- x2 <- x3 <- cormat(rho, length(times), corstr)
      for(i in 1:nrow(m)) {
        x1[m[i, 1], m[i, 2]] <- 
          generate.cond.cor(a1 = 2 * a1ind - 1, r = 0, a2R = 1, a2NR = 1,
                            t1 = times[m[i, 1]], t2 = times[m[i, 2]], spltime = spltime,
                            design = 1, rprob = get(paste0("r", a1ind)),
                            sigma.t1 = sigma[m[i, 1]], sigma.t2 = sigma[m[i, 2]], 
                            sigma.t1.ref = get(paste0("sigma.r", a1ind, 1), envir = varEnv)[m[i, 1]],
                            sigma.t2.ref = get(paste0("sigma.r", a1ind, 1), envir = varEnv)[m[i, 2]],
                            rho = ifelse(corstr == "ar1", rho^abs(m[i, 1] - m[i, 2]), rho),
                            rho.ref = get(paste0("cormat.r", a1ind, 1), envir = varEnv)[m[i, 1], m[i, 2]],
                            gammas = gammas, lambdas = lambdas)
      }
      assign(paste0("cormat.nr", a1ind, "1"), x1, envir = varEnv)
      
      for(i in 1:nrow(m)) {
        x2[m[i, 1], m[i, 2]] <- 
          generate.cond.cor(a1 = 2 * a1ind - 1, r = 1, a2R = -1, a2NR = 1,
                            t1 = times[m[i, 1]], t2 = times[m[i, 2]], spltime = spltime,
                            design = 1, rprob = get(paste0("r", a1ind)),
                            sigma.t1 = sigma[m[i, 1]], sigma.t2 = sigma[m[i, 2]], 
                            sigma.t1.ref = get(paste0("sigma.nr", a1ind, 1), envir = varEnv)[m[i, 1]],
                            sigma.t2.ref = get(paste0("sigma.nr", a1ind, 1), envir = varEnv)[m[i, 2]],
                            rho = ifelse(corstr == "ar1", rho^abs(m[i, 1] - m[i, 2]), rho),
                            rho.ref = get(paste0("cormat.nr", a1ind, 1), envir = varEnv)[m[i, 1], m[i, 2]],
                            gammas = gammas, lambdas = lambdas)
      }
      assign(paste0("cormat.r", a1ind, "0"), x2, envir = varEnv)
      
      for(i in 1:nrow(m)) {
        x3[m[i, 1], m[i, 2]] <- 
          generate.cond.cor(a1 = 2 * a1ind - 1, r = 0, a2R = -1, a2NR = -1,
                            t1 = times[m[i, 1]], t2 = times[m[i, 2]], spltime = spltime,
                            design = 1, rprob = get(paste0("r", a1ind)),
                            sigma.t1 = sigma[m[i, 1]], sigma.t2 = sigma[m[i, 2]], 
                            sigma.t1.ref = get(paste0("sigma.r", a1ind, 0), envir = varEnv)[m[i, 1]],
                            sigma.t2.ref = get(paste0("sigma.r", a1ind, 0), envir = varEnv)[m[i, 2]],
                            rho = ifelse(corstr == "ar1", rho^abs(m[i, 1] - m[i, 2]), rho),
                            rho.ref = get(paste0("cormat.r", a1ind, 0), envir = varEnv)[m[i, 1], m[i, 2]],
                            gammas = gammas, lambdas = lambdas)
      }
      assign(paste0("cormat.nr", a1ind, "0"), x3, envir = varEnv)
    }
  } else if (design == 2) {
    assign("sigma.r10", sigma.r1, envir = varEnv)
    assign("sigma.r00", sigma.r0, envir = varEnv)
    
    assign("cormat.r10", cormat.r1, env = varEnv)
    assign("cormat.r00", cormat.r0, env = varEnv)
    
    if ((!is.null(uneqsdDTR) & is.null(uneqsd)) | (is.null(uneqsdDTR) & !is.null(uneqsd)))
      stop("For design 2, you must provide either both uneqsdDTR and uneqsd or neither.")
    invisible(sapply(1:nDTR, function(i) {
      if (!is.null(uneqsdDTR)) {
        # Check if the current embedded DTR (x) is in uneqsdDTR; if it is, store the index
        # (If it's not, uneqsdDTRindex has length 0)
        uneqsdDTRindex <- which(sapply(uneqsdDTR, function(dtr) {
          sum(c(dtrs$a1[i], dtrs$a2R[i], dtrs$a2NR[i]) == dtr) == 3
        }))
        if (length(uneqsdDTRindex) == 0) {
          sigmaStar <- sigma
        } else {
          # If the current DTR is in uneqsdDTR, extract the appropriate SD and vectorize it as usual
          # i.e., create a length(times)-vector where each entry is the marginal SD at that time
          sigmaStar <- uneqsd[uneqsdDTRindex]
          if (length(sigmaStar) == 1 & sum(times > spltime) > 1) {
            sigmaStar <- rep(sigmaStar, sum(times > spltime))
          }
          if (length(sigmaStar) != length(times)) {
            sigmaStar <- c(sigma[times <= spltime], sigmaStar)
          }
        }
      } else {
        # If no uneqsdDTR, use the common SD
        sigmaStar <- sigma[i, ]
      }
      
      # Generate conditional SDs for non-responders and assign them to appropriate variables.
      assign(paste0("sigma.nr", (dtrs$a1[i] + 1) / 2, (dtrs$a2NR[i] + 1) / 2),
             unname(sapply(1:length(times), function(j) {
               # Generate conditional standard deviation for non-responders
               generate.cond.sd(a1 = dtrs$a1[i], r = 0, a2R = dtrs$a2R[i], a2NR = dtrs$a2NR[i],
                                obsTime = times[j], spltime, design = 2, 
                                rprob = get(paste0("r", (dtrs$a1[i] + 1) / 2)),
                                sigma = sigmaStar[j],
                                sigma.cond = get(paste0("sigma.r", (dtrs$a1[i] + 1) / 2, 0), envir = varEnv)[j],
                                gammas = gammas, lambdas = lambdas)
             })), envir = varEnv)
      
      x <- cormat(rho, length(times), corstr)
      
      for (k in 1:nrow(m)) {
        x[m[k, 1], m[k, 2]] <- 
          generate.cond.cor(a1 = dtrs$a1[i], r = 0, a2R = dtrs$a2R[i], a2NR = dtrs$a2NR[i],
                            t1 = times[m[k, 1]], t2 = times[m[k, 2]], spltime = spltime,
                            design = 2, rprob = get(paste0("r", (dtrs$a1[i] + 1) / 2)),
                            sigma.t1 = sigmaStar[m[k, 1]], sigma.t2 = sigmaStar[m[k, 2]],
                            sigma.t1.ref = get(paste0("sigma.r", (dtrs$a1[i] + 1) / 2, 0), envir = varEnv)[m[k, 1]],
                            sigma.t2.ref = get(paste0("sigma.r", (dtrs$a1[i] + 1) / 2, 0), envir = varEnv)[m[k, 2]],
                            rho = ifelse(corstr == "ar1", rho^abs(m[k, 1] - m[k, 2]), rho),
                            rho.ref = get(paste0("cormat.r", (dtrs$a1[i] + 1) / 2, 0), envir = varEnv)[m[k, 1], m[k, 2]],
                            gammas = gammas, lambdas = lambdas)
      }
      assign(paste0("cormat.nr", (dtrs$a1[i] + 1) / 2, (dtrs$a2NR[i] + 1) / 2), x, envir = varEnv)
    }))
    
    # NRgrid <- expand.grid("a1ind" = c(0, 1), "a2NRind" = c(0, 1))
    # 
    # for (j in 1:nrow(NRgrid)) {
    #   a1ind <- NRgrid$a1ind[j]
    #   a2NRind <- NRgrid$a2NRind[j]
    #   
    #   x <- cormat(rho, length(times), corstr)
    #   
    #   for (i in 1:nrow(m)) {
    #     x[m[i, 1], m[i, 2]] <- 
    #       generate.cond.cor(a1 = 2 * a1ind - 1, r = 0, a2R = 0, a2NR = 2 * a2NRind - 1,
    #                         t1 = times[m[i, 1]], t2 = times[m[i, 2]], spltime = spltime,
    #                         design = 2, rprob = get(paste0("r", a1ind)),
    #                         sigma.t1 = sigmaStar[m[i, 1]], sigma.t2 = sigmaStar[m[i, 2]],
    #                         sigma.t1.ref = get(paste0("sigma.r", a1ind, 0), envir = varEnv)[m[i, 1]],
    #                         sigma.t2.ref = get(paste0("sigma.r", a1ind, 0), envir = varEnv)[m[i, 2]],
    #                         rho = ifelse(corstr == "ar1", rho^abs(m[i, 1] - m[i, 2]), rho),
    #                         rho.ref = get(paste0("cormat.r", a1ind, 0), envir = varEnv)[m[i, 1], m[i, 2]],
    #                         gammas = gammas, lambdas = lambdas)
    #   }
    #   assign(paste0("cormat.nr", a1ind, a2NRind), x, envir = varEnv)
    # }
  } else if (design == 3) {
    assign("sigma.r10", sigma.r1, envir = varEnv)
    assign("sigma.r00", sigma.r0, envir = varEnv)
    
    assign("cormat.r10", cormat.r1, env = varEnv)
    assign("cormat.r00", cormat.r0, env = varEnv)
    
    invisible(sapply(1:nDTR, function(i) {
      assign(paste0("sigma.nr", (dtrs$a1[i] + 1) / 2, abs(dtrs$a2NR[i]) * (dtrs$a2NR[i] + 1) / 2),
             unname(sapply(1:length(times), function(j) {
               generate.cond.sd(a1 = dtrs$a1[i], r = 0, a2R = dtrs$a2R[i], a2NR = dtrs$a2NR[i],
                                obsTime = times[j], spltime, design = 3, 
                                rprob = get(paste0("r", (dtrs$a1[j] + 1) / 2)),
                                sigma = sigma[j],
                                sigma.cond = get(paste0("sigma.r", (dtrs$a1[i] + 1) / 2, 0), envir = varEnv)[j],
                                gammas = gammas, lambdas = lambdas)
             })), envir = varEnv)
    }))
    
    NRgrid <- expand.grid("a1ind" = c(0, 1), "a2NRind" = c(0, 1))
    NRgrid <- subset(NRgrid, !(a1ind == 0 & a2NRind == 1))
    
    for (j in 1:nrow(NRgrid)) {
      a1ind <- NRgrid$a1ind[j]
      a2NRind <- NRgrid$a2NRind[j]
      
      x <- cormat(rho, length(times), corstr)
      
      for (i in 1:nrow(m)) {
        x[m[i, 1], m[i, 2]] <- 
          generate.cond.cor(a1 = 2 * a1ind - 1, r = 0, a2R = 0, a2NR = 2 * a2NRind - 1,
                            t1 = times[m[i, 1]], t2 = times[m[i, 2]], spltime = spltime,
                            design = 3, rprob = get(paste0("r", a1ind)),
                            sigma.t1 = sigma[m[i, 1]], sigma.t2 = sigma[m[i, 2]],
                            sigma.t1.ref = get(paste0("sigma.r", a1ind, 0), envir = varEnv)[m[i, 1]],
                            sigma.t2.ref = get(paste0("sigma.r", a1ind, 0), envir = varEnv)[m[i, 2]],
                            rho = ifelse(corstr == "ar1", rho^abs(m[i, 1] - m[i, 2]), rho),
                            rho.ref = get(paste0("cormat.r", a1ind, 0), envir = varEnv)[m[i, 1], m[i, 2]],
                            gammas = gammas, lambdas = lambdas)
      }
      assign(paste0("cormat.nr", a1ind, a2NRind), x, envir = varEnv)
    }
    
  } else {
    warning("Design must be 1-3")
  }
  
  condVarmats <- lapply(1:nrow(A), function(dtr) {
    a1ind   <- (A$a1[dtr] + 1) / 2
    R       <- A$r[dtr]
    a2Rind  <- ifelse(A$a2R[dtr] == 0, 0, (A$a2R[dtr] + 1) / 2)
    a2NRind <- abs(A$a2NR[dtr]) * (A$a2NR[dtr] + 1) / 2
    xsd  <- diag(get(paste0("sigma.", ifelse(R == 1, "r", "nr"),
                            a1ind, ifelse(R == 1, a2Rind, a2NRind)),
                     envir = varEnv))
    xcor <- get(paste0("cormat.", ifelse(R == 1, "r", "nr"),
                       a1ind, ifelse(R == 1, a2Rind, a2NRind)),
                envir = varEnv)
    xsd %*% xcor %*% xsd
  })
  names(condVarmats) <- apply(A, 1, function(x) paste(x, collapse = " "))
  condVarmats
}

### Construct a t-by-t correlation matrix
cormat <- function(rho, t, corstr = c("identity", "exchangeable", "ar1", "unstructured")) {
  corstr <- match.arg(corstr)
  if (corstr == "identity") {
    diag(rep(1, t))
  } else if (corstr == "exchangeable") {
    m <- matrix(rep(rho, t^2), nrow = t)
    diag(m) <- rep(1,t)
    m
  } else if (corstr == "ar1") {
    m <- diag(t)
    rho^(abs(row(m) - col(m)))
  } else if (corstr == "unstructured") {
    if (length(rho) != t) stop("for unstructured corstr, must have rho of length t")
    m <- diag(3)
    tpairs <- combn(1:t, 2)
    for (j in 1:ncol(tpairs)) {
      m[tpairs[1, j], tpairs[2, j]] <- m[tpairs[2, j], tpairs[1, j]] <- rho[j]
    }
    m
  }
}

dtrIndex <- function(design) {
  if (design == 1) {
    a1   <- rep(c(1, -1), each = 4)
    a2R  <- rep(rep(c(1, -1), each = 2), 2)
    a2NR <- rep(c(1, -1), 4)
  } else if (design == 2) {
    a1   <- rep(c(1, -1), each = 2)
    a2R  <- rep(0, 4)
    a2NR <- rep(c(1, -1), 2)
  } else if (design == 3) {
    a1   <- c(1, 1, -1)
    a2R  <- c(0, 0, 0)
    a2NR <- c(1, -1, 0)
  }
  return(list("a1" = a1, "a2R" = a2R, "a2NR" = a2NR))
}

### Compute value of the estimating equations
esteqn.compute <- function(d, V, times, spltime, design, gammas) {
  # d is an unreplicated data set in LONG format
  # V is a working covariance matrix
  
  n <- length(unique(d$id))
  
  nDTR <- switch(design, 8, 4, 3)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR)
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  
  mvec <- meanvec(times, spltime, design, gammas)
  deriv <- mod.derivs(times, spltime, design)
  
  resids <- matrix(rep(d$Y, nDTR), ncol = nDTR) -
    matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE)
  resids <- cbind("id" = d$id, resids * d[, grep("dtr", names(d))] *
                    matrix(rep(d$weight, nDTR), ncol = nDTR))
  
  DVinv <- lapply(1:nDTR, function(i) crossprod(deriv[[i]], solve(V[[i]])))
  
  Reduce("+", lapply(split.data.frame(resids, resids$id), function(x) {
    Reduce("+", lapply(1:nDTR, function(dtr) {
      DVinv[[dtr]] %*% matrix(x[[paste0("dtr", dtr)]], ncol = 1)
    }))
  })
  ) / n
}

esteqn.jacobian <- function(d, V, times, spltime, design) {
  nDTR <- switch(design, 8, 4, 3)
  deriv <- mod.derivs(times, spltime, design)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR) {
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  }
  
  Reduce("+", lapply(split.data.frame(d, d$id), function(x) {
    Reduce("+", lapply(1:nDTR, function(dtr) {
      unique(x[[paste0("dtr", dtr)]]) * unique(x$weight) *
        crossprod(deriv[[dtr]],solve(V[[dtr]])) %*% deriv[[dtr]]
    }))
  })) / -length(unique(d$id))
}

estimate.condCov <- function(d, times, spltime, design, pool.dtr = F) {
  nDTR <- switch(design, 8, 4, 3)
  
  stage1times <- times[times <= spltime]
  stage2times <- times[times > spltime]
  m <- expand.grid("stage1" = stage1times, "stage2" = stage2times)
  x <- lapply(1:nDTR, function(j) {
    condCov <- matrix(rep(NA, length = nrow(m) * 2), nrow = 2,
                      dimnames = list(apply(m, 1, function(x) paste(x, collapse = "")), c("NR", "R")))
    for (i in 1:nrow(m)) {
      condCov[i, 1] <- with(d[d$R == 0 & d[[paste0("dtr", j)]] == 1, ],
                            cov(Y[time == m$stage1[i]], Y[time == m$stage2[i]]))
      condCov[i, 2] <- with(d[d$R == 1 & d[[paste0("dtr", j)]] == 1, ],
                            cov(Y[time == m$stage1[i]], Y[time == m$stage2[i]]))
    }
    condCov
  })
  if (pool.dtr) 
    Reduce("+", x) / nDTR
  else 
    x
}

estimate.params <- function(d, V, times, spltime, design, start, maxiter.solver, tol) {
  nDTR <- switch(design, 8, 4, 3)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR)
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  
  params.hat <- start
  j <- 0
  epsilon <- 10
  J <- esteqn.jacobian(d, V, times, spltime, design)
  
  for (j in 1:maxiter.solver) {
    params.hat.new <-
      params.hat - solve(J) %*% esteqn.compute(d, V, times, spltime, design, gammas = params.hat)
    epsilon <- norm(params.hat.new - params.hat, type = "F")
    params.hat <- params.hat.new
    if (epsilon <= tol) break
  }
  params.hat
}

estimate.paramvar <- function(d, V, times, spltime, design, gammas) {
  nDTR <- switch(design, 8, 4, 3)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR)
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  
  n <- length(unique(d$id))
  J <- esteqn.jacobian(d, V, times, spltime, design)
  solve(-J) %*%  meat.compute(d, V, times, spltime, design, gammas) %*% solve(-J) / n
}

estimate.respCor <- function(d, design, times, spltime, gammas, causal = F) {
  # Input is WIDE DATA
  if (!("SMART" %in% class(d))) stop("'d' must be of class 'SMART'")
  
  # Get all distinct pairs of times in stage 1
  corPairs <- combn(times[times <= spltime], 2)
  # Add "squared" times (i.e., non-distinct pairs like 0,0)
  corPairs <- cbind(corPairs,
                    do.call(cbind, lapply(times[times <= spltime],
                                          function(x) rep(x, 2))))
  # sort corPairs by time
  corPairs <- corPairs[, order(corPairs[1, ], corPairs[2, ])]
  
  x <- lapply(unique(d$A1), function(a1) {
    # subset data to only look at those consistent with one DTR at a time
    with(d[d$A1 == a1, ],
         # loop over columns of corPairs
         sapply(1:ncol(corPairs), function(j) {
           cor(R, (get(paste0("Y", corPairs[1, j])) - 
                     marginal.model(a1, 0, 0, corPairs[1, j], spltime, design, gammas)) * 
                 (get(paste0("Y", corPairs[2, j])) - 
                    marginal.model(a1, 0, 0, corPairs[2, j], spltime, design, gammas)))
         })
         )
  })
  x <- do.call(cbind, x)
  colnames(x) <- unique(d$A1)
  rownames(x) <- apply(corPairs, 2, function(x) paste0(x, collapse = ""))
  x
}

estimate.rho <- function(d, times, spltime, design, sigma, gammas,
                         corstr = c("exchangeable", "ar1", "unstructured"),
                         pool.dtr = T) {
  ## FIXME: Allow for different rhos across DTRs
  
  corstr <- match.arg(corstr)
  if (any(is(d$Y) == "NULL")) stop("d has to be in long format.")
  
  n <- length(unique(d$id))
  nDTR <- switch(design, 8, 4, 3)
  mvec <- meanvec(times, spltime, design, gammas)
  Dmat <- mod.derivs(times, spltime, design)
  
  sigma <- reshapeSigma(sigma, times, design)
  
  # if (sum(rownames(sigma) == times) == length(times) & ncol(sigma) != nDTR) {
  #   sigma <- matrix(rep(sigma, nDTR), ncol = nDTR)
  # } else if (sum(grepl("dtr", rownames(sigma))) != 0) {
  #   sigma <- t(matrix(rep(sigma, length(times)), ncol = length(times)))
  # } else if (nrow(sigma) == ncol(sigma) & ncol(sigma) == 1) {
  #   sigma <- matrix(rep(sigma, nDTR * length(times)), ncol = nDTR)
  # } else if (nrow(sigma) == length(times) & ncol(sigma) == 1) {
  #   sigma <- matrix(rep(sigma, nDTR), ncol = nDTR)
  # }
  
  colnames(sigma) <- paste0("dtr", 1:nDTR)
  rownames(sigma) <- paste0("time", times)
  
  # Compute residuals (Y_{it} - mu_{t}(a1, a2))
  resids <- matrix(rep(d$Y, nDTR), ncol = nDTR) -
    matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE)
  
  # Create data.frame from resids, indexed by time and ID
  # Multiply residuals by DTR indicators (dtr1[i] = 1 iff i is consistent with DTR 1)
  resids <- cbind("id" = d$id, "time" = d$time,
                  resids * d[, grep("dtr", names(d))])
  
  weights <- aggregate(weight ~ id, d, unique)
  # Create matrix of weights per person-time per DTR
  weightmat <- cbind("id" = d$id, "time" = d$time, d$weight * d[, grep("dtr", names(d))])
  
  # Sum weights over individuals 
  # (but only use one weight per person -- weightmat has duplicated rows: 1 per time)
  sumweights.time <- Reduce(function(...) merge(..., by = "time"), 
                            lapply(1:sum(grepl("dtr", names(d))), function(dtr) {
                              x <- list(weightmat[[paste0("dtr", dtr)]])
                              sumwts <- aggregate(x = setNames(x, paste0("dtr", dtr)),
                                                  by = list("time" = resids$time), sum)
                            }))[, -1]
  rownames(sumweights.time) <- times
  sumweights <- apply(sumweights.time, 2, unique)
  
  if (corstr == "exchangeable") {
    # For every id and every dtr, compute (\sum_{s<t} r_{s} * r_{t} / sigma_{s} sigma_{t}
    r <- do.call(rbind,
                 lapply(split.data.frame(resids, resids$id),
                        function(subdat) {
                          sapply(1:sum(grepl("dtr", names(subdat))), function(dtr) {
                            sum(sapply(2:length(subdat$time), function(ti) {
                              sum((subdat[ti, paste0("dtr", dtr)] / sigma[ti, dtr]) *
                                    (subdat[1:(ti - 1), paste0("dtr", dtr)] / sigma[1:(ti - 1), dtr]))
                            }))
                          }) * weights$weight[weights$id == unique(subdat$id)]
                        }))
    colnames(r) <- grep("dtr", names(resids), value = T)
    denominator <- sumweights * (length(times) * (length(times) - 1) / 2) - length(gammas)
  } else if (corstr == "ar1") {
    r <- do.call(rbind,
                 lapply(split.data.frame(resids, resids$id),
                        function(subdat) {
                          sapply(1:sum(grepl("dtr", names(subdat))), function(dtr) {
                            sum(sapply(2:length(subdat$time), function(time) {
                              sum((subdat[time, paste0("dtr", dtr)] / sigma[time, dtr]) *
                                    (subdat[time - 1, paste0("dtr", dtr)] / sigma[time - 1, dtr]))
                            }))
                          }) * weights$weight[weights$id == unique(subdat$id)]
                        }))
    colnames(r) <- grep("dtr", names(resids), value = T)
    denominator <- sumweights * (length(times) - 1) - length(gammas)
  } else if (corstr == "unstructured") {
    m <- combn(times, 2)
    r <- lapply(split.data.frame(resids, resids$id),
                function(subdat) {
                  x <- do.call(rbind, lapply(1:ncol(m), function(tpair) {
                    # data.frame("tpair" = paste(m[, tpair], collapse = ""), 
                               sapply(1:sum(grepl("dtr", names(subdat))), function(dtr) {
                                 prod(subdat[subdat$time %in% m[, tpair], paste0("dtr", dtr)]) /
                                   prod(sigma[rownames(sigma) %in% m[, tpair], paste0("dtr", dtr)])
                               }) * weights$weight[weights$id == unique(subdat$id)]
                  }))
                  rownames(x) <- sapply(1:ncol(m), function(j) paste(m[, j], collapse = ""))
                  x
                })
    r <- lapply(1:ncol(m), function(tpair) {
      do.call(rbind,
              lapply(r, function(x) x[paste(m[, tpair], collapse = ""), ]))
    })
    names(r) <- sapply(1:ncol(m), function(j) paste(m[, j], collapse = ""))
  }
  
  if (is.list(r)) {
    numerator <- lapply(r, colSums)
    denominator <- sumweights - length(gammas)
    rhoHat <- do.call(rbind, lapply(numerator, function(x) x / denominator))
    if (pool.dtr) 
      rhoHat <- apply(rhoHat, 1, mean)
  } else {
    numerator <- colSums(r)
    rhoHat <- numerator / denominator
    if (pool.dtr)
      rhoHat <- mean(rhoHat)
  }
  
  return(rhoHat)
}


estimate.sigma2 <- function(d, times, spltime, design, gammas, pool.time = F, pool.dtr = F) {
  n <- length(unique(d$id))
  nDTR <- switch(design, 8, 4, 3)
  mvec <- meanvec(times, spltime, design, gammas)
  Dmat <- mod.derivs(times, spltime, design)
  
  if (any(is(d$Y) == "NULL")) stop("d has to be in long format")
  
  resids <- (matrix(rep(d$Y, nDTR), ncol = nDTR) -
               matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE))^2
  resids <- cbind("id" = d$id, "time" = d$time,
                  resids * d[, grep("dtr", names(d))] *
                    matrix(rep(d$weight, nDTR), ncol = nDTR))
  
  weightmat <- cbind("id" = d$id, "time" = d$time, d$weight * d[, grep("dtr", names(d))])
  
  if (pool.time & pool.dtr) {
    numerator   <- sum(subset(resids, select = grep("dtr", names(resids))))
    denominator <- (sum(subset(weightmat, select = grep("dtr", names(weightmat)))) - length(gammas))
  } else if (pool.time & !pool.dtr) {
    numerator   <- apply(subset(resids, select = grep("dtr", names(resids))), 2, sum)
    denominator <- apply(subset(weightmat, select = grep("dtr", names(weightmat))), 2, sum) - length(gammas)
  } else if (!pool.time & pool.dtr) {
    numerator <- do.call("c", lapply(times, function(x) {
      sum(subset(resids, time == x, select = grep("dtr", names(resids))))
    }))
    denominator <- do.call("c", lapply(times, function(x) {
      sum(subset(weightmat, time == x, select = grep("dtr", names(weightmat))))
    }))
    names(numerator) <- names(denominator) <- paste0("time", times)
  } else {
    numerator <- Reduce(function(...) merge(..., by = "time"),
                        lapply(1:nDTR, function(dtr) {
                          x <- list(resids[[paste0("dtr", dtr)]])
                          sumsqrs <- aggregate(x = setNames(x, paste0("dtr", dtr)),
                                               by = list("time" = resids$time), sum)
                          x <- list(weightmat[[paste0("dtr", dtr)]])
                          sumwts <- aggregate(x = setNames(x, paste0("dtr", dtr)),
                                              by = list("time" = resids$time), sum)
                          sumsqrs[, 2] <- sumsqrs[, 2] / sumwts[, 2]
                          sumsqrs
                        }))
    rownames(numerator) <- numerator$time 
    numerator <- as.matrix(numerator[, -1])
    denominator <- 1
  }
  as.matrix(numerator / denominator)
}

# Clean up output of combine.results, to be used in the foreach loop after
# combining everything
finalize.results <- function(x) {
  param.hat            <- apply(as.matrix(x$param.hat), 2, mean, na.rm = T)
  sigma2.hat           <- x$sigma2.hat / x$valid
  param.var            <- x$param.var / x$valid
  param.var.est        <- apply(as.matrix(x$param.hat), 2, var, na.rm = T)
  rho.hat              <- mean(x$rho.hat, na.rm = T)
  coverage             <- as.numeric(x$coverage / x$valid)
  condVars             <- lapply(x$condVars, function(v) v / x$valid)
  respCor              <- x$respCor / x$valid
  condCov              <- x$condCov / x$valid
  assumptionViolations <- lapply(x$assumptionViolations, function(v) v / x$valid)
  
  cat("Finishing...\n")
  
  result <- list("n" = x$n, "alpha" = x$alpha, "power.target" = x$power.target, "delta" = x$delta,
                 "niter" = x$niter, "corstr" = x$corstr, "pval" = x$pval, "param.hat" = param.hat,
                 "sigma2.hat" = sigma2.hat, "iter" = x$iter, "param.var" = param.var,
                 "rho.hat" = rho.hat, "valid" = x$valid, "coverage" = coverage, "condVars" = condVars,
                 "condVars" = condVars, "respCor" = respCor, "condCov" = condCov,
                 "assumptionViolations" = assumptionViolations)
  
  if ("data" %in% names(x))
    result[["data"]] <- x$data
  
  return(result)
}

generate.sd <- function(sigma, a1, a2R, a2NR, obsTime, spltime, design, rprob, gammas, lambdas) {
  dtrs <- dtrIndex(design)
  mu <- sapply(1:length(dtrs$a1), function(i) {
    c(conditional.model(dtrs$a1[i], 1, dtrs$a2R[i], 0, obsTime, spltime, design, rprob, gammas, lambdas),
      conditional.model(dtrs$a1[i], 0, 0, dtrs$a2NR[i], obsTime, spltime, design, rprob, gammas, lambdas))
  })
  stage1 <- (a1 + 1) / 2
  sqrt(sigma^2 + stage1 * a2R * a2NR *
         rprob * (1 - rprob) * ((mu[1, 1] - mu[2, 1])^2 + (mu[1, 3] - mu[2, 2])^2 
                                - (mu[1, 1] - mu[2, 2])^2 - (mu[1, 3] - mu[2, 1])^2) +
         (1 - stage1) * a2R * a2NR *
         rprob * (1 - rprob) * ((mu[1, 5] - mu[2, 5])^2 + (mu[1, 7] - mu[2, 6])^2 
                                - (mu[1, 5] - mu[2, 6])^2 - (mu[1, 7] - mu[2, 5])^2))
}

generate.cond.cor <- function(a1, r, a2R, a2NR, t1, t2, spltime, design, rprob,
                              sigma.t1, sigma.t2, sigma.t1.ref, sigma.t2.ref,
                              rho, rho.ref, gammas, lambdas) {
  if (t1 == t2) {
    rho <- rho.r <- 1
  }
  # rmean.t1  <- conditional.model(a1, 1, a2R, a2NR, t1, spltime, design, rprob, gammas, lambdas)
  # rmean.t2  <- conditional.model(a1, 1, a2R, a2NR, t2, spltime, design, rprob, gammas, lambdas)
  # nrmean.t1 <- conditional.model(a1, 0, a2R, a2NR, t1, spltime, design, rprob, gammas, lambdas)
  # nrmean.t2 <- conditional.model(a1, 0, a2R, a2NR, t2, spltime, design, rprob, gammas, lambdas)
  
  sigma.t1.cond <- generate.cond.sd(a1, r, a2R, a2NR, t1, spltime, design, rprob, sigma.t1, sigma.t1.ref, gammas, lambdas)
  sigma.t2.cond <- generate.cond.sd(a1, r, a2R, a2NR, t2, spltime, design, rprob, sigma.t2, sigma.t2.ref, gammas, lambdas)
  
  (1/(r * rprob + (1 - r) * (1 - rprob))) * # upweight by appropriate response probability
    ((rho * sigma.t1 * sigma.t2) - # compute marginal covariance
       (r * (1 - rprob) + (1 - r) * rprob) * # weight by appropriate response probability
       (rho.ref * sigma.t1.ref * sigma.t2.ref)) / # compute conditional reference covariance
    (sigma.t1.cond * sigma.t2.cond) # divide by product of standard deviations in conditional target group
  
  # (-1/(sigma.t1.nr * sigma.t2.nr)) * (1 / (1 - rprob)) * 
  #   (rprob * rho.ref * sigma.t1.r * sigma.t2.r -
  #      rho * sigma.t1 * sigma.t2 +
  #      rprob * (1 - rprob) * (rmean.t1 - nrmean.t1) * (rmean.t2 - nrmean.t2))
}

generate.cond.sd <- function(a1, r, a2R, a2NR, obsTime, spltime, design, rprob, 
                             sigma, sigma.cond, gammas, lambdas) {
  rmean  <- conditional.model(a1, 1, a2R, a2NR, obsTime, spltime, design, rprob, gammas, lambdas)
  nrmean <- conditional.model(a1, 0, a2R, a2NR, obsTime, spltime, design, rprob, gammas, lambdas)
  
  if (obsTime <= spltime) {
    sigma
  } else {
    sqrt((1 / (r * rprob + (1 - r) * (1 - rprob))) *
           (
             sigma^2 - ((1 - r) * rprob + r * (1 - rprob)) * sigma.cond ^ 2 -
               rprob * (1 - rprob) * (rmean - nrmean)^2
           ))
  }
}

marginal.model <- function(a1, a2R, a2NR, obsTime, spltime, design, gammas) {
  if (obsTime <= spltime) {
    mu <- gammas[1] + obsTime * (gammas[2] + gammas[3] * a1)
  } else {
    if (design == 1) {
      if (length(gammas) != 9) stop("for design 1 gammas must be length 9.")
      mu <- gammas[4] + gammas[5] * a1 + gammas[6] * a2R + gammas[7] * a2NR +
        gammas[8] * a1 * a2R + gammas[9] * a1 * a2NR
    } else if (design == 2) {
      if (length(gammas) != 7) stop("for design 1 gammas must be length 7.")
      mu <- gammas[4] + gammas[5] * a1 + gammas[6] * a2NR + gammas[7] * a1 * a2NR
    } else if (design == 3) {
      if (length(gammas) != 6) stop("for design 3 gammas must be length 6.")
      mu <- gammas[4] + gammas[5] * a1 + gammas[6] * (a1 == 1) * a2NR
    }
    mu <- gammas[1] + spltime * (gammas[2] + gammas[3] * a1) + (obsTime - spltime) * mu
  }
  mu
}

meanvec <- function(times, spltime, design, gammas) {
  A <- dtrIndex(design)
  do.call(cbind, lapply(1:length(A$a1), function(dtr) {
    do.call(rbind, lapply(times, function(obsTime) {
      marginal.model(A$a1[dtr], A$a2R[dtr], A$a2NR[dtr], obsTime, spltime, design, gammas)
    }))
  }))
}

meat.compute <- function(d, V, times, spltime, design, gammas) {
  # d should be LONG
  n <- length(unique(d$id))
  nDTR <- switch(design, 8, 4, 3)
  mvec <- meanvec(times, spltime, design, gammas)
  dmat <- mod.derivs(times, spltime, design)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR)
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  
  resids <- matrix(rep(d$Y, nDTR), ncol = nDTR) -
    matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE)
  resids <- cbind("id" = d$id, resids * d[, grep("dtr", names(d))] *
                    matrix(rep(d$weight, nDTR), ncol = nDTR))
  
  # sum over IDs
  Reduce("+", lapply(split.data.frame(resids, resids$id), function(obs) {
    # compute weighted residual matrix for individual
    # resid <- obs$weight * obs[, grepl("dtr", names(obs))] * 
    # (obs$Y - meanvec(times, spltime, design, gammas))
    # Sum over DTRs
    m <- Reduce("+", lapply(1:nDTR, function(dtr) {
      t(dmat[[dtr]]) %*% solve(V[[dtr]]) %*% as.matrix(obs[[paste0("dtr", dtr)]], ncol = 1)
    }))
    m %*% t(m)
  })) / n
}

mod.derivs <- function(times, spltime, design) {
  nDTR <- switch(design, 8, 4, 3)
  A <- dtrIndex(design)
  
  dmat <- lapply(1:nDTR, function(dtr) {
    do.call(rbind, lapply(times, function(ti) {
      matrix(c(
        1,
        (ti <= spltime) * matrix(c(ti, ti * A$a1[dtr], rep(0, 6)), nrow = 1) +
          (ti > spltime) * matrix(
            c(
              spltime,
              spltime * A$a1[dtr],
              (ti - spltime),
              (ti - spltime) * A$a1[dtr],
              (ti - spltime) * A$a2R[dtr],
              (ti - spltime) * A$a2NR[dtr],
              (ti - spltime) * A$a1[dtr] * A$a2R[dtr],
              (ti - spltime) * A$a1[dtr] * A$a2NR[dtr]
            ),
            nrow = 1
          )
      ),
      nrow = 1)
    }))
  })
  
  # Since we constructed dmat using the model for design 1, and the models for
  # designs 2 and 3 are nested inside the design 1 model (see paper supplement),
  # we can remove the columns of dmat that are all zeros, as these correspond to
  # parameters that aren't in the nested model.
  
  # Find columns where all entries are zero for all DTRs
  removeCols <- vector("integer")
  if (design == 2)
    removeCols <- c(6, 8)
  else if (design == 3)
    removeCols <- c(6, 8, 9)
  
  # Remove those columns (if necessary)
  if (length(removeCols) != 0) {
    lapply(dmat, function(x) x[, -removeCols])
  } else
    dmat
}

print.simResult <- function(sim, quote = FALSE) {
  # NOTE: 'quote=FALSE' is provided ONLY for compatibility with slackr. 
  # The argument can be safely ignored.
  
  cat("\tSimulation Results\n\n")
  # cat("Call:\n")
  # print(sim$call)
  # cat("\n\n")
  cat(paste("Input summary:\nnumber of participants =", sim$n, "number of trials = ", sim$niter,
            "number of valid trials = ", sim$valid, "\n"))
  cat(paste("effect size =", sim$delta, "\nr0 =", sim$r0, "r1 =", sim$r1, "rho = ", sim$rho, "\n"))
  cat("\nPower results:\n")
  fp <- format.pval(sim$pval.power, digits = max(1L, getOption('digits') - 3L))
  cat(paste("target power =", sim$power.target,
            "estimated power =", round(sim$power, max(1L, getOption('digits') - 3L)),
            "p-value", ifelse(substr(fp, 1L, 1L) == "<", fp, paste("=", fp)), "\n"))
  cat(paste("\nParameter estimates:\nmodel parameters =", paste0("(", paste(round(sim$param.hat, 5), collapse = ", "), ")"),
            "coverage =", round(sim$coverage, max(1L, getOption('digits') - 3L)), "\n"))
  cat(paste("\nMarginal estimates:\nvariances =", paste0("(", paste(round(sim$sigma2.hat, 5), collapse = ", "), ")"),
            "correlation =", round(sim$rho.hat, 5), "\n"))
  cat("\nEstimate of correlation between response status and products of first-stage residuals:\n")
  print(sim$respCor)
  cat("\n\nEstimate of covariance between t_i>=t* and t_j<t* among responders and non-responders:\n")
  print(sim$condCov)
}

reshapeSigma <- function(sigma, times, design) {
  nDTR <- switch(design, 8, 4, 3)
  
  sigma <- as.matrix(sigma)
  
  # Restructure sigma to a length(times)-by-nDTR matrix
  if (sum(dim(sigma) == c(1, 1)) == 2) {
    # scalar case
    sigma <- matrix(rep(sigma, length(times) * nDTR), nrow = length(times))
  } else if (nrow(sigma) == length(times) & 
             length(grep("time", rownames(sigma))) == length(times) & 
             ncol(sigma) == 1) {
    # case in which there is one row per time (i.e., pool.time = F, pool.dtr = T)
    sigma <- matrix(rep(sigma, nDTR), ncol = nDTR, byrow = F)
  } else if (nrow(sigma) == nDTR & 
             length(grep("dtr", rownames(sigma))) == nDTR & 
             ncol(sigma) == 1) {
    # case in which there is one row per DTR (i.e., pool.time = T, pool.dtr = F)
    sigma <- t(matrix(rep(sigma, length(times)), ncol = length(times), byrow = F))
  }
  
  sigma
}

response.beta <- function(d, gammas, r1, r0, respDirection = NULL, sigma, causal = F) {
  shape1 <- r1 / (1 - r1)
  shape0 <- r0 / (1 - r0)
  
  if (causal) {
    x1 <- pnorm(d$Y1.1, mean = gammas[1] + gammas[2] + gammas[3], sd = sigma)
    x0 <- pnorm(d$Y1.0, mean = gammas[1] + gammas[2] - gammas[3], sd = sigma)
    
    d$respProb.1 <- qbeta(x1, shape1 = shape1, shape2 = 1)
    d$respProb.0 <- qbeta(x0, shape1 = shape0, shape2 = 1)
    
    d$R.1 <- sapply(1:nrow(d), function(i) rbinom(1, 1, d$respProb.1[i]))
    d$R.0 <- sapply(1:nrow(d), function(i) rbinom(1, 1, d$respProb.0[i]))
  } else {
    x <- pnorm(d$Y1, mean = gammas[1] + gammas[2] + gammas[3]*d$A1, sd = sigma)
    respProb <- vector("numeric", nrow(d))
    
    respProb[d$A1 ==  1] <- qbeta(x[d$A1 ==  1], shape1 = shape1, shape2 = 1)
    respProb[d$A1 == -1] <- qbeta(x[d$A1 == -1], shape1 = shape0, shape2 = 1)
    d$respProb <- respProb
    d$R <- sapply(1:nrow(d), function(i) rbinom(1, 1, respProb[i]))
  }
  
  return(list("data" = d, "r1" = r1, "r0" = r0))
}

response.indep <- function(d, gammas, r1, r0, respDirection = NULL, causal = F, ...) {
  if (causal) {
    d$R.1 <- rbinom(nrow(d), 1, r1)
    d$R.0 <- rbinom(nrow(d), 1, r0)
  } else{
    d$R <- NA
    d$R[d$A1 ==  1] <- rbinom(sum(d$A1 ==  1), 1, r1)
    d$R[d$A1 == -1] <- rbinom(sum(d$A1 == -1), 1, r0)
  }
  return(list("data" = d, "r1" = r1, "r0" = r0))
}

response.oneT <- function(d, gammas, r1, r0, respDirection = c("high", "low"), sigma, causal = F) {
  respDirection <- match.arg(respDirection)
  tail <- switch(respDirection, "high" = F, "low" = T)
  
  upsilon <- qnorm(r1, as.numeric(sum(gammas[1:3])), sigma, lower.tail = tail)
  if (causal) {
    d$R.0 <- as.numeric(d$Y1.0 >= upsilon)
    d$R.1 <- as.numeric(d$Y1.1 >= upsilon)
  } else {
    d$R <- as.numeric(d$Y1 >= upsilon)
  }
  r0temp <- pnorm(upsilon, sum(gammas[1:2]) - gammas[3], sigma, lower.tail = FALSE)
  if (r0temp != r0) {
    warning(paste("Overwriting the provided value of r0 to accomodate the oneThreshold respModel.",
                  "The provided value is", round(r0, 3), "and the new value is", round(r0temp, 3)))
    r0 <- r0temp
  }
  return(list("data" = d, "r1" = r1, "r0" = r0))
}

response.sq <- function(d, gammas, r1, r0, respDirection = c("high", "low"), sigma, causal = F) {
  respDirection <- match.arg(respDirection)
  
  upsilon1.up   <- qnorm(r1/2, sum(gammas[1:3]),             sigma, lower.tail = F)
  upsilon0.up   <- qnorm(r0/2, sum(gammas[1:2]) - gammas[3], sigma, lower.tail = F)
  upsilon1.down <- qnorm(r1/2, sum(gammas[1:3]),             sigma)
  upsilon0.down <- qnorm(r0/2, sum(gammas[1:2]) - gammas[3], sigma)
  
  if (causal) {
    d$R.1 <- d$R.0 <- NA
    if (respDirection == "high") {
      d$R.1 <- as.numeric(d$Y1.1 >= upsilon1.up | d$Y1.1 <= upsilon1.down)
      d$R.0 <- as.numeric(d$Y1.0 >= upsilon0.up | d$Y1.0 <= upsilon0.down)
    } else {
      d$R.1 <- as.numeric(d$Y1.1 >= -upsilon1.down & d$Y1.1 <= upsilon1.up)
      d$R.0 <- as.numeric(d$Y1.0 >= -upsilon0.down & d$Y1.0 <= upsilon0.up)
    }
  } else {
    d$R <- NA
    if (respDirection == "high") {
      d$R[d$A1 ==  1] <- as.numeric(d$Y1[d$A1 ==  1] >= upsilon1.up | d$Y1[d$A1 ==  1] <= upsilon1.down)
      d$R[d$A1 == -1] <- as.numeric(d$Y1[d$A1 == -1] >= upsilon0.up | d$Y1[d$A1 == -1] <= upsilon0.down)
    } else {
      d$R[d$A1 ==  1] <- as.numeric(d$Y1[d$A1 ==  1] >= -upsilon1.down & d$Y1[d$A1 ==  1] <= upsilon1.up)
      d$R[d$A1 == -1] <- as.numeric(d$Y1[d$A1 == -1] >= -upsilon0.down & d$Y1[d$A1 == -1] <= upsilon0.up)
    }
  }
  return(list("data" = d, "r1" = r1, "r0" = r0))
}

response.twoT <- function(d, gammas, r1, r0, respDirection = c("high", "low"), sigma, causal = F) {
  respDirection <- match.arg(respDirection)
  tail <- switch(respDirection, "high" = F, "low" = T)
  
  upsilon1 <- qnorm(r1, sum(gammas[1:3]),             sigma, lower.tail = tail)
  upsilon0 <- qnorm(r0, sum(gammas[1:2]) - gammas[3], sigma, lower.tail = tail)
  
  if (causal) {
    d$R.1 <- d$R.0 <- NA
    d$R.1 <- as.numeric(d$Y1.1 >= upsilon1)
    d$R.0 <- as.numeric(d$Y1.0 >= upsilon0)
  } else {
    d$R <- NA
    d$R[d$A1 ==  1] <- as.numeric(d$Y1[d$A1 ==  1] >= upsilon1)
    d$R[d$A1 == -1] <- as.numeric(d$Y1[d$A1 == -1] >= upsilon0)
  }
  
  return(list("data" = d, "r1" = r1, "r0" = r0))
}


#' Compile sim results into a LaTeX-ready table
#'
#' @param results A list of objects of class \code{simResult} to be tabulated
#'
#' @return
#' @export
#'
#' @examples
resultTable <- function(results, alternative = c("two.sided", "less", "greater"), paper = FALSE) {
  alternative <- match.arg(alternative)
  
  d <- do.call("rbind", lapply(results, function(l) {
    if (alternative == "two.sided") {
      pval.power <- l$pval.power
    } else {
      pval.power <- binom.test(x = sum(l$pval <= l$alpha / 2, na.rm = T) + sum(l$pval >= 1 - l$alpha / 2, na.rm = T),
                               n = l$valid, p = l$power.target, alternative = alternative)$p.value
    }
    
    respCor.mean <- apply(l$respCor, 1, function(x) x[which.max(abs(x))])
    
    x <- data.frame("delta" = l$delta,
                    "rprob" = min(l$r0, l$r1),
                    "rho" = l$rho,
                    "rho.size" = ifelse(is.null(l$rho.size), l$rho, l$rho.size),
                    "sharp" = ifelse(l$sharp, "sharp", "cons."),
                    "n" = l$n,
                    "estPwr" = l$power,
                    "pval.power" = pval.power,
                    "coverage" = l$coverage,
                    "pval.coverage" = tryCatch(binom.test(x = l$coverage * l$valid,
                                                          n = l$valid, p = .95, 
                                                          alternative = "two.sided")$p.value,
                                               error = function(e) NULL),
                    "corY0sq" = as.numeric(respCor.mean[1]),
                    "corY1sq" = as.numeric(respCor.mean[3]),
                    "corY0Y1" = as.numeric(respCor.mean[2]),
                    "nTrial" = l$niter,
                    "nValid" = l$valid
    )
    # if (paper) cbind(d, l$design)
  }))
  
  # d <- d[order(d$delta, d$rprob, d$rho, d$sharp), ]
  
  oldColnames <- colnames(d)
  
  colnames(d) <- c("$\\delta$", "$\\min\\left\\{r_{-1},r_{1}\\right\\}$", "$\\rho$",
                   "$\\rho_{\\text{size}}$",
                   "Formula", "$n$", "$1-\\hat{\\beta}$", "$p$ value", "Coverage",
                   "$p$ value", 
                   "$\\cor(R,Y_0^2)$", "$\\cor(R,Y_1^2)$", "$\\cor(R,Y_0Y_1)$",
                   "Num. Runs", "Valid Runs")
  
  if (paper) {
    
  }
  
  print(xtable(d, digits = c(1,1,1,1,1,1,0,3,3,3,3,3,3,3,0,0)),
        sanitize.text.function = identity, booktabs = TRUE,
        include.rownames = F,
        include.colnames = T)
  
  colnames(d) <- oldColnames
  return(invisible(d))
}

sample.size <- function(delta, r, r1 = r, r0 = r, rho, alpha = 0.05, power = .8,
                        design = 2, rounding = c("up", "down"), conservative = TRUE) {
  rounding <- match.arg(rounding)
  
  if (is.null(r)) r <- mean(c(r1, r0))
  
  # Input checks
  if (!is.null(rounding) & !(rounding %in% c("up", "down"))) stop("rounding must be either up or down")
  nid <- (4 * (qnorm(1 - alpha / 2) + qnorm(power))^2) / delta^2
  correction <- 1 - rho^2
  if (design == 1) {
    designEffect <- 2
    if (!conservative) {
      correction <- (1-rho^2) - (1-rho)*rho^2 / (2*(1+rho))
    }
  } else if (design == 2) {
    designEffect <- ((2 - r1) + (2 - r0)) / 2
    if (!conservative) {
      correction <- ((1 - rho) * (2 * rho + 1)) / (1 + rho) + (1 - rho)/(1 + rho) * rho^2 / (2 - r)
    }
  } else if (design == 3) {
    designEffect <- (3 - r) / 2
    if (!conservative) {
      correction <- (1-rho)*(3-r+6*rho-2*r*rho+2*rho^2)/((3-r)*(1+rho))
    }
  }
  else stop("Not a valid design indicator.")
  
  if (rounding == "up") {
    n <- ceiling(nid * designEffect * correction)
  } else if (rounding == "down") {
    n <- floor(nid * designEffect * correction)
  }
  message('generating sample size')
  return(n)
}

### Wrapper for all the estimation functions
#' Estimation for repeated-measures SMART data
#'
#' @param data A `data.frame` with the appropriate structure (see `generateSMART`)
#' @param corstr A string indicating the type of working correlation structure. One of `"identity"`/`"id"`, `"exchangeable"`/`"exch"`. 
#' @param times A vector of times at which the repeated-measures outcome has been collected
#' @param spltime The time point (in `times`) immediately after which response/non-response is assessed and participants are re-randomized
#' @param design An integer, one of `1`, `2`, or `3`, indicating the design type of the SMART
#' @param start A vector of the same length as the regression parameter vector indicating the start point for the solver
#' @param maxiter.solver The maximum number of iterations the solver is allowed to perform until convergence. Defaults to `1000`
#' @param tol The tolerance required to achieve covergence. Defaults to `10e-8`
#' @param pool.time Logical. Should estimates of the variance be pooled over time?
#' @param pool.dtr Logical. Should estimates of the variance be pooled over DTR?
#'
#' @return
#' @export
#'
#' @examples
SMART.estimate <- function(data, corstr = c("identity", "exchangeable", "ar1"),
                           times, spltime, design,
                           start, maxiter.solver = 1000, tol = 10e-8,
                           pool.time = FALSE, pool.dtr = FALSE) {
  corstr <- match.arg(corstr)
  
  d1 <- reshape(data, varying = list(grep("Y", names(data))), ids = data$id, 
                times = times, direction = "long", v.names = "Y")
  d1 <- d1[order(d1$id, d1$time), ]
  
  param.hat <- estimate.params(d1, diag(rep(1, length(times))), times, spltime,
                               design, start, maxiter.solver, tol)
  
  sigma2.hat <- estimate.sigma2(d1, times, spltime, design, param.hat,
                                pool.time = constant.var.time, pool.dtr = constant.var.dtr)
  if (corstr != "id") 
    rho.hat <- estimate.rho(d1, times, spltime, design, sqrt(sigma2.hat), param.hat)
  
  # Compute variance matrices for all conditional cells
  condVars <- lapply(split.SMART(d$data), function(x) {
    var(subset(x, select = grep("Y", names(x), value = TRUE)))
  })
  
  # Iterate parameter estimation
  if (corstr == "identity" | (corstr == "exchangeable" & rho == 0)) {
    outcome.var <- varmat(sigma2.hat, 0, times, design, "exch")
    param.var <- estimate.paramvar(d1, diag(rep(1, length(times))), times, spltime, design, gammas = param.hat)
    iter <- 1
  } else if (corstr == "exchangeable" & rho != 0) {
    outcome.var <- varmat(sigma2.hat, rho.hat, times, design, "exch")
    param.hat <- estimate.params(d1, outcome.var, times, spltime, design, param.hat, maxiter.solver, tol)
    # Iterate until estimates of gammas and rho converge
    for (i in 1:maxiter.solver) {
      sigma2.new <- estimate.sigma2(d1, times, spltime, design, param.hat,
                                    pool.time = constant.var.time, pool.dtr = constant.var.dtr)
      rho.new <- estimate.rho(d1, times, spltime, design, sqrt(sigma2.hat), param.hat)
      outcomeVar.new <- varmat(sigma2.new, rho.new, times, design, "exch")
      param.new <- estimate.params(d1, outcomeVar.new, times, spltime, design, start = param.hat, maxiter.solver, tol)
      if (norm(param.new - param.hat, type = "F") <= tol & norm(as.matrix(sigma2.new) - as.matrix(sigma2.hat), type = "F") <= tol &
          (rho.new - rho.hat)^2 <= tol) {
        param.hat <- param.new
        sigma2.hat <- sigma2.new
        rho.hat <- rho.new
        iter <- i
        break
      } else {
        param.hat <- param.new
        sigma2.hat <- sigma2.new
        rho.hat <- rho.new
        iter <- i
      }
    }
    param.var <- estimate.paramvar(d1, cormat(rho.hat, length(times), corstr), times, spltime, design, gammas = param.hat)
  }
  
  
}

### split.data.frame doesn't quite work to split data from a SMART into conditional cells.
### This function does that.
split.SMART <- function(d, marginal = FALSE) {
  if (marginal) {
    l <- lapply(grep("dtr", names(d), value = TRUE), function(s) {
      subset(d, get(s) == 1)
    })
    names(l) <- grep("dtr", names(d), value = TRUE)
  } else {
    l <- split.data.frame(d, list(d$A1, d$R, d$A2R, d$A2NR))
    x <- apply(unique(subset(d, select = c("A1", "R", "A2R", "A2NR"))), 1, function(v) {
      paste(v, collapse = ".")
    })
    l <- lapply(1:length(x), function(i) l[[x[i]]])
    names(l) <- x
  }
  l
}

### Check if generated trial is 'valid'; i.e., there is at least one participant per cell
validTrial <- function(d, design) {
  l <- split.SMART(d)
  if (design == 1) {
    length(l) == 8
  } else if (design == 2) {
    length(l) == 6
  } else if (design == 3) {
    length(l) == 5
  } else stop("design must be in 1-3.")
}

### Construct V matrix 
varmat <- function(sigma2, rho, times, design, corstr = c("identity", "exchangeable", "ar1", "unstructured")) {
  nDTR <- switch(design, 8, 4, 3)
  corstr <- match.arg(corstr)
  
  sigma2 <- reshapeSigma(sigma2, times, design)
  
  if (length(rho) == 1 & corstr == "unstructured") {
    warning("rho is length 1 with unstructured corstr. Setting corstr to exchangeable.")
    corstr <- "exchangeable"
  } else if (corstr == "unstructured") {
    if (!is.matrix(rho) | (is.matrix(rho) & !all.equal(dim(rho), c(length(times), nDTR)))) {
      stop("rho is not of proper dimension: must be T by nDTR")
    } else {
      lapply(1:nDTR, function(dtr) {
        diag(sqrt(sigma2[, dtr])) %*% cormat(rho[, dtr], length(times),corstr) %*% diag(sqrt(sigma2[, dtr]))
      })
    }
  } else {
    lapply(1:nDTR, function(dtr) {
      diag(sqrt(sigma2[, dtr])) %*% cormat(rho, length(times),corstr) %*% diag(sqrt(sigma2[, dtr]))
    })
  }
}
