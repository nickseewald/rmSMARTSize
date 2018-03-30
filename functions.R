### Functions for simulations of SMARTs with continuous repeated-measures outcomes
### Nick Seewald

### Custom combine function for the foreach %dopar% loop in sim()
combine.results <- function(list1, list2) {
  pval       <- c(list1$pval, list2$pval)
  param.hat  <- rbind(list1$param.hat, list2$param.hat)
  sigma2.hat <- rbind(list1$sigma2.hat, list2$sigma2.hat)
  param.var  <- Reduce(function(x, y) {x[is.na(x)] <- 0; y[is.na(y)] <- 0; x + y},
                       list(list1$param.var, list2$param.var))
  rho.hat    <- c(list1$rho.hat, list2$rho.hat)
  coverage   <- list1$coverage + list2$coverage
  iter       <- c(list1$iter, list2$iter)
  valid      <- list1$valid + list2$valid
  condVars   <- Map(function(x, y) {x[is.na(x)] <- 0; y[is.na(y)] <- 0; x + y},
                    list1$condVars, list2$condVars)

  result <- list("pval" = pval, "param.hat" = param.hat, "sigma2.hat" = sigma2.hat, "iter" = iter,
                 "param.var" = param.var, "rho.hat" = rho.hat, "valid" = valid, "coverage" = coverage,
                 "condVars" = condVars)

  if ("data" %in% names(list1)) {
    dataList <- c(list1$data, list2$data)
    result[["data"]] <- dataList
  }
  
  return(result)
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
conditional.model <- function(a1, r, a2R, a2NR, t, spltime, design, rprob, gammas, lambdas) {
  if (t <= spltime) {
    mu <- gammas[1] + t * (gammas[2] + gammas[3] * a1)
  } else {
    if (design == 1) {
      if (length(gammas) != 9) stop("for design 1 gammas must be length 9.")
      mu <- (t - spltime) * (gammas[4] + gammas[5] * a1 + r * (gammas[6] * a2R + gammas[8] * a1 * a2R) / rprob + 
                               (1 - r) * (gammas[7] * a2NR + gammas[9] * a1 * a2NR) / (1 - rprob) +
                               (r - rprob) * (lambdas[1] + lambdas[2] * a1))
    } else if (design == 2) {
      if (length(gammas) != 7) stop("for design 2, gammas must be length 7.")
      mu <- (t - spltime) * (gammas[4] + gammas[5] * a1 + gammas[6] * (1 - r) * a2NR + gammas[7] * a1 * (1 - r) * a2NR) +
        (t - spltime) * ((rprob / (1 - rprob)) * (gammas[6] * (1 - r) * a2NR + gammas[7] * (1 - r) *a1 * a2NR)) + 
        (t - spltime) * (r - rprob) * (lambdas[1] + lambdas[2] * a1)
    } else if (design == 3) {
      if (length(gammas) != 6) stop("for design 3, gammas must be length 6.")
      mu <- (t - spltime) * (gammas[4] + gammas[5] * a1 + gammas[6] * (a1 == 1) * a2NR * (1 - r) / (1 - rprob)) +
        (t - spltime) * (r - rprob) * (lambdas[1] + lambdas[2] * a1)
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
  
  if (length(sigma) == 1) {
    sigma <- rep(sigma, length(times))
  } else if (length(sigma) != length(times)) {
    stop("sigma must either be length-1 or the same length as times.")
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
    sigma.r1 <- c(sigma[times <= spltime], sigma.r1)
  }
  if (length(sigma.r0) != length(times)) {
    sigma.r0 <- c(sigma[times <= spltime], sigma.r0)
  }
  
  varEnv <- new.env()
  corstr <- match.arg(corstr)
  A <- conditionalIndex(design)
  
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
  if (corstr == "identity") {
    cormat.r1 <- cormat.r0 <- diag(rep(1, length(times)))
  } else if (corstr == "exchangeable") {
    cormat.r1 <- cormat.r0 <- cormat.exch(rho, length(times))
  } else if (corstr == "ar1") {
    cormat.r1 <- cormat.r0 <- cormat.ar1(rho, length(times))
  }
  # Create responders' correlation matrices by looping over m and computing appropriate entries
  for (i in 1:dim(m)[1]) {
    cormat.r1[m[i, 1], m[i, 2]] <- rho.r1
    cormat.r0[m[i, 1], m[i, 2]] <- rho.r0
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
        sapply(times, function(t) {
          generate.sd(sigma[which(times == t)], a1 = dtr[1], a2R = dtr[2], a2NR = dtr[3], t = t, spltime = spltime, 
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
               generate.cond.sd(a1 = 2 * a1ind - 1, r = 0, a2R = 1, a2NR = 1, t = times[i], spltime = spltime,
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
               generate.cond.sd(a1 = 2 * a1ind - 1, r = 1, a2R = -1, a2NR = 1, t = times[i], spltime = spltime,
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
               generate.cond.sd(a1 = 2 * a1ind - 1, r = 0, a2R = -1, a2NR = -1, t = times[i], spltime = spltime,
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
      if (corstr %in% c("exchangeable", "identity")) {
        x1 <- x2 <- x3 <- cormat.exch(rho, length(times))
      } else if (corstr == "ar1") {
        x1 <- x2 <- x3 <- cormat.ar1(rho, length(times))
      }
      for(i in 1:nrow(m)) {
        x1[m[i, 1], m[i, 2]] <- 
          generate.cond.cor(a1 = 2 * a1ind - 1, r = 0, a2R = 1, a2NR = 1,
                            t1 = times[m[i, 1]], t2 = times[m[i, 2]], spltime = spltime,
                            design = 1, rprob = get(paste0("r", a1ind)),
                            sigma.t1 = sigma[m[i, 1]], sigma.t2 = sigma[m[i, 2]], 
                            sigma.t1.ref = get(paste0("sigma.r", a1ind, 1), envir = varEnv)[m[i, 1]],
                            sigma.t2.ref = get(paste0("sigma.r", a1ind, 1), envir = varEnv)[m[i, 2]],
                            rho = rho,
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
                            rho = rho,
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
                            rho = rho,
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
    
    dtrs <- dtrIndex(2)
    if ((!is.null(uneqsdDTR) & is.null(uneqsd)) | (is.null(uneqsdDTR) & !is.null(uneqsd)))
      stop("For design 2, you must provide either both uneqsdDTR and uneqsd or neither.")
    invisible(sapply(1:4, function(i) {
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
        sigmaStar <- sigma
      }
      
      # Generate conditional SDs for non-responders and assign them to appropriate variables.
      assign(paste0("sigma.nr", (dtrs$a1[i] + 1) / 2, (dtrs$a2NR[i] + 1) / 2),
             unname(sapply(1:length(times), function(j) {
               generate.cond.sd(a1 = dtrs$a1[i], r = 0, a2R = dtrs$a2R[i], a2NR = dtrs$a2NR[i],
                                t = times[j], spltime, design = 2, 
                                rprob = get(paste0("r", (dtrs$a1[j] + 1) / 2)),
                                sigma = sigmaStar[j],
                                sigma.cond = get(paste0("sigma.r", (dtrs$a1[i] + 1) / 2, 0), envir = varEnv)[j],
                                gammas = gammas, lambdas = lambdas)
             })), envir = varEnv)
    }))
    
    NRgrid <- expand.grid("a1ind" = c(0, 1), "a2NRind" = c(0, 1))
    
    for (j in 1:nrow(NRgrid)) {
      a1ind <- NRgrid$a1ind[j]
      a2NRind <- NRgrid$a2NRind[j]
      
      if (corstr %in% c("exchangeable", "identity")) {
        x <- cormat.exch(rho, length(times))
      } else if (corstr == "ar1") {
        x <- cormat.ar1(rho, length(times))
      }
      
      for (i in 1:nrow(m)) {
        x[m[i, 1], m[i, 2]] <- 
          generate.cond.cor(a1 = 2 * a1ind - 1, r = 0, a2R = 0, a2NR = 2 * a2NRind - 1,
                            t1 = times[m[i, 1]], t2 = times[m[i, 2]], spltime = spltime,
                            design = 2, rprob = get(paste0("r", a1ind)),
                            sigma.t1 = sigma[m[i, 1]], sigma.t2 = sigma[m[i, 2]],
                            sigma.t1.ref = get(paste0("sigma.r", a1ind, 0), envir = varEnv)[m[i, 1]],
                            sigma.t2.ref = get(paste0("sigma.r", a1ind, 0), envir = varEnv)[m[i, 2]],
                            rho = rho,
                            rho.ref = get(paste0("cormat.r", a1ind, 0), envir = varEnv)[m[i, 1], m[i, 2]],
                            gammas = gammas, lambdas = lambdas)
      }
      assign(paste0("cormat.nr", a1ind, a2NRind), x, envir = varEnv)
    }
  } else 
    stop("design 3 hasn't been implemented in conditionalVarmat")
  
  
  condVarmats <- lapply(1:nrow(A), function(dtr) {
    a1ind   <- (A$a1[dtr] + 1) / 2
    R       <- A$r[dtr]
    a2Rind  <- ifelse(A$a2R[dtr] == 0, 0, (A$a2R[dtr] + 1) / 2)
    a2NRind <- (A$a2NR[dtr] + 1) / 2
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

### Construct a t-by-t AR(1) correlation matrix with parameter rho
cormat.ar1 <- function(rho, t) {
  m <- diag(t)
  rho^(abs(row(m) - col(m)))
}

### Construct a t-by-t exchangeable correlation matrix with correlation rho
cormat.exch <- function(rho, t) {
  m <- matrix(rep(rho, t^2), nrow = t)
  diag(m) <- rep(1,t)
  m
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
  
  DVinv <- lapply(1:nDTR, function(i) t(deriv[[i]]) %*% solve(V[[i]]))
  
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
  
  Reduce("+", lapply(split.data.frame(d, d$id), function(x) {
    Reduce("+", lapply(1:nDTR, function(dtr) {
      unique(x[[paste0("dtr", dtr)]]) * unique(x$weight) *
        t(deriv[[dtr]]) %*% solve(V) %*% deriv[[dtr]]
    }))
  })) / -length(unique(d$id))
}

  estimate.params <- function(d, V, times, spltime, design, start, maxiter.solver, tol) {
  params.hat <- start
  j <- 0
  epsilon <- 10
  J <- esteqn.jacobian(d, V, times, spltime, design)
  
  for (j in 1:maxiter.solver) {
    params.hat.new <- params.hat - solve(J) %*%
      esteqn.compute(d, V, times, spltime, design, gammas = params.hat)
    epsilon <- norm(params.hat.new - params.hat, type = "F")
    params.hat <- params.hat.new
    if (epsilon <= tol) break
  }
  params.hat
}

estimate.paramvar <- function(d, V, times, spltime, design, gammas) {
  n <- length(unique(d$id))
  J <- esteqn.jacobian(d, V, times, spltime, design)
  solve(-J) %*%  meat.compute(d, V, times, spltime, design, gammas) %*% solve(-J) / n
}

estimate.rho <- function(d, times, spltime, design, sigma, gammas, corstr = "exchangeable") {
  n <- length(unique(d$id))
  nDTR <- switch(design, 8, 4, 3)
  mvec <- meanvec(times, spltime, design, gammas)
  Dmat <- mod.derivs(times, spltime, design)
  
  if (any(is(d$Y) == "NULL")) stop("d has to be in long format.")
  if (length(sigma) != length(times) & length(sigma) != 1)
    stop("sigma must be length one OR a vector with the same length as times.")
  
  # If a time-invariant sigma is provided, replicate it so there's one copy per time
  if (length(sigma) == 1)
    sigma <- rep(sigma, length(times))
  
  # Compute residuals (Y_{it} - mu_{t}(a1, a2))
  resids <- matrix(rep(d$Y, nDTR), ncol = nDTR) -
    matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE)
  
  # Create data.frame from resids, indexed by time and ID
  # Multiply residuals by DTR indicators (dtr1[i] = 1 iff i is consistent with DTR 1)
  resids <- cbind("id" = d$id, "time" = d$time,
                  resids * d[, grep("dtr", names(d))])
  
  if (corstr == "exchangeable") {
    # For every id and every dtr, compute (\sum_{s<t} r_{s} * r_{t} / sigma_{s} sigma_t()
    r <- do.call(rbind,
                 lapply(split.data.frame(resids, resids$id),
                        function(d) {
                          sapply(1:sum(grepl("dtr", names(d))), function(dtr) {
                            sum(sapply(2:length(d$time), function(t) {
                              sum((d[t, paste0("dtr", dtr)] / sigma[t]) *
                                    (d[1:(t - 1), paste0("dtr", dtr)] / sigma[1:(t - 1)]))
                            }))
                          })
                        })
    )
    
    # Extract weights from data frame
    weights <- aggregate(weight ~ id, d, unique)
    
    # Weight residual sums
    r <- weights$weight * r
    
    # Construct numerator for estimator of rho in exchangeable corstr
    numerator <- apply(r, 2, sum)
    
    # Create matrix of weights per person-time per DTR
    weightmat <- cbind("id" = d$id, "time" = d$time, d$weight * d[, grep("dtr", names(d))])
    
    # Sum weights over individuals 
    # (but only use one weight per person -- weightmat has duplicated rows: 1 per time)
    sumweights <- apply(Reduce(function(...) merge(..., by = "time"), 
                               lapply(1:sum(grepl("dtr", names(d))), function(dtr) {
                                 x <- list(weightmat[[paste0("dtr", dtr)]])
                                 sumwts <- aggregate(x = setNames(x, paste0("dtr", dtr)),
                                                     by = list("time" = resids$time), sum)
                               }))[, -1],
                        2, unique)
    
    # Construct denominator for estimator of alpha in exchangeable corstr
    denominator <- sumweights * (length(times) * (length(times) - 1) / 2) - length(gammas)
    
    # Average over DTRs
    return(mean(numerator/denominator))
  } else NULL
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
  numerator / denominator
}


# Clean up output of combine.results, to be used in the foreach loop after
# combining everything
finalize.results <- function(x) {
  param.hat     <- apply(as.matrix(x$param.hat), 2, mean, na.rm = T)
  sigma2.hat    <- apply(as.matrix(x$sigma2.hat), 2, mean, na.rm = T)
  param.var     <- x$param.var / x$valid
  param.var.est <- apply(as.matrix(x$param.hat), 2, var, na.rm = T)
  rho.hat       <- mean(x$rho.hat, na.rm = T)
  coverage      <- as.numeric(x$coverage / x$valid)
  condVars      <- lapply(x$condVars, function(v) v / x$valid)

  cat("Finishing...\n")

  result <- list("n" = x$n, "alpha" = x$alpha, "power.target" = x$power.target, "delta" = x$delta,
                 "niter" = x$niter, "corstr" = x$corstr, "pval" = x$pval, "param.hat" = param.hat,
                 "sigma2.hat" = sigma2.hat, "iter" = x$iter, "param.var" = param.var,
                 "rho.hat" = rho.hat, "valid" = x$valid, "coverage" = coverage, "condVars" = condVars)
  
  if ("data" %in% names(x))
    result[["data"]] <- x$data
  
  return(result)
}

generate.sd <- function(sigma, a1, a2R, a2NR, t, spltime, design, rprob, gammas, lambdas) {
  dtrs <- dtrIndex(design)
  mu <- sapply(1:length(dtrs$a1), function(i) {
    c(conditional.model(dtrs$a1[i], 1, dtrs$a2R[i], 0, t, spltime, design, rprob, gammas, lambdas),
    conditional.model(dtrs$a1[i], 0, 0, dtrs$a2NR[i], t, spltime, design, rprob, gammas, lambdas))
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

generate.cond.sd <- function(a1, r, a2R, a2NR, t, spltime, design, rprob, 
                           sigma, sigma.cond, gammas, lambdas) {
  rmean  <- conditional.model(a1, 1, a2R, a2NR, t, spltime, design, rprob, gammas, lambdas)
  nrmean <- conditional.model(a1, 0, a2R, a2NR, t, spltime, design, rprob, gammas, lambdas)
  
  if (t <= spltime) {
    sigma
  } else {
    sqrt((1 / (r * rprob + (1 - r) * (1 - rprob))) *
           (
             sigma^2 - ((1 - r) * rprob + r * (1 - rprob)) * sigma.cond ^ 2 -
               rprob * (1 - rprob) * (rmean - nrmean)^2
           ))
  }
}

marginal.model <- function(a1, a2R, a2NR, t, spltime, design, gammas) {
  if (t <= spltime) {
    mu <- gammas[1] + t * (gammas[2] + gammas[3] * a1)
  } else {
    if (design == 1) {
      if (length(gammas) != 9) stop("for design 1 gammas must be length 9.")
      mu <- gammas[4] + gammas[5] * a1 + gammas[6] * a2R + gammas[7] * a2NR + gammas[8] * a1 * a2R + gammas[9] * a1 * a2NR
    } else if (design == 2) {
      if (length(gammas) != 7) stop("for design 1 gammas must be length 7.")
      mu <- gammas[4] + gammas[5] * a1 + gammas[6] * a2NR + gammas[7] * a1 * a2NR
    } else if (design == 3) {
      if (length(gammas) != 6) stop("for design 3 gammas must be length 6.")
      mu <- gammas[4] + gammas[5] * a1 + gammas[6] * (a1 == 1) * a2NR
    }
    mu <- gammas[1] + spltime * (gammas[2] + gammas[3] * a1) + (t - spltime) * mu
  }
  mu
}

meanvec <- function(times, spltime, design, gammas) {
  A <- dtrIndex(design)
  do.call(cbind, lapply(1:length(A$a1), function(dtr) {
    do.call(rbind, lapply(times, function(t) {
      marginal.model(A$a1[dtr], A$a2R[dtr], A$a2NR[dtr], t, spltime, design, gammas)
    }))
  }))
}

meat.compute <- function(d, V, times, spltime, design, gammas) {
  # d should be LONG
  
  n <- length(unique(d$id))
  nDTR <- switch(design, 8, 4, 3)
  mvec <- meanvec(times, spltime, design, gammas)
  dmat <- mod.derivs(times, spltime, design)
  
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
      t(dmat[[dtr]]) %*% solve(V) %*% as.matrix(obs[[paste0("dtr", dtr)]], ncol = 1)
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

print.simResult <- function(sim) {
  cat("\tSimulation Results\n")
  cat(paste("input summary:\nnumber of participants =", sim$n, "number of trials = ", sim$niter, "number of valid trials = ", sim$valid, "\n"))
  cat(paste("effect size =", sim$delta, "\nr0 =", sim$r0, "r1 =", sim$r1, "rho = ", sim$rho, "\n"))
  cat("power results:\n")
  fp <- format.pval(sim$pval.power, digits = max(1L, getOption('digits') - 3L))
  cat(paste("target power =", sim$power.target, "estimated power =", sim$power, "p-value",
            ifelse(substr(fp, 1L, 1L) == "<", fp, paste("=", fp)), "\n"))
  cat(paste("parameter estimates:\nmodel parameters =", paste0("(", paste(round(sim$param.hat, 5), collapse = ", "), ")"),
            "coverage =", sim$coverage, "\n"))
  cat(paste("marginal estimates:\nvariances =", paste0("(", paste(round(sim$sigma2.hat, 5), collapse = ", "), ")"),
            "correlation =", round(sim$rho.hat, 5), "\n"))
}

#' Compile sim results into a LaTeX-ready table
#'
#' @param results A list of objects of class \code{simResult} to be tabulated
#'
#' @return
#' @export
#'
#' @examples
resultTable <- function(results) {
  d <- do.call("rbind", lapply(results, function(l) {
    data.frame("delta" = l$delta,
               "rprob" = min(l$r0, l$r1),
               "rho" = l$rho,
               "sharp" = ifelse(l$sharp, "sharp", "cons."),
               "n" = l$n,
               "estPwr" = l$power,
               "pval.power" = l$pval.power,
               "coverage" = l$coverage,
               "pval.coverage" = tryCatch(binom.test(x = l$coverage * l$valid,
                                            n = l$valid, p = .95, 
                                            alternative = "two.sided")$p.value,
                                          error = function(e) NULL),
               "nTrial" = l$niter,
               "nValid" = l$valid
               )
  }))
  
  colnames(d) <- c("$\\delta$", "$\\min\\left\\{r_{-1},r_{1}\\right\\}$", "$\\rho$",
                   "Formula", "$n$", "$1-\\hat{\\beta}$", "$p$ value", "Coverage",
                   "$p$ value", "Num. Runs", "Valid Runs")
  
  print(xtable(d, digits = c(1,1,1,1,1,0,3,3,3,3,0,0)),
        sanitize.text.function = identity, booktabs = TRUE,
        include.rownames = F,
        include.colnames = T)
}

sample.size <- function(delta, r, rho, alpha = 0.05, power = .8,
                        design = 2, round = "up", conservative = TRUE) {
  # Input checks
  if (!(round %in% c("up", "down"))) stop("round must be either up or down")
  nid <- (4 * (qnorm(1 - alpha / 2) + qnorm(power))^2) / delta^2
  correction <- 1 - rho^2
  if (design == 1) {
    designEffect <- 2
    if (!conservative) {
      correction <- (1-rho^2) - (1-rho)*rho^2 / (2*(1+rho))
    }
  } else if (design == 2) {
    designEffect <- 2 - r
    if (!conservative) {
      correction <- ((1 - rho) * (2 * rho + 1)) / (1 + rho) + (1 - rho)/(1 + rho) * rho^2 / (2 - r)
    }
  } else if (design == 3) {
    designEffect <- (3 - r) / 2
  }
  else stop("Not a valid design indicator.")
  
  if (round == "up") {
    ceiling(nid * designEffect * correction)
  } else if (round == "down") {
    floor(nid * designEffect * correction)
  }
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
    outcome.var <- varmat(sigma2.hat, 0, times, "exch")
    param.var <- estimate.paramvar(d1, diag(rep(1, length(times))), times, spltime, design, gammas = param.hat)
    iter <- 1
  } else if (corstr == "exchangeable" & rho != 0) {
    outcome.var <- varmat(sigma2.hat, rho.hat, times, "exch")
    param.hat <- estimate.params(d1, outcome.var, times, spltime, design, param.hat, maxiter.solver, tol)
    # Iterate until estimates of gammas and rho converge
    for (i in 1:maxiter.solver) {
      sigma2.new <- estimate.sigma2(d1, times, spltime, design, param.hat,
                                    pool.time = constant.var.time, pool.dtr = constant.var.dtr)
      rho.new <- estimate.rho(d1, times, spltime, design, sqrt(sigma2.hat), param.hat)
      outcomeVar.new <- varmat(sigma2.new, rho.new, times, "exch")
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
    param.var <- estimate.paramvar(d1, cormat.exch(rho.hat, length(times)), times, spltime, design, gammas = param.hat)
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
# NB: ONLY WORKS WHEN POOLING OVER DTRs RIGHT NOW!!
varmat <- function(sigma2, alpha, times, corstr = c("identity", "exchangeable", "ar1")) {
  corstr <- match.arg(corstr)
  if (length(sigma2) != 1 & length(sigma2) != length(times)) {
    stop("sigma must be either length 1 (assuming constant variance over time) or have the same length as times.")
  } else if (length(sigma2) == 1) {
    sigma2 <- rep(sigma2, length(times))
  }
  
  sigma2 <- as.vector(sigma2)
  
  if (corstr == "exchangeable") {
    diag(sqrt(sigma2)) %*% cormat.exch(alpha, length(times)) %*% diag(sqrt(sigma2))
  }
}
