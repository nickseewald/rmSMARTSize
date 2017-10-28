### Functions for simulations of SMARTs with continuous repeated-measures outcomes
### Nick Seewald

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

### Custom combine function for the foreach %dopar% loop in sim()
combine.results <- function(list1, list2) {
  pval       <- c(list1$pval, list2$pval)
  param.hat  <- rbind(list1$param.hat, list2$param.hat)
  sigma2.hat <- rbind(list1$sigma2.hat, list2$sigma2.hat)
  param.var  <- Reduce("+", list(list1$param.var, list2$param.var))
  rho.hat    <- c(list1$rho.hat, list2$rho.hat)
  coverage   <- list1$coverage + list2$coverage
  iter       <- c(list1$iter, list2$iter)
  valid      <- list1$valid + list2$valid
  condVars   <- Map("+", list1$condVars, list2$condVars)
  
  return(list("pval" = pval, "param.hat" = param.hat, "sigma2.hat" = sigma2.hat, "iter" = iter,
              "param.var" = param.var, "rho.hat" = rho.hat, "valid" = valid, "coverage" = coverage,
              "condVars" = condVars))
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
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:4, function(i) V)
  } else if (length(V) != 4)
    stop("V must be a single matrix or a list of 4 matrices.")
  
  mvec <- meanvec(times, spltime, design, gammas)
  deriv <- mod.derivs(times, spltime)
  
  resids <- matrix(rep(d$Y, 4), ncol = 4) -
    matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE)
  resids <- cbind("id" = d$id, resids * d[, grep("dtr", names(d))] *
                    matrix(rep(d$weight, 4), ncol = 4))
  
  DVinv <- lapply(1:4, function(i) t(deriv[[i]]) %*% solve(V[[i]]))
  
  Reduce("+", lapply(split.data.frame(resids, resids$id), function(x) {
    Reduce("+", lapply(1:4, function(dtr) {
      DVinv[[dtr]] %*% matrix(x[[paste0("dtr", dtr)]], ncol = 1)
    }))
  })
  ) / n
}

esteqn.jacobian <- function(d, V, times, spltime) {
  
  deriv <- mod.derivs(times, spltime)
  
  Reduce("+", lapply(split.data.frame(d, d$id), function(x) {
    Reduce("+", lapply(1:4, function(dtr) {
      unique(x[[paste0("dtr", dtr)]]) * unique(x$weight) *
        t(deriv[[dtr]]) %*% solve(V) %*% deriv[[dtr]]
    }))
  })) / -length(unique(d$id))
}

estimate.params <- function(d, V, times, spltime, design, start, maxiter.solver, tol) {
  params.hat <- start
  j <- 0
  epsilon <- 10
  J <- esteqn.jacobian(d, V, times, spltime)
  
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
  J <- esteqn.jacobian(d, V, times, spltime)
  solve(-J) %*%  meat.compute(d, V, times, spltime, design, gammas) %*% solve(-J) / n
}

estimate.rho <- function(d, times, spltime, design, sigma, gammas, corstr = "exch") {
  n <- length(unique(d$id))
  mvec <- meanvec(times, spltime, design, gammas)
  Dmat <- mod.derivs(times, spltime)
  
  if (any(is(d$Y) == "NULL")) stop("d has to be in long format.")
  if (length(sigma) != length(times) & length(sigma) != 1)
    stop("sigma must be length one OR a vector with the same length as times.")
  
  # If a time-invariant sigma is provided, replicate it so there's one copy per time
  if (length(sigma) == 1)
    sigma <- rep(sigma, length(times))
  
  # Compute residuals (Y_{it} - mu_{t}(a1, a2))
  resids <- matrix(rep(d$Y, 4), ncol = 4) -
    matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE)
  
  # Create data.frame from resids, indexed by time and ID
  # Multiply residuals by DTR indicators (dtr1[i] = 1 iff i is consistent with DTR 1)
  resids <- cbind("id" = d$id, "time" = d$time,
                  resids * d[, grep("dtr", names(d))])
  
  if (corstr == "exch") {
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
    
    # Construct numerator for estimator of alpha in exchangeable corstr
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
  mvec <- meanvec(times, spltime, design, gammas)
  Dmat <- mod.derivs(times, spltime)
  
  if (any(is(d$Y) == "NULL")) stop("d has to be in long format")
  
  resids <- (matrix(rep(d$Y, 4), ncol = 4) -
               matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE))^2
  resids <- cbind("id" = d$id, "time" = d$time,
                  resids * d[, grep("dtr", names(d))] *
                    matrix(rep(d$weight, 4), ncol = 4))
  
  weightmat <- cbind("id" = d$id, "time" = d$time, d$weight * d[, grep("dtr", names(d))])
  
  sigma2.t.dtr <- Reduce(function(...) merge(..., by = "time"), 
                         lapply(1:4, function(dtr) {
                           x <- list(resids[[paste0("dtr", dtr)]])
                           sumsqrs <- aggregate(x = setNames(x, paste0("dtr", dtr)),
                                                by = list("time" = resids$time), sum)
                           x <- list(weightmat[[paste0("dtr", dtr)]])
                           sumwts <- aggregate(x = setNames(x, paste0("dtr", dtr)),
                                               by = list("time" = resids$time), sum)
                           sumsqrs[, 2] <- sumsqrs[, 2] / (sumwts[, 2] - length(gammas))
                           sumsqrs
                         }))
  
  if (pool.time & pool.dtr) {
    mean(apply(sigma2.t.dtr[, -1], 2, mean))
  } else if (pool.time & !pool.dtr) {
    matrix(apply(sigma2.t.dtr[, -1], 2, mean), nrow = 1, 
           dimnames = list(NULL, sapply(1:4, function(dtr) paste0('dtr', dtr))))
  } else if (!pool.time & pool.dtr) {
    matrix(apply(sigma2.t.dtr[, -1], 1, mean), nrow = 1,
           dimnames = list(NULL, times))
  } else {
    sigma2.t.dtr
  }
}

# Clean up output of combine.results, to be used in the foreach loop after combining everything
finalize.results <- function(x) {
  param.hat     <- apply(as.matrix(x$param.hat), 2, mean)
  sigma2.hat    <- apply(as.matrix(x$sigma2.hat), 2, mean)
  param.var     <- x$param.var / x$valid
  param.var.est <- apply(as.matrix(x$param.hat), 2, var)
  rho.hat       <- mean(x$rho.hat)
  coverage      <- as.numeric(x$coverage / x$valid)
  condVars      <- lapply(x$condVars, function(v) v / x$valid)
  
  cat("Finishing...\n")
  
  return(list("n" = x$n, "alpha" = x$alpha, "power.target" = x$power.target, "delta" = x$delta,
              "niter" = x$niter, "corstr" = x$corstr, "pval" = x$pval, "param.hat" = param.hat,
              "sigma2.hat" = sigma2.hat, "iter" = x$iter, "param.var" = param.var, 
              "rho.hat" = rho.hat, "valid" = x$valid, "coverage" = coverage, "condVars" = condVars))
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

generate.nr.cor <- function(a1, a2R, a2NR, t1, t2, spltime, design, rprob,
                            sigma.t1, sigma.t2, sigma.t1.r, sigma.t2.r,
                            rho, rho.r, gammas, lambdas) {
  if (t1 == t2) {
    rho <- rho.r <- 1
  }
  rmean.t1  <- conditional.model(a1, 1, a2R, a2NR, t1, spltime, design, rprob, gammas, lambdas)
  rmean.t2  <- conditional.model(a1, 1, a2R, a2NR, t2, spltime, design, rprob, gammas, lambdas)
  nrmean.t1 <- conditional.model(a1, 0, a2R, a2NR, t1, spltime, design, rprob, gammas, lambdas)
  nrmean.t2 <- conditional.model(a1, 0, a2R, a2NR, t2, spltime, design, rprob, gammas, lambdas)
  
  sigma.t1.nr <- generate.nr.sd(a1, a2R, a2NR, t1, spltime, design, rprob, sigma.t1, sigma.t1.r, gammas, lambdas)
  sigma.t2.nr <- generate.nr.sd(a1, a2R, a2NR, t2, spltime, design, rprob, sigma.t2, sigma.t2.r, gammas, lambdas)
  
  (-1/(sigma.t1.nr * sigma.t2.nr)) * (1 / (1 - rprob)) * 
    (rprob * rho.r * sigma.t1.r * sigma.t2.r -
       rho * sigma.t1 * sigma.t2 +
       rprob * (1 - rprob) * (rmean.t1 - nrmean.t1) * (rmean.t2 - nrmean.t2))
}

generate.cond.sd <- function(a1, r, a2R, a2NR, t, spltime, design, rprob, 
                           sigma, sigma.cond, gammas, lambdas) {
  rmean  <- conditional.model(a1, 1, a2R, a2NR, t, spltime, design, rprob, gammas, lambdas)
  nrmean <- conditional.model(a1, 0, a2R, a2NR, t, spltime, design, rprob, gammas, lambdas)
  
  if (t <= spltime) {
    sigma
  } else {
    sqrt((1 / (r * rprob + (1 - r) * (1 - rprob))) * (sigma^2 - ((1 - r) * rprob + r * (1 - rprob)) * sigma.cond^2 - 
                                                       rprob * (1 - rprob) * (rmean - nrmean)^2)
    )
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
  mvec <- meanvec(times, spltime, design, gammas)
  dmat <- mod.derivs(times, spltime)
  
  resids <- matrix(rep(d$Y, 4), ncol = 4) -
    matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE)
  resids <- cbind("id" = d$id, resids * d[, grep("dtr", names(d))] *
                    matrix(rep(d$weight, 4), ncol = 4))
  
  # sum over IDs
  Reduce("+", lapply(split.data.frame(resids, resids$id), function(obs) {
    # compute weighted residual matrix for individual
    # resid <- obs$weight * obs[, grepl("dtr", names(obs))] * 
    # (obs$Y - meanvec(times, spltime, design, gammas))
    # Sum over DTRs
    m <- Reduce("+", lapply(1:4, function(dtr) {
      t(dmat[[dtr]]) %*% solve(V) %*% as.matrix(obs[[paste0("dtr", dtr)]], ncol = 1)
    }))
    m %*% t(m)
  })) / length(unique(d$id))
}

mod.derivs <- function(times, spltime) {
  lapply(1:4, function(dtr) {
    a1 <- switch(dtr, 1, 1, -1, -1)
    a2 <- switch(dtr, 1, -1, 1, -1)
    do.call(rbind, lapply(times, function(t) {
      matrix(c(1, (t <= spltime) * matrix(c(t, t * a1, 0, 0, 0, 0), nrow = 1) + 
                 (t > spltime) * matrix(c(spltime, spltime * a1, (t - spltime),
                                          (t - spltime) * a1, (t - spltime) * a2,
                                          (t - spltime) * a1 * a2), nrow = 1)),
             nrow = 1)
    }))
  })
}

print.simResult <- function(sim) {
  test <- binom.test(x = sum(sim$pval <= sim$alpha / 2) + sum(sim$pval >= 1 - sim$alpha / 2),
                     n = sim$valid, p = sim$power.target, alternative = "two.sided")
  cat("\tSimulation Results\n")
  cat(paste("input summary:\nnumber of participants =", sim$n, "number of trials = ", sim$niter, "number of valid trials = ", sim$valid, "\n"))
  cat("power results:\n")
  fp <- format.pval(test$p.value, digits = max(1L, getOption('digits') - 3L))
  cat(paste("target power =", sim$power.target, "estimated power =", test$estimate, "p-value",
            ifelse(substr(fp, 1L, 1L) == "<", fp, paste("=", fp)), "\n"))
  cat(paste("parameter estimates:\nmodel parameters =", paste0("(", paste(round(sim$param.hat, 5), collapse = ", "), ")"),
            "coverage =", sim$coverage, "\n"))
  cat(paste("marginal estimates:\nvariances =", paste0("(", paste(round(sim$sigma2.hat, 5), collapse = ", "), ")"),
            "correlation =", round(sim$rho.hat, 5), "\n"))
}

sample.size <- function(delta, r, rho, alpha = 0.05, power = .8, design = 2, round = "up", conservative = TRUE) {
  # Input checks
  if (!(round %in% c("up", "down"))) stop("round must be either up or down")
  if (design == 2) {
    if (conservative) {
      ceiling(2 * (qnorm(1 - alpha / 2) + qnorm(power))^2 * 2 * (2 - r) * (1 - rho^2) / delta^2)
    } else {
      nid <- (2 * (qnorm(1 - alpha / 2) + qnorm(power))^2 * 2 * (2 - r)) / delta^2
      correction <- ((1 - rho) * (2 * rho + 1)) / (1 + rho) + (1 - rho)/(1 + rho) * rho^2 / (2 - r)
      if (round == "up") {
        ceiling(nid * correction)
      } else if (round == "down") {
        floor(nid * correction)
      }
    }
  } else stop("Not a valid design indicator.")
}

sim <- function(n = NULL, gammas, lambdas, times, spltime,
                alpha = .05, power = .8, delta, design = 2, round = "up", conservative = TRUE,
                r = NULL, r1 = r, r0 = r,
                uneqvarDTR = NULL, sigma, sigma.r1 = sigma, sigma.r0 = sigma,
                corstr = "identity",
                rho = NULL, rho.r1 = rho, rho.r0 = rho,
                L.eos = NULL, L.auc = NULL,
                constant.var.time = TRUE, constant.var.dtr = TRUE,
                uneqVarDTR = NULL, perfect = FALSE,
                niter = 5000, tol = 1e-8, maxiter.solver = 1000,
                notify = FALSE, pbDevice = NULL, pbIdentifier = NULL) {
  
  # n:        number of participants in trial
  # gammas:   Parameters from marginal mean model
  # lambdas:  Parameters for response "offset" in conditional model
  # r:        Probability of response to first-stage treatment (assuming equal for both txts)
  # r1:       Probability of response to A1 = 1  (assuming not equal to r0; r must be NULL)
  # r0:       Probability of response to A1 = -1 (assuming not equal to r1; r must be NULL)
  # times:    Vector of times at which measurements are collected (currently limited to length three)
  # spltime:  Time (contained in times) at which re-randomization occurs
  # uneqvarDTR: list of vectors of (a1, a2R, a2NR) for DTR(s) which do not have the same variance as the others
  # sigma:    Marginal variance of Y (assumed constant over time and DTR)
  # sigma.r1: Conditional variance of Y for responders to treatment A1 = 1
  # sigma.r0: Conditional variance of Y for responders to treatment A1 = -1
  # corstr:   Character string, one of "identity", "exch"/"exchangeable", "unstr"/"unstructured", or "custom".
  #            This is the TRUE form of the within-person correlation structure
  #           Form of V to use in estimating equations.
  # rho:      If corstr is either "exch" or "exchangeable", the within-person correlation used to generate data 
  #            (assumed constant across DTRs)
  
  # if (!is.null(r) & (r < 0 | r > 1)) stop("r must be between 0 and 1.")
  if (is.null(r) & is.null(r1) & is.null(r0)) stop("You must provide either r or both r1 and r0.")
  ## TODO: Finish input handling
  
  # Handle uneqvarDTR
  if (design == 1 & is.null(uneqvarDTR)) 
    stop("For design I, you must provide a list of vectors of (a1, a2R, a2NR) for DTRs which do not have the same variance as others.")
  if (design == 1 & !is.list(uneqvarDTR))
    stop("uneqvarDTR must be a list.")
  
  ## Handle correlation structure
  if (corstr %in% c("id", "identity")) rho <- rho.r1 <- rho.r0 <- 0
  
  ## If design == 1, constant.var.dtr is impossible to satisfy in a generative model. Set it to false.'
  if (design == 1) constant.var.dtr <- FALSE
  
  # If n is not provided, compute it from the other inputs
  if (is.null(n)) 
    n <- sample.size(delta, ifelse(is.null(r), min(r0, r1), r), rho, alpha, power, design, round, conservative)
  
  ## Simulate
  cat(paste0("********************\nStarting simulation!\n", ifelse(is.null(pbIdentifier), NULL, pbIdentifier), "\nn = ", n, "\n",
             niter, " iterations.\n********************\n"))
  results <- foreach(i = 1:niter, .combine = combine.results, .final = finalize.results,
                     .verbose = FALSE, .errorhandling = "stop", .multicombine = FALSE) %dorng% { 
                       
                       d <- SMART.generate(n, times, spltime, r1, r0, gammas, lambdas, design = design, sigma, sigma.r1, sigma.r0, corstr = corstr,
                                           rho, rho.r1, rho.r0, perfect = perfect)
                       if (d$valid == FALSE) {
                         return(list("pval" = NA, "param.hat" = rep(NA, length(gammas)), 
                              "param.var" = matrix(0, ncol = length(gammas), nrow = length(gammas)),
                              "sigma2.hat" = matrix(ncol = (length(times) * (1 - constant.var.time) * constant.var.dtr) +
                                                      4 * (1 - constant.var.dtr) +
                                                      (constant.var.time * constant.var.dtr),
                                                    nrow = (length(times) * (1 - constant.var.dtr) +
                                                              (4 * (1 - constant.var.dtr) * constant.var.time) +
                                                              (constant.var.time * constant.var.dtr))),
                              "rho.hat" = NA, "valid" = 0, "coverage" = 0, "iter" = NULL,
                              "condVars" = lapply(1:length(dtrIndex(design)$a1), function(x) matrix(0, ncol = length(times), nrow = length(times))),
                              "alpha" = alpha, "power.target" = power))
                       } else {
                         d1 <- reshape(d$data, varying = list(grep("Y", names(d$data))), ids = d$data$id, 
                                       times = times, direction = "long", v.names = "Y")
                         d1 <- d1[order(d1$id, d1$time), ]
                         
                         param.hat <- estimate.params(d1, diag(rep(1, length(times))), times, spltime, design, rep(0, length(gammas)), maxiter.solver, tol)
                         
                         sigma2.hat <- estimate.sigma2(d1, times, spltime, design, param.hat,
                                                       pool.time = constant.var.time, pool.dtr = constant.var.dtr)
                         
                         rho.hat <- estimate.rho(d1, times, spltime, design, sqrt(sigma2.hat), param.hat)
                         
                         # Compute variance matrices for all conditional cells
                         condVars <- lapply(split.SMART(d$data), function(x) {
                           var(subset(x, select = grep("Y", names(x), value = TRUE)))
                         })
                         
                         # Iterate parameter estimation
                         if (corstr %in% c("id", "identity") | (corstr == "exch" & rho == 0)) {
                           outcome.var <- varmat(sigma2.hat, 0, times, "exch")
                           param.var <- estimate.paramvar(d1, diag(rep(1, length(times))), times, spltime, design, gammas = param.hat)
                           iter <- 1
                         } else if (corstr == "exch" & rho != 0) {
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
                         
                         confLB <- L.eos %*% param.hat - sqrt(L.eos %*% param.var %*% L.eos) * qnorm(.975)
                         confUB <- L.eos %*% param.hat + sqrt(L.eos %*% param.var %*% L.eos) * qnorm(.975)
                         
                         coverage <- ifelse(confLB <= L.eos %*% gammas & L.eos %*% gammas <= confUB, 1, 0)
                         
                         pval <- 1 - pnorm(as.numeric((L.eos %*% param.hat) / sqrt(L.eos %*% param.var %*% L.eos)))
                         
                         list("pval" = pval, "param.hat" = t(param.hat), "param.var" = param.var,
                              "sigma2.hat" = sigma2.hat, "rho.hat" = rho.hat, "valid" = 1, "coverage" = coverage,
                              "iter" = iter, "condVars" = condVars)
                       }
                     }
  pbMessageText <- paste0(ifelse(!is.null(pbIdentifier), paste0(pbIdentifier, "\n"), ""), "n = ", n,
                          "\neffect size = ", delta, "\npower = ", power,
                          "\ntrue corr structure = ", 
                          switch(corstr, "id" =, "identity" = "identity", 
                                 "exch" =, "exchangeable" = paste0("exchangeable(", rho, ")")))
  if (notify) {
    pbPost("note", "Simulation Complete!", body = pbMessageText, recipients = pbDevice)
  }
  
  results <- c(list("n" = n, "alpha" = alpha, "power.target" = power, "delta" = delta,
                    "niter" = niter, "corstr" = corstr), results)
  
  class(results) <- c("simResult", class(results))
  return(results)
}


SMART.generate <- function(n, times, spltime, r1, r0, gammas, lambdas, design,
                           sigma, sigma.r1, sigma.r0, corstr = "identity", uneqvarDTR = NULL,
                           rho = NULL, rho.r1 = rho, rho.r0 = rho, perfect = FALSE) {
  
  # n:        number of participants in trial
  # gammas:   Parameters from marginal mean model
  # mean.r1:  Vector of means for responders to A1 = 1 in the second stage (IF DESIGN = 1, mean for A2R = 1)
  # mean.r0:  Vector of means for responders to A1 = -1 in the second stage (IF DESIGN = 1, mean for A2R = 1)
  # r:        Probability of response to first-stage treatment (assuming equal for both txts)
  # r1:       Probability of response to A1 = 1  (assuming not equal to r0; r must be NULL)
  # r0:       Probability of response to A1 = -1 (assuming not equal to r1; r must be NULL)
  # times:    Vector of times at which measurements are collected (currently limited to length three)
  # spltime:  Time (contained in times) at which re-randomization occurs
  # sigma:    Marginal variance of Y (assumed constant over time and DTR)
  # sigma.r1: Vector of conditional variances of Y for responders to treatment A1 = 1 in the second stage 
  # sigma.r0: Vector of conditional variances of Y for responders to treatment A1 = -1 in the second stage
  # corstr:   True correlation structure
  # rho:      If corstr is either "exch" or "exchangeable", the within-person correlation used to compute V 
  #           (assumed constant across DTRs)
  # rho.r1:   Vector of conditional correlations between Ys for responders to treatment A1 = 1 in the second stage
  # rho.r0:   Exchangeable correlation for responders to A1 == -1
  
  ## Input Checks
  if (sum(times > spltime) == 0) stop("Currently spltime must be less than max(times)")
  
  ## Handle correlation structure
  if (corstr %in% c("id", "identity")) {
    rho <- rho.r1 <- rho.r0 <- 0
    corstr <- "id"
  }
  if (corstr %in% c("exch", "exchangeable")) corstr <- "exch"
  
  ## Make sure design input is valid
  if (!(design %in% 1:3)) stop("Invalid design choice. Must be one of 1-3.")
  
  ## If sigma is not a vector of length length(times), make it that way (or throw an error)
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
  
  ## Generate treatment allocations and response
  d <- data.frame("id"   = 1:n,
                  "A1"   = 2 * rbinom(n, 1, .5) - 1,
                  "R"    = NA,
                  "A2R"  = 0,
                  "A2NR" = 0)
  
  d$R[d$A1 ==  1] <- rbinom(sum(d$A1 == 1), 1, r1)
  d$R[d$A1 == -1] <- rbinom(sum(d$A1 == -1), 1, r0)
  
  if (design == 1) {
    d$A2R[d$A1 == 1  & d$R == 1] <- 2 * rbinom(sum(d$A1 == 1  & d$R == 1), 1, .5) - 1
    d$A2R[d$A1 == -1 & d$R == 1] <- 2 * rbinom(sum(d$A1 == -1 & d$R == 1), 1, .5) - 1
    d$A2NR[d$A1 == 1  & d$R == 0] <- 2 * rbinom(sum(d$A1 == 1  & d$R == 0), 1, .5) - 1
    d$A2NR[d$A1 == -1 & d$R == 0] <- 2 * rbinom(sum(d$A1 == -1 & d$R == 0), 1, .5) - 1
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
    d$dtr1 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR == -1))))
    d$dtr1 <- as.numeric(with(d, (A1 == -1)))
  }
  
  ## Check validity of treatment allocations
  if (!validTrial(d, design)) {
    return(list("data" = d, "valid" = F))
  }
    
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
  
  ## Compute conditional standard deviations for non-responders
  if (design == 1) {
    dtrs <- do.call("rbind", dtrIndex(design))
    uneqvarDTRindex <- sapply(1:2, function(i) 
      which((dtrs[1, ] == uneqvarDTR[[i]][1]) + 
              (dtrs[2, ] == uneqvarDTR[[i]][2]) +
              (dtrs[3, ] == uneqvarDTR[[i]][3]) == 3)
      )
    uneqvarDTR.sd <- lapply(uneqvarDTR, function(dtr) {
      sapply(times, function(t) {
        generate.sd(sigma[which(times == t)], a1 = dtr[1], a2R = dtr[2], a2NR = dtr[3], t = t, spltime = spltime, 
                    design, rprob = get(paste0("r", (dtr[1] + 1) / 2)), gammas, lambdas)
      })
      })
    sigma.r11 <- sigma.r1
    sigma.r01 <- sigma.r0
    
    #Loop over indices for first-stage treatment
    invisible(sapply(c(0, 1), function(a1ind) {
      # Use DTR (a1ind, 1, 1) to 
      assign(paste0("sigma.nr", a1ind, "1"),
             sapply(1:length(times), function(i) {
               generate.cond.sd(a1 = 2 * a1ind - 1, r = 0, a2R = 1, a2NR = 1, t = times[i], spltime = spltime,
                                design = 1, rprob = get(paste0("r", a1ind)),
                                sigma = ifelse(sum(sapply(uneqvarDTR, function(dtr) sum(dtr == c(2 * a1ind - 1, 1, 1)) == 3)) > 0,
                                               uneqvarDTR.sd[[which(sapply(uneqvarDTR, function(dtr) 
                                                 sum(dtr == c(2 * a1ind -1, 1, 1)) == 3))]][i],
                                               sigma[i]),
                                sigma.cond = get(paste0("sigma.r", a1ind, "1"))[i], gammas, lambdas)
               }), envir = .GlobalEnv)
      assign(paste0("sigma.r", a1ind, "0"),
             sapply(1:length(times), function(i) {
               generate.cond.sd(a1 = 2 * a1ind - 1, r = 1, a2R = -1, a2NR = 1, t = times[i], spltime = spltime,
                                design = 1, rprob = get(paste0("r", a1ind)),
                                sigma = ifelse(sum(sapply(uneqvarDTR, function(dtr) sum(dtr == c(2 * a1ind - 1, -1, 1)) == 3)) > 0,
                                               uneqvarDTR.sd[[which(sapply(uneqvarDTR, function(dtr) 
                                                 sum(dtr == c(2 * a1ind -1, -1, 1)) == 3))]][i],
                                               sigma[i]),
                                sigma.cond = get(paste0("sigma.nr", a1ind, "1"))[i], gammas, lambdas)
             }), envir = .GlobalEnv)
      assign(paste0("sigma.nr", a1ind, "0"),
             sapply(1:length(times), function(i) {
               generate.cond.sd(a1 = 2 * a1ind - 1, r = 0, a2R = -1, a2NR = -1, t = times[i], spltime = spltime,
                                design = 1, rprob = get(paste0("r", a1ind)),
                                sigma = ifelse(sum(sapply(uneqvarDTR, function(dtr) sum(dtr == c(2 * a1ind - 1, -1, -1)) == 3)) > 0,
                                               uneqvarDTR.sd[[which(sapply(uneqvarDTR, function(dtr) 
                                                 sum(dtr == c(2 * a1ind -1, -1, -1)) == 3))]][i],
                                               sigma[i]),
                                sigma.cond = get(paste0("sigma.r", a1ind, "0"))[i], gammas, lambdas)
             }), envir = .GlobalEnv)
    }))
    
      invisible(apply(dtrs, 2, function(dtr) {
      a1ind   <- (dtr[1] + 1) / 2
      a2Rind  <- (dtr[2] + 1) / 2 
      a2NRind <- (dtr[3] + 1) / 2
      dtrNum <- which((dtrs[1, ] == dtr[1]) + 
                          (dtrs[2, ] == dtr[2]) +
                          (dtrs[3, ] == dtr[3]) == 3)
      sapply(c(0, 1), function(r) {
        assign(paste0("sigma.", switch(r + 1, "nr", "r"), a1ind, ifelse(r == 1, a2Rind, a2NRind)),
               unname(sapply(times, function(t) {
                 generate.cond.sd(a1 = dtr[1], r = r, a2R = dtr[2], a2NR = dtr[3], t = t, spltime, design = 1, rprob = rprob,
                                  sigma = ifelse(dtrNum %in% uneqvarDTRindex, 
                                                 uneqvarDTR.sd[[which(dtrNum == uneqvarDTRindex)]][which(times == t)], 
                                                 sigma[which(times == t)]), 
                                  sigma.cond = get(paste0("sigma.", ifelse(r == 1, "nr", "r"), a1ind, ifelse(r == 1, a2NRind, a2Rind))),
                                  gammas = gammas, lambdas = lambdas)
               })),
               envir = .GlobalEnv)
      })
    }))
  } else if (design == 2) {
    invisible(apply(unique(subset(d, R == 0, select = c("A1", "A2R", "A2NR"))), 1, function(x) {
      assign(paste0("sigma.nr", (x[1] + 1) / 2, (x[3] + 1) / 2),
             unname(sapply(1:length(times), function(i) {
               generate.cond.sd(a1 = x[1], r = 0, a2R = x[2], a2NR = x[3], t = times[i], spltime, design = 2, 
                              rprob = get(paste0("r", (x[1] + 1) / 2)),
                              sigma = sigma[i],
                              sigma.cond = get(paste0("sigma.r", (x[1] + 1) / 2))[i],
                              gammas = gammas, lambdas = lambdas)
             })),
             envir = .GlobalEnv)
    }))
  }
  
  if (corstr == "id") {
    cormat.r1 <- cormat.r0 <- diag(rep(1, length(times)))
  } else if (corstr == "exch") {
    cormat.r1 <- cormat.r0 <- cormat.exch(rho, length(times))
  } else if (corstr == "ar1") {
    cormat.r1 <- cormat.r0 <- cormat.ar1(rho, length(times))
  }
  
  m <- expand.grid(which(times > spltime), which(times <= spltime))
  m <- rbind(as.matrix(m), as.matrix(m[, 2:1]))
  m <- rbind(m, expand.grid(which(times > spltime), which(times > spltime)))
  m <- subset(m, Var1 != Var2)
  for (i in 1:dim(m)[1]) {
    if (min(m[i, ]) > which(times == spltime)) {
      cormat.r1[m[i, 1], m[i, 2]] <- rho.r1
      cormat.r0[m[i, 1], m[i, 2]] <- rho.r0
    } else {
      cormat.r1[m[i, 1], m[i, 2]] <- cormat.r1[m[i, 1], m[i, 2]] * sigma[max(m[i, ])] / sigma.r1[max(m[i, ])]
      cormat.r0[m[i, 1], m[i, 2]] <- cormat.r0[m[i, 1], m[i, 2]] * sigma[max(m[i, ])] / sigma.r0[max(m[i, ])]
    }
  }
  
  d <-
    do.call("rbind", lapply(split.SMART(d), function(x) {
      a1ind   <- as.character((unique(x$A1) + 1) / 2)
      a2Rind  <- as.character((unique(x$A2R) + 1) / 2)
      a2NRind <- as.character((unique(x$A2NR) + 1) / 2)
      if (unique(x$R) == 1) {
        x[, grepl("Y", names(x))] <- x[, grepl("Y", names(x))] +
          mvrnorm(
            dim(x)[1],
            mu = rep(0, length(times)),
            Sigma = diag(get(paste0("sigma.r", a1ind))) %*%
              get(paste0("cormat.r", (unique(x$A1 + 1) / 2))) %*%
              diag(get(paste0("sigma.r", a1ind)))
          )
      } else {
        cormat.nr <- cormat.exch(rho, length(times))
        for (i in 1:dim(m)[1]) {
          cormat.nr[m[i, 1], m[i, 2]] <- 
            generate.nr.cor(a1 = unique(x$A1), a2R = unique(x$A2R), a2NR = unique(x$A2NR),
                            t1 = times[m[i, 1]], t2 = times[m[i, 2]], spltime = spltime,
                            design = design, rprob = get(paste0("r", a1ind)),
                            sigma.t1 = sigma[m[i, 1]], sigma.t2 = sigma[m[i, 2]], 
                            sigma.t1.r = get(paste0("sigma.r", a1ind))[m[i, 1]],
                            sigma.t2.r = get(paste0('sigma.r', a1ind))[m[i, 2]],
                            rho = cormat.nr[m[i, 1], m[i, 2]],
                            rho.r = get(paste0("cormat.r", a1ind))[m[i, 1], m[i, 2]],
                            gammas = gammas, lambdas = lambdas)
        }
        x[, grepl("Y", names(x))] <- x[, grepl("Y", names(x))] +
          mvrnorm(
            dim(x)[1],
            mu = rep(0, length(times)),
            Sigma = diag(get(paste0("sigma.nr", a1ind, a2Rind, a2NRind))) %*%
              get(paste0("cormat.", switch(unique(x$R) + 1, "nr", paste0("r", a1ind)))) %*%
              diag(get(paste0("sigma.nr", a1ind, a2NRind)))
          )
      }
      x
    }))
  
  # ## Add error to Y's
  # d[d$A1 == 1 & d$R == 1, grepl("Y", names(d))] <-
  #   d[d$A1 == 1 & d$R == 1, grepl("Y", names(d))] +
  #   mvrnorm(sum(d$A1 == 1 & d$R == 1), mu = rep(0, 3),
  #           Sigma = diag(c(sigma, sigma, sigma.r1)) %*%
  #             cormat.r1 %*% 
  #             diag(c(sigma, sigma, sigma.r1))
  #           )
  # 
  # d[d$A1 == 1 & d$R == 0 & d$A2 == 1, grepl("Y", names(d))] <-
  #   d[d$A1 == 1 & d$R == 0 & d$A2 == 1, grepl("Y", names(d))] +
  #   mvrnorm(sum(d$A1 == 1 & d$R == 0 & d$A2 == 1), mu = rep(0, 3), 
  #           # Sigma = diag(c(sigma, sigma, sigma.r1)) %*%
  #           # cormat.exch(rho.r1, length(times)) %*%
  #           # diag(c(sigma, sigma, sigma.r1))
  #           Sigma = diag(c(sigma, sigma, sigma.nr11)) %*%
  #             sapply(1:length(times), function(i) {
  #               sapply(1:length(times), function(j) {
  #                 generate.nr.cor(r1, 1, 1, times[i], times[j], spltime, sigma, sigma, sigma.r1, sigma.r1, rho, rho.r1, lambdas, gammas)
  #               })
  #             }) %*%
  #             diag(c(sigma, sigma, sigma.nr11))
  #           )
  # 
  # d[d$A1 == 1 & d$R == 0 & d$A2 == -1, grepl("Y", names(d))] <-
  #   d[d$A1 == 1 & d$R == 0 & d$A2 == -1, grepl("Y", names(d))] +
  #   mvrnorm(sum(d$A1 == 1 & d$R == 0 & d$A2 == -1), mu = rep(0, 3), 
  #           # Sigma = diag(c(sigma, sigma, sigma.r1)) %*%
  #           # cormat.exch(rho.r1, length(times)) %*%
  #           # diag(c(sigma, sigma, sigma.r1))
  #           Sigma = diag(c(sigma, sigma, sigma.nr10)) %*%
  #             sapply(1:length(times), function(i) {
  #               sapply(1:length(times), function(j) {
  #                 generate.nr.cor(r1, 1, -1, times[i], times[j], spltime, sigma, sigma, sigma.r1, sigma.r1, rho, rho.r1, lambdas, gammas)
  #               })
  #             }) %*%
  #             diag(c(sigma, sigma, sigma.nr10))
  #           )
  # 
  # d[d$A1 == -1 & d$R == 1, grepl("Y", names(d))] <-
  #   d[d$A1 == -1 & d$R == 1, grepl("Y", names(d))] +
  #   mvrnorm(sum(d$A1 == -1 & d$R == 1), mu = rep(0, 3), 
  #           # Sigma = diag(c(sigma, sigma, sigma.r1)) %*%
  #           #   cormat.exch(rho.r1, length(times)) %*%
  #           #   diag(c(sigma, sigma, sigma.r1))
  #           Sigma = diag(c(sigma, sigma, sigma.r0)) %*%
  #             cormat.exch(rho.r0, 3) %*%
  #             diag(c(sigma, sigma, sigma.r0))
  #           )
  # 
  # d[d$A1 == -1 & d$R == 0 & d$A2 == 1, grepl("Y", names(d))] <-
  #   d[d$A1 == -1 & d$R == 0 & d$A2 == 1, grepl("Y", names(d))] +
  #   mvrnorm(sum(d$A1 == -1 & d$R == 0 & d$A2 == 1), mu = rep(0, 3), 
  #           # Sigma = diag(c(sigma, sigma, sigma.r1)) %*%
  #           # cormat.exch(rho.r1, length(times)) %*%
  #           # diag(c(sigma, sigma, sigma.r1))
  #           Sigma = diag(c(sigma, sigma, sigma.nr01)) %*%
  #             sapply(1:length(times), function(i) {
  #               sapply(1:length(times), function(j) {
  #                 generate.nr.cor(r0, -1, 1, times[i], times[j], spltime, sigma, sigma, sigma.r0, sigma.r0, rho, rho.r0, lambdas, gammas)
  #               })
  #             }) %*%
  #             diag(c(sigma, sigma, sigma.nr01))
  #           )
  # 
  # d[d$A1 == -1 & d$R == 0 & d$A2 == -1, grepl("Y", names(d))] <-
  #   d[d$A1 == -1 & d$R == 0 & d$A2 == -1, grepl("Y", names(d))] +
  #   mvrnorm(sum(d$A1 == -1 & d$R == 0 & d$A2 == -1), mu = rep(0, 3), 
  #           # Sigma = diag(c(sigma, sigma, sigma.r1)) %*%
  #           # cormat.exch(rho.r1, length(times)) %*%
  #           # diag(c(sigma, sigma, sigma.r1))
  #           Sigma = diag(c(sigma, sigma, sigma.nr00)) %*%
  #             sapply(1:length(times), function(i) {
  #               sapply(1:length(times), function(j) {
  #                 generate.nr.cor(r0, -1, -1, times[i], times[j], spltime, sigma, sigma, sigma.r0, sigma.r0, rho, rho.r0, lambdas, gammas)
  #               })
  #             }) %*%
  #             diag(c(sigma, sigma, sigma.nr00))
  #           )
  d <- d[order(d$id), ]
  rownames(d) <- 1:n
  return(list("data" = d, "valid" = T))
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
varmat <- function(sigma2, alpha, times, corstr) {
  if (length(sigma2) != 1 & length(sigma2) != length(times)) {
    stop("sigma must be either length 1 (assuming constant variance over time) or have the same length as times.")
  } else if (length(sigma2) == 1) {
    sigma2 <- rep(sigma2, length(times))
  }
  
  sigma2 <- as.vector(sigma2)
  
  if (corstr == "exch") {
    diag(sqrt(sigma2)) %*% cormat.exch(alpha, length(times)) %*% diag(sqrt(sigma2))
  }
}
