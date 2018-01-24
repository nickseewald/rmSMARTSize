### Functions for simulations of SMARTs with continuous repeated-measures outcomes
### Nick Seewald

library(xtable)

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
  mvec <- meanvec(times, spltime, design, gammas)
  Dmat <- mod.derivs(times, spltime)
  
  if (any(is(d$Y) == "NULL")) stop("d has to be in long format")
  
  resids <- (matrix(rep(d$Y, 4), ncol = 4) -
               matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE))^2
  resids <- cbind("id" = d$id, "time" = d$time,
                  resids * d[, grep("dtr", names(d))] *
                    matrix(rep(d$weight, 4), ncol = 4))
  
  weightmat <- cbind("id" = d$id, "time" = d$time, d$weight * d[, grep("dtr", names(d))])
  
  if (pool.time & pool.dtr) {
    numerator <- apply(subset(resids, select = grep("dtr", names(resids))), 2, sum)
    denominator <- apply(subset(weightmat, select = grep("dtr", names(weightmat))), 2, sum) - length(gamm)
    numerator / denominator
  }
  
  sum(subset(resids, select = grep("dtr", names(resids)))) / (sum(subset(weightmat, select = grep("dtr", names(weightmat)))) - length(gammas))
  
  sigma2.t.dtr <- Reduce(function(...) merge(..., by = "time"), 
                         lapply(1:4, function(dtr) {
                           x <- list(resids[[paste0("dtr", dtr)]])
                           sumsqrs <- aggregate(x = setNames(x, paste0("dtr", dtr)),
                                                by = list("time" = resids$time), sum)
                           x <- list(weightmat[[paste0("dtr", dtr)]])
                           sumwts <- aggregate(x = setNames(x, paste0("dtr", dtr)),
                                               by = list("time" = resids$time), sum)
                           sumsqrs[, 2] <- sumsqrs[, 2] / sumwts[, 2]
                           sumsqrs
                         }))
  
  if (pool.time & pool.dtr) {
    sum(sigma2.t.dtr) / (length(times) * sum(grepl("dtr", names(sigma2.t.dtr))))
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
               "pval.coverage" = binom.test(x = l$coverage * l$valid,
                                            n = l$valid, p = .95, 
                                            alternative = "two.sided")$p.value,
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
                uneqsdDTR = NULL, uneqsd = NULL, 
                sigma, sigma.r1 = sigma, sigma.r0 = sigma,
                corstr = "identity",
                rho = NULL, rho.r1 = rho, rho.r0 = rho,
                L.eos = NULL, L.auc = NULL,
                constant.var.time = TRUE, constant.var.dtr = TRUE, perfect = FALSE,
                niter = 5000, tol = 1e-8, maxiter.solver = 1000,
                save.data = FALSE,
                notify = FALSE, pbDevice = NULL, postIdentifier = NULL) {
  
  # n:        number of participants in trial
  # gammas:   Parameters from marginal mean model
  # lambdas:  Parameters for response "offset" in conditional model
  # r:        Probability of response to first-stage treatment (assuming equal for both txts)
  # r1:       Probability of response to A1 = 1  (assuming not equal to r0; r must be NULL)
  # r0:       Probability of response to A1 = -1 (assuming not equal to r1; r must be NULL)
  # times:    Vector of times at which measurements are collected (currently limited to length three)
  # spltime:  Time (contained in times) at which re-randomization occurs
  # uneqsdDTR: list of vectors of (a1, a2R, a2NR) for DTR(s) which do not have the same variance as the others
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
  
  # Handle uneqsdDTR
  if (design == 1 & is.null(uneqsdDTR)) {
    if (length(unique(marginal.model(dtrIndex(1)$a1, dtrIndex(1)$a2R, dtrIndex(1)$a2NR, 2, spltime, 1, gammas))) == 8)
      stop(paste("For design I, you must provide a list of vectors of (a1, a2R, a2NR) for DTRs which do not have the same variance as others.\n",
               "Alternatively, you may specify gammas such that there are only 4 unique DTR means (2 per first-stage treatment)."))
  }
  if (design == 1 & !is.null(uneqsdDTR) & !is.list(uneqsdDTR))
    stop("uneqsdDTR must be either NULL or a list.")
  
  ## Handle correlation structure
  if (corstr %in% c("id", "identity")) rho <- rho.r1 <- rho.r0 <- 0
  
  ## If design == 1, constant.var.dtr is impossible to satisfy in a generative model. Set it to false.'
  if (design == 1) constant.var.dtr <- FALSE
  
  # If n is not provided, compute it from the other inputs
  if (is.null(n)) 
    n <- sample.size(delta, ifelse(is.null(r), min(r0, r1), r), rho, alpha, power, design, round, conservative)
  
  ## Simulate
  cat(paste0("********************\n",
             "Starting simulation!\n",
             "delta = ", delta, "\n",
             "corstr = exch(", rho, ")\n",
             ifelse(is.null(postIdentifier), NULL, postIdentifier),
             "\nn = ", n, "\n",
             niter, " iterations.\n",
             "********************\n"))
  results <- foreach(i = 1:niter, .combine = combine.results, .final = finalize.results,
                     .verbose = FALSE, .errorhandling = "stop", .multicombine = FALSE, .inorder = FALSE) %dorng% { 
                       
                       d <- SMART.generate(n, times, spltime, r1, r0, gammas, lambdas, design = design, sigma, sigma.r1, sigma.r0, corstr = corstr,
                                           rho, rho.r1, rho.r0, uneqsd = NULL, uneqsdDTR = NULL)
                       if (d$valid == FALSE) {
                         result <- list("pval" = NA, "param.hat" = rep(NA, length(gammas)), 
                                        "param.var" = matrix(0, ncol = length(gammas), nrow = length(gammas)),
                                        "sigma2.hat" = matrix(ncol = (length(times) * (1 - constant.var.time) * constant.var.dtr) +
                                                                4 * (1 - constant.var.dtr) +
                                                                (constant.var.time * constant.var.dtr),
                                                              nrow = (length(times) * (1 - constant.var.dtr) +
                                                                        (4 * (1 - constant.var.dtr) * constant.var.time) +
                                                                        (constant.var.time * constant.var.dtr))),
                                        "rho.hat" = NA, "valid" = 0, "coverage" = 0, "iter" = NULL,
                                        "condVars" = lapply(1:length(dtrIndex(design)$a1), function(x) matrix(0, ncol = length(times), nrow = length(times))),
                                        "alpha" = alpha, "power.target" = power)
                         # if (save.data) {
                         #   result[["data"]] <- list(d$data)
                         # }
                         return(result)
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
                         
                         result <- list("pval" = pval, "param.hat" = t(param.hat), "param.var" = param.var,
                                        "sigma2.hat" = sigma2.hat, "rho.hat" = rho.hat, "valid" = 1, "coverage" = coverage,
                                        "iter" = iter, "condVars" = condVars)
                         if (save.data) {
                           result[["data"]] <- list(d$data)
                         }
                         
                         return(result)
                       }
                     }
  
  results <- c(list("n" = n, "alpha" = alpha, "power.target" = power, "delta" = delta,
                    "r0" = r0, "r1" = r1, "niter" = niter, "corstr" = corstr, "sharp" = !conservative), results)
  
  test <- binom.test(x = sum(results$pval <= results$alpha / 2, na.rm = T) + sum(results$pval >= 1 - results$alpha / 2, na.rm = T),
                     n = results$valid, p = results$power.target, alternative = "two.sided")

  results <- c(results, "power" = unname(test$estimate), "pval.power" = unname(test$p.value))

  postMessageText <- paste0(ifelse(!is.null(postIdentifier), paste0(postIdentifier, "\n"), ""),
                          "\ntrue corr structure = ",
                          switch(corstr, "id" =, "identity" = "identity",
                                 "exch" =, "exchangeable" = paste0("exchangeable(", rho, ")")),
                          "\nr0 = ", r0, ", r1 = ", r1,
                          "\neffect size = ", delta,
                          "\nn = ", n,
                           "\ntarget power = ", power,
                          "\nempirical power = ", results$power,
                          " (p = ", round(results$pval.power, 3), ")"
                          )
  
  class(results) <- c("simResult", class(results))
  
  if (notify) {
    try(slackr_bot(postMessageText))
    try(slackr_bot(print(results)))
    # pbPost("note", "Simulation Complete!", body = postMessageText, recipients = pbDevice)
  }
  
  return(results)
}

## Function to generate SMART data
SMART.generate <- function(n, times, spltime, r1, r0, gammas, lambdas, design,
                           sigma, sigma.r1, sigma.r0, corstr = "identity", uneqsdDTR = NULL, uneqsd = NULL,
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
  
  # uneqvar:  Vector of same length as uneqsdDTR. If not all DTRs have the same variance (i.e., sigma^2),
  #           this is the sigma
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
    sigma.r11 <- sigma.r1
    sigma.r01 <- sigma.r0
    # rm("sigma.r1", "sigma.r0")
    
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
                                sigma.cond = get(paste0("sigma.r", a1ind, "1"))[i], gammas, lambdas)
               }), envir = .GlobalEnv)
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
                                sigma.cond = get(paste0("sigma.nr", a1ind, "1"))[i], gammas, lambdas)
             }), envir = .GlobalEnv)
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
                                sigma.cond = get(paste0("sigma.r", a1ind, "0"))[i], gammas, lambdas)
             }), envir = .GlobalEnv)
    }))
    rm(sigma.r1, sigma.r0)
  } else if (design == 2) { ## Generate conditional standard deviations in Design 2
    if ((!is.null(uneqsdDTR) & is.null(uneqsd)) | (is.null(uneqsdDTR) & !is.null(uneqsd)))
      stop("For design 2, you must provide either both uneqsdDTR and uneqsd or neither.")
    sigma.r10 <- sigma.r1
    sigma.r00 <- sigma.r0
    # rm("sigma.r1", "sigma.r0")
    invisible(apply(unique(subset(d, R == 0, select = c("A1", "A2R", "A2NR"))), 1, function(x) {
      if (!is.null(uneqsdDTR)) {
        # Check if the current embedded DTR (x) is in uneqsdDTR; if it is, store the index
        # (If it's not, uneqsdDTRindex has length 0)
        uneqsdDTRindex <- which(sapply(uneqsdDTR, function(dtr) {
          sum(c(x[1], x[2], x[3]) == dtr) == 3
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
      assign(paste0("sigma.nr", (x[1] + 1) / 2, (x[3] + 1) / 2),
             unname(sapply(1:length(times), function(i) {
               generate.cond.sd(a1 = x[1], r = 0, a2R = x[2], a2NR = x[3], t = times[i], spltime, design = 2, 
                                rprob = get(paste0("r", (x[1] + 1) / 2)),
                                sigma = sigmaStar[i],
                                sigma.cond = get(paste0("sigma.r", (x[1] + 1) / 2))[i],
                                gammas = gammas, lambdas = lambdas)
             })), envir = parent.frame(n = 2))
    }))
  } else {
    
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
  if (corstr == "id") {
    cormat.r1 <- cormat.r0 <- diag(rep(1, length(times)))
  } else if (corstr == "exch") {
    cormat.r1 <- cormat.r0 <- cormat.exch(rho, length(times))
  } else if (corstr == "ar1") {
    cormat.r1 <- cormat.r0 <- cormat.ar1(rho, length(times))
  }
  # Create responders' correlation matrices by looping over m and computing appropriate entries
  for (i in 1:dim(m)[1]) {
    # If both times are after the second randomization, correlation between times is the provided value
    # if (min(m[i, ]) > which(times == spltime)) {
      cormat.r1[m[i, 1], m[i, 2]] <- rho.r1
      cormat.r0[m[i, 1], m[i, 2]] <- rho.r0
    # } else {
      # If both times aren't after the second randomization, compute the correlation required to 
      # achieve the proper marginalized covariance
      # cormat.r1[m[i, 1], m[i, 2]] <- cormat.r1[m[i, 1], m[i, 2]] * sigma[max(m[i, ])] / sigma.r1[max(m[i, ])]
      # cormat.r0[m[i, 1], m[i, 2]] <- cormat.r0[m[i, 1], m[i, 2]] * sigma[max(m[i, ])] / sigma.r0[max(m[i, ])]
    # }
  }
  
  if (design == 1) {
    cormat.r11 <- cormat.r1
    cormat.r01 <- cormat.r0
    
    for (a1ind in c(0, 1)) {
      if (corstr %in% c("exch", "id")) {
        x <- cormat.exch(rho, length(times))
      } else if (corstr == "ar1") {
        x <- cormat.ar1(rho, length(times))
      }
      for(i in 1:dim(m)[1]) {
        x[m[i, 1], m[i, 2]] <- 
          generate.cond.cor(a1 = 2 * a1ind - 1, r = 0, a2R = 1, a2NR = 1,
                            t1 = times[m[i, 1]], t2 = times[m[i, 2]], spltime = spltime,
                            design = 1, rprob = r1,
                            sigma.t1 = sigma[m[i, 1]], sigma.t2 = sigma[m[i, 2]], 
                            sigma.t1.ref = sigma.r11[m[i, 1]],
                            sigma.t2.ref = sigma.r11[m[i, 2]],
                            rho = rho,
                            rho.ref = cormat.r11[m[i, 1], m[i, 2]],
                            gammas = gammas, lambdas = lambdas)
      }
      assign(paste0("cormat.nr", a1ind, "1"), x)
    }
    
    cormat.nr11 <- cormat.exch(rho, length(times))
    for (i in 1:dim(m)[1]) {
      cormat.nr11[m[i, 1], m[i, 2]] <- 
        generate.cond.cor(a1 = 1, r = 0, a2R = 1, a2NR = 1,
                          t1 = times[m[i, 1]], t2 = times[m[i, 2]], spltime = spltime,
                          design = 1, rprob = r1,
                          sigma.t1 = sigma[m[i, 1]], sigma.t2 = sigma[m[i, 2]], 
                          sigma.t1.ref = sigma.r11[m[i, 1]],
                          sigma.t2.ref = sigma.r11[m[i, 2]],
                          rho = rho,
                          rho.ref = cormat.r11[m[i, 1], m[i, 2]],
                          gammas = gammas, lambdas = lambdas)
    }
    
    # FIXME: Finish this!
    
  } else if (design == 2) {
    cormat.r10 <- cormat.r1
    cormat.r00 <- cormat.r0
  }
  
  d <- do.call("rbind",
               lapply(split.SMART(d), function(x) {
                 a1ind   <- as.character((unique(x$A1) + 1) / 2)
                 a2Rind  <- ifelse(unique(x$A2R) == 0,  "0", as.character((unique(x$A2R) + 1) / 2))
                 a2NRind <- ifelse(unique(x$A2NR) == 0, "0", as.character((unique(x$A2NR) + 1) / 2))
                 if (unique(x$R) == 1 & (design != 1 | (design == 1 & a2Rind == 0))) {
                   x[, grepl("Y", names(x))] <- x[, grepl("Y", names(x))] +
                     mvrnorm(
                       dim(x)[1],
                       mu = rep(0, length(times)),
                       Sigma =
                         diag(get(paste0("sigma.r", a1ind, a2Rind))) %*%
                         get(paste0("cormat.r", a1ind, a2Rind)) %*%
                         diag(get(paste0("sigma.r", a1ind, a2Rind)))
                     )
                 } else if (unique(x$R) == 1 & design == 1 & a2Rind == 1) {
                   cormat.r <- cormat.exch(rho, length(times))
                   for (i in 1:dim(m)[1]) {
                     cormat.nr[m[i, 1], m[i, 2]] <- 
                       generate.cond.cor(a1 = unique(x$A1), r = unique(x$R), a2R = unique(x$A2R), a2NR = unique(x$A2NR),
                                         t1 = times[m[i, 1]], t2 = times[m[i, 2]], spltime = spltime,
                                         design = design, rprob = get(paste0("r", a1ind)),
                                         sigma.t1 = sigma[m[i, 1]], sigma.t2 = sigma[m[i, 2]], 
                                         sigma.t1.ref = get(paste0("sigma.nr", a1ind))[m[i, 1]],
                                         sigma.t2.ref = get(paste0('sigma.nr', a1ind))[m[i, 2]],
                                         rho = cormat.nr[m[i, 1], m[i, 2]],
                                         rho.r = get(paste0("cormat.nr", a1ind))[m[i, 1], m[i, 2]],
                                         gammas = gammas, lambdas = lambdas)
                   }
                   x[, grepl("Y", names(x))] <- x[, grepl("Y", names(x))] +
                     mvrnorm(
                       dim(x)[1],
                       mu = rep(0, length(times)),
                       Sigma = diag(get(paste0("sigma.r", a1ind, a2Rind))) %*%
                         get(paste0("cormat.r", a1ind, a2Rind)) %*%
                         diag(get(paste0("sigma.r", a1ind, a2Rind)))
                     )
                 } else if (design == 2) {
                   if (unique(x$R) == 0) {
                     cormat.nr <- cormat.exch(rho, length(times))
                     for (i in 1:dim(m)[1]) {
                       cormat.nr[m[i, 1], m[i, 2]] <- 
                         generate.cond.cor(a1 = unique(x$A1), r = 0, a2R = unique(x$A2R), a2NR = unique(x$A2NR),
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
                         Sigma = diag(get(paste0("sigma.nr", a1ind, a2NRind))) %*%
                           cormat.nr %*%
                           diag(get(paste0("sigma.nr", a1ind, a2NRind)))
                       )
                   }
                   else {
                     x[, grepl("Y", names(x))] <- x[, grepl("Y", names(x))] +
                       mvrnorm(
                         dim(x)[1],
                         mu = rep(0, length(times)),
                         Sigma = diag(get(paste0("sigma.r", a1ind, a2Rind))) %*%
                           get(paste0("cormat.r", a1ind, a2Rind)) %*%
                           diag(get(paste0("sigma.r", a1ind, a2Rind)))
                       )
                   }
                 }
                 x
               })
  )
  
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
