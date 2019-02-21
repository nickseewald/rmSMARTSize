# simulateSMART.R
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

#' Simulation wrapper for SMARTs with repeated-measures outcomes
#'
#' @param n number of participants in trial
#' @param gammas Parameters from marginal mean model
#' @param lambdas Parameters for response "offset" in conditional model
#' @param times Vector of times at which measurements are collected (currently 
#' limited to length three)
#' @param spltime Time (contained in times) at which re-randomization occurs
#' @param alpha Type-I error rate for hypothesis tests
#' @param power Target power
#' @param delta Deprecated. Effect size (determined by gammas and lambdas; this
#' exists only for identification purposes)
#' @param design One of 1, 2, or 3. Identifies which type of SMART to simulate.
#' @param rounding The direction to round computed decimal sample sizes. 
#' Passed to sample.size; either "up" or "down".
#' @param conservative Logical. Whether to use conservative or "sharp" sample 
#' size formulae; passed to sample.size
#' @param r  Probability of response to first-stage treatment (assuming equal 
#' for both txts)
#' @param r1 Probability of response to A1 = 1  (assuming not equal to r0; 
#' r must be NULL)
#' @param r0  Probability of response to A1 = -1 (assuming not equal to r1; 
#' r must be NULL)
#' @param uneqsdDTR list of vectors of (a1, a2R, a2NR) for DTR(s) which do not
#' have the same variance as the others
#' @param uneqsd 
#' @param sigma Marginal variance of Y (assumed constant over time and DTR)
#' @param sigma.r1 Conditional variance of Y for responders to treatment A1 = 1
#' @param sigma.r0 Conditional variance of Y for responders to treatment A1 = -1
#' @param corstr Character string, one of "identity", "exch"/"exchangeable", 
#' "unstr"/"unstructured". This is the TRUE form of the within-person 
#' correlation structure
#' @param rho If corstr is either "exch" or "exchangeable", the within-person 
#' correlation used to generate data (assumed constant across DTRs)
#' @param rho.r1 
#' @param rho.r0 
#' @param rho.size 
#' @param L 
#' @param varmats 
#' @param respFunction 
#' @param respDirection 
#' @param pool.time 
#' @param pool.dtr 
#' @param niter 
#' @param tol 
#' @param maxiter.solver 
#' @param save.data 
#' @param empirical 
#' @param balanceRand 
#' @param notify 
#' @param pbDevice 
#' @param postIdentifier 
#'
#' @return
#' @export
#'
#' @examples


simulateSMART <- function(n = NULL,
                          gammas,
                          lambdas,
                          times,
                          spltime,
                          alpha = .05,
                          power = .8,
                          delta,
                          design = 2,
                          rounding = "up",
                          conservative = TRUE,
                          r = NULL,
                          r1 = r,
                          r0 = r,
                          uneqsdDTR = NULL,
                          uneqsd = NULL,
                          sigma,
                          sigma.r1 = sigma,
                          sigma.r0 = sigma,
                          corstr = c("identity", "exchangeable", "ar1"),
                          rho = NULL,
                          rho.r1 = rho,
                          rho.r0 = rho,
                          rho.size = rho,
                          L = NULL,
                          varmats = NULL,
                          variances = NULL,
                          respFunction = response.oneT,
                          respDirection = c("high", "low"),
                          pool.time = TRUE,
                          pool.dtr = TRUE,
                          niter = 5000,
                          tol = 1e-8,
                          maxiter.solver = 1000,
                          save.data = FALSE,
                          empirical = FALSE,
                          balanceRand = FALSE,
                          notify = FALSE,
                          pbDevice = NULL,
                          postIdentifier = NULL,
                          ...
) {
  
  call <- match.call()
  respDirection <- match.arg(respDirection)
  
  # if (!is.null(r) & (r < 0 | r > 1)) stop("r must be between 0 and 1.")
  if (is.null(r) & is.null(r1) & is.null(r0)) 
    stop("You must provide either r or both r1 and r0.")
  ## TODO: Finish input handling
  
  # Handle uneqsdDTR
  if (design == 1 & is.null(uneqsdDTR)) {
    if (length(unique(marginal.model(dtrIndex(1)$a1, dtrIndex(1)$a2R, 
                                     dtrIndex(1)$a2NR, 2, spltime, 1, gammas)))
        == 8)
      stop(paste("For design I, you must provide a list of vectors of", 
                 "(a1, a2R, a2NR) for DTRs which do not have the same variance", 
                 "as others.\n Alternatively, you may specify gammas such that",
                 "there are only 4 unique DTR means (2 per first-stage", 
                 "treatment)."))
  }
  if (design == 1 & !is.null(uneqsdDTR) & !is.list(uneqsdDTR))
    stop("uneqsdDTR must be either NULL or a list.")
  
  nDTR <- switch(design, 8, 4, 3)
  
  ## Handle correlation structure
  corstr <- match.arg(corstr)
  if (corstr == "identity") {
    rho <- rho.r1 <- rho.r0 <- 0
    if (is.null(rho.size)) rho.size <- rho
  }
  
  ## If design == 1, pool.dtr is impossible to satisfy in a generative model. 
  ## Set it to false.
  # if (design == 1) pool.dtr <- FALSE
  
  # If n is not provided, compute it from the other inputs
  if (is.null(n)) {
    n <- sample.size(delta = delta, r = r, r1 = r1, r0 = r0,
                     rho = rho.size, alpha = alpha, power = power,
                     design = design, rounding = rounding,
                     conservative = conservative)
  }
  
  # Compute conditional variances
  if (old) {
    varmats <- conditionalVarmat(times, spltime, design, r1, r0,
                                 corstr = corstr,
                                 sigma, sigma.r1 = sigma.r1, sigma.r0 = sigma.r0,
                                 uneqsd = NULL, uneqsdDTR = NULL,
                                 rho, rho.r1, rho.r0,
                                 gammas, lambdas)
  }
  
  ## Construct string describing simulation parameters
  designText <- paste0("Design ", design, "\n",
                       ifelse(is.null(postIdentifier), "", postIdentifier),
                       "\ndelta = ", delta, "\n",
                       "true corstr = ", corstr, "(", rho, ")\n",
                       "sized for exchangeable(", rho.size, ")\n",
                       "r0 = ", round(r0, 3), ", r1 = ", round(r1, 3),
                       "\nn = ", n, "\n",
                       niter, " iterations.")
  
  ## Print message in console describing the simulation that's currently running
  message(
    paste0("********************\n",
           "Starting simulation!\n",
           designText,
           "\n********************\n")
  )
  if (notify)
    slackr_bot(designText)
  
  results <-
    foreach(i = 1:niter, .combine = combine.results, .final = finalize.results,
            .export = ls(envir = .GlobalEnv),
            .packages = c("MASS", "xtable", "slackr"),
            .verbose = FALSE, .errorhandling = "remove", .multicombine = FALSE,
            .inorder = FALSE) %dorng% {
              
              d <- generateSMART(n, times, spltime, r1, r0, gammas, lambdas, 
                                 design = design, sigma, sigma.r1, sigma.r0,
                                 corstr = corstr, rho, rho.r1, rho.r0, 
                                 uneqsd = NULL, uneqsdDTR = NULL,
                                 varmats = varmats, variances = variances,
                                 balanceRand = balanceRand, 
                                 empirical = empirical,
                                 respFunction = respFunction,
                                 respDirection = respDirection)#, ...)
              
              if (d$valid == FALSE) {
                ## If a non-valid trial has been generated (i.e., fewer than one
                ## observation per cell), return a "blank" result
                
                sigma2.blank <-
                  matrix(
                    ncol = (length(times) * (1 - pool.time) * pool.dtr) +
                      4 * (1 - pool.dtr) + (pool.time * pool.dtr),
                    nrow = (
                      length(times) * (1 - pool.dtr) +
                        (4 * (1 - pool.dtr) * pool.time) +
                        (pool.time * pool.dtr)
                    )
                  )
                
                result <- list(
                  "pval" = NA,
                  "param.hat" = rep(NA, length(gammas)),
                  "respCor" = matrix(0,
                                     ncol = sum(times <= spltime),
                                     nrow = nDTR),
                  "param.var" = matrix(0, ncol = length(gammas),
                                       nrow = length(gammas)),
                  "sigma2.hat" = sigma2.blank,
                  "rho.hat" = NA,
                  "valid" = 0,
                  "varFailCount" = varFailCount,
                  "coverage" = 0,
                  "iter" = NULL,
                  "sigma.r0" = NULL,
                  "sigma.r1" = NULL,
                  "condVars" = lapply(1:length(dtrIndex(design)$a1),
                                      function(x)
                                        matrix(0, ncol = length(times),
                                               nrow = length(times))),
                  "respCor" = matrix(rep(0, 
                                         ncol(combn(times[times <= spltime], 2))
                                         * 2),
                                     ncol = 2),
                  "condCov" = matrix(rep(0, nrow(
                    expand.grid(times[times <= spltime], 
                                times[times > spltime])) * 2),
                    ncol = 2)
                  # "assumptionViolations" = d$assumptions
                )
                
                if (save.data) {
                  result[["data"]] <- list(d$data)
                }
                return(result)
                
              } else {
                ### FOR VALID TRIALS ###
                
                d1 <- reshape(d$data, varying = list(grep("Y", names(d$data))),
                              ids = d$data$id, 
                              times = times, direction = "long", v.names = "Y")
                d1 <- d1[order(d1$id, d1$time), ]
                
                # Check working assumption re: correlation between response and
                # products of residuals
                respCor <- estimate.respCor(d$data, design, times, spltime, 
                                            gammas)
                
                condCov <- estimate.condCov(d1, times, spltime,
                                            design, pool.dtr = T)
                
                param.hat <-
                  try(estimate.params(d1,
                                      diag(rep(1, length(times))),
                                      times,
                                      spltime,
                                      design,
                                      rep(0, length(gammas)),
                                      maxiter.solver,
                                      tol))
                
                # param.hat2 <-
                #   multiroot(
                #     esteqn.compute,
                #     start = rep(0, 7),
                #     jacfunc = esteqn.jacobian,
                #     d = d1,
                #     V = diag(rep(1, length(times))),
                #     times = times,
                #     spltime = spltime,
                #     design = 2
                #   )
                
                sigma2.hat <- estimate.sigma2(d1,
                                              times,
                                              spltime,
                                              design,
                                              param.hat,
                                              pool.time = pool.time,
                                              pool.dtr = pool.dtr)
                
                rho.hat <-
                  estimate.rho(d1, times, spltime, design, sqrt(sigma2.hat),
                               param.hat, corstr = "exchangeable")
                
                # Compute variance matrices for all conditional cells
                condVars <- lapply(split.SMART(d$data), function(x) {
                  var(subset(x, select = grep("Y", names(x), value = TRUE)))
                })
                
                # Iterate parameter estimation
                outcome.var <-
                  varmat(sigma2.hat, rho.hat, times, design, 
                         corstr = "exchangeable")
                
                param.hat <-
                  estimate.params(d1,
                                  outcome.var,
                                  times,
                                  spltime,
                                  design,
                                  param.hat,
                                  maxiter.solver,
                                  tol)
                
                # Iterate until estimates of gammas and rho converge
                for (j in 1:maxiter.solver) {
                  sigma2.new <- estimate.sigma2(d1,
                                                times,
                                                spltime,
                                                design,
                                                param.hat,
                                                pool.time = pool.time,
                                                pool.dtr = pool.dtr)
                  rho.new <-
                    estimate.rho(d1,
                                 times,
                                 spltime,
                                 design,
                                 sqrt(sigma2.hat),
                                 param.hat,
                                 corstr = "exchangeable")
                  outcomeVar.new <-
                    varmat(sigma2.new, rho.new, times, design, 
                           corstr = "exchangeable")
                  param.new <-
                    estimate.params(d1,
                                    outcomeVar.new,
                                    times,
                                    spltime,
                                    design,
                                    start = param.hat,
                                    maxiter.solver,
                                    tol)
                  if (norm(param.new - param.hat, type = "F") <= tol &
                      norm(as.matrix(sigma2.new) - as.matrix(sigma2.hat), 
                           type = "F") <= tol &
                      (rho.new - rho.hat) ^ 2 <= tol) {
                    param.hat <- param.new
                    sigma2.hat <- sigma2.new
                    rho.hat <- rho.new
                    iter <- j
                    break
                  } else {
                    param.hat <- param.new
                    sigma2.hat <- sigma2.new
                    rho.hat <- rho.new
                    iter <- j
                  }
                }
                param.var <-
                  estimate.paramvar(d1,
                                    cormat(rho.hat, length(times), 
                                           corstr = "exchangeable"),
                                    times,
                                    spltime,
                                    design,
                                    gammas = param.hat)
                
                confLB <- 
                  L %*% param.hat - sqrt(L %*% param.var %*% L) * qnorm(.975)
                confUB <-
                  L %*% param.hat + sqrt(L %*% param.var %*% L) * qnorm(.975)
                
                coverage <- 
                  ifelse(confLB <= L %*% gammas & L %*% gammas <= confUB, 1, 0)
                
                pval <-
                  1 - pnorm(as.numeric((L %*% param.hat) / 
                                         sqrt(L %*% param.var %*% L)))
                cat("\npval")
                result <-
                  list(
                    "pval" = pval,
                    "param.hat" = t(param.hat),
                    "param.var" = param.var,
                    "sigma2.hat" = sigma2.hat,
                    "rho.hat" = rho.hat,
                    "valid" = 1,
                    "varFailCount" = d$varFailCount,
                    "coverage" = coverage,
                    "respCor" = respCor,
                    "condCov" = condCov,
                    "iter" = iter,
                    "sigma.r0" = d$params$sigma.r0,
                    "sigma.r1" = d$params$sigma.r1,
                    "condVars" = condVars
                    # "assumptionViolations" = d$assumptions
                  )
                
                if (save.data) {
                  result[["data"]] <- list(d$data)
                }
                return(result)
              }
            }
  
  seedList <- attr(results, 'rng')
  
  results <-
    c(
      list(
        "call" = call,
        "n" = n,
        "alpha" = alpha,
        "power.target" = power,
        "delta" = delta,
        "corstr" = corstr,
        "rho" = rho,
        "rho.r0" = rho.r0,
        "rho.r1" = rho.r1,
        "rho.size" = rho.size,
        "sigma" = sigma,
        "sigma.r0" = sigma.r0,
        "sigma.r1" = sigma.r1,
        "r0" = r0,
        "r1" = r1,
        "niter" = niter,
        "sharp" = !conservative
      ),
      results
    )
  
  test <-
    binom.test(
      x = sum(results$pval <= results$alpha / 2, na.rm = T) + 
        sum(results$pval >= 1 - results$alpha / 2, na.rm = T),
      n = results$valid,
      p = results$power.target,
      alternative = "two.sided"
    )
  
  results <-
    c(
      results,
      "power" = unname(test$estimate),
      "pval.power" = unname(test$p.value)
    )
  
  class(results) <- c("simResult", class(results))
  
  attr(results, 'rng') <- seedList
  
  if (notify) {
    try(slackr_bot(print(results)))
  }
  
  return(results)
}
