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


#' Title
#'
#' @param n 
#' @param gammas 
#' @param lambdas 
#' @param times 
#' @param spltime 
#' @param alpha 
#' @param power 
#' @param delta 
#' @param design 
#' @param round 
#' @param conservative 
#' @param r 
#' @param r1 
#' @param r0 
#' @param uneqsdDTR 
#' @param uneqsd 
#' @param sigma 
#' @param sigma.r1 
#' @param sigma.r0 
#' @param corstr 
#' @param rho 
#' @param rho.r1 
#' @param rho.r0 
#' @param L 
#' @param constant.var.time 
#' @param constant.var.dtr 
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
simulateSMART <- function(n = NULL, gammas, lambdas, times, spltime,
                alpha = .05, power = .8, delta, design = 2, round = "up", conservative = TRUE,
                r = NULL, r1 = r, r0 = r,
                uneqsdDTR = NULL, uneqsd = NULL, 
                sigma, sigma.r1 = sigma, sigma.r0 = sigma,
                corstr = c("identity", "exchangeable", "ar1"),
                rho = NULL, rho.r1 = rho, rho.r0 = rho,
                L = NULL, varmats = NULL,
                constant.var.time = TRUE, constant.var.dtr = TRUE,
                niter = 5000, tol = 1e-8, maxiter.solver = 1000,
                save.data = FALSE, empirical = FALSE, balanceRand = FALSE,
                notify = FALSE, pbDevice = NULL, postIdentifier = NULL) {

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
  corstr <- match.arg(corstr)
  if (corstr == "identity") rho <- rho.r1 <- rho.r0 <- 0
  
  ## If design == 1, constant.var.dtr is impossible to satisfy in a generative model. Set it to false.'
  # if (design == 1) constant.var.dtr <- FALSE
  
  # If n is not provided, compute it from the other inputs
  if (is.null(n)) 
    n <- sample.size(delta = delta, r = ifelse(is.null(r), min(r0, r1), r),
                     rho = rho, alpha = alpha, power = power,
                     design = design, round = round,
                     conservative = conservative)
  
  # Compute conditional variances
  varmats <- conditionalVarmat(times, spltime, design, r1, r0, 
                               corstr = corstr,
                               sigma, sigma.r1, sigma.r0,
                               uneqsd = NULL, uneqsdDTR = NULL,
                               rho, rho.r1, rho.r0,
                               gammas, lambdas)
  
  ## Simulate
  cat(paste0("********************\n",
             "Starting simulation!\n",
             "delta = ", delta, "\n",
             "corstr = exch(", rho, ")\n",
             "r0 = ", r0, ", r1 = ", r1, "\n",
             ifelse(is.null(postIdentifier), NULL, postIdentifier),
             "\nn = ", n, "\n",
             niter, " iterations.\n",
             "********************\n"))
  results <- foreach(i = 1:niter, .combine = combine.results, .final = finalize.results,
                     .verbose = FALSE, .errorhandling = "stop", .multicombine = FALSE, .inorder = FALSE) %dorng% { 
                       
                       d <- generateSMART(n, times, spltime, r1, r0, gammas, lambdas, design = design, sigma, sigma.r1, sigma.r0, corstr = corstr,
                                          rho, rho.r1, rho.r0, uneqsd = NULL, uneqsdDTR = NULL, varmats = varmats,
                                          balanceRand = balanceRand, empirical = empirical)
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
                         if (corstr == "identity" | (corstr == "exchangeable" & rho == 0)) {
                           outcome.var <- varmat(sigma2.hat, 0, times, "exchangeable")
                           param.var <- estimate.paramvar(d1, diag(rep(1, length(times))), times, spltime, design, gammas = param.hat)
                           iter <- 1
                         } else if (corstr == "exchangeable" & rho != 0) {
                           outcome.var <- varmat(sigma2.hat, rho.hat, times, "exchangeable")
                           param.hat <- estimate.params(d1, outcome.var, times, spltime, design, param.hat, maxiter.solver, tol)
                           # Iterate until estimates of gammas and rho converge
                           for (i in 1:maxiter.solver) {
                             sigma2.new <- estimate.sigma2(d1, times, spltime, design, param.hat,
                                                           pool.time = constant.var.time, pool.dtr = constant.var.dtr)
                             rho.new <- estimate.rho(d1, times, spltime, design, sqrt(sigma2.hat), param.hat)
                             outcomeVar.new <- varmat(sigma2.new, rho.new, times, "exchangeable")
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
                         
                         confLB <- L %*% param.hat - sqrt(L %*% param.var %*% L) * qnorm(.975)
                         confUB <- L %*% param.hat + sqrt(L %*% param.var %*% L) * qnorm(.975)
                         
                         coverage <- ifelse(confLB <= L %*% gammas & L %*% gammas <= confUB, 1, 0)
                         
                         pval <- 1 - pnorm(as.numeric((L %*% param.hat) / sqrt(L %*% param.var %*% L)))
                         
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
                    "corstr" = corstr, "rho" = rho, "rho.r0" = rho.r0, "rho.r1" = rho.r1,
                    "sigma" = sigma, "sigma.r0" = sigma.r0, "sigma.r1" = sigma.r1,
                    "r0" = r0, "r1" = r1, "niter" = niter, "sharp" = !conservative), results)
  
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