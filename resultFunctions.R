# resultFunctions.R
# Functions for compiling simulation results into useful tables
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

#' Compile sim results into a LaTeX-ready table
#'
#' @param results A list of objects of class \code{simResult} to be tabulated
#'
#' @return
#' @export
#'
#' @examples
resultTable <- function(results, 
                        alternative = c("two.sided", "less", "greater"),
                        paper = FALSE, alpha = 0.05) {
  alternative <- match.arg(alternative)
  
  if (alternative == "two.sided")
    
    if (is.null(names(results)))
      stop("results must be a named list.")
  
  d <- do.call("rbind", lapply(results, function(l) {
    if (alternative == "two.sided") {
      pval.power <- l$pval.power
    } else {
      pval.power <- binom.test(x = sum(l$pval <= l$alpha / 2, na.rm = T) +
                                 sum(l$pval >= 1 - l$alpha / 2, na.rm = T),
                               n = l$valid, p = l$power.target,
                               alternative = alternative)$p.value
    }
    
    respCor.mean <- apply(l$respCor, 1, function(x) x[which.max(abs(x))])
    
    condCov.diff <- apply(l$condCov, 1, 
                          function(x) ((x[2] >= x[1]) - (x[2] < x[1])) *
                            abs(x[2] - x[1]) / abs(x[1]))
    
    x <- data.frame(
      "delta" = l$delta,
      # "rprob" = min(l$r0, l$r1),
      "r0" = l$r0,
      "r1" = l$r1,
      "rho" = l$rho,
      "rho.size" = ifelse(is.null(l$rho.size), l$rho, l$rho.size),
      "sharp" = ifelse(l$sharp, "sharp", "cons."),
      "n" = l$n,
      "estPwr" = l$power,
      "pval.power" = pval.power,
      "coverage" = l$coverage,
      "pval.coverage" = tryCatch(
        binom.test(
          x = l$coverage * l$valid,
          n = l$valid,
          p = .95,
          alternative = "two.sided"
        )$p.value,
        error = function(e)
          NULL
      ),
      "corY0sq" = as.numeric(respCor.mean[1]),
      "corY1sq" = as.numeric(respCor.mean[3]),
      "corY0Y1" = as.numeric(respCor.mean[2]),
      "condCov02" = as.numeric(condCov.diff[1]),
      "condCov12" = as.numeric(condCov.diff[2]),
      "nTrial" = l$niter,
      "nValid" = l$valid
    )
    # if (paper) cbind(d, l$design)
  }))
  
  resNames <- strsplit(names(results), "\\.")
  resNames <- do.call(rbind,
                      lapply(resNames, function(rn) {
                        x <- rn[length(rn)]
                        if (x %in% c("sharp", "old")) {
                          x <- rn[length(rn) - 1]
                          if (x %in% c("sharp", "old"))
                            x <- rn[length(rn) - 2]
                        }
                        x
                      })
  )
  
  d$respFunc <- resNames
  d <- d[, c("delta", "respFunc", names(d)[2:(ncol(d) - 1)])]
  
  d <- d[order(d$delta, d$respFunc,
               d$r0, d$r1, d$rho, d$sharp), ]
  
  oldColnames <- colnames(d)
  
  colnames(d) <- c("$\\delta$", 
                   "Resp. Func.",
                   # "$\\min\\left\\{r_{-1},r_{1}\\right\\}$",
                   "$r_{-1}$",
                   "$r_{1}$",
                   "$\\rho$",
                   "$\\rho_{\\text{size}}$",
                   "Formula", "$n$", "$1-\\hat{\\beta}$", "$p$ value",
                   "Coverage",
                   "$p$ value", 
                   "$\\cor(R,Y_0^2)$", "$\\cor(R,Y_1^2)$", "$\\cor(R,Y_0Y_1)$",
                   "$\\cov(Y_0, Y_2)$", "$\\cov(Y_1, Y_2)$",
                   "Num. Runs", "Valid Runs")
  
  if (paper) {
    
  }
  
  print(
    xtable(d, digits = c(1, 1, 0, 1, 1, 1, 1, 0, 0,
                         3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0)),
    sanitize.text.function = identity,
    booktabs = TRUE,
    include.rownames = F,
    include.colnames = T,
    floating = F,
    tabular.environment = "longtable"
  )
  
  colnames(d) <- oldColnames
  return(invisible(d))
}