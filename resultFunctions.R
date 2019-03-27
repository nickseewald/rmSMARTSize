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
                        paper = FALSE, alpha = 0.05, caption = NULL) {
  
  alternative <- match.arg(alternative)
  
  if (length(caption) == 1) {
    caption <- c(caption, "")
  } else if (length(caption) > 2) {
    warning(paste("caption argument can be at most length 2; ignoring all but",
                  "the first two elements."))
  }
    
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
    
    pval.coverage <- tryCatch(
      binom.test(
        x = l$coverage * l$valid,
        n = l$valid,
        p = .95,
        alternative = "two.sided"
      )$p.value,
      error = function(e)
        NULL
    )
    
    x <- data.frame(
      "delta" = l$delta,
      "sharp" = ifelse(l$sharp, "sharp", "cons."),
      # "rprob" = min(l$r0, l$r1),
      "r0" = l$r0,
      "r1" = l$r1,
      "rho" = l$rho,
      "rho.size" = ifelse(is.null(l$rho.size), l$rho, l$rho.size),
      "n" = l$n,
      "estPwr" = l$power,
      "pval.power" = formatOutput(pval.power, pval.power, 0.05, 1, 1),
      "coverage" = l$coverage,
      "pval.coverage" = formatOutput(pval.coverage, pval.coverage, 0.05, 1, 1),
      "corY0sq" = as.numeric(respCor.mean[1]),
      "corY1sq" = as.numeric(respCor.mean[3]),
      "corY0Y1" = as.numeric(respCor.mean[2]),
      "condCov02" = as.numeric(condCov.diff[1]),
      "condCov12" = as.numeric(condCov.diff[2]),
      "propValid" = 100 * l$valid / l$niter
    )
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
  
  d <- d[order(d$delta, d$respFunc, d$sharp,
               d$r0, d$r1, d$rho), ]
  
  if (paper) {
    d$respFunc[d$respFunc == "indepHigh"]    <- "$R_{\\perp\\!\\!\\!\\perp}$"
    d$respFunc[d$respFunc == "betaHigh"]     <- "$R_{+}$"
    d$respFunc[d$respFunc == "betaLow"]      <- "$R_{-}$"
    d$respFunc[d$respFunc == "sqHigh"]       <- "$R_{+}^{2}$"
    d$respFunc[d$respFunc == "sqLow"]        <- "$R_{-}^{2}$"
    d$respFunc[d$respFunc == "oneTHigh"]     <- "$R_{+}^{1\\kappa}$"
    d$respFunc[d$respFunc == "oneTLow"]      <- "$R_{-}^{1\\kappa}$"
    d$respFunc[d$respFunc == "twoTHigh"]     <- "$R_{+}^{2\\kappa}$"
    d$respFunc[d$respFunc == "twoTLow"]      <- "$R_{-}^{2\\kappa}$"
    d$respFunc[d$respFunc == "nrCorrDown50"] <- "$R_{\\perp\\!\\!\\!\\perp}$"
  }
  
  oldColnames <- colnames(d)
  
  colnames(d) <- c("$\\delta$", 
                   "Response",
                   "Formula",
                   "$r_{-1}$",
                   "$r_{1}$",
                   "$\\rho$",
                   "$\\rho_{\\text{size}}$", "$n$", "$1-\\hat{\\beta}$", "$p$ value",
                   "Coverage",
                   "$p$ value", 
                   "$\\cor(R,Y_0^2)$", "$\\cor(R,Y_1^2)$", "$\\cor(R,Y_0Y_1)$",
                   "$\\cov(Y_0, Y_2)$", "$\\cov(Y_1, Y_2)$",
                   "\\% Valid")
  
  if (paper) {
    xtab <- xtable(d, digits = c(0, 1, 0, 0, 1, 1, 1, 1, 0,
                                 3, 3, 3, 3, 3, 3, 3, 3, 3, 1))
    
    ## Inspired by https://stackoverflow.com/a/33166479: 
    cat(paste0("\n",
               "\\begingroup\\footnotesize \n",
              "\\begin{longtable}{",
              paste(rep("l", ncol(xtab)), collapse = ""),
              "}\n",
              "\\caption[", caption[2], "]{", caption[1], "} \\\\ \\toprule \n", 
              "\\endfirsthead",
              "\\caption*{", caption[1], "} \\\\ \\toprule \n",
              paste(colnames(xtab), collapse = " & "), 
              " \\\\ \\midrule \\endhead",
              "\\midrule \n",
              "\\multicolumn{4}{l}{\\footnotesize{Continued on next page}} \n",
              "\\endfoot \\endlastfoot"))
    
    print(xtab,
          include.colnames = TRUE,
          include.rownames = FALSE,
          only.contents = TRUE,
          sanitize.text.function = identity
    )
    
    cat("\\end{longtable} \n \\endgroup")
  } else {
    print(
      xtable(d, digits = c(0, 1, 0, 0, 1, 1, 1, 1, 0,
                           3, 3, 3, 3, 3, 3, 3, 3, 3, 1), caption = caption),
      sanitize.text.function = identity,
      booktabs = F,
      include.rownames = F,
      include.colnames = T,
      floating = F,
      tabular.environment = "longtable",
      size = "footnotesize",
      caption.placement = "top",
      caption.width = "\\textwidth",
      hline.after = c(-1, nrow(d)),
      add.to.row = list(pos = list(0),
                        command = paste("\\hline \\endhead "))
    )
  }
  
  colnames(d) <- oldColnames
  return(invisible(d))
}

parseResultFileName <- function(filename) {
  x <- unlist(strsplit(filename, "-"), use.names = F)
  x[length(x)] <- unlist(strsplit(x[length(x)], "\\."), use.names = F)[1]
  designString <- switch(x[1], 
                         "simsDesign1" = "\\Cref{sub@fig:design-allrerand}",
                         "simsDesign2" = "\\Cref{sub@fig:design-prototypical",
                         "simsDesign3" = "\\Cref{sub@fig:design-autism}")
  condition <- switch(x[3], 
                      "ar1" = "True correlation structure is AR(1)",
                      "basic" = "All working assumptions satisfied",
                      "sLBdown25" = paste0("\\var\\left( Y_{2}^{(d)}", 
                                           "\\mid R^{(a_1)} \\right) ",
                                           "25\\% below lower bound"))
  paste0(designString, ": ", condition)
}

formatOutput <- function(x, p, alpha, nvalid, niter) {
  if (p <= alpha & nvalid < niter) {
    out <- paste0(formatC(round(x, 3), format = "f", digits = 3), "$^{*\\S}$")
  } else if (p <= alpha & nvalid == niter) {
    out <- paste0(formatC(round(x, 3), format = "f", digits = 3), "$^*$")
  } else if (p > alpha & nvalid < niter) {
    out <- paste0(formatC(round(x, 3), format = "f", digits = 3), "$^{\\S}$")
  } else {
    out <- formatC(round(x, 3), format = "f", digits = 3)
  }
  out
}
