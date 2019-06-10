# sampleSize.R
# Function to compute sample size for a SMART with a longitudinal outcome
# in which the primary aim is to compare two embedded DTRs.
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

#' Sample Size for End-of-Study Comparisons of Embedded DTRs in
#' SMARTs with Longitudinal Outcomes
#'
#' @param delta The target standardized effect size
#' @param r (optional) The common probability of response to first-stage
#' intervention options
#' @param r1 Probability of response to the first-stage intervention option
#' denoted a1 = 1
#' @param r0 Probability of response to the first-stage intervention option
#' denoted a1 = -1
#' @param rho Within-person correlation
#' @param alpha Permissible type-I error of the test
#' @param power Desired power of the test
#' @param design 1, 2, or 3 (see details)
#' @param corstr 
#' @param rounding 
#' @param conservative 
#'
#' @return
#' @export
#'
#' @examples
sample.size <- function(delta, r = NULL, r1 = r, r0 = r, rho,
                        alpha = 0.05, power = .8,
                        design = 2, corstr = c("exchangeable", "ar1"),
                        path = c("separate", "shared"),
                        rounding = c("up", "down"), sigma = NULL,
                        conservative = TRUE) {
  
  rounding <- match.arg(rounding)
  corstr   <- match.arg(corstr)
  path     <- match.arg(path)
  
  if (is.null(r)) r <- mean(c(r1, r0))
  
  # Input checks
  if (!is.null(rounding) & !(rounding %in% c("up", "down"))) 
    stop("rounding must be either up or down")
  
  # standard sample size formula
  nid <- (4 * (qnorm(1 - alpha / 2) + qnorm(power))^2) / delta^2
  
  # longitudinal deflation
  correction <- 1 - rho^2
  
  # design effects (and modifications to deflation factor for sharp formula)
  if (design == 1) {
    if (!conservative)
      correction <- (1-rho^2) - (1-rho)*rho^2 / (2*(1+rho))
    if (path == "separate") {
      if (corstr == "exchangeable")
        designEffect <- 2
      else if (corstr == "ar1")
        designEffect <- 2 + rho^2
    } else
      # designEffect <
      stop("shared path not ready for design 1 yet")
    
  } else if (design == 2) {
    if (path == "shared") {
        if (conservative) {
          designEffect <- 2 / (1 - r)
        } else {
          designEffect <- 1
          correction <- (1 - r) * (1 + rho - 2*rho^2) / (1 + rho) - r / (2 * sigma^2)
        }
    } else {
      designEffect <- ((2 - r1) + (2 - r0)) / 2
      if (corstr == "ar1") 
        designEffect <- designEffect + rho^2
      if (!conservative) {
        correction <- ((1 - rho) * (2 * rho + 1)) / (1 + rho) + 
          (1 - rho)/(1 + rho) * rho^2 / (2 - r)
      }
    }
  
    } else if (design == 3) {
    designEffect <- (3 - r1) / 2
    if (!conservative) {
      correction <- ((3-r1)*(1+2*rho) + 2*rho^2) * (1-rho)/((3-r1)*(1+rho))
    }
  }
  else stop("Not a valid design indicator.")
  
  if (rounding == "up") {
    n <- ceiling(nid * designEffect * correction)
  } else if (rounding == "down") {
    n <- floor(nid * designEffect * correction)
  }
  return(n)
}
