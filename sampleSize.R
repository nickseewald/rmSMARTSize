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

sample.size <- function(delta, r = NULL, r1 = r, r0 = r, rho, alpha = 0.05, power = .8,
                        design = 2, rounding = c("up", "down"),
                        conservative = TRUE) {
  rounding <- match.arg(rounding)
  
  if (is.null(r)) r <- mean(c(r1, r0))
  
  # Input checks
  if (!is.null(rounding) & !(rounding %in% c("up", "down"))) 
    stop("rounding must be either up or down")
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
      correction <- ((1 - rho) * (2 * rho + 1)) / (1 + rho) + 
        (1 - rho)/(1 + rho) * rho^2 / (2 - r)
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
  message('generating sample size')
  return(n)
}
