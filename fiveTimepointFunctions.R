conditional.model <- function(a1, r, a2, t, spltime, rprob, lambdas, gammas) {
  gammas[1] * (t == 0) + 
    (a1 == 1) * (gammas[2] * (t == 1) + gammas[3] * (t == 2) + 
                   (a2 %in% c(0, 1)) * (gammas[4] * (t == 3) + gammas[5] * (t == 4) + 
                                          (t > spltime) * (t - spltime) * (r - rprob) * lambdas[1] ) +
                   (a2 == -1) * (gammas[6] * (t == 3) + gammas[7] * (t == 4) +
                                   (t > spltime) * (t - spltime) * (r - rprob) * lambdas[2])) +
    (a1 == -1) * (gammas[8] * (t == 1) + gammas[9] * (t == 2) +
                    
                    (a2 %in% c(0, 1)) * (gammas[10] * (t == 3) + gammas[11] * (t == 4) + 
                                           (t > spltime) * (t - spltime) * (r - rprob) * lambdas[3]) +
                    (a2 == -1) * (gammas[12] * (t == 3) + gammas[13] * (t == 4) + 
                                    (t > spltime) * (t - spltime) * (r - rprob) * lambdas[4])) 
  }

marginal.model <- function(gammas, t, spltime, a1, a2) {
  gammas[1] * (t == 0) + 
    (a1 == 1) * (gammas[2] * (t == 1) + gammas[3] * (t == 2) + 
                   (a2 == 1) * (gammas[4] * (t == 3) + gammas[5] * (t == 4)) +
                   (a2 == -1) * (gammas[6] * (t == 3) + gammas[7] * (t == 4))) +
    (a1 == -1) * (gammas[8] * (t == 1) + gammas[9] * (t == 2) + 
                    (a2 == 1) * (gammas[10] * (t == 3) + gammas[11] * (t == 4)) +
                    (a2 == -1) * (gammas[12] * (t == 3) + gammas[13] * (t == 4)))
}

mod.derivs <- function(times, spltime) {
  lapply(1:4, function(dtr) {
    a1 <- switch(dtr, 1, 1, -1, -1)
    a2 <- switch(dtr, 1, -1, 1, -1)
    if (a1 == 1 & a2 == 1) {
      t(matrix(c(1, rep(0, 12), 1, 1, rep(0, 11), 1, 0, 1, rep(0, 10), 1, 0, 0, 1, rep(0, 9), 1, 0, 0, 0, 1, rep(0, 8)), ncol = 5))
    } else if (a1 == 1 & a2 == -1) {
      t(matrix(c(1, rep(0, 12), 1, 1, rep(0, 11), 1, 0, 1, rep(0, 10), 1, rep(0, 4), 1, rep(0, 7), 1, rep(0, 5), 1, rep(0, 6)), ncol = 5))
    } else if (a1 == -1 & a2 == 1) {
      t(matrix(c(1, rep(0, 12), 1, rep(0, 6), 1, rep(0, 5), 1, rep(0, 7), 1, rep(0, 4), 1, rep(0, 8), 1, rep(0, 3), 1, rep(0, 9), 1, rep(0, 2)), ncol = 5))
    } else if (a1 == -1 & a2 == -1) {
      t(matrix(c(1, rep(0, 12), 1, rep(0, 6), 1, rep(0, 5), 1, rep(0, 7), 1, rep(0, 4), 1, rep(0, 10), 1, rep(0, 1), 1, rep(0, 11), 1), ncol = 5))
    }
  })
}