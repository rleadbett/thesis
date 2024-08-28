# This R script generates the prior predictive simulations for Figure 3.1 in
# the manuscript Leadbetter, R.L, Gonzalez Caceres, G. J. and Phatak, A. (2024)
# `Fitting noisy gamma processes using RStan'

library(tidyverse)
library(tidybayes)
library(cowplot)
library(bayesplot)

theme_set(theme_tidybayes() + panel_border())

# Generate the simulations for the three alternative priors
set.seed(509480983)

## Ga(1, 0.001) for shape and rate
time <- 0:11
I <- length(time)
N <- 100
beta <- 1; xi <- 0.001
shape_rate <- matrix(rgamma(2 * N, shape = beta, rate = xi), ncol = 2)
jumps <- apply(shape_rate, 1, function(x) rgamma(I - 1, shape = x[1], rate = x[2]))
trace1 <- rbind(0, apply(jumps, 2, cumsum))

## Ga(0.001, 0.001) for shape and rate
time <- 0:11
I <- length(time)
N <- 100
beta <- 0.001; xi <- 0.001
shape_rate <- matrix(rgamma(2 * 100, shape = beta, rate = xi), ncol = 2)
jumps <- apply(shape_rate, 1, function(x) rgamma(I - 1, shape = x[1], rate = x[2]))
trace2 <- rbind(0, apply(jumps, 2, cumsum))

## mu ~ N(1, 0.5); nu ~ trunc_T(3)
time <- 0:11
I <- length(time)
N <- 100
mu <- rnorm(N, mean = 1, sd = 0.5)
nu <- abs(rt(N, df = 3))
beta <- 1/(nu^2); xi <- 1/(mu * nu^2)
shape_rate <- cbind(beta, xi)
jumps <- apply(shape_rate, 1, function(x) rgamma(I - 1, shape = x[1], rate = x[2]))
trace3 <- rbind(0, apply(jumps, 2, cumsum))

# generate plot
pdf(
  file.path("..", "..", "figures", "ch-4", "PPCs.pdf"),
  height = 4, width = 10
)
par(mfrow = c(1, 3))
par(cex = 1.1)
matplot(time, trace1, type = "l", lty = 1, col = scales::alpha("black", 0.2),
        xlab = "time (arbitrary)", ylab = "degradation (arbitrary)",
        ylim = c(0, 35))
text(0.1, 32.5, "(a)", adj = 0)
matplot(time, trace2, type = "l", lty = 1, col = scales::alpha("black", 0.2),
        xlab = "time (arbitrary)", ylab = "degradation (arbitrary)",
        ylim = c(0, 35))
text(0.1, 32.5, "(b)", adj = 0)
matplot(time, trace3, type = "l", lty = 1, col = scales::alpha("black", 0.2),
        xlab = "time (arbitrary)", ylab = "degradation (arbitrary)",
        ylim = c(0, 35))
text(0.1, 32.5, "(c)", adj = 0)
dev.off()