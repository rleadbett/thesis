library(rstan)
library(stringr)
library(posterior)
library(tidybayes)
library(ggplot2)
library(ggdist)
library(tidyr)

library(multidplyr)
library(dplyr, warn.conflicts = FALSE)

# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Functions to simulate FT distribution
ftDraw <- function(y_filt_current, t_current, nu, mu, phi, limit) {
  # starting from most recent spline coef
  y <- y_filt_current
  t <- t_current
  delta_t <- 0.001

  # average over variability along length of belt
  z <- lapply(
    1:20,
    function(itr) {
      y_noisy <- rTruncT(
        n = length(y), 
        mean = y, 
        nu = 50, 
        sd = phi * y, 
        lb = 0
      )
      # calculate values of spline at measurement locations
      z <- B %*% y_noisy
      return(z)
    }
  ) %>%
    unlist()
  
  

  if(max(z) < limit) {
    while(max(z) < limit){
      t <- t + delta_t
      delta_y <- rgamma(
        n = length(y),
        shape = (delta_t / nu^2),
        rate = 1 / (mu * nu^2)
      )
      y <- y + delta_y
      z <- lapply(
        1:20,
        function(itr) {
          y_noisy <- rTruncT(
            n = length(y), 
            mean = y, 
            nu = 10, 
            sd = phi * y, 
            lb = 0
          )
          # calculate values of spline at measurement locations
          z <- B %*% y_noisy
          return(z)
        }
      ) %>%
        unlist()
    }
    return(t)
  } else {
    return(t)
  }
}

ftCDFDraw <- function(n = 500, ...) {
  emp_cdf <- lapply(
    1:n,
    function(i) ftDraw(...)
  ) %>% unlist()

  return(
    emp_cdf
  )
}

rTruncT <- function(n, nu, mean, sd, lb) {  
  lb_norm <- (lb - mean) / sd
  unif_rv <- runif(n, pt(lb_norm, df = nu), 1)
  t_rv <- (qt(unif_rv, df = nu) * sd) + mean
  
  return(t_rv)
}

# Load data and model
stan_data <- readRDS(
  file = "stan-data.rds"
)
stan_model_gp <- readRDS(
  file = "compiled-stan-model-gp.rds"
)

# FT for eight observations

stan_data_8 <- list(
  I = 8,
  N = stan_data$N,
  M = stan_data$M,
  t = stan_data$t[1:8],
  z = stan_data$z[1:8, ],
  B = stan_data$B,
  a_hat = stan_data$a_hat,
  b_hat = stan_data$b_hat
)

stan_fit_8_gp <- sampling(
  stan_model_gp,
  stan_data_8,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(
    adapt_delta = 0.999, 
    max_treedepth = 14
  )
)

## set up cluster
cores <- parallel::detectCores()
cluster <- new_cluster(cores - 2)
cluster

cluster_library(cluster, "dplyr")
cluster_assign(
  cluster,
  ftDraw = ftDraw,
  ftCDFDraw = ftCDFDraw,
  rTruncT = rTruncT,
  B = stan_data$B
)

## thin number of samples to reduce the computation
subset_draws <- sample(1:(1000 * 4), size = 500)#size = 1000)

## Calculate FT
ft_dist_8_gp <- stan_fit_8_gp %>% 
  as_draws_rvars() %>%
  spread_draws(
    mu[m],
    nu[m],
    phi,
    y[i, m]
  ) %>%
  filter((i == (8 - 1)) & (`.draw` %in% subset_draws)) %>%
  mutate(
    t = stan_data$t[8],
    y = y * mu
  ) %>%
  group_by(`.draw`) %>%
  partition(cluster) %>%
  summarise(
    ft = ftCDFDraw(
      y_filt_current = y,
      t_current = t,
      mu = mu, 
      nu = nu,
      phi = unique(phi),
      limit = 25
    )
  ) %>%
  collect()

saveRDS(ft_dist_8_gp, "ft-8-obs-gp.rds")


## Linear general path model

stan_model_lm <- readRDS(
  file = "compiled-stan-model-lm.rds"
)

# FT for eight observations
stan_fit_8_lm <- sampling(
  stan_model_lm,
  stan_data_8,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(
    adapt_delta = 0.999, 
    max_treedepth = 14
  )
)

# Functions to simulate FT distribution
ftDraw_lm <- function(t_current, mu, phi, limit) {
  # starting from most recent spline coef
  t <- t_current
  y <- t * mu
  delta_t <- 0.001

  # average over variability along length of belt
  z <- lapply(
    1:20,
    function(itr) {
      y_noisy <- rTruncT(
        n = length(y), 
        mean = y, 
        nu = 10, 
        sd = phi * y, 
        lb = 0
      )
      # calculate values of spline at measurement locations
      z <- B %*% y_noisy
      return(z)
    }
  ) %>%
    unlist()
  
  if(max(z) < limit) {
    while(max(z) < limit){
      t <- t + delta_t
      y <- t * mu
      z <- lapply(
        1:20,
        function(itr) {
          y_noisy <- rTruncT(
            n = length(y), 
            mean = y, 
            nu = 10, 
            sd = phi * y, 
            lb = 0
          )
          # calculate values of spline at measurement locations
          z <- B %*% y_noisy
          return(z)
        }
      ) %>%
        unlist()
    }
    return(t)
  } else {
    return(t)
  }
}

ftCDFDraw_lm <- function(n = 500, ...) {
  emp_cdf <- lapply(
    1:n,
    function(i) ftDraw_lm(...)
  ) %>% unlist()

  return(emp_cdf)
}

cluster_assign(
  cluster,
  ftDraw_lm = ftDraw_lm,
  ftCDFDraw_lm = ftCDFDraw_lm
)

ft_dist_8_lm <- stan_fit_8_lm %>% 
  as_draws_rvars() %>%
  spread_draws(
    mu[m],
    phi
  ) %>%
  filter(`.draw` %in% subset_draws) %>%
  mutate(t = stan_data$t[8]) %>%
  group_by(`.draw`) %>%
  partition(cluster) %>%
  summarise(
    ft = ftCDFDraw_lm(
      t_current = t,
      mu = mu, 
      phi = unique(phi),
      limit = 25
    )
  ) %>%
  collect()

saveRDS(ft_dist_8_lm, "ft-8-obs-lm.rds")

# FT for 7 obse

stan_data_7 <- list(
  I = 7,
  N = stan_data$N,
  M = stan_data$M,
  t = stan_data$t[1:7],
  z = stan_data$z[1:7, ],
  B = stan_data$B,
  a_hat = stan_data$a_hat,
  b_hat = stan_data$b_hat
)

stan_fit_7_gp <- sampling(
  stan_model_gp,
  stan_data_7,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(
    adapt_delta = 0.999, 
    max_treedepth = 14
  )
)

stan_fit_7_lm <- sampling(
  stan_model_lm,
  stan_data_7,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(
    adapt_delta = 0.999, 
    max_treedepth = 14
  )
)

ft_dist_7_gp <- stan_fit_7_gp %>% 
  as_draws_rvars() %>%
  spread_draws(
    mu[m],
    nu[m],
    phi,
    y[i, m]
  ) %>%
  filter((i == (7 - 1)) & (`.draw` %in% subset_draws)) %>%
  mutate(
    t = stan_data$t[7],
    y = y * mu
  ) %>%
  group_by(`.draw`) %>%
  partition(cluster) %>%
  summarise(
    ft = ftCDFDraw(
      y_filt_current = y,
      t_current = t,
      mu = mu, 
      nu = nu,
      phi = unique(phi),
      limit = 25
    )
  ) %>%
  collect()

ft_dist_7_lm <- stan_fit_7_lm %>% 
  as_draws_rvars() %>%
  spread_draws(
    mu[m],
    phi
  ) %>%
  filter(`.draw` %in% subset_draws) %>%
  mutate(t = stan_data$t[7]) %>%
  group_by(`.draw`) %>%
  partition(cluster) %>%
  summarise(
    ft = ftCDFDraw_lm(
      t_current = t,
      mu = mu, 
      phi = unique(phi),
      limit = 25
    )
  ) %>%
  collect()

saveRDS(ft_dist_7_gp, "ft-7-obs-gp.rds")
saveRDS(ft_dist_7_lm, "ft-7-obs-lm.rds")


# FT for nine observations

stan_fit_9_gp <- sampling(
  stan_model_gp,
  stan_data,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(
    adapt_delta = 0.999, 
    max_treedepth = 14
  )
)


ft_dist_9_gp <- stan_fit_9_gp %>%
  as_draws_rvars() %>%
  spread_draws(
    mu[m],
    nu[m],
    phi,
    y[i, m]
  ) %>%
  filter((i == (9 - 1)) & (`.draw` %in% subset_draws)) %>%
  mutate(
    t = stan_data$t[9],
    y = y * mu
  ) %>%
  group_by(`.draw`) %>%
  partition(cluster) %>%
  summarise(
    ft = ftCDFDraw(
      y_filt_current = y,
      t_current = t,
      mu = mu, 
      nu = nu,
      phi = unique(phi),
      limit = 20
    )
  ) %>%
  collect()

saveRDS(ft_dist_9_gp, "ft-9-obs-gp.rds")


# All 9 observations
stan_fit_9_lm <- sampling(
  stan_model_lm,
  stan_data,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(
    adapt_delta = 0.999, 
    max_treedepth = 14
  )
)

ft_dist_9_lm <- stan_fit_9_lm %>% 
  as_draws_rvars() %>%
  spread_draws(
    mu[m],
    phi
  ) %>%
  filter(`.draw` %in% subset_draws) %>%
  mutate(t = stan_data$t[9]) %>%
  group_by(`.draw`) %>%
  partition(cluster) %>%
  summarise(
    ft = ftCDFDraw_lm(
      t_current = t,
      mu = mu, 
      phi = unique(phi),
      limit = 20
    )
  ) %>%
  collect()

saveRDS(ft_dist_9_lm, "ft-9-obs-lm.rds")