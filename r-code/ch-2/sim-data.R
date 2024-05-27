library(dplyr)
library(ggplot2)
library(rstan)
library(bayesplot)
library(posterior)
library(tidybayes)
library(kableExtra)
library(cowplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Function to simulate data
sim_data <- function(
  beta,
  eta,
  n_units,
  t_start,
  t_end
){
  weibull_05_quant <- qweibull(0.05, shape = beta, scale = eta)

  n_lifetimes_per_unit <- ceiling(t_end / weibull_05_quant) * 2

  small_sim <- data.frame(
      unit = factor(
        rep(1:n_units, each = n_lifetimes_per_unit),
        1:n_units
      ),
      lifetime = rweibull(
        n_units * n_lifetimes_per_unit,
        shape = beta,
        scale = eta
      )
    ) %>% 
    group_by(unit) %>%
    mutate(
      failure_time = cumsum(lifetime),
      install_time = lag(failure_time),
      lifetime_id = 1:n()
    ) %>%
    ungroup() %>%
    replace(is.na(.), 0) %>%
    filter(
      between(install_time, t_start, t_end) |
      between(failure_time, t_start, t_end)
    ) %>%
    mutate(
      int_censored = !between(install_time, t_start, t_end),
      right_censored = !between(failure_time, t_start, t_end)
    )

  small_sim$install_time[small_sim$int_censored] <- t_start
  small_sim$failure_time[small_sim$right_censored] <- t_end
  small_sim$observed_lifetime <- small_sim$failure_time - small_sim$install_time

  return(small_sim)
}

# Function for the likelihood of the censored weibull data
ll_weibull <- function(
  params,
  observed_lifetimes,
  right_cens_lifetimes,
  int_cens_lifetimes_lower,
  int_cens_lifetimes_upper
) {
  beta = params[1]
  eta = params[2]

  ll <- sum(log(dweibull(observed_lifetimes, shape = beta, scale = eta))) +
    sum(log(1 - pweibull(right_cens_lifetimes, shape = beta, scale = eta))) +
    sum(
      log(
        pweibull(int_cens_lifetimes_upper, shape = beta, scale = eta) -
          pweibull(int_cens_lifetimes_lower, shape = beta, scale = eta)
      )
    )
  
  return(-ll)
}
# Function for the likelihood of the censored weibull data
ll_weibull_rcense <- function(
  params,
  observed_lifetimes,
  right_cens_lifetimes,
  int_cens_lifetimes_lower,
  int_cens_lifetimes_upper
) {
  beta = params[1]
  eta = params[2]

  ll <- sum(log(dweibull(observed_lifetimes, shape = beta, scale = eta))) +
    sum(log(1 - pweibull(right_cens_lifetimes, shape = beta, scale = eta))) +
    sum(log(1 - pweibull(int_cens_lifetimes_lower, shape = beta, scale = eta)))
  
  return(-ll)
}
# Load Stan models
stan_model_non_informative <- stan_model(
  file = file.path("..", "..", "stan_models", "non_inf_weibull.stan")
)

stan_model_non_informative_alt <- stan_model(
  file = file.path("..", "..", "stan_models", "non_inf_weibull_alt.stan")
)

# Simulate and plot data
set.seed(154)
beta_truth <- 1.05
eta_truth <- 1000
sim_data_df <- sim_data(
  beta = beta_truth,
  eta = eta_truth,
  n_units = 100,
  t_start = 1500,
  t_end = 2800
)

p_sorted <- sim_data_df %>%
  mutate(censored = right_censored | int_censored) %>%
  arrange(censored, observed_lifetime) %>%
  mutate(index = (1:n()) / n()) %>%
  ggplot() +
  geom_point(aes(x = observed_lifetime, y = index, colour = censored)) +
  theme_minimal()

p_true_cdf <- ggplot() +
  geom_step(
    data = sim_data_df %>%
      arrange(lifetime) %>%
      mutate(CDF = (1:n())/n()),
    aes(x = lifetime, y = CDF),
    colour = "black",
    linetype = 2
  ) +
  geom_step(
    data = sim_data_df %>%
      arrange(observed_lifetime) %>%
      mutate(index = n():1) %>%
      filter(!(right_censored | int_censored)) %>%
      mutate(
        one_minus_q_hat = 1 - (1 / index),
        S_hat = cumprod(one_minus_q_hat),
        F_hat = 1 - S_hat
      ),
    aes(x = observed_lifetime, y = F_hat)
  ) +
  geom_function(
    fun = pweibull,
    args = list(
      shape = beta_truth,
      scale = eta_truth
    ),
    colour = "red"
  ) +
  theme_minimal()

p_sim_lifetime_overview <- cowplot::plot_grid(p_sorted, p_true_cdf, nrow = 1)
pdf(
  file.path("..", "..", "figures", "sim_data_desc.pdf")
)
p_sim_lifetime_overview
dev.off()

# Write first 10 rows of simulated data to latex style table
sim_data_df %>%
  mutate(
    censoring = ifelse(int_censored, "interval", NA),
    censoring = ifelse(right_censored, "right", censoring),
    lifetime = round(lifetime, 1),
    install_time = round(install_time, 1),
    failure_time = round(failure_time, 1),
    observed_lifetime = round(observed_lifetime, 1)
  ) %>%
  select(
    unit,
    `true lifetime` = lifetime,
    `install time` = install_time,
    `failure time` = failure_time,
    censoring,
    `observed lifetime` = observed_lifetime
  ) %>%
  head(10) %>%
  kbl(
    booktabs = T,
    format = "latex",
    caption = "Example of simulated censored lifetime data.",
    label = "sim_censored_units"
  ) %>%
  kable_styling(latex_options = "striped") %>% 
  save_kable(
    file = file.path("..", "..", "tables", "sim_data_head.tex"),
    keep_tex = TRUE
  )


# MLE estimate
fully_obs_lifetimes <- sim_data_df %>%
  filter(!(right_censored | int_censored)) %>%
  pull(observed_lifetime)
right_cens_lifetimes <- sim_data_df %>%
  filter(right_censored) %>%
  pull(observed_lifetime)
int_cens_lifetimes <- sim_data_df %>%
  filter(int_censored) %>%
  pull(observed_lifetime)

ll_fit <- optim(
  c(0.5, 1000),
  fn = function(params){
    ll_weibull(
      params,
      observed_lifetimes = fully_obs_lifetimes,
      right_cens_lifetimes = right_cens_lifetimes,
      int_cens_lifetimes_lower = int_cens_lifetimes,
      int_cens_lifetimes_upper = int_cens_lifetimes + 1500
    )
  }
)

ll_fit_rcense <- optim(
  c(0.5, 1000),
  fn = function(params){
    ll_weibull_rcense(
      params,
      observed_lifetimes = fully_obs_lifetimes,
      right_cens_lifetimes = right_cens_lifetimes,
      int_cens_lifetimes_lower = int_cens_lifetimes,
      int_cens_lifetimes_upper = int_cens_lifetimes + 1500
    )
  }
)


p_mle_v_emperical <- ggplot() +
  geom_step(
    data = sim_data_df %>%
      arrange(observed_lifetime) %>%
      mutate(index = n():1) %>%
      filter(!(right_censored | int_censored)) %>%
      mutate(
        one_minus_q_hat = 1 - (1 / index),
        S_hat = cumprod(one_minus_q_hat),
        F_hat = 1 - S_hat
      ),
    aes(x = observed_lifetime, y = F_hat)
  ) +
  geom_step(
    data = sim_data_df %>%
      arrange(observed_lifetime) %>%
      mutate(index = n():1) %>%
      mutate(
        one_minus_q_hat = 1 - (1 / index),
        S_hat = cumprod(one_minus_q_hat),
        F_hat = 1 - S_hat
      ),
    aes(x = observed_lifetime, y = F_hat),
    linetype = 2
  ) +
  geom_step(
    data = sim_data_df %>%
      mutate(
        observed_lifetime = ifelse(
          (right_censored | int_censored),
          2500,
          observed_lifetime
        )
      ) %>%
      arrange(observed_lifetime) %>%
      mutate(index = n():1) %>%
      mutate(
        one_minus_q_hat = 1 - (1 / index),
        S_hat = cumprod(one_minus_q_hat),
        F_hat = 1 - S_hat
      ),
    aes(x = observed_lifetime, y = F_hat),
    linetype = 2
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = beta_truth, scale = eta_truth),
    colour = "red"
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = ll_fit$par[1], scale = ll_fit$par[2]),
    colour = "blue"
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = ll_fit_rcense$par[1], scale = ll_fit_rcense$par[2]),
    colour = "blue",
    linetype = 2
  ) +
  theme_minimal()

p_mle_v_emperical

# Repeat simulations

sim_data_df <- lapply(
  1:1000,
  function(nsim) {
    data <- sim_data(
      beta = beta_truth,
      eta = eta_truth,
      n_units = 100,
      t_start = 1500,
      t_end = 2800
    )
    # prep data
    fully_obs_lifetimes <- data %>%
      filter(!(right_censored | int_censored)) %>%
      pull(observed_lifetime)
    right_cens_lifetimes <- data %>%
      filter(right_censored) %>%
      pull(observed_lifetime)
    int_cens_lifetimes <- data %>%
      filter(int_censored) %>%
      pull(observed_lifetime)
    # right and interval
    ll_fit <- optim(
      c(0.5, 1000),
      fn = function(params){
        ll_weibull(
          params,
          observed_lifetimes = fully_obs_lifetimes,
          right_cens_lifetimes = right_cens_lifetimes,
          int_cens_lifetimes_lower = int_cens_lifetimes,
          int_cens_lifetimes_upper = int_cens_lifetimes + 1500
        )
      }
    )
    # right
    ll_fit_rcense <- optim(
      c(0.5, 1000),
      fn = function(params){
        ll_weibull_rcense(
          params,
          observed_lifetimes = fully_obs_lifetimes,
          right_cens_lifetimes = right_cens_lifetimes,
          int_cens_lifetimes_lower = int_cens_lifetimes,
          int_cens_lifetimes_upper = int_cens_lifetimes + 1500
        )
      }
    )
    grid <- seq(0.001, 0.999, length.out = 50)
    quant_grid <- qweibull(grid, beta_truth, eta_truth)
    df <- data.frame(
      nsim = nsim,
      quantile = grid,
      beta_est = rep(
        c(ll_fit_rcense$par[1], ll_fit$par[1]),
        each = length(quant_grid)
      ),
      eta_est = rep(
        c(ll_fit_rcense$par[2], ll_fit$par[2]),
        each = length(quant_grid)
      ),
      cdf_est = c(
        pweibull(
          quant_grid,
          ll_fit_rcense$par[1],
          ll_fit_rcense$par[2]
        ),
        pweibull(
          quant_grid,
          ll_fit$par[1],
          ll_fit$par[2]
        ) 
      ),
      model = rep(
        c("right", "both"),
        each = length(quant_grid)
      )
    ) %>%
    mutate(
      cdf_dev = cdf_est - quantile
    )
  }
)

p_mle_bias_cdf <- sim_data_df %>%
  bind_rows() %>%
  mutate(
    nsim = factor(nsim),
    model = factor(model)
  ) %>%
  ggplot(
    aes(
      x = quantile,
      y = cdf_dev,
      group = interaction(nsim, model),
      colour = model
    )
  ) +
  geom_line(
    alpha = 0.5
  ) +
  geom_hline(yintercept = 0) +
  ylim(-0.45, 0.15) +
  theme_minimal()

p_mle_bias_pars <- sim_data_df %>%
  bind_rows() %>%
  select(nsim, beta_est, eta_est, model) %>%
  unique() %>%
  ggplot() +
  geom_point(aes(x = beta_est, y = eta_est, colour = model), alpha = 0.4)+
  geom_point(x = beta_truth, y = eta_truth, colour = "black") +
  xlim(0.7, 2.2) +
  ylim(750, 3750) +
  theme_minimal()

plot_grid(p_mle_bias_pars, p_mle_bias_cdf, nrow = 1)

# Add Bayesian estimates
stan_data <- list(
  N_obs = length(fully_obs_lifetimes),
  N_Icens = length(right_cens_lifetimes),
  N_Rcens = length(int_cens_lifetimes),
  lifetime_obs = fully_obs_lifetimes,
  lifetime_Rcens = right_cens_lifetimes,
  lifetime_Icens_Upper = int_cens_lifetimes + 1500,
  lifetime_Icens_Lower = int_cens_lifetimes
)

stan_fit_non_informative <- sampling(
  stan_model_non_informative,
  stan_data,
  cores = 4,
  iter = 1000,
  warmup = 300
)

p_nonif_post <- mcmc_scatter(
  stan_fit_non_informative,
  pars = c("beta", "eta")
) + theme_minimal()

stan_fit_non_informative_alt <- sampling(
  stan_model_non_informative_alt,
  stan_data,
  cores = 4,
  iter = 1000,
  warmup = 300
)

p_nonif_post_alt <- mcmc_scatter(
  stan_fit_non_informative_alt,
  pars = c("beta", "eta")
) + theme_minimal()

cowplot::plot_grid(p_nonif_post, p_nonif_post_alt, nrow = 1)

# Comparison of parameter estimates
## Parameter joint dist
p_joint_fits <- stan_fit_non_informative_alt %>%
  as_draws_df() %>%
  spread_draws(beta, eta) %>%
  ggplot() +
  stat_density_2d(aes(x = beta, y = eta, fill = ..level..), geom = "polygon") +
  geom_point(
    data = data.frame(
      beta = ll_fit$par[1],
      eta = ll_fit$par[2]
    ),
    aes(x = beta, y = eta),
    colour = "green",
    shape = 17,
    size = 3
  ) +
  geom_point(
    data = data.frame(
      beta = beta_truth,
      eta = eta_truth
    ),
    aes(x = beta, y = eta),
    colour = "red",
    shape = 17,
    size = 3
  ) +
  theme_minimal() +
  theme(legend.position = "none")
## CDFs
p_cdf_fits <- ggplot() +
  geom_step(
    data = sim_data_df %>%
      arrange(lifetime) %>%
      mutate(CDF = (1:n())/n()),
    aes(x = lifetime, y = CDF),
    colour = "black",
    linetype = 2
  ) +
  geom_function(
    fun = pweibull,
    args = list(
      shape = beta_truth,
      scale = eta_truth
    ),
    colour = "red"
  ) +
  geom_function(
    fun = pweibull,
    args = list(
      shape = summary(stan_fit_non_informative)$summary["beta", "50%"],
      scale = summary(stan_fit_non_informative)$summary["eta", "50%"]
    ),
    colour = "blue"
  ) +
  geom_function(
    fun = pweibull,
    args = list(
      shape = ll_fit$par[1],
      scale = ll_fit$par[2]
    ),
    colour = "green"
  ) +
  theme_minimal()

cowplot::plot_grid(p_joint_fits, p_cdf_fits, nrow = 1)

p_mle_v_emperical +
  geom_function(
    fun = pweibull,
    args = list(
      shape = summary(stan_fit_non_informative)$summary["beta", "50%"],
      scale = summary(stan_fit_non_informative)$summary["eta", "50%"]
    ),
    colour = "green"
  ) +
  geom_function(
    fun = pweibull,
    args = list(
      shape = summary(stan_fit_non_informative_alt)$summary["beta", "50%"],
      scale = summary(stan_fit_non_informative_alt)$summary["eta", "50%"]
    ),
    colour = "green",
    linetype = 2
  )

mcmc_intervals(
  stan_fit_non_informative_alt,
  regex_pars = c("Y_Rcens")
)

mcmc_intervals(
  stan_fit_non_informative_alt,
  regex_pars = c("Y_Icens")
)

# Next plot the true ordered plot with the censoring intervals and Bayesian predictions for the missing data.

# -> plot the true lifetime values with censored observations in Red
# -> plot the censoring intervals as horizontal black lines
# -> plot the inferred values for the "missing" data as interval plots

stan_fit_non_informative_alt %>%
  as_draws_rvars() %>%
  gather_rvars(Y_Rcens[i], Y_Icens[i]) %>%
  mutate(t_hat = E(.value)) %>%
  select(
    observed_lifetime = t_hat,
    dist = .value
  ) %>%
  mutate(
    type = "pred"
  ) %>%
  rbind(
    sim_data_df %>%
      filter((!right_censored) & (!int_censored)) %>%
      select(observed_lifetime) %>%
      mutate(
        type = "observed",
        dist = rvar(
          array(
            rep(observed_lifetime, each = (700 * 4)), 
            dim = c(700, 4, n())
          ),
          with_chains = TRUE
        )
      )
  ) %>%
  arrange(observed_lifetime) %>%
  mutate(index = n():1) %>%
  mutate(
    one_minus_q_hat = 1 - (1 / index),
    S_hat = cumprod(one_minus_q_hat),
    F_hat = 1 - S_hat,
    ecdf = (1:n()) / n()
  ) %>%
  ggplot(aes()) +
  #geom_step(aes(y = F_hat)) +
  stat_pointinterval(
    aes(x = observed_lifetime, y = ecdf, xdist = dist),
    point_interval = "mean_qi"
  ) +
  geom_point(
    aes(x = observed_lifetime, y = ecdf, colour = type)
  ) +
  geom_step(
    data = sim_data_df %>%
      arrange(lifetime) %>%
      mutate(ecdf = 1:n() / n()),
    aes(x = lifetime, y = ecdf),
    colour = "black"
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = beta_truth, scale = eta_truth),
    colour = "red"
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = ll_fit$par[1], scale = ll_fit$par[2]),
    colour = "blue"
  ) +
  theme_minimal()


stan_fit_non_informative_alt %>%
  as_draws_df() %>%
  gather_draws(Y_Rcens[i], Y_Icens[i]) %>%
  select(
    draw = .draw,
    observed_lifetime = .value
  ) %>%
  rbind(
    data.frame(
      observed_lifetime = sim_data_df %>%
        filter((!right_censored) & (!int_censored)) %>%
        pull(observed_lifetime) %>%
        rep(., each = (700 * 4)),
      draw = 1:(700 * 4)
    )
  ) %>% 
  group_by(draw) %>%
  arrange(observed_lifetime) %>%
  mutate(ecdf = 1:n() / n()) %>%
  ungroup() %>%
  filter(draw %in% sample(1:(700 * 4), 200)) %>%
  ggplot() +
  geom_step(
    aes(x = observed_lifetime, y = ecdf, group = draw),
    colour = "gray", alpha = 0.4
  ) +
  geom_step(
    data = sim_data_df %>%
      arrange(lifetime) %>%
      mutate(ecdf = 1:n() / n()),
    aes(x = lifetime, y = ecdf),
    colour = "black"
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = beta_truth, scale = eta_truth),
    colour = "red"
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = ll_fit$par[1], scale = ll_fit$par[2]),
    colour = "blue"
  ) +
  theme_minimal()

# Small sim for bayesian model

bayes_sim_data_df <- lapply(
  1:100,
  function(nsim) {
    # sim data
    data <- sim_data(
      beta = beta_truth,
      eta = eta_truth,
      n_units = 100,
      t_start = 1500,
      t_end = 2800
    )
    # prep data
    fully_obs_lifetimes <- data %>%
      filter(!(right_censored | int_censored)) %>%
      pull(observed_lifetime)
    right_cens_lifetimes <- data %>%
      filter(right_censored) %>%
      pull(observed_lifetime)
    int_cens_lifetimes <- data %>%
      filter(int_censored) %>%
      pull(observed_lifetime)
    stan_data <- list(
      N_obs = length(fully_obs_lifetimes),
      N_Icens = length(right_cens_lifetimes),
      N_Rcens = length(int_cens_lifetimes),
      lifetime_obs = fully_obs_lifetimes,
      lifetime_Rcens = right_cens_lifetimes,
      lifetime_Icens_Upper = int_cens_lifetimes + 1500,
      lifetime_Icens_Lower = int_cens_lifetimes
    )
    # fit stan model
    repeat {
      #print(nsim)
      stanfit <- tryCatch({
        sampling(
          stan_model_non_informative_alt,
          stan_data,
          cores = 4,
          iter = 700,
          warmup = 300,
          refresh = 0
        )
      },
      error = function(e) {
        NULL
      })
      if (!is.null(stanfit)) {
        break
      }
      Sys.sleep(1)
    }

    beta_est <- summary(stanfit)$summary["beta", "50%"]
    eta_est <- summary(stanfit)$summary["eta", "50%"]
    # calculate difference from true cdf
    grid <- seq(0.001, 0.999, length.out = 50)
    quant_grid <- qweibull(grid, beta_truth, eta_truth)
    df <- data.frame(
      nsim = nsim,
      quantile = grid,
      beta_est = beta_est,
      eta_est = eta_est,
      cdf_est = pweibull(
        quant_grid,
        beta_est,
        eta_est
      )
    ) %>%
    mutate(
      cdf_dev = cdf_est - quantile
    )
    
    return(df)
  }
)

p_bayes_bias_cdf <- bayes_sim_data_df  %>%
  bind_rows() %>%
  mutate(
    nsim = factor(nsim)
  ) %>%
  ggplot(
    aes(
      x = quantile,
      y = cdf_dev,
      group = nsim
    )
  ) +
  geom_line(
    alpha = 0.5,
    colour = "green"
  ) +
  geom_hline(yintercept = 0) +
  ylim(-0.45, 0.15) +
  theme_minimal()

p_bayes_bias_pars <- bayes_sim_data_df %>%
  bind_rows() %>%
  select(nsim, beta_est, eta_est) %>%
  unique() %>%
  ggplot() +
  geom_point(
    aes(x = beta_est, y = eta_est),
    alpha = 0.4,
    colour = "green"
  ) +
  geom_point(x = beta_truth, y = eta_truth, colour = "black") +
  xlim(0.7, 2.2) +
  ylim(750, 3750) +
  theme_minimal()

plot_grid(p_bayes_bias_pars, p_bayes_bias_cdf, nrow = 1)

# Informative Bayesian model
## Call stan model
stan_model_joint_alt <- stan_model(
  file = file.path("..", "..", "stan_models", "joint_weibull_alt.stan")
)
# Support function for PlotJointPrior used when calculating Weibull params
fn <- function(tCDF) log(-log(1 - tCDF))
# Calculate the beta parameters from mean and variance
GetBetaPars <- function(t_mean, t_var){
  a <- ((t_mean ^ 2) * (1 - t_mean) / t_var) - t_mean;
  b <- (a / t_mean) - a;

  return(c(a, b)) 
}
# Generate samples from lower bound truncated beta distribution
rTruncBeta <- function(n, a, b, lb) {
  lb_trans <- pbeta(lb, a, b)
  unif_samples <- runif(
    n = n,
    min = lb_trans,
    # Note: Numerical issues if unif close to 1 and beta has large mass at zero
    max = 0.99
  )
  trunc_beta_samples <- qbeta(unif_samples, a, b)

  return(trunc_beta_samples)
}
# pweibull that supports rvar objects
pweibullRvar <- rfun(pweibull)

PlotJointPrior <- function(
  t1,
  t2,
  t1_mean,
  t1_sd,
  t2_mean,
  t2_sd,
  n = 1000
) {
  if (t2 < t1) stop("t2 must be greater than t1")
  if (t2_mean < t1_mean) stop("CDF must be monotonic increasing")
  
  # Calculate the parameters of the two beta distributions
  pars_t1 <- GetBetaPars(t1_mean, t1_sd ^ 2)
  pars_t2 <- GetBetaPars(t2_mean, t2_sd ^ 2)

  # Sample pairs of points along CDF
  t1CDF <- rbeta(n, pars_t1[1], pars_t1[2])
  t2CDF <- rTruncBeta(
    n = n,
    a = pars_t2[1],
    b = pars_t2[2],
    lb = t1CDF
  )

  # Calculate Weibull params from draws
  beta <- (fn(t2CDF) - fn(t1CDF)) / log(t2 / t1)
  eta <- exp(log(t1) - (fn(t1CDF) / beta))

  rvar_beta <- rvar(beta)
  rvar_eta <- rvar(eta)
  
  grid <- qweibull(
    seq(0.0001, 0.9999, 0.01),
    shape = median(beta),
    scale = median(eta)
  )
  
  p <- data.frame(
    exposure = grid,
    cdf = pweibullRvar(
      grid,
      shape = rvar_beta,
      scale = rvar_eta
    )
  ) %>%
  ggplot() +
  stat_dist_lineribbon(aes(x = exposure, dist = cdf)) +
  scale_fill_brewer() +
  ylim(0, 1) +
  theme_minimal()

  return(p)
}

PlotJointPrior(
  t1 = 2000,
  t2 = 5000,
  t1_mean = 0.875,
  t1_sd = 0.01,
  t2_mean = 0.995,
  t2_sd = 0.0025
) +
geom_step(
  data = sim_data_df %>%
    arrange(lifetime) %>%
    mutate(ecdf = 1:n() / n()),
  aes(x = lifetime, y = ecdf),
  colour = "black"
) +
geom_function(
  fun = pweibull,
  args = list(shape = beta_truth, scale = eta_truth),
  colour = "red"
)

stan_data_joint_comp <- list(
  t_1 = 2000,
  t_2 = 5000,
  t1_mean = 0.875,
  t1_var = 0.01 ^ 2,
  t2_mean = 0.995,
  t2_var = 0.0025 ^ 2
)

stan_fit_joint_alt <- sampling(
  stan_model_joint_alt,
  c(
    stan_data,
    stan_data_joint_comp
  ),
  cores = 4,
  iter = 2000,
  warmup = 300,
  control = list(
      adapt_delta = 0.99,
      max_treedepth = 14
  )
)

stan_fit_joint_alt %>%
  as_draws_df() %>%
  gather_draws(Y_Rcens[i], Y_Icens[i]) %>%
  select(
    draw = .draw,
    observed_lifetime = .value
  ) %>%
  rbind(
    data.frame(
      observed_lifetime = sim_data_df %>%
        filter((!right_censored) & (!int_censored)) %>%
        pull(observed_lifetime) %>%
        rep(., each = ((2000 - 300) * 4)),
      draw = 1:((2000 - 300) * 4)
    )
  ) %>%
  group_by(draw) %>%
  arrange(observed_lifetime) %>%
  mutate(ecdf = 1:n() / n()) %>%
  ungroup() %>%
  filter(draw %in% sample(1:((2000 - 300) * 4), 200)) %>%
  ggplot() +
  geom_step(
    aes(x = observed_lifetime, y = ecdf, group = draw),
    colour = "gray", alpha = 0.4
  ) +
  geom_step(
    data = sim_data_df %>%
      arrange(lifetime) %>%
      mutate(ecdf = 1:n() / n()),
    aes(x = lifetime, y = ecdf),
    colour = "black"
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = beta_truth, scale = eta_truth),
    colour = "red"
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = ll_fit$par[1], scale = ll_fit$par[2]),
    colour = "blue"
  ) +
  geom_function(
    fun = pweibull,
    args = list(
      shape = summary(stan_fit_joint_alt)$summary["beta", "50%"],
      scale = summary(stan_fit_joint_alt)$summary["eta", "50%"]
    ),
    colour = "green"
  ) +
  theme_minimal()


stan_fit_joint_alt %>%
  as_draws_rvars() %>%
  gather_rvars(Y_Rcens[i], Y_Icens[i]) %>%
  mutate(t_hat = E(.value)) %>%
  select(
    observed_lifetime = t_hat,
    dist = .value
  ) %>%
  mutate(
    type = "pred"
  ) %>%
  rbind(
    sim_data_df %>%
      filter((!right_censored) & (!int_censored)) %>%
      select(observed_lifetime) %>%
      mutate(
        type = "observed",
        dist = rvar(
          array(
            rep(observed_lifetime, each = ((2000 - 300) * 4)), 
            dim = c((2000 - 300), 4, n())
          ),
          with_chains = TRUE
        )
      )
  ) %>%
  arrange(observed_lifetime) %>%
  mutate(index = n():1) %>%
  mutate(
    one_minus_q_hat = 1 - (1 / index),
    S_hat = cumprod(one_minus_q_hat),
    F_hat = 1 - S_hat,
    ecdf = (1:n()) / n()
  ) %>%
  ggplot(aes()) +
  #geom_step(aes(y = F_hat)) +
  stat_pointinterval(
    aes(x = observed_lifetime, y = ecdf, xdist = dist),
    point_interval = "mean_qi"
  ) +
  geom_point(
    aes(x = observed_lifetime, y = ecdf, colour = type)
  ) +
  geom_step(
    data = sim_data_df %>%
      arrange(lifetime) %>%
      mutate(ecdf = 1:n() / n()),
    aes(x = lifetime, y = ecdf),
    colour = "black"
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = beta_truth, scale = eta_truth),
    colour = "red"
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = ll_fit$par[1], scale = ll_fit$par[2]),
    colour = "blue"
  ) +
  geom_function(
    fun = pweibull,
    args = list(
      shape = summary(stan_fit_joint_alt)$summary["beta", "50%"],
      scale = summary(stan_fit_joint_alt)$summary["eta", "50%"]
    ),
    colour = "green"
  ) +
  theme_minimal()

mcmc_scatter(
  stan_fit_joint_alt,
  pars = c("beta", "eta")
) +
geom_point(
  data = data.frame(
    beta = ll_fit$par[1],
    eta = ll_fit$par[2]
  ),
  aes(x = beta, y = eta),
  colour = "green",
  shape = 17,
  size = 3
) +
geom_point(
  data = data.frame(
    beta = beta_truth,
    eta = eta_truth
  ),
  aes(x = beta, y = eta),
  colour = "red",
  shape = 17,
  size = 3
) +
theme_minimal()


stan_fit_joint_alt %>%
  as_draws_df() %>%
  spread_draws(beta, eta) %>%
  ggplot() +
  stat_density_2d(aes(x = beta, y = eta, fill = ..level..), geom = "polygon") +
  geom_point(
    data = data.frame(
      beta = ll_fit$par[1],
      eta = ll_fit$par[2]
    ),
    aes(x = beta, y = eta),
    colour = "green",
    shape = 17,
    size = 3
  ) +
  geom_point(
    data = data.frame(
      beta = beta_truth,
      eta = eta_truth
    ),
    aes(x = beta, y = eta),
    colour = "red",
    shape = 17,
    size = 3
  ) +
  theme_minimal() +
  theme(legend.position = "none")
