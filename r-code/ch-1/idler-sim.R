# A function to generate a fictitious simulated data set.
library(rstan)
library(posterior)
library(tidyr)
library(splines2)
library(tidyr)
library(ggplot2)
library(ggdist)
library(cowplot)
library(tidybayes)
library(stringr)
library(ggpubr)
library(wesanderson)
library(dplyr)
library(splines2)
library(RColorBrewer)
library(animation)
library(stringr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Plot theme
background_col <- "#fcfbf9"
my_theme <- function() {
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = background_col), #transparent panel bg
    plot.background = element_rect(fill = background_col, color = NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill = background_col), #transparent legend bg
    legend.box.background = element_rect(fill = background_col) #transparent legend panel
  )
}

GenerateCensoredData <- function(
  t_end,
  t_start = 0,
  replacement = TRUE,
  n_units = 100,
  n_samples_per_unit = 99,
  beta,
  eta
) {

  if (!replacement) {
    if ((t_start != 0) || !(n_samples_per_unit %in% c(1, 99))) {
      message("
        !For replacement = FALSE the start time must be zero and number 
        of samples per unit must be one. Changing t_start to zero and 
        n_samples_per_unit to one.
      ")
    }
    t_start  <- 0
    n_samples_per_unit <- 1
  }

  # Sample lifetimes from Weibull distribution and arrange in a matrix.
  sample_size <- n_units * n_samples_per_unit

  samples <- rweibull(sample_size, shape = beta, scale = eta)

  unit_lifetimes <- matrix(
    samples,
    nrow = n_units
  )

  # Calculate the replacement times.
  unit_replacement_times <- apply(
    unit_lifetimes,
    1,
    cumsum
  ) %>%
    t()

  # Censor the observations according to the specified parameters.
  if (!replacement) {
    censored_logical <- unit_replacement_times > t_end
    unit_replacement_times[censored_logical] <- t_end

    replacements_df <- data.frame(
      install_time = 0,
      replacement_time = unit_replacement_times,
      censoring = ifelse(censored_logical, "right", NA),
      unit = 1:n_units
    )
  } else {

    if (
      sum(unit_replacement_times[, ncol(unit_replacement_times)] < t_end) >= 1
    ) {
      message("
        !For the inputs you have selected the repeated replacements of some of 
        the units stops before t_end.
      ")
    }

    replacements_df <- lapply(
      1:n_units,
      function(i) {
        CensorUnitReplacements(
          unit_replacement_times[i, ],
          t_start,
          t_end
        ) %>% mutate(unit = i)
      }
    ) %>%
      bind_rows()
  }

  # Return a data frame with columns for; the install time of the unit, the
  # replacement times, if the lifetime is right or interval censored, and which
  # unit.
  return(replacements_df)

}

CensorUnitReplacements <- function(
  unit_replacements_vec,
  t_start,
  t_end
) {

  in_window_logical <- between(
    unit_replacements_vec,
    left = t_start,
    right = t_end
  )

  install_times  <- lag(unit_replacements_vec[in_window_logical])
  replacement_times  <- unit_replacements_vec[in_window_logical]
  install_times[1] <- t_start

  if (length(replacement_times) == 0) {
    last_life_included <- FALSE
  } else {
    last_life_included <- tail(replacement_times, 1) ==
      tail(unit_replacements_vec, 1)
  }

  if (!last_life_included) {
    install_times <- append(
      install_times,
      replacement_times[length(replacement_times)]
    )
    replacement_times <- append(
      replacement_times,
      t_end
    )
  }

  if (length(replacement_times) == 1) {
    df <- data.frame(
      install_time = install_times,
      replacement_time = replacement_times,
      censoring = "right"
    )
  } else {
    df <- data.frame(
      install_time = install_times,
      replacement_time = replacement_times,
      censoring = c(
        ifelse(t_start == 0, NA, "int"),
        rep(NA, (length(install_times) - 2)),
        ifelse(last_life_included, NA, "right")
      )
    )
  }

  return(df)
}

set.seed(2468)
#set.seed(12)
beta_truth <-1.1
eta_truth <- 1100
fake_data <- GenerateCensoredData(
  t_end = 4000,
  t_start = 2200,
  replacement = TRUE,
  n_units = 140,
  n_samples_per_unit = 99,
  beta = beta_truth,
  eta = eta_truth
) %>%
  mutate(lifetime = replacement_time - install_time)

fake_data %>%
  mutate(
    censored = ifelse(is.na(censoring), FALSE, TRUE)
  ) %>%
  group_by(censored) %>%
  arrange(censored, lifetime) %>%
  ungroup() %>%
  mutate(ID = 1:n()) %>% 
  ggplot(aes(x = lifetime, y = ID, col = censored)) +
  geom_point() +
  my_theme()

# KM plot
ecdf_upper <- fake_data %>%
  arrange(lifetime) %>%
  mutate(
    id = n():1,
    q_hat = 1 / id,
    s_hat = cumprod(1 - q_hat),
    ecdf = 1 - s_hat
  )

ecdf_est <- fake_data %>%
  arrange(lifetime) %>%
  mutate(
    id = n():1
  ) %>%
  filter(is.na(censoring)) %>%
  mutate(
    q_hat = 1 / id,
    s_hat = cumprod(1 - q_hat),
    ecdf = 1 - s_hat
  )

ecdf_lower <- fake_data
ecdf_lower$censoring[is.na(ecdf_lower$censoring)] <- ""
ecdf_lower$lifetime[ecdf_lower$censoring == "right"] <- NA
ecdf_lower$lifetime[ecdf_lower$censoring == "int"] <- NA

ecdf_lower <- ecdf_lower %>%
  arrange(lifetime) %>%
  mutate(
    id = n():1
  ) %>%
  mutate(
    q_hat = 1 / id,
    s_hat = cumprod(1 - q_hat),
    ecdf = 1 - s_hat
  )

mc_sim <- data.frame(
  lifetime = rweibull(374 * 1000, shape = beta_truth, scale = eta_truth),
  itr = rep(1:1000, each = 374)
) %>%
  arrange(itr, lifetime) %>%
  group_by(itr) %>%
  mutate(
    id = n():1,
    q_hat = 1 / id,
    s_hat = cumprod(1 - q_hat),
    ecdf = 1 - s_hat
  )

ll_weibull <- function(params) {
  beta = params[1]
  eta = params[2]
  obs_lifetimes <- fake_data$lifetime[is.na(fake_data$censoring)]
  cense_lifetimes <- fake_data$lifetime[!is.na(fake_data$censoring)]

  ll <- sum(log(dweibull(obs_lifetimes, shape = beta, scale = eta))) +
    sum(log(1 - pweibull(cense_lifetimes, shape = beta, scale = eta)))
  
  return(-ll)
}

ll_fit <- optim(c(0.5, 1000), ll_weibull)

ecdf_est %>%
  ggplot(aes(x = lifetime, y = ecdf)) +
  geom_step(data = mc_sim, aes(group = itr), colour = "gray", alpha = 0.2) +
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
    args = list(shape = 1.42, scale = 1253.32),
    colour = "green"
  ) +
  geom_step() +
  geom_step(data = ecdf_upper, linetype = 2) +
  geom_step(data = ecdf_lower, linetype = 2) +
  ylim(0, 1) +
  xlim(0, 1800) +
  my_theme()

naive_stan_model <- stan_model(
  model_code = "
  data {
    int N_obs;
    int N_Rcens;
    int N_Icens;
    real lifetime_obs[N_obs];
    real lifetime_Rcens[N_Rcens];
    real lifetime_Icens_Upper[N_Icens];
    real lifetime_Icens_Lower[N_Icens];
  }
  parameters {
    real<lower = 0> beta;
    real<lower = 0> eta;
  }
  model{

  // Likelihood
  // non-censored portion
  for(i in 1:N_obs){
    target += weibull_lpdf(lifetime_obs[i]|beta, eta);
  }
  // right censored portion
  for(j in 1:N_Rcens){
    target += weibull_lccdf(lifetime_Rcens[j]|beta, eta);
  }
  // interval portion
  for(k in 1:N_Icens){
    target += log_diff_exp(
      weibull_lcdf(lifetime_Icens_Upper[k]|beta, eta),
      weibull_lcdf(lifetime_Icens_Lower[k]|beta, eta)
    );
  }
    
  // Prior models
  eta ~ normal(1100, 100);
  beta ~ normal(1.1, 0.08);
  }
  ")

fake_data_obs <- fake_data %>%
  filter(is.na(censoring))
fake_data_Rcens <- fake_data %>%
  filter(censoring == "right")
fake_data_Icens <- fake_data %>%
  filter(censoring == "int")

stan_data <- list(
  N_obs = nrow(fake_data_obs),
  N_Rcens = nrow(fake_data_Rcens),
  N_Icens = nrow(fake_data_Icens),
  lifetime_obs = fake_data_obs$lifetime,
  lifetime_Rcens = fake_data_Rcens$lifetime,
  lifetime_Icens_Upper = fake_data_Icens$replacement_time,
  lifetime_Icens_Lower = fake_data_Icens$lifetime
)

naive_stan_fit <- sampling(
  naive_stan_model,
  stan_data,
  chains = 4,
  iter = 2000,
  warmup = 500,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 14
  )
)

marginal_post <- naive_stan_fit %>%
  as_draws_rvars() %>%
  spread_rvars(beta, eta)

res <- 100
data.frame(
  t = seq(1, 3000, length.out = res),
  beta = rep(marginal_post$beta, res),
  eta = rep(marginal_post$eta, res),
  beta_true = beta_truth,
  eta_true = eta_truth
) %>%
  mutate(
    hazard = (beta / eta) * ((t / eta)^(beta - 1)),
    hazard_true = (beta_true / eta_true) * ((t / eta_true)^(beta_true - 1))
  ) %>%
  ggplot(aes(x = t, dist = hazard)) +
  stat_dist_lineribbon() +
  scale_fill_brewer() +
  geom_line(aes(y = hazard_true), color = "red") +
  ylab("hazard") +
  my_theme()

rvarPweibull <- rfun(stats::pweibull)

data.frame(
  t = seq(1, 6000, length.out = res),
  beta = rep(marginal_post$beta, res),
  eta = rep(marginal_post$eta, res),
  beta_true = beta_truth,
  eta_true = eta_truth
) %>%
  mutate(
    cdf = rvarPweibull(t, shape = beta, scale = eta),
    cdf_true = stats::pweibull(t, shape = beta_true, scale = eta_true)
  ) %>%
  ggplot(aes(x = t, dist = cdf)) +
  stat_dist_lineribbon() +
  scale_fill_brewer() +
  geom_line(aes(y = cdf_true), color = "red") +
  ylab("cdf") +
  my_theme() 


####
data.frame(
  lifetime = seq(1, 6000, length.out = res),
  beta = rep(marginal_post$beta, res),
  eta = rep(marginal_post$eta, res),
  beta_true = beta_truth,
  eta_true = eta_truth
) %>%
  mutate(
    ecdf = rvarPweibull(lifetime, shape = beta, scale = eta),
    cdf_true = stats::pweibull(lifetime, shape = beta_true, scale = eta_true)
  ) %>%
  ggplot(aes(x = lifetime, dist = ecdf)) +
  geom_step(data = mc_sim, aes(y = ecdf, group = itr), colour = "gray", alpha = 0.5) +
  geom_function(
    fun = pweibull,
    args = list(shape = beta_truth, scale = eta_truth),
    colour = "red"
  ) +
  stat_dist_lineribbon(alpha = 0.2) +
  scale_fill_brewer() +
  geom_function(
    fun = pweibull,
    args = list(shape = ll_fit$par[1], scale = ll_fit$par[2]),
    colour = "blue"
  ) +
  geom_step(data = ecdf_est, aes(y = ecdf)) +
  geom_step(data = ecdf_upper, linetype = 2, aes(y = ecdf)) +
  geom_step(data = ecdf_lower, linetype = 2, aes(y = ecdf)) +
  ylim(0, 1) +
  xlim(0, 1800) +
  my_theme()

####
data.frame(
  beta = rnorm(5000, 1.1, 0.08),
  eta = rnorm(5000, 1100, 100)
) %>%
  ggplot(aes(x = beta, y = eta)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  geom_point(x = beta_truth, y = eta_truth, col = "red") +
  xlim(0.5, 1.5) +
  ylim(800, 1500) +
  my_theme() 

naive_stan_fit %>%
  as_draws_df() %>%
  spread_draws(beta, eta) %>%
  ggplot(aes(x = beta, y = eta)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  geom_point(x = beta_truth, y = eta_truth, col = "red") +
  xlim(0.5, 1.5) +
  ylim(800, 1500) +
  my_theme() 

informative_stan_model <- stan_model(
  model_code = "
functions {
  real fn(real tCDF) {
    return log(-log1m(tCDF));
  }
}

data {
  int N_obs;
  int N_Icens;
  int N_Rcens;
  real lifetime_obs[N_obs];
  real lifetime_Rcens[N_Rcens];
  real lifetime_Icens_Upper[N_Icens];
  real lifetime_Icens_Lower[N_Icens];
  real t1;   // should be 3.82
  real t2;   // should be 15
  real t1cdf_a;
  real t1cdf_b;
  real t2cdf_a;
  real t2cdf_b;
}

parameters {
  real<lower = 0> t1CDF;
  real<lower = t1CDF> t2CDF;  // the CDF at t2 must be greater than at t1
}

transformed parameters {
  real<lower = 0> beta;
  real<lower = 0> eta;

  // calculate Weibull paramaters based on the
  // draws from the CDF at t1 and t2.
  beta = (fn(t2CDF) - fn(t1CDF)) / log(t2 / t1);
  eta = exp(log(t1) - (fn(t1CDF) / beta));
}

model {
  // Likelihood
  // non-censored portion
  for(i in 1:N_obs){
    target += weibull_lpdf(lifetime_obs[i]|beta, eta);
  }
  // censored portion
  for(j in 1:N_Rcens){
    target += weibull_lccdf(lifetime_Rcens[j]|beta, eta);
  }
  // interval portion
  for(k in 1:N_Icens){
    target += log_diff_exp(
      weibull_lcdf(lifetime_Icens_Upper[k]|beta, eta),
      weibull_lcdf(lifetime_Icens_Lower[k]|beta, eta)
    );
  }
  
  // Prior models
  // The prior was constructed by simulateing 100 datasets of size 
  // n = 100 from the true Weibull distribution and estimating the 
  // paramaters via MLE and calculating to value of the estimated 
  // CDF at t1 and t2 to get a distribution.
  t1CDF ~ beta(t1cdf_a, (t1cdf_b - t1cdf_a));
  t2CDF ~ beta(t2cdf_a, (t2cdf_b - t2cdf_a));
}
")

t_1 <- 800
cdf_est_t1 <- 0.5
cdf_sd_t1 <- 0.05
t1cdf_a <- ((cdf_est_t1^2) * (1 - cdf_est_t1) / (cdf_sd_t1^2)) - cdf_est_t1
t1cdf_b <- t1cdf_a / cdf_est_t1

t_2 <- 9000
cdf_est_t2 <- 0.98
cdf_sd_t2 <- 0.01
t2cdf_a <- ((cdf_est_t2^2) * (1 - cdf_est_t2) / (cdf_sd_t2^2)) - cdf_est_t2
t2cdf_b <- t2cdf_a / cdf_est_t2

# Plot joint prior
cdf1 <- rbeta(5000, t1cdf_a, (t1cdf_b - t1cdf_a))
cdf2 <- rbeta(5000, t2cdf_a, (t2cdf_b - t2cdf_a))
sample_check <- cdf2 > cdf1

cdf1_sta <- log(-log(1 - cdf1[sample_check]))
cdf2_sta <- log(-log(1 - cdf2[sample_check]))

beta_draws <- (cdf2_sta - cdf1_sta) / log(t_2/t_1)
eta_draws <- exp(log(t_1) - (cdf1_sta / beta_draws)) 

data.frame(
  beta = beta_draws,
  eta = eta_draws
) %>%
  ggplot(aes(x = beta, y = eta)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  xlim(0.5, 1.5) +
  ylim(800, 1500) +
  geom_point(x = beta_truth, y = eta_truth, col = "red") +
  my_theme()

# Draw from posterior

stan_data <- list(
  N_obs = nrow(fake_data_obs),
  N_Rcens = nrow(fake_data_Rcens),
  N_Icens = nrow(fake_data_Icens),
  lifetime_obs = fake_data_obs$lifetime,
  lifetime_Rcens = fake_data_Rcens$lifetime,
  lifetime_Icens_Upper = fake_data_Icens$replacement_time,
  lifetime_Icens_Lower = fake_data_Icens$lifetime,
  t1 = t_1,
  t2 = t_2,
  t1cdf_a = t1cdf_a,
  t1cdf_b = t1cdf_b,
  t2cdf_a = t2cdf_a,
  t2cdf_b = t2cdf_b
)

informative_stan_fit <- sampling(
  informative_stan_model,
  stan_data,
  chains = 4,
  iter = 2000,
  warmup = 500,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 14
  )
)

informative_stan_fit %>%
  as_draws_df() %>%
  spread_draws(beta, eta) %>%
  ggplot(aes(x = beta, y = eta)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  xlim(0.5, 1.5) +
  ylim(800, 1500) +
  geom_point(x = beta_truth, y = eta_truth, col = "red") +
  my_theme()

inf_marginal_post <- informative_stan_fit %>%
  as_draws_rvars() %>%
  spread_rvars(beta, eta)

res <- 100

data.frame(
  t = seq(1, 3000, length.out = res),
  beta = rep(inf_marginal_post$beta, res),
  eta = rep(inf_marginal_post$eta, res),
  beta_true = beta_truth,
  eta_true = eta_truth
) %>%
  mutate(
    hazard = (beta / eta) * ((t / eta)^(beta - 1)),
    hazard_true = (beta_true / eta_true) * ((t / eta_true)^(beta_true - 1))
  ) %>%
  ggplot(aes(x = t, dist = hazard)) +
  stat_dist_lineribbon() +
  scale_fill_brewer() +
  geom_line(aes(y = hazard_true), color = "red") +
  ylab("hazard") +
  my_theme()

data.frame(
  t = seq(1, 6000, length.out = res),
  beta = rep(inf_marginal_post$beta, res),
  eta = rep(inf_marginal_post$eta, res),
  beta_true = beta_truth,
  eta_true = eta_truth
) %>%
  mutate(
    cdf = rvarPweibull(t, shape = beta, scale = eta),
    cdf_true = stats::pweibull(t, shape = beta_true, scale = eta_true)
  ) %>%
  ggplot(aes(x = t, dist = cdf)) +
  stat_dist_lineribbon() +
  scale_fill_brewer() +
  geom_line(aes(y = cdf_true), color = "red") +
  ylab("cdf") +
  my_theme() 






informative_stan_model <- stan_model(
  model_code = "
functions {
  real fn(real tCDF) {
    return log(-log1m(tCDF));
  }
}

data {
  int N_obs;
  int N_Icens;
  int N_Rcens;
  vector<lower=0>[N_obs] lifetime_obs;
  vector<lower=0>[N_Rcens] lifetime_Rcens;
  vector<lower=0>[N_Icens] lifetime_Icens_Upper;
  vector<lower=0>[N_Icens] lifetime_Icens_Lower;
}

parameters {
  real<lower = 0> t1CDF;
  real<lower = t1CDF> t2CDF; 
  vector<lower = lifetime_Rcens>[N_Rcens] Y_Rcens;
  vector<lower=lifetime_Icens_Lower, upper=lifetime_Icens_Upper>[N_Icens] Y_Icens;
}

transformed parameters {
  real<lower = 0> beta;
  real<lower = 0> eta;

  // calculate Weibull paramaters based on the
  // draws from the CDF at t1 and t2.
  beta = (fn(t2CDF) - fn(t1CDF)) / log(3000 / 800);
  eta = exp(log(800) - (fn(t1CDF) / beta));
}

model {
  // Likelihood
  lifetime_obs ~ weibull(beta, eta);
  Y_Rcens ~ weibull(beta, eta);
  Y_Icens ~ weibull(beta, eta);
  
  // Prior models
  t1CDF ~ beta(49.5, (99 - 49.5));
  t2CDF ~ beta(191.1, (195 - 191.1));
}
")

naive_stan_fit <- sampling(
  informative_stan_model,
  stan_data,
  chains = 4,
  iter = 800,
  warmup = 400#,
#  control = list(
#    adapt_delta = 0.99,
#    max_treedepth = 14
# )
)
