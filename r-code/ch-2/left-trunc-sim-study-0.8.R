library(dplyr)
library(ggplot2)
library(ggdist)
library(cowplot)
library(rstan)
library(posterior)
library(tidybayes)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
rvar_dweibull <- rfun(dweibull)
GetFactorValue <- function(x) as.numeric(as.character(x))

# Define true parameters
beta <- 0.8
eta <- 1

# Define data gen function
SimData <- function(
  beta,
  eta,
  n_units,
  t_start,
  t_end
){
  # Calculate how many lifetimes to sample
  weibull_05_quant <- qweibull(0.05, shape = beta, scale = eta)
  n_lifetimes_per_unit <- ceiling(t_end / weibull_05_quant) * 2

  # Create simulation data frame
  small_sim <- data.frame(
      # Define units
      unit = factor(
        rep(1:n_units, each = n_lifetimes_per_unit),
        1:n_units
      ),
      # Sample for Weibull distribution
      lifetime = rweibull(
        n_units * n_lifetimes_per_unit,
        shape = beta,
        scale = eta
      )
    ) %>% 
    group_by(unit) %>%
    # Calculate the failure times and install times
    mutate(
      failure_time = cumsum(lifetime),
      install_time = lag(failure_time),
      lifetime_id = 1:n()
    ) %>%
    ungroup() %>%
    # Replace NAs created by lag() with t = 0
    replace(is.na(.), 0) %>%
    # Discard any lifetimes that didn't fail or within the observation period
    filter(
      between(install_time, t_start, t_end) |
      between(failure_time, t_start, t_end) |
      ((install_time < t_start) & (failure_time > t_end))
    ) %>%
    # Create right and interval censoring indicator variables
    mutate(
      int_censored = !between(install_time, t_start, t_end),
      right_censored = !between(failure_time, t_start, t_end),
      install_time_obs = ifelse(int_censored, t_start, install_time),
      failure_time_obs = ifelse(right_censored, t_end, failure_time)
    )

  return(small_sim)
}

# Define Stan models
## left trunc observed
### Weak prior
stan_model_lt_r_imputed <- stan_model(
  model_code = "
data {
int N_obs;
int N_rc;
int N_lt;
int N_lt_rc;
vector<lower=0>[N_obs] y_obs;
vector<lower=0>[N_rc] y_rc;
vector<lower=0>[N_lt] y_lt;
vector<lower=0>[N_lt] t_lt;
vector<lower=0>[N_lt_rc] y_lt_rc;
vector<lower=0>[N_lt_rc] t_lt_rc;
}
parameters {
real<lower = 0> beta;
real<lower = 0> eta;
vector<lower = y_rc>[N_rc] y_rc_hat;
vector<lower = y_lt_rc>[N_lt_rc] y_lt_rc_hat;
}
model{
// Data model
// fully observed lifetimes
y_obs ~ weibull(beta, eta);
// right censored lifetimes
y_rc_hat ~ weibull(beta, eta);
// left truncated lifetimes
for (i in 1:N_lt) {
  y_lt[i] ~ weibull(beta, eta) T[t_lt[i], ]; 
}
// left truncated and right censored lifetimes
for (j in 1:N_lt_rc) {
  y_lt_rc_hat[j] ~ weibull(beta, eta) T[t_lt_rc[j], ]; 
}

// Prior model
eta ~ normal(1, 1);
beta ~ normal(0.8, 1);
}

"
)
### Strong prior
stan_model_lt_r_imputed_inf <- stan_model(
  model_code = "
functions {
// function to simplify the calculation of eta and beta
real fn(real tCDF) {
  return log(-log1m(tCDF));
}
}
data {
int N_obs;
int N_rc;
int N_lt;
int N_lt_rc;
vector<lower=0>[N_obs] y_obs;
vector<lower=0>[N_rc] y_rc;
vector<lower=0>[N_lt] y_lt;
vector<lower=0>[N_lt] t_lt;
vector<lower=0>[N_lt_rc] y_lt_rc;
vector<lower=0>[N_lt_rc] t_lt_rc;
// Define the prior
real t_1;
real t_2;
real t1_mean;
real t1_var;
real t2_mean;
real t2_var;
}
parameters {
real<lower = 0, upper = 1> t1CDF;
real<lower = t1CDF, upper = 1> t2CDF;
vector<lower = y_rc>[N_rc] y_rc_hat;
vector<lower = y_lt_rc>[N_lt_rc] y_lt_rc_hat;
}
transformed parameters {
real<lower = 0> beta;
real<lower = 0> eta;

// calculate Weibull paramaters based on the
// draws from the CDF at t1 and t2.
beta = (fn(t2CDF) - fn(t1CDF)) / log(t_2 / t_1);
eta = exp(log(t_1) - (fn(t1CDF) / beta));
}
model{
// Data model
// fully observed lifetimes
y_obs ~ weibull(beta, eta);
// right censored lifetimes
y_rc_hat ~ weibull(beta, eta);
// left truncated lifetimes
for (i in 1:N_lt) {
  y_lt[i] ~ weibull(beta, eta) T[t_lt[i], ]; 
}
// left truncated and right censored lifetimes
for (j in 1:N_lt_rc) {
  y_lt_rc_hat[j] ~ weibull(beta, eta) T[t_lt_rc[j], ]; 
}

// Prior model
t1CDF ~ normal(t1_mean, t1_var);
t2CDF ~ normal(t2_mean, t2_var);
}

"
)
## left trunc discarded
### Weak prior
stan_model_no_lt <- stan_model(
  model_code = "
data {
int N_obs;
int N_rc;
vector<lower=0>[N_obs] y_obs;
vector<lower=0>[N_rc] y_rc;
}
parameters {
real<lower = 0> beta;
real<lower = 0> eta;
}
model{
// Data model
// fully observed lifetimes
target += weibull_lpdf(y_obs|beta, eta);
// right censored lifetimes
target += weibull_lccdf(y_rc|beta, eta);

// Prior model
eta ~ normal(1, 1);
beta ~ normal(0.8, 1);
}

"
)
### Strong prior
stan_model_no_lt_inf <- stan_model(
  model_code = "
functions {
// function to simplify the calculation of eta and beta
real fn(real tCDF) {
  return log(-log1m(tCDF));
}
}
data {
int N_obs;
int N_rc;
vector<lower=0>[N_obs] y_obs;
vector<lower=0>[N_rc] y_rc;
// Define the prior
real t_1;
real t_2;
real t1_mean;
real t1_var;
real t2_mean;
real t2_var;
}
parameters {
real<lower = 0, upper = 1> t1CDF;
real<lower = t1CDF, upper = 1> t2CDF;
}
transformed parameters {
real<lower = 0> beta;
real<lower = 0> eta;

// calculate Weibull paramaters based on the
// draws from the CDF at t1 and t2.
beta = (fn(t2CDF) - fn(t1CDF)) / log(t_2 / t_1);
eta = exp(log(t_1) - (fn(t1CDF) / beta));
}
model{
// Data model
// fully observed lifetimes
target += weibull_lpdf(y_obs|beta, eta);
// right censored lifetimes
target += weibull_lccdf(y_rc|beta, eta);

// Prior model
t1CDF ~ normal(t1_mean, t1_var);
t2CDF ~ normal(t2_mean, t2_var);
}

"
)
## left trunc imputed
### Weak prior
stan_model_unknown_lt_rc <- stan_model(
  model_code = "
data {
int N_obs;                             # N fully observed lives
int N_rc;                              # N right censored only lives
int N_lt;                              # N left truncated only lives
int N_lt_rc;                           # N right cens and left trunc lives
array[N_obs] real<lower=0> y_obs;      # Fully observed lifetimes
array[N_rc] real<lower=0> y_rc;        # Right censored lifetimes
array[N_lt] real<lower=0> y_lt;        # Left trunc lifetimes
array[N_lt_rc] real<lower=0> y_lt_rc;  # right cens and left trunc lifetimes
real<lower=0> t_start;                 # start of the observation window
}
transformed data{
array[N_lt] real<lower=0> y_lt_upper;  # The upper bound of the left trunc lives

for (m in 1:N_lt){
  y_lt_upper[m] = y_lt[m] + t_start;   # Upper bound = lower bound + start of observation
}

}
parameters {
real<lower= 0> beta;     # weibull shape
real<lower= 0> eta;      # weibull scale
array[N_rc] real<lower=y_rc> y_rc_hat;   # imputed right censored values
array[N_lt] real<lower=y_lt, upper=y_lt_upper> y_lt_hat;  # imputed left trunc values
array[N_lt_rc] real<lower=y_lt_rc> y_lt_rc_hat;   # imputed left trunc and right cens values
array[N_lt_rc] real<lower=0, upper=1> t_lt_rc_st; # imputed left truncation times for left trunc and right cens values (standardised)
}
transformed parameters{
array[N_lt] real t_lt;  # imputed left trunc times for left trunc values
array[N_lt_rc] real<lower=0, upper=t_start> t_lt_rc_upper;
array[N_lt_rc] real<lower=0, upper=t_lt_rc_upper> t_lt_rc;  # imputed left trunc times for left trunc and right cens values

for (i in 1:N_lt) {
  t_lt[i] = y_lt_hat[i] - y_lt[i];
}

for (k in 1:N_lt_rc){
  if ((y_lt_rc_hat[k] - y_lt_rc[k]) < t_start)
    t_lt_rc_upper[k] = y_lt_rc_hat[k] - y_lt_rc[k];
  else
    t_lt_rc_upper[k] = t_start;

  t_lt_rc[k] = t_lt_rc_st[k] * t_lt_rc_upper[k];
}
}
model{
// Data model
// fully observed lifetimes
y_obs ~ weibull(beta, eta);
// right censored lifetimes
y_rc_hat ~ weibull(beta, eta);
// left truncated lifetimes
for (i in 1:N_lt) {
  y_lt_hat[i] ~ weibull(beta, eta) T[t_lt[i], ]; 
}
// left truncated and right censored lifetimes
for (j in 1:N_lt_rc) {
  y_lt_rc_hat[j] ~ weibull(beta, eta) T[t_lt_rc[j], ]; 
}

// Prior model
eta ~ normal(1, 1);
beta ~ normal(0.8, 1);
t_lt_rc_st ~ uniform(0, 1);
}

"
)
### Strong prior
stan_model_unknown_lt_rc_inf <- stan_model(
  model_code = "
functions {
// function to simplify the calculation of eta and beta
real fn(real tCDF) {
  return log(-log1m(tCDF));
}
}
data {
int N_obs;                             # N fully observed lives
int N_rc;                              # N right censored only lives
int N_lt;                              # N left truncated only lives
int N_lt_rc;                           # N right cens and left trunc lives
array[N_obs] real<lower=0> y_obs;      # Fully observed lifetimes
array[N_rc] real<lower=0> y_rc;        # Right censored lifetimes
array[N_lt] real<lower=0> y_lt;        # Left trunc lifetimes
array[N_lt_rc] real<lower=0> y_lt_rc;  # right cens and left trunc lifetimes
real<lower=0> t_start;                 # start of the observation window
// Define the prior
real t_1;
real t_2;
real t1_mean;
real t1_var;
real t2_mean;
real t2_var;
}
transformed data{
array[N_lt] real<lower=0> y_lt_upper;  # The upper bound of the left trunc lives

for (m in 1:N_lt){
  y_lt_upper[m] = y_lt[m] + t_start;   # Upper bound = lower bound + start of observation
}

}
parameters {
real<lower = 0, upper = 1> t1CDF;
real<lower = t1CDF, upper = 1> t2CDF;
array[N_rc] real<lower=y_rc> y_rc_hat;   # imputed right censored values
array[N_lt] real<lower=y_lt, upper=y_lt_upper> y_lt_hat;  # imputed left trunc values
array[N_lt_rc] real<lower=y_lt_rc> y_lt_rc_hat;   # imputed left trunc and right cens values
array[N_lt_rc] real<lower=0, upper=1> t_lt_rc_st; # imputed left truncation times for left trunc and right cens values (standardised)
}
transformed parameters{
real<lower = 0> beta;
real<lower = 0> eta;
array[N_lt] real t_lt;  # imputed left trunc times for left trunc values
array[N_lt_rc] real<lower=0, upper=t_start> t_lt_rc_upper;
array[N_lt_rc] real<lower=0, upper=t_lt_rc_upper> t_lt_rc;  # imputed left trunc times for left trunc and right cens values

// calculate Weibull paramaters based on the
// draws from the CDF at t1 and t2.
beta = (fn(t2CDF) - fn(t1CDF)) / log(t_2 / t_1);
eta = exp(log(t_1) - (fn(t1CDF) / beta));

for (i in 1:N_lt) {
  t_lt[i] = y_lt_hat[i] - y_lt[i];
}

for (k in 1:N_lt_rc){
  if ((y_lt_rc_hat[k] - y_lt_rc[k]) < t_start)
    t_lt_rc_upper[k] = y_lt_rc_hat[k] - y_lt_rc[k];
  else
    t_lt_rc_upper[k] = t_start;

  t_lt_rc[k] = t_lt_rc_st[k] * t_lt_rc_upper[k];
}
}
model{
// Data model
// fully observed lifetimes
y_obs ~ weibull(beta, eta);
// right censored lifetimes
y_rc_hat ~ weibull(beta, eta);
// left truncated lifetimes
for (i in 1:N_lt) {
  y_lt_hat[i] ~ weibull(beta, eta) T[t_lt[i], ]; 
}
// left truncated and right censored lifetimes
for (j in 1:N_lt_rc) {
  y_lt_rc_hat[j] ~ weibull(beta, eta) T[t_lt_rc[j], ]; 
}

// Prior model
t1CDF ~ normal(t1_mean, t1_var);
t2CDF ~ normal(t2_mean, t2_var);
t_lt_rc_st ~ uniform(0, 1);
}

"
)

# Define the factor levels for experiment
expected_lifetime <- eta * gamma(1 + (1 / beta))
M_itteration <- 50
## Dataset sim factors
N_levels <- c(
  10,
  100,
  500
)
t_start_levels <- c(
  round(expected_lifetime * 1, 2),
  round(expected_lifetime * 5, 2),
  round(expected_lifetime * 15, 2)
)
t_window_levels <- c(
  round(expected_lifetime, 2),
  round(expected_lifetime * 3, 2),
  round(expected_lifetime * 6, 2)
)
## Model factors
left_trunc_treatment_levels <- c(
  "known",
  "discarded",
  "imputed"
)
prior_levels <- c(
  "week",
  "strong"
)

# Get grid of factor combinations
factor_sets_df <- expand.grid(
  N = factor(N_levels),
  t_start = factor(t_start_levels),
  t_window = factor(t_window_levels)
)

# Create test set for model validation
test_sets <- lapply(
  1:100,
  function(i) rweibull(100, beta, eta)
)

# Simulations
experiment_results_df <- parallel::mclapply(
  1:nrow(factor_sets_df),
  function(set_id) {
    factor_set <- factor_sets_df[set_id, ]
    factor_set_results_df <- lapply(
      1:100,
      function(itr) {
        #   Sim data   #
        sim_data <- SimData(
          beta = beta,
          eta = eta,
          n_units = GetFactorValue(factor_set$N),
          t_start = GetFactorValue(factor_set$t_start),
          t_end = GetFactorValue(factor_set$t_start) + GetFactorValue(factor_set$t_window)
        )
        # separate different censoring-truncation types
        sim_data_obs <- sim_data %>%
          filter(!(right_censored | int_censored))
        sim_data_rc <- sim_data %>%
          filter(right_censored & !int_censored)
        sim_data_lt <- sim_data %>%
          filter(!right_censored & int_censored)
        sim_data_lt_rc <- sim_data %>%
          filter(right_censored & int_censored)
        # prepare for stan
        stan_data_full <- list(
          N_obs = nrow(sim_data_obs),
          N_rc = nrow(sim_data_rc),
          N_lt = nrow(sim_data_lt),
          N_lt_rc = nrow(sim_data_lt_rc),
          y_obs = as.array(
            sim_data_obs$failure_time - sim_data_obs$install_time
          ),
          y_rc = as.array(
            sim_data_rc$failure_time_obs - sim_data_rc$install_time
          ),
          y_lt = as.array(
            sim_data_lt$failure_time - sim_data_lt$install_time
          ),
          t_lt = as.array(
            sim_data_lt$install_time_obs - sim_data_lt$install_time
          ),
          y_lt_rc = as.array(
            sim_data_lt_rc$failure_time_obs - sim_data_lt_rc$install_time
          ),
          t_lt_rc = as.array(
            sim_data_lt_rc$install_time_obs - sim_data_lt_rc$install_time
          )
        )
        stan_data_no_lt <- list(
          N_obs = nrow(sim_data_obs),
          N_rc = nrow(sim_data_rc),
          y_obs = as.array(
            sim_data_obs$failure_time - sim_data_obs$install_time
          ),
          y_rc = as.array(
            sim_data_rc$failure_time_obs - sim_data_rc$install_time
          )
        )
        stan_data_unkown_lt <- list(
          N_obs = nrow(sim_data_obs),
          N_rc = nrow(sim_data_rc),
          N_lt = nrow(sim_data_lt),
          N_lt_rc = nrow(sim_data_lt_rc),
          y_obs = as.array(
            sim_data_obs$failure_time - sim_data_obs$install_time
          ),
          y_rc = as.array(
            sim_data_rc$failure_time_obs - sim_data_rc$install_time
          ),
          y_lt = as.array(
            sim_data_lt$failure_time_obs - sim_data_lt$install_time_obs
          ),
          y_lt_rc = as.array(
            sim_data_lt_rc$failure_time_obs - sim_data_lt_rc$install_time_obs
          ),
          t_start = GetFactorValue(factor_set$t_start)
        )
        #  Fit models  #
        # weak prior
        stan_fit_full_weak <- sampling(
          stan_model_lt_r_imputed,
          stan_data_full,
          chains = 3,
          cores = 1,
          iter = 500,
          warmup = 300,
          refresh = 0
        )
        stan_fit_no_lt_weak <- sampling(
          stan_model_no_lt,
          stan_data_no_lt,
          chains = 3,
          cores = 1,
          iter = 500,
          warmup = 300,
          refresh = 0
        )
        stan_fit_imp_lt_weak <- sampling(
          stan_model_unknown_lt_rc,
          stan_data_unkown_lt,
          chains = 3,
          cores = 1,
          iter = 500,
          warmup = 300,
          refresh = 0
        )
        # strong prior
        stan_fit_full_strong <- sampling(
          stan_model_lt_r_imputed_inf,
          c(
            stan_data_full,
            t_1 = qweibull(0.8, beta, eta),
            t_2 = qweibull(0.98, beta, eta),
            t1_mean = 0.8,
            t1_var = 0.1,
            t2_mean = 0.98,
            t2_var = 0.05
          ),
          chains = 3,
          cores = 1,
          iter = 500,
          warmup = 300,
          refresh = 0
        )
        stan_fit_no_lt_strong <- sampling(
          stan_model_no_lt_inf,
          c(
            stan_data_no_lt,
            t_1 = qweibull(0.8, beta, eta),
            t_2 = qweibull(0.98, beta, eta),
            t1_mean = 0.8,
            t1_var = 0.1,
            t2_mean = 0.98,
            t2_var = 0.05
          ),
          chains = 3,
          cores = 1,
          iter = 500,
          warmup = 300,
          refresh = 0
        )
        stan_fit_imp_lt_strong <- sampling(
          stan_model_unknown_lt_rc_inf,
          c(
            stan_data_unkown_lt,
            t_1 = qweibull(0.8, beta, eta),
            t_2 = qweibull(0.98, beta, eta),
            t1_mean = 0.8,
            t1_var = 0.1,
            t2_mean = 0.98,
            t2_var = 0.05
          ),
          chains = 3,
          cores = 1,
          iter = 500,
          warmup = 300,
          refresh = 0
        )
        # Calc metrics #
        sim_results_df <- lapply(
          c(
            "stan_fit_full_weak", "stan_fit_full_strong",
            "stan_fit_no_lt_weak", "stan_fit_no_lt_strong",
            "stan_fit_imp_lt_weak", "stan_fit_imp_lt_strong"
          ),
          function(stan_fit_name) {
            stan_fit_obj <- get(stan_fit_name)
            joint_post <- stan_fit_obj %>%
              as_draws_rvars()
            model_results_df <- factor_set %>%
              mutate(
                itr = itr,
                model = stan_fit_name,
                p_value_beta = Pr(joint_post$beta > beta),
                p_value_eta = Pr(joint_post$eta > eta),
                E_log_score = lapply(
                  test_sets,
                  function(test_set) {
                    elppd <- rvar_dweibull(
                      test_set,
                      joint_post$beta,
                      joint_post$eta
                    ) %>%
                      log() %>%
                      rvar_sum() %>%
                      E()
                    if (is.infinite(elppd)) elppd <- NA
                    return(elppd)
                  }
                ) %>%
                  unlist() %>%
                  mean(na.rm = TRUE)
              )
            return(model_results_df)
          }
        ) %>%
          bind_rows()
        return(sim_results_df)
      }
    ) %>%
      bind_rows()
    return(factor_set_results_df)
  },
  mc.cores = (parallel::detectCores() - 1)
) %>%
  bind_rows()

saveRDS(
  experiment_results_df,
  file = "LT_experiment_results_df_0.8.rds"
)
