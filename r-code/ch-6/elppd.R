library(rstan)
library(dplyr)
library(posterior)
library(tidybayes)
library(tidyr)
library(kableExtra)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load data and model
stan_data <- readRDS(
  file = "stan-data.rds"
)
stan_model_gp <- readRDS(
  file = "compiled-stan-model-gp.rds"
)
stan_model_lm <- readRDS(
  file = "compiled-stan-model-lm.rds"
)

sa_stan_fit <- function(i, stan_model) {
  stan_data_training <- list(
    I = i,
    N = stan_data$N,
    M = stan_data$M,
    t = stan_data$t[1:i],
    z = stan_data$z[1:i, ],
    B = stan_data$B,
    a_hat = stan_data$a_hat,
    b_hat = stan_data$b_hat
  )
  
  stan_fit <- sampling(
    stan_model,
    stan_data_training,
    chains = 4,
    iter = 1500,
    warmup = 800,
    refresh = 0,
    control = list(
      adapt_delta = 0.99, 
      max_treedepth = 14
    )
  )

  return(stan_fit)
}

rTruncT <- function(n, nu, mean, sd, lb) {  
  lb_norm <- (lb - mean) / sd
  unif_rv <- runif(n, pt(lb_norm, df = nu), 1)
  t_rv <- (qt(unif_rv, df = nu) * sd) + mean
  
  return(t_rv)
}
rvar_rTruncT  <- rfun(rTruncT)
rvar_rgamma <- rfun(rgamma)
rvar_dt <- rfun(dt)

max.LL <- function(par, z){
  z_mean <- stan_data$B %*% matrix(par, ncol = 1)
  ll <- sum(log(dnorm(z, z_mean, 0.39)))

  return(-ll)
}

elppd_gp <- lapply(
  5:8,
  function(I) lapply(
    (I+1):stan_data$I,
    function(I_p1) {
      forecast_step <- stan_data$t[I_p1] - stan_data$t[I]
      # Fit model to training data
      stan_model <- sa_stan_fit(I, stan_model_gp)
      # Observed coefs
      y_forecast_obs <- optim(
        par = rep(1, 8), 
        fn = max.LL, 
        z = stan_data$z[I_p1, ],
        lower = 0
      )[["par"]]
      # Generate ppd from posterior
      ppd_forecast <- stan_model %>%
        as_draws_rvars() %>%
        spread_rvars(y[i, m], mu[m], nu[m], phi) %>%
        filter(i == (I - 1)) %>%
        mutate(
          y = y * mu,
          y_jump = rvar_rgamma(
            8, 
            shape = (forecast_step / nu^2), 
            rate = (1 / (mu * nu^2))
          ),
          y_forecast = y + y_jump
        )
      y_forecast_obs_scaled <- (y_forecast_obs - ppd_forecast$y_forecast) / 
        (unique(ppd_forecast$phi) * ppd_forecast$y_forecast)
      elpd <- rvar_dt(y_forecast_obs_scaled, df = 10) %>%
        rvar_prod() %>%
        E() %>%
        log()
      return(elpd)
    }
  )
)

elppd_lm <- lapply(
  5:8,
  function(I) lapply(
    (I+1):stan_data$I,
    function(I_p1) {
      # Fit model to training data
      stan_model <- sa_stan_fit(I, stan_model_lm)
      # Observed coefs
      y_forecast_obs <- optim(
        par = rep(1, 8), 
        fn = max.LL, 
        z = stan_data$z[I_p1, ],
        lower = 0
      )[["par"]]
      # Generate ppd from posterior
      ppd_forecast <- stan_model %>%
        as_draws_rvars() %>%
        spread_rvars(mu[m], phi) %>%
        mutate(
          y_forecast = mu * stan_data$t[I_p1]
        )
      y_forecast_obs_scaled <- (y_forecast_obs - ppd_forecast$y_forecast) / 
        (unique(ppd_forecast$phi) * ppd_forecast$y_forecast)
      elpd <- rvar_dt(y_forecast_obs_scaled, df = 10) %>%
        rvar_prod() %>%
        E() %>%
        log()
      return(elpd)
    }
  ) 
)

options(knitr.kable.NA = '')
data.frame(
  `$I$` = c(
    rep(5, 4),
    rep(6, 3),
    rep(7, 2),
    rep(8, 1)
  ),
  `$I + 1$` = c(
    6:9,
    7:9,
    8:9,
    9
  ),
  `$\\mbox{ELS}_{\\textit{gamma process}}$` = round(unlist(elppd_gp), 3),
  `$\\mbox{ELS}_{\\textit{linear path}}$` = round(unlist(elppd_lm), 3),
  check.names = FALSE
) %>% 
bind_rows(
  .,
  tibble(
    `$\\mbox{ELS}_{\\textit{gamma process}}$` = sum(.[["$\\mbox{ELS}_{\\textit{gamma process}}$"]]),
    `$\\mbox{ELS}_{\\textit{linear path}}$` = sum(.[["$\\mbox{ELS}_{\\textit{linear path}}$"]]),
  )
) %>% 
  kbl(
    booktabs = T,
    format = "latex",
    caption = "The expected log score ($ELS$) for each model when fitting to a portion of the data and predicting n-steps ahead. $I$ is the maximum observation that the model was fit to and $I + 1$ is the withheld observation that the forecast is generated for. The summation of the elppd scores are displayed at the bottom of the table.",
    escape = FALSE,
    label = "elppd-beltwear"
  ) %>%
  row_spec(11, bold = TRUE) %>%
  kable_styling(latex_options = "striped") %>%
  save_kable(
    file = file.path(
      "..",
      "..",
      "tables",
      "ch-6",
      "stan_fit_cp_table.tex"
    ),
    keep_tex = TRUE
  )
