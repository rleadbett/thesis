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

#cols <- wes_palette("Zissou1", cross_sections*2, type = "continuous")
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

load(
  file.path("..", "data", "Example_beltwear_data.RData")
)
belt_ut_data <- example_belts[["belt_A"]]

# Define true design matrix
numb_knots <- 8
B <- bSpline(seq(0, 15, length.out = 22),
             knots = seq(1, 14, length.out = numb_knots),
             degree = 3,
             intercept = T)[, 3:10]

# Define fine desing matrix for smooth plotting
fine_grid <- seq(0, 15, 0.1)
fine_basis_mat <- bSpline(fine_grid,
                          knots = seq(1, 14, 
                                      length.out = numb_knots),
                          degree = 3,
                          intercept = T)[, 3:10]


# Function to fit B-spline coeficients
min.RSS <- function(par, y){
  par_mat <- matrix(par, ncol = 1)
  res <- (B %*% par)[2:21] - y
  return(sum(res^2))
}

# Fit splines to UT data
cumulative_tonnes <- unique(belt_ut_data$tonnes)
fitted_spline_coefs <- lapply(
  cumulative_tonnes,
  function(tt) {
    ut_measurements <- belt_ut_data %>%
      filter(tonnes == tt) %>%
      pull(wear)

    spline_coefs <- optim(
      par = rep(0.1, 8), 
      fn = min.RSS, 
      y = ut_measurements,
      lower = 0
    )[["par"]]

    return(
      spline_coefs
    )
  }
)

# Plot smoothed wear profile
PlotFittedBsplines <- function(obs) {
  n_grid <- seq(0, 22, length.out = length(fine_grid))
  coef_mat <- matrix(fitted_spline_coefs[[obs]], ncol = 1)
  functional_obs <- fine_basis_mat %*% coef_mat
  p_spline <- data.frame(
      x = n_grid,
      y = functional_obs
    ) %>%
    ggplot(aes(x = x, y = y)) +
    geom_line(
      color = brewer.pal(n = 9, "YlOrRd")[obs],
      size = 1.2
    ) +
      geom_point(
      data = belt_ut_data %>%
        mutate(
          tonnes = factor(tonnes),
          x = rep(seq(0, 22, length.out = 22)[2:21], 9)
        ) %>%
        rename(
          y = wear
        ),
      aes(col = tonnes),,
      alpha = 0.35,
      size = 3
    ) +
    scale_color_brewer(palette = "YlOrRd") +
    geom_point(
      data = data.frame(
        x = seq(0, 22, length.out = 22)[2:21],
        y = belt_ut_data %>%
          filter(tonnes == cumulative_tonnes[obs]) %>%
          pull(wear)
      ),
      color = brewer.pal(n = 9, "YlOrRd")[obs],
      size = 3
    ) +
    ylim(0, 16) +
    xlab("") +
    ylab("wear (mm)") +
    my_theme() +
    theme(
      legend.position = "none"
    )

  # Plot of weighted basis funcitons
  p_weighted_basis <- data.frame(
      t(t(fine_basis_mat) * fitted_spline_coefs[[obs]])
    ) %>%
    mutate(n = (fine_grid / 15) * 21) %>%
    pivot_longer(-c("n")) %>%
    ggplot(aes(x = n, y = value, group = name, col = name)) +
    geom_line() +
    ylab("weighted value") +
    xlab("") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    ylim(0, 11) +
    my_theme() +
    theme(
      legend.position = "none"
    )

  # Plot of unweighted basis functions
  p_raw_basis <- data.frame(fine_basis_mat) %>%
    mutate(n = (fine_grid / 15) * 21) %>%
    pivot_longer(-c("n")) %>%
    ggplot(aes(x = n, y = value, group = name, col = name)) +
    geom_line() +
    ylab("unweighted value") +
    xlab("n") +
    my_theme() +
    theme(
      legend.position = "none"
    )

  # Join_plots
  p <- plot_grid(p_spline, p_weighted_basis, p_raw_basis, ncol = 1)
  return(p)
}

saveGIF(
  {
    for (i in 1:9) {
      p <- PlotFittedBsplines(i)
      print(p)
    }
  },
  movie.name = "spline-fits.gif",
  ani.width = 2500,
  ani.height = 2700,
  ani.res = 400
)

p_coef_pathways <- data.frame(
  spline_coef = factor(str_c("coef ", rep(1:9, 8))),
  tonnes = rep(cumulative_tonnes, each = 8),
  coef_value = unlist(fitted_spline_coefs)
) %>%
  ggplot(aes(x = tonnes, y = coef_value, col = spline_coef)) +
  geom_line() +
  facet_wrap(vars(spline_coef), nrow = 1) +
  my_theme() +
  theme(
    legend.position = "none"
  ) +
  scale_x_continuous(breaks=seq(0, 0.5, 1))

fitSpline <- function(Y){
  optim(par = rep(0.1, 8), 
        fn = min.RSS, 
        y = Y,
        lower = 0)[["par"]] %>%
    return()
}

historic_spline_coefs <- rbind(
  example_belts[["belt_B"]] %>%
    mutate(belt = "belt_B"),
  example_belts[["belt_C"]] %>%
    mutate(belt = "belt_C")
) %>% 
  group_by(belt, tonnes) %>%
  summarise(
    spline_coef = fitSpline(wear),
    coef = 1:8
  )

coef_lm_fits <- lapply(
  1:8,
  function(m) {
    fit <- lm(
      data = historic_spline_coefs %>%
        filter(coef == m),
      formula = spline_coef ~ tonnes + 0
    )
    param_est <- summary(fit)
    return(
      c(
        a_hat = param_est$coefficients[1, 1],
        b_hat = param_est$coefficients[1, 2]
      )
    )
  }
)

set.seed(64467364)

t <- unique(belt_ut_data$tonnes)
I <- length(t)
N <- unique(belt_ut_data$measurement_pos)
M <- ncol(B)

# Get informative hyper parameters for mu
mu_hyper_par <- coef_lm_fits %>% 
  bind_rows() %>%
  mutate(b_hat = b_hat * 5)

# Sample values of the parameters
## mu
rTruncNorm <- function(n, mean, sd, lb) {
  lb_norm <- (lb - mean) / sd
  unif_rv <- runif(n, pnorm(lb_norm), 1)
  norm_rv <- (qnorm(unif_rv) * sd) + mean
  
  return(norm_rv)
}

mu <- lapply(
  1:M,
  function(m){
    rTruncNorm(
      1000, 
      mean = mu_hyper_par$a_hat[m], 
      sd = mu_hyper_par$b_hat[m],
      lb = 0
    ) %>%
    rvar() %>%
    return()
  }
)

## nu
rTruncT <- function(n, nu, mean, sd, lb) {  
  lb_norm <- (lb - mean) / sd
  unif_rv <- runif(n, pt(lb_norm, df = nu), 1)
  t_rv <- (qt(unif_rv, df = nu) * sd) + mean
  
  return(t_rv)
}

nu_mean <- rTruncT(    # sample from t_3
  1000, 
  nu = 3, 
  mean = 0, 
  sd = 0.5,
  lb = 0
) %>%
  rvar()

nu_sd <- rTruncT(       # sample from cauchy
  1000, 
  nu = 1, 
  mean = 0, 
  sd = 0.25, 
  lb = 0
) %>%
  rvar()

rvar_rTruncNorm <- rfun(rTruncNorm)
rvar_rTruncT <- rfun(rTruncT)

nu <- rvar_rTruncNorm(
  8,
  mean = nu_mean,
  sd = nu_sd,
  lb = 0
)

# Sample the jumps of the spline coefficients
rvar_rgamma <- rfun(rgamma)

jumps <- lapply(
  1:M,
  function(m){
    shape_par <- diff(t) / (nu[[m]]^2)
    rate_par <- 1 / (mu[[m]] * nu[[m]]^2)
    rvar_rgamma(
      (I - 1), 
      shape = shape_par, 
      rate = rate_par
    ) %>%
      return()
  }
)

# Calculate the value of the coefs and arange in matrix
spline_params <- lapply(
  jumps,
  cumsum
)

coef_mat <- lapply(
  spline_params,
  function(param){
    param %>%
      as.matrix() %>%
      return()
  }
) %>%
  do.call(cbind, .)

# Map to data space
ut_mat <- t(fine_basis_mat %**% t(coef_mat))
ut_df <- ut_mat %>%
  as.data.frame() %>%
  mutate(i = 1:(I - 1)) %>%
  pivot_longer(cols = -c("i")) %>%
  mutate(n = rep(1:nrow(fine_basis_mat), (I - 1)))

# Functions to generate PPC
ppc <- function() {
  draw_id <- sample(1:1000, 1)
  nu_mean_draw <- ut_df %>%
    mutate(i = as.factor(i),
          draw = lapply(value, 
                        function(val) {
                          draws_of(val)[draw_id]
                        }) %>%
                  unlist()) %>%
    ggplot(aes(x = n, y = draw, col = i)) +
    geom_line() +
    scale_color_brewer(palette = "YlOrRd") +
    theme(legend.position = "none") +
    annotate(
      "text", 
      x = 50, y = 75, 
      label = as.expression(
        bquote(
          "("~
          nu["mean"]~
          " ="~
          .(round(draws_of(nu_mean)[draw_id], digits = 2))~
          ","~
          nu["sd"]~
          " ="~
          .(round(draws_of(nu_sd)[draw_id], digits = 2))~
          ")"
        )
      ), size = 2
    )
}

# Plot 16 PPC in a 4x4 grid
ppc_list <- lapply(
  1:16, 
  function(x) {
    ppc() +
      my_theme() +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
      ) +
      ylim(0, 80)
  }
)

p_prior_predictive_checks <- plot_grid(plotlist = ppc_list, ncol = 4)

stan_model <- readRDS(
  file.path(".", "compiled-stan-model.rds")
)
t <- unique(belt_ut_data$tonnes)
Ut_matrix <- matrix(belt_ut_data$wear,
                    nrow = length(t),
                    byrow = TRUE)

stan_data_full <- list(
  I = length(t) - 1,
  N = ncol(Ut_matrix),
  M = ncol(B),
  t = t[1:8],
  z = Ut_matrix[1:8, ],
  B = B[2:21, ],
  a_hat = lapply(coef_lm_fits, function(m) m["a_hat"]) %>% unlist(),
  b_hat = lapply(coef_lm_fits, function(m) m["b_hat"]) %>% unlist() * 5
)

stan_fit_full <- rstan::sampling(
  stan_model,
  stan_data_full,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  control = list(
    adapt_delta = 0.99, 
    max_treedepth = 14
  )
)

rvar_rTruncNorm <- rfun(rTruncNorm)

wear_post <- stan_fit_full %>%
  as_draws_rvars() %>%
  spread_rvars(y_noisy[i, m], sigma) %>%
  as.data.frame()

wear_ut <- wear_post %>%
  rbind(
    data.frame(
      y_noisy = rep(rvar(rep(0, 3000*4), nchains = 4) , 8),
      i = 0,
      m = 1:8,
      sigma = rep(unique(wear_post$sigma), 8)
    )
  ) %>%
  arrange(m, i) %>%
  group_by(i) %>%
  summarise(z_smoothed = y_noisy %**% t(fine_basis_mat))


p_posterior <- wear_ut$z_smoothed %>%
  as.data.frame() %>%
  pivot_longer(cols = everything()) %>%
  mutate(
    n = rep(seq(0, 21, length.out = nrow(fine_basis_mat)), 8),
    t = factor(rep(round(t[1:8], 2), each = nrow(fine_basis_mat))),
    sd_ut = unique(wear_post$sigma),
    z_noisy = rvar_rTruncNorm(1, 
                              mean = value,
                              sd = sd_ut,
                              lb = -5)
  ) %>%
  ggplot(aes(x = n, dist = z_noisy, col = t, fill = t)) +
  stat_dist_lineribbon(.width = c(0.95)) +
  scale_colour_brewer(palette = "YlOrRd") +
  scale_fill_manual(
    values = alpha(
      brewer.pal(n = 8, name = "YlOrRd"), 
      0.3
    )
  ) +
  geom_point(
    data = belt_ut_data %>% 
    filter(tonnes < 1) %>%
    mutate(
      n = rep(1:20, 8),
      t = as.factor(round(tonnes, 2)),
      z_noisy = wear
    ),
    aes(y = z_noisy)
  ) +
  ylim(-2, 30) +
  xlab("measurement location") +
  ylab("wear (mm)") +
  labs(col = "cumulative tonnage") +
  my_theme() +
  theme(legend.position = "bottom")

SimulateSpline <- function(draw, tt) {
  sim_gp <- stan_fit_full %>%
    as_draws_df() %>%
    spread_draws(mu[m], nu[m], phi) %>%
    filter(.draw == draw) %>%
    mutate(
      jump = rgamma(
        n(),
        shape = (tt / nu^2),
        rate = 1 / (mu * nu^2)
      ),
      noisy_jump = rTruncT(
        n(),
        nu = 10,
        mean = jump,
        sd = jump * phi,
        lb = 0
      )
    )
  function_vals <- sim_gp$noisy_jump %*% t(fine_basis_mat)
  return(
    data.frame(
      tonnes = tt,
      n = seq(0, 21, length.out = length(function_vals)),
      wear = as.numeric(function_vals)
    )
  )
}

PosteriorPredictiveCheck <- function(draw) {
  p <- lapply(
    t[2:9],
    function(tt) SimulateSpline(draw, tt)
  ) %>%
    bind_rows() %>%
    mutate(tonnes = factor(round(tonnes, 2))) %>%
    ggplot(aes(x = n, y = wear, col = tonnes)) +
    geom_line() +
    scale_color_brewer(palette = "YlOrRd") +
    ylim(0, 30) +
    my_theme() +
    theme(legend.position = "none")

  return(p)
}

draws <- sample(1:12000, 3)

p_post_pc_list <- lapply(
  draws,
  function(draw) PosteriorPredictiveCheck(draw)
)

p_post_pc <- plot_grid(plotlist = p_post_pc_list, nrow = 3)

p_forecast <- stan_fit_full %>%
  as_draws_rvars() %>%
  spread_rvars(y[i, m], mu[m], nu[m], phi) %>%
  filter(i == max(i)) %>%
  mutate(
    t_step = t[9] - t[8],
    y_jump = rvar_rgamma(
      8, 
      shape = (t_step / nu^2), 
      rate = (1 / nu^2)
    ),
    y_forecast_std = y + y_jump,
    y_forecast = y_forecast_std * mu ,
    y_forecast_noisy = rvar_rTruncT(
      8,
      nu = 10,
      mean = y_forecast,
      sd = phi * y_forecast,
      lb = 0
    )
  ) %>% 
  summarise(z_forecast = y_forecast_noisy %**% t(fine_basis_mat)) %>% 
  .$z_forecast %>%
  as.data.frame() %>%
  pivot_longer(everything()) %>%
  mutate(
    n = seq(0, 21, length.out = nrow(fine_basis_mat))
  ) %>%
  ggplot(aes(x = n, dist = value)) +
  stat_dist_lineribbon() +
  scale_fill_brewer() +
  geom_point(
    data = data.frame(
      n = 1:20,
      y_true = Ut_matrix[9, ],
      value = NA
    ),
    aes(y = y_true)
  ) +
  ylim(0, 30)  +
  xlab("measurement location") +
  ylab("wear (mm)") +
  my_theme()

ft_dist_9 <- readRDS(
  file.path(".", "images", "ft-9-obs.rds")
)

p_failure_time <- ft_dist_9 %>% 
  arrange(.draw, ft) %>%
  mutate(
    CDF = seq(0, 1, length.out = n())
  ) %>%
  group_by(CDF) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  ggplot(aes(x = ft, y = CDF, xmin = .lower, xmax = .upper)) +
  geom_lineribbon() +
  scale_fill_brewer() +
  xlim(0, 2.6) +
  xlab("failure time") +
  theme(legend.position = "bottom") +
  my_theme()

path <- file.path(".", "images")

saveRDS(
  p_coef_pathways,
  file.path(path, "p_coef_pathways.rds")
)
saveRDS(
  p_prior_predictive_checks,
  file.path(path, "p_prior_predictive_checks.rds")
)
saveRDS(
  p_posterior,
  file.path(path, "p_posterior.rds")
)
saveRDS(
  p_post_pc,
  file.path(path, "p_post_pc.rds")
)
saveRDS(
  p_forecast,
  file.path(path, "p_forecast.rds")
)
saveRDS(
  p_failure_time,
  file.path(path, "p_failure_time.rds")
)

jpeg("p_posterior.jpg", width = 1500, height = 1000)
p_posterior +
theme(legend.position = "none")
dev.off()

jpeg("p_post_pc.jpg", width = 1500, height = 1000)
p_post_pc
dev.off()

jpeg("p_forecast.jpg", width = 750, height = 1000)
p_forecast
dev.off()

jpeg("p_failure_time.jpg", width = 750, height = 1000)
p_failure_time
dev.off()