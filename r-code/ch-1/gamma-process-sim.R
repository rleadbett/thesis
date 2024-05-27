library(ggplot2)
library(cowplot)
library(dplyr)
library(posterior)
library(tidybayes)
library(ggdist)
library(bayesplot)
library(rstan)

base_dir <- file.path(".")

# Stan options
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

# Set random seed
set.seed(7319)

# Define values of parameters 
mu_true <- 1
nu_true <- 0.3
sigma_true <- 0.05

# Define time steps
delta_t <- 0.01
N <- 50

# Calculate shape and rate parameter of gamma distribution
shape_true <- delta_t / (nu_true^2)
rate_true <- 1 / (mu_true * nu_true^2)

# Sample degradation jumps
delta_z <- rgamma(
  n = N,
  shape = shape_true,
  rate = rate_true
)

# Calculate the degradation pathway from jumps
z <- cumsum(delta_z)

# Add noise
y <- rnorm(
  n = length(z),
  mean = z,
  sd = sigma_true
)

sim_data_df <- data.frame(
  time = cumsum(rep(delta_t, N)),
  true_degradation = z,
  observed_degradation = y
)

jpeg(
  filename = file.path(base_dir, "figures", "sim_GP_data.jpg"),
  width = (1950 * 2), height = (1080 * 2),
  res = 500,
  pointsize = 1000,
  quality = 2000,
  bg="transparent"
)

sim_data_df %>%
  ggplot(aes(x = time)) +
  geom_line(aes(y = true_degradation)) +
  ylab("degradation") +
  my_theme() +
  ylim(0, 0.9)

dev.off()

jpeg(
  filename = file.path(base_dir, "figures", "sim_noisy_GP_data.jpg"),
  width = (1950 * 2), height = (1080 * 2),
  res = 500,
  pointsize = 1000,
  quality = 2000,
  bg="transparent"
)

sim_data_df %>%
  ggplot(aes(x = time)) +
  geom_line(aes(y = true_degradation)) +
  geom_point(aes(y = observed_degradation)) +
  geom_line(aes(y = observed_degradation), linetype = 2) +
  ylab("degradation") +
  my_theme() +
  ylim(0, 0.9)

dev.off()

noisy_gp_stan_model <- stan_model(
  model_code = "
data{
// data
int I;         // number of observations
vector[I] t;   // the set of observation times
vector[I] y;   // the set of noisy degradation measurements

// hyper parameters for mu
real a_hat;    // our estimate of the mean wear rate
real b_hat;    // our uncertanty of this estimate
}
transformed data{
vector[I] t_diff;              // The time steps between each observation

// calculate the time steps
t_diff[1] = t[1];
t_diff[2:I] = t[2:I] - t[1:I-1];
}
parameters{
real<lower = 0> sigma;         // the standard deviation of the measurement error
real<lower = 0> mu;            // the mean wear rate of the gamma process
real<lower = 0> tau_cov;       // the reparameterisation for the CoV
real<lower = 0> alpha_cov;     // ...
vector<lower = 0>[I] z_diff;   // the degradation jumps
}
transformed parameters{
real<lower = 0> CoV;           // coefficient of variation of GP
vector<lower = 0>[I] z;        // the filtered degradation measurements

// reparameterise for CoV
CoV = (alpha_cov / sqrt(tau_cov));

// calculate the degradation measurements from the stochastic jumps
z = cumulative_sum(z_diff);
}
model{
real half_nu_cov = 3*0.5;

// priors //
// for process model
tau_cov ~ gamma(half_nu_cov, half_nu_cov);
alpha_cov ~ normal(0, 1);
mu ~ normal(a_hat, b_hat);
// for data model
sigma ~ uniform(0, 100);

// process model //
z_diff ~ gamma(t_diff / pow(CoV, 2), 1 / (mu * pow(CoV, 2)));

// data model //
y ~ normal(z, sigma);

}
generated quantities{
vector[I] y_pred;
real log_CoV;

// generate posterior predictive samples within Stan (much quicker than in R)
for(i in 1:I){
  y_pred[i] = normal_rng(z[i], sigma);
}

// used in diagnostics
log_CoV = log(CoV);
}
"
)

stan_data <- list(
  I = nrow(sim_data_df),
  t = sim_data_df$time,
  y = sim_data_df$observed_degradation,
  a_hat = 1,
  b_hat = 0.2
)

noisy_gp_stan_fit <- sampling(
    noisy_gp_stan_model,
    stan_data,
    chains = 4,
    iter = 2000,
    warmup = 500,
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 14
    )
)

p_fit <- noisy_gp_stan_fit %>%
  as_draws_rvars() %>%
  spread_rvars(z[i]) %>%
  mutate(
    time = sim_data_df$time,
    true_path = sim_data_df$true_degradation,
    noisy_observations = sim_data_df$observed_degradation
  ) %>%
  ggplot(aes(x = time, dist = z)) +
  stat_dist_lineribbon(alpha = 0.75, col = "#08306B", n = 5000) +
  scale_fill_brewer() +
  geom_line(aes(y = true_path), col = "red", size = 1) +
  geom_point(aes(y = noisy_observations)) +
  my_theme() +
  theme(legend.position = "bottom") +
  ylab("degradation")

jpeg(
  filename = file.path(base_dir, "figures", "noisy_GP_background.jpg"),
  width = (1950 * 2), height = (1080 * 2),
  res = 500,
  pointsize = 1000,
  quality = 2000,
  bg="transparent"
)
p_fit +
  geom_line(aes(y = true_path), size = 1) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "#fcfbf9"),
    plot.background = element_rect(fill = "#fcfbf9", color = NA)
  )
dev.off()

p_post_orthogonal <- noisy_gp_stan_fit %>%
  as_draws_df() %>%
  spread_draws(mu, CoV) %>%
  ggplot(aes(x = mu, y = CoV)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  my_theme() +
  theme(legend.position = "none") +
  ylab("coefficient of variation") +
  xlab("mean wear rate per unit time")

p_post_nonorthogonal <- noisy_gp_stan_fit %>%
  as_draws_df() %>%
  spread_draws(mu, CoV) %>%
  mutate(
    shape = 1 / CoV^2,
    rate = 1 / (mu * CoV^2)
  ) %>%
  ggplot(aes(x = shape, y = rate)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  my_theme() +
  theme(legend.position = "none") +
  ylab("shape") +
  xlab("rate")

jpeg(
  filename = file.path(base_dir, "figures", "noisy_GP_post_shape_rate.jpg"),
  width = 900, height = 800,
  res = 250,
  pointsize = 500,
  quality = 1000,
  bg="transparent"
)
p_post_nonorthogonal
dev.off()

jpeg(
  filename = file.path(base_dir, "figures", "noisy_GP_post_mean_cov.jpg"),
  width = 900, height = 800,
  res = 250,
  pointsize = 500,
  quality = 1000,
  bg="transparent"
)
p_post_orthogonal
dev.off()

p_post_mu_nu <- p_post_orthogonal +
  geom_point(
    x = mu_true,
    y = nu_true,
    col = "red",
    shape = 3,
    size = 4
  )

p_post_sigma <- noisy_gp_stan_fit %>%
  as_draws_df() %>%
  spread_draws(sigma) %>%
  ggplot(aes(x = sigma)) +
  geom_density(fill = "#08306B", alpha = 0.8) +
  my_theme() +
  geom_vline(xintercept = sigma_true, col = "red", size = 1) +
  xlim(0, 0.1)


p_posts <- plot_grid(
  p_post_mu_nu,
  p_post_sigma,
  ncol = 1,
  rel_heights = c(1, 0.75)
)

jpeg(
  filename = file.path(base_dir, "figures", "noisy_GP_fit.jpg"),
  width = (1950 * 2), height = (1080 * 2),
  res = 500,
  pointsize = 1000,
  quality = 2000,
  bg="transparent"
)
p <- plot_grid(p_fit, p_posts, nrow = 1, rel_widths = c(1.3, 1))
p
dev.off()
saveRDS(
  p,
  file.path(base_dir, "figures", "noisy_GP_fit.RDS")
)


