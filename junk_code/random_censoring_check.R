library(ggplot2)
library(dplyr)
library(cowplot)

# Sample lifetimes
n_obs <- 100
beta <- 1.15
eta <- 1100

lifetimes <- rweibull(
  n_obs,
  shape = beta,
  scale = eta
) %>% round()

# Create censoring times
censoring_times <- rexp(10, 1/1000)#runif(n_obs, 0, 1500)
censoring_times[censoring_times > 1500] <- 1500

# DF of true values
df_true <- data.frame(
  lifetimes = lifetimes,
  ecdf = rank(lifetimes) / n_obs
)

# DF for ecdf
df_ecdf <- data.frame(
  lifetimes = ifelse(
    lifetimes < censoring_times,
    lifetimes,
    censoring_times
  ),
  censoring_ind = lifetimes > censoring_times
) %>%
  arrange(lifetimes, censoring_ind) %>%
  mutate(
    ind = n():1
  ) %>%
  filter(!censoring_ind) %>%
  mutate(
    q = 1 - (1 / ind),
    ecdf = 1 - cumprod(q)
  )

# Worst case with censoring
df_worst <- data.frame(
  lifetimes = ifelse(
    lifetimes < censoring_times,
    lifetimes,
    censoring_times
  )
) %>%
  mutate(
    ecdf = rank(lifetimes) / n_obs
  )

# best case with censoring
df_best <- data.frame(
  lifetimes = ifelse(
    lifetimes < censoring_times,
    lifetimes,
    6000
  )
) %>%
  mutate(
    ecdf = rank(lifetimes) / n_obs
  )

# ML estimated parameters
censored_lifetimes <- ifelse(
  lifetimes < censoring_times,
  lifetimes,
  censoring_times
)
censoring_indicator <- lifetimes > censoring_times

ll_weibull <- function(params) {
  beta = params[1]
  eta = params[2]
  obs_lifetimes <- censored_lifetimes[!censoring_indicator]
  cense_lifetimes <- censored_lifetimes[censoring_indicator]

  ll <- sum(log(dweibull(obs_lifetimes, shape = beta, scale = eta))) +
    sum(log(1 - pweibull(cense_lifetimes, shape = beta, scale = eta)))
  
  return(-ll)
}

ll_fit <- optim(c(0.5, 1000), ll_weibull)

# Plot the ECDFs
ggplot() +
  geom_step(
    data = df_ecdf,
    aes(x = lifetimes, y = ecdf)
  ) +
  geom_step(
    data = df_true,
    aes(x = lifetimes, y = ecdf),
    colour = "red"
  ) +
  geom_step(
    data = df_worst,
    aes(x = lifetimes, y = ecdf),
    linetype = 2
  ) +
  geom_step(
    data = df_best,
    aes(x = lifetimes, y = ecdf),
    linetype = 2
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = beta, scale = eta),
    colour = "green"
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = ll_fit$par[1], scale = ll_fit$par[2]),
    colour = "blue"
  ) +
  xlim(0, 5500) +
  theme_minimal()


idler_data <- readRDS("./data/idler_frame_life_example.RDS")
idler_data %>% head()
idler_data %>%
  filter(censored_at_start | censored_at_end) %>%
  ggplot(aes(lifetime)) +
  geom_histogram()

df_ecdf <- data.frame(
  lifetimes = idler_data$lifetime,
  censoring_ind = (idler_data$censored_at_start | idler_data$censored_at_end)
) %>%
  arrange(lifetimes, censoring_ind) %>%
  mutate(
    ind = n():1
  ) %>%
  filter(!censoring_ind) %>%
  mutate(
    q = 1 - (1 / ind),
    ecdf = 1 - cumprod(q)
  )

# Worst case with censoring
df_worst <- data.frame(
  lifetimes = idler_data$lifetime
) %>%
  mutate(
    ecdf = rank(lifetimes) / nrow(idler_data)
  )

# best case with censoring
df_best <- data.frame(
  lifetimes = ifelse(
    (idler_data$censored_at_start | idler_data$censored_at_end),
    1500,
    idler_data$lifetime
  )
) %>%
  mutate(
    ecdf = rank(lifetimes) / nrow(idler_data)
  )


censored_lifetimes <- as.numeric(idler_data$lifetime)
censoring_indicator <- (idler_data$censored_at_start | idler_data$censored_at_end)
ll_fit <- optim(c(0.5, 1000), ll_weibull)

ggplot() +
  geom_step(
    data = df_ecdf,
    aes(x = lifetimes, y = ecdf)
  ) +
  geom_function(
    fun = pweibull,
    args = list(shape = ll_fit$par[1], scale = ll_fit$par[2]),
    colour = "blue"
  ) +
  geom_step(
    data = df_worst,
    aes(x = lifetimes, y = ecdf),
    linetype = 2
  ) +
  geom_step(
    data = df_best,
    aes(x = lifetimes, y = ecdf),
    linetype = 2
  ) +
  theme_minimal()


censored_lifetimes <- as.numeric(df_best$lifetimes)
censoring_indicator <- rep(FALSE, nrow(idler_data))
ll_fit <- optim(c(0.5, 1000), ll_weibull)
ll_fit$par

censored_lifetimes <- as.numeric(df_worst$lifetimes)
censoring_indicator <- rep(FALSE, nrow(idler_data))
ll_fit <- optim(c(0.5, 1000), ll_weibull)
ll_fit$par

ll_weibull(params = c(1.1, 1000))


beta_min <- 0.1
beta_max <- 5
beta_int <- 0.01

eta_min <- 100
eta_max <- 2500
eta_int <- 10

beta_seq <- seq(from = beta_min, to = beta_max, by = beta_int)
eta_seq <- seq(from = eta_min, to = eta_max, by = eta_int)

grid <- expand.grid(x = beta_seq, y = eta_seq)

ll_grid <- lapply(
  1:nrow(grid),
  function(i) {
    -ll_weibull(params = c(grid[i, 1], grid[i, 2]))
  }
) %>% unlist()

data.frame(
  x = grid[, 1],
  y = grid[, 2],
  ll = ll_grid
) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = ll)) +
  scale_fill_distiller(
    palette = "Spectral",
    direction = -1,
    limits = c(-500, -350)
  )


beta <- 1.05
eta <- 1200

n <- 500
m <- 100

lower_bound <- 3500
upper_bound <- 5000
junk_df <- data.frame(
  unit = rep(1:n, each = m),
  lifetimes = rweibull(n * m, shape = beta, scale = eta)
) %>% 
group_by(unit) %>%
mutate(
  failure_times = cumsum(lifetimes),
  install_times = lag(failure_times),
  install_times = ifelse(is.na(install_times), 0, install_times)
) %>% 
ungroup() %>%
filter(
  between(failure_times, lower_bound, upper_bound) | between(install_times, lower_bound, upper_bound)
) %>%
mutate(
  censored = ifelse(
    !(between(failure_times, lower_bound, upper_bound) & between(install_times, lower_bound, upper_bound)),
    TRUE,
    FALSE
  ),
  failure_times = ifelse(
    (failure_times > upper_bound),
    upper_bound,
    failure_times
  ),
  install_times = ifelse(
    (install_times < lower_bound),
    lower_bound,
    install_times
  ),
  obs_lifetime = failure_times - install_times
)

censored_lifetimes <- junk_df$obs_lifetime
censoring_indicator <- junk_df$censored
table(censoring_indicator)

ll_fit <- optim(c(1, 1000), ll_weibull)
ll_fit$par


set.seed(2468)
small_sim <- data.frame(
    unit = factor(
      rep(1:3, each = 10),
      c(1, 2, 3)
    ),
    lifetime = rweibull(3 * 10, shape = 1.1, scale = 1200)
  ) %>% 
  group_by(unit) %>%
  mutate(
    failure_time = cumsum(lifetime),
    install_time = lag(failure_time),
    lifetime_id = 1:n()
  ) %>%
  ungroup() %>%
  replace(is.na(.), 0) %>%
  filter(failure_time < 10000) %>%
  mutate(
    observed_failure = between(failure_time, 5000, 7500)
  )

p_2 <- small_sim %>%
  ggplot() +
  geom_point(
    aes(x = failure_time, y = unit, colour = observed_failure),
    shape = 4, size = 4
  ) +
  scale_colour_manual(values = c("grey", "black")) +
  geom_segment(
    data = data.frame(
      x = rep(0, 3), y = c(3, 2, 1),
      xend = rep(5000, 3), yend = c(3, 2, 1)
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "gray",
    linetype = 2
  ) +
  geom_segment(
    data = data.frame(
      x = rep(5000, 3), y = c(3, 2, 1),
      xend = rep(7500, 3), yend = c(3, 2, 1)
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "black"
  ) +
  geom_segment(
    data = data.frame(
      x = rep(7500, 3), y = c(3, 2, 1),
      xend = rep(10000, 3), yend = c(3, 2, 1)
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "gray",
    linetype = 2
  ) +
  geom_vline(
    xintercept = c(5000, 7500),
    colour = "red"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  annotate(
    "text",
    x = 5000, y = 3.2,
    label = "start of observation",
    hjust = -0.1,
    colour = "red"
  ) +
  annotate(
    "text",
    x = 7500, y = 3.2,
    label = "end of observation",
    hjust = -0.1,
    colour = "red"
  ) +
  xlab("operational time")

p_1 <- small_sim %>%
  filter(
    ((unit == 3) & (lifetime_id == 5)) |
    ((unit == 2) & (lifetime_id == 2)) |
    ((unit == 1) & (lifetime_id == 8))
  ) %>%
  mutate(
    observed_install = install_time > 5000
  ) %>%
  ggplot() +
  geom_point(
    aes(x = failure_time, y = unit, colour = observed_failure),
    shape = 4, size = 4
  ) +
  geom_point(
    aes(x = install_time, y = unit, colour = observed_install),
    shape = 19, size = 4
  ) +
  scale_colour_manual(values = c("grey", "black")) +
  geom_vline(
    xintercept = c(5000, 7500),
    colour = "red"
  ) +
  geom_segment(
    data = data.frame(
      x = c(5923, 5000, 6499), y = c(3, 2, 1),
      xend = c(7416, 7430, 7500), yend = c(3, 2, 1)
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "black",
    linetype = 1
  ) +
  geom_segment(
    data = data.frame(
      x = c(3212, 7500), y = c(2, 1),
      xend = c(5000, 8771), yend = c(2, 1)
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "grey",
    linetype = 2
  ) +
  geom_point(
    data = data.frame(
      x = c(5000, 7500),
      y = c(2, 1)
    ),
    aes(x = x, y = y),
    shape = 1, size  = 4
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  annotate(
    "text",
    x = 5000, y = 3.2,
    label = "start of observation",
    hjust = -0.1,
    colour = "red"
  ) +
  annotate(
    "text",
    x = 7500, y = 3.2,
    label = "end of observation",
    hjust = -0.1,
    colour = "red"
  ) +
  xlab("operational time")

pdf("")
plot_grid(p_1, p_2, nrow = 1, labels = "AUTO")
