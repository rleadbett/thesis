library(dplyr)
library(ggplot2)

# Censoring example figure
set.seed(123)
y_obs <- rweibull(3, 1.1, 1)

pdf(
  file.path("..", "..", "figures", "censoring_example.pdf"),
  width = 10,
  height = 6
)

data.frame(
  failure_times = y_obs,
  install_times = 0
) %>%
arrange(y_obs) %>%
mutate(
  units = factor(
    c("left censoring", "interval censoring", "right censoring"),
    levels = c("left censoring", "interval censoring", "right censoring"),
    ordered = TRUE
  ),
  observation_end = c(0, 0.5, 1)
) %>%
ggplot() +
geom_point(
  aes(x = failure_times, y = units),
  shape = 4
) +
geom_vline(xintercept = 0.5, colour = "red", linetype = "dashed") +
geom_vline(xintercept = 1, colour = "red", linetype = "dashed") +
geom_segment(
  aes(x = install_times, y = units, xend = observation_end, yend = units),
  colour = "black",
  linetype = 1
) +
geom_segment(
  aes(x = observation_end, y = units, xend = failure_times, yend = units),
  colour = "grey",
  linetype = 2
) +
geom_point(
  aes(x = install_times, y = units),
  shape = 16
) +
theme_minimal() +
xlab("time") +
ylab("") +
annotate(
  "text",
  x = 0.51, y = 0.5,
  label = "t_1",
  hjust = -0.01,
  size = 2.5,
  colour = "red"
) +
annotate(
  "text",
  x = 1.01, y = 0.5,
  label = "t_2",
  hjust = -0.01,
  size = 2.5,
  colour = "red"
)

dev.off()

# Left-truncation example
set.seed(112)

left_trunc_df <- data.frame(
  y = rweibull(10, 1.1, 1),
  trunc_time = 0.5
) %>%
  mutate(
    observed = y > trunc_time,
    origin = 0,
    partial_obs = ifelse(observed, trunc_time, y),
    observed_start = ifelse(observed, trunc_time, NA),
    observed_end = ifelse(observed, y, NA)
  )

p_left_cens <- left_trunc_df %>%
mutate(n = factor(1:n(), levels = 1:10, ordered = T)) %>%
ggplot() +
geom_vline(xintercept = 0.5, colour = "red", linetype = "dashed") +
geom_segment(
  aes(x = observed_start, y = n, xend = observed_end, yend = n),
  colour = "black",
  linetype = 1
) +
geom_segment(
  aes(x = origin, y = n, xend = partial_obs, yend = n),
  colour = "grey",
  linetype = 2
) +
geom_point(
  data = left_trunc_df %>%
    mutate(n = factor(1:n(), levels = 1:10, ordered = T)) %>%
    filter(observed),
  aes(x = y, y = n),
  shape = 4
) +
geom_point(
  data = left_trunc_df %>%
    mutate(n = factor(1:n(), levels = 1:10, ordered = T)) %>%
    filter(!observed),
  aes(x = y, y = n),
  shape = 4,
  colour = "gray"
) +
geom_point(
  aes(x = origin, y = n),
  shape = 16,
  colour = "gray"
) +
theme_minimal() +
xlab("age") +
ylab("")

p_left_trunc <- left_trunc_df %>%
filter(observed) %>%
mutate(n = factor(1:n(), levels = 1:10, ordered = T)) %>%
ggplot() +
geom_vline(xintercept = 0.5, colour = "red", linetype = "dashed") +
geom_segment(
  aes(x = observed_start, y = n, xend = observed_end, yend = n),
  colour = "black",
  linetype = 1
) +
geom_segment(
  aes(x = origin, y = n, xend = partial_obs, yend = n),
  colour = "grey",
  linetype = 2
) +
geom_point(
  data = left_trunc_df %>%
    filter(observed) %>%
    mutate(n = factor(1:n(), levels = 1:10, ordered = T)),
  aes(x = y, y = n),
  shape = 4
) +
geom_point(
  aes(x = origin, y = n),
  shape = 16,
  colour = "gray"
) +
theme_minimal() +
xlab("age") +
ylab("")


pdf(
  file.path("..", "..", "figures", "left_truncation_example.pdf"),
  width = 10,
  height = 6
)
plot_grid(
  p_left_cens,
  p_left_trunc,
  ncol = 2, nrow = 1,
  labels = c(
    "(a)", "(b)"
  ),
  label_fontfamily = "Times",
  label_face = "plain"
)
dev.off()
data.frame(
  y_obs = y_obs,
  install_times = runif(10, 0, 2),
  unit = factor(
    as.character(1:10), levels = as.character(1:10), ordered = TRUE
  ),
  trunc_time = 0.5
) %>%
mutate(
  failure_times = install_times + y_obs,
  trunc_time = ifelse(install_times > trunc_time, install_times, trunc_time),
  trunc_time = ifelse(failure_times < trunc_time, failure_times, trunc_time)
) %>%
ggplot() +
geom_vline(xintercept = 0.5, colour = "red", linetype = "dashed") +
geom_segment(
  aes(x = trunc_time, y = unit, xend = failure_times, yend = unit),
  colour = "black",
  linetype = 1
) +
geom_segment(
  aes(x = install_times, y = unit, xend = trunc_time, yend = unit),
  colour = "grey",
  linetype = 2
) +
geom_point(
  aes(x = failure_times, y = unit),
  shape = 4
) +
geom_point(
  aes(x = install_times, y = unit),
  shape = 16
) +
annotate(
  "text",
  x = 0.51, y = 0.5,
  label = "t_start",
  hjust = -0.01,
  size = 2.5,
  colour = "red"
) +
theme_minimal() +
xlab("time")
# Left-truncation and right censoring plot

set.seed(111)
y_obs <- rweibull(3 * 10, 1.1, 1)
t_start <- 2
t_end <- 4
example_data <- data.frame(
  unit = factor(
    rep(c("unit1", "unit2", "unit3"), 10),
    levels = c("unit1", "unit2", "unit3"),
    ordered = TRUE
  ),
  lifetimes = y_obs
) %>%
group_by(unit) %>%
mutate(
  failure_times = cumsum(lifetimes),
  install_times = lag(failure_times),
  install_times = ifelse(is.na(install_times), 0, install_times),
  censoring_times = ifelse(failure_times < t_end, failure_times, t_end),
  missing = ifelse(failure_times < t_start, TRUE, FALSE)
) %>%
filter(install_times < t_end)

pdf(
  file.path("..", "..", "figures", "left_truncation_w_right_censoring_example.pdf"),
  width = 10,
  height = 6
)

ggplot() +
geom_segment(
  data = example_data %>%
    filter(missing == TRUE),
  aes(
    x = install_times, y = unit,
    xend = censoring_times, yend = unit
  ),
  linetype = "dashed",
  colour = "grey"
) +
geom_segment(
  data = example_data %>%
    filter(missing == FALSE),
  aes(
    x = install_times, y = unit,
    xend = censoring_times, yend = unit
  ),
  colour = "black"
) +
geom_segment(
  data = example_data,
  aes(x = censoring_times, y = unit, xend = failure_times, yend = unit),
  colour = "grey",
  linetype = "dashed"
) +
geom_point(
  data = example_data %>%
    filter((missing == TRUE) | failure_times > t_end),
  aes(x = failure_times, y = unit),
  shape = 4,
  colour = "grey"
) +
geom_point(
  data = example_data %>%
    filter((missing == FALSE) & (failure_times < t_end)),
  aes(x = failure_times, y = unit),
  shape = 4
) +
geom_point(
  data = example_data %>%
    filter(missing == TRUE),
  aes(x = install_times, y = unit),
  shape = 16,
  colour = "grey"
) +
geom_point(
  data = example_data %>%
    filter(missing == FALSE),
  aes(x = install_times, y = unit),
  shape = 16
) +
geom_vline(xintercept = t_start, colour = "red", linetype = "dashed") +
geom_vline(xintercept = t_end, colour = "red", linetype = "dashed") +
annotate(
  "text",
  x = 2.05, y = 0.5,
  label = "t_start",
  hjust = -0.01,
  size = 2.5,
  colour = "red"
) +
annotate(
  "text",
  x = 4.05, y = 0.5,
  label = "t_end",
  hjust = -0.01,
  size = 2.5,
  colour = "red"
) +
theme_minimal() +
xlab("time") +
xlim(0, max(example_data$failure_times))

dev.off()