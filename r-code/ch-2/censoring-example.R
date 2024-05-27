library(ggplot2)
library(dplyr)
library(cowplot)

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
  geom_vline(
    xintercept = c(5000, 7500),
    colour = "red"
  ) +
  geom_point(
    aes(x = failure_time, y = unit, colour = observed_failure),
    shape = 4, size = 3
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
  theme_classic() +
  theme(legend.position = "none") +
  annotate(
    "text",
    x = 5000, y = 3.2,
    label = "start of observation",
    hjust = -0.01,
    size = 2.5,
    colour = "red"
  ) +
  annotate(
    "text",
    x = 7500, y = 3.2,
    label = "end of observation",
    hjust = -0.01,
    size = 2.5,
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
  geom_vline(
    xintercept = c(5000, 7500),
    colour = "red"
  ) +
  geom_point(
    aes(x = failure_time, y = unit, colour = observed_failure),
    shape = 4, size = 3
  ) +
  geom_point(
    aes(x = install_time, y = unit, colour = observed_install),
    shape = 19, size = 3
  ) +
  scale_colour_manual(values = c("grey", "black")) +
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
    shape = 1, size  = 3
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  annotate(
    "text",
    x = 5000, y = 3.2,
    label = "start of observation",
    hjust = -0.01,
    size = 2.5,
    colour = "red"
  ) +
  annotate(
    "text",
    x = 7500, y = 3.2,
    label = "end of observation",
    hjust = -0.01,
    size = 2.5,
    colour = "red"
  ) +
  xlab("operational time")

pdf(
  file.path("..", "..", "figures", "censoring_example.pdf"),
  width = 10,
  height = 6
)
plot_grid(p_1, p_2, nrow = 1, labels = "AUTO")
dev.off()