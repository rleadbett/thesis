library(dplyr)
library(ggplot2)
library(cowplot)

idler_lifetime_data <- readRDS(
  file.path("..", "data", "idler_frame_life_example.RDS")
)

# Plot the idler lifetimes by frame
pdf(
  file.path("..", "figures", "idler_frame_lifetimes.pdf"),
  height = 5,
  width = 7
)

p_idler_lifetimes <- idler_lifetime_data %>%
  mutate(
    frame_number = as.numeric(frame_number),
    censored = ((censored_at_start + censored_at_end) > 0),
    log_lifetime = log(as.numeric(lifetime))
  ) %>%
  ggplot(aes(x = frame_number, y = log_lifetime, col = censored, shape = censored)) +
  geom_point() +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  xlab("frame number") +
  ylab("lifetime (log(days))")

p_idler_lifetimes

dev.off()

saveRDS(
  p_idler_lifetimes,
  file.path("..", "figures", "idler_frame_lifetimes.RDS")
)

# Plot the empirical CDF of the idler frames
pdf(
  file.path("..", "figures", "idler_eCDF.pdf"),
  height = 5.2,
  width = 7
)

p_idler_eCDF <- idler_lifetime_data %>%
  mutate(
    censored = ((censored_at_start + censored_at_end) > 0),
    lifetime = as.numeric(lifetime)
  ) %>%
  arrange(
    censored, lifetime
  ) %>%
  mutate(
    eCDF = (1:n()) / n()
  ) %>%
  ggplot(aes(x = lifetime, y = eCDF, shape = censored, col = censored)) +
  geom_point() +
  theme_minimal() +
  ylab("density") +
  xlab("lifetime (days)") +
  labs("censored lifetime") +
  theme(
    legend.position = "bottom",
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

p_idler_eCDF

dev.off()

saveRDS(
  p_idler_eCDF,
  file.path("..", "figures", "idler_eCDF.RDS")
)

# Plot together

pdf(
  file.path("..", "figures", "idler_data_desc.pdf"),
  height = 8,
  width = 7
)

p_idler_data_desc <- plot_grid(
  p_idler_lifetimes,
  p_idler_eCDF,
  nrow = 2,
  rel_heights = c(5, 5.2),
  labels = c("(a)", "(b)"),
  greedy = FALSE
)

p_idler_data_desc

dev.off()

saveRDS(
  p_idler_data_desc,
  file.path("..", "figures", "idler_data_desc.RDS")
)
