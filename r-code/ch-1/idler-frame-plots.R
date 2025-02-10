library(dplyr)
library(ggplot2)
library(cowplot)

theme_update(axis.title = element_text(size = 16))

idler_lifetime_data <- readRDS(
  file.path("..", "..", "data", "idler_frame_life_example.RDS")
)

# Plot the idler lifetimes by frame
p_idler_lifetimes <- idler_lifetime_data %>%
  mutate(
    frame_number = as.numeric(frame_number),
    censored = ((censored_at_start + censored_at_end) > 0),
    lifetime = as.numeric(lifetime)
  ) %>%
  ggplot(aes(x = frame_number, y = lifetime, col = censored, shape = censored)) +
  geom_point() +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  scale_y_log10() +
  xlab("frame number") +
  ylab("lifetime (days)")

# Plot the empirical CDF of the idler frames
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
  ylab("empirical cumulative distribution") +
  xlab("lifetime (days)") +
  labs("censored lifetime") +
  theme(
    legend.position = "bottom",
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

# Plot together

pdf(
  file.path("..", "..", "figures", "ch-1", "idler_data_desc.pdf"),
  height = 8,
  width = 7
)

plot_grid(
  p_idler_lifetimes,
  p_idler_eCDF,
  ncol = 1,
  nrow = 2,
  rel_heights = c(5, 5.2),
  labels = c("(a)", "(b)"),
  greedy = FALSE,
  label_fontfamily = "Times",
  label_face = "plain"
)

dev.off()
