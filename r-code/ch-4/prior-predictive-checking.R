library(dplyr)
library(ggplot2)
library(cowplot)

set.seed(111)
N_itr <- 100
ppd_df <- data.frame(
  t = rep(1:20, each = N_itr),
  itr = rep(1:N_itr, 20),
  rate_1 = rep(rgamma(N_itr, shape = 1, rate = 0.001), 20),
  shape_1 = rep(rgamma(N_itr, shape = 1, rate = 0.001), 20),
  rate_2 = rep(rgamma(N_itr, shape = 0.001, rate = 0.001), 20),
  shape_2 = rep(rgamma(N_itr, shape = 0.001, rate = 0.001), 20),
  mu = rep(rnorm(N_itr, 1, 0.5), 20),
  nu = rep(abs(rt(N_itr, df = 3)), 20)
) %>%
  mutate(
    jumps_1 = rgamma(n(), shape = shape_1, rate = rate_1),
    jumps_2 = rgamma(n(), shape = shape_2, rate = rate_2),
    jumps_3 = rgamma(n(), shape = 1 / (nu^2), rate = 1 / (mu * nu^2))
  ) %>%
  group_by(itr) %>%
  mutate(
    GP_1 = cumsum(jumps_1),
    GP_2 = cumsum(jumps_2),
    GP_3 = cumsum(jumps_3)
  )

p_ppd_1 <- ppd_df %>%
  ggplot() +
  geom_line(
    aes(x = t, y = GP_1, group = itr),
    alpha = 0.2
  ) +
  xlim(0, 20) +
  ylim(0, 50) +
  xlab("time (arbitrary)") +
  ylab("degradation (arbitrary)") +
  theme_minimal()
  
p_ppd_2 <- ppd_df %>%
  ggplot() +
  geom_line(
    aes(x = t, y = GP_2, group = itr),
    alpha = 0.2
  ) +
  xlim(0, 20) +
  ylim(0, 50) +
  xlab("time (arbitrary)") +
  ylab("degradation (arbitrary)") +
  theme_minimal()

p_ppd_3 <- ppd_df %>%
  ggplot() +
  geom_line(
    aes(x = t, y = GP_3, group = itr),
    alpha = 0.2
  ) +
  xlim(0, 20) +
  ylim(0, 50) +
  xlab("time (arbitrary)") +
  ylab("degradation (arbitrary)") +
  theme_minimal()

pdf(
  file.path("..", "..", "figures", "ch-4", "PPCs.pdf"),
  height = 4, width = 10
)
plot_grid(
  p_ppd_1, p_ppd_2, p_ppd_3,
  labels = c("(a)", "(b)", "(c)"),
  label_fontfamily = "Times",
  label_face = "plain",
  nrow = 1
)
dev.off()