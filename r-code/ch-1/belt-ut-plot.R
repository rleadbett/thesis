library(dplyr)
library(ggplot2)

load(
  file.path("..", "..", "data", "Example_beltwear_data.RData")
)

belt_ut_data <- example_belts[["belt_A"]]

pdf(
  file.path("..", "..", "figures", "ch-1", "belt_wear_ut.pdf"),
  height = 6,
  width = 8
)

p_belt_ut_data <- belt_ut_data %>%
  mutate(
    tonnes = as.factor(round(tonnes, 2))
  ) %>%
  ggplot(
    aes(
      x = measurement_pos,
      y = wear,
      col = tonnes,
      group = tonnes
    )
  ) +
  geom_line() +
  scale_color_brewer(palette = "YlOrRd") +
  theme_minimal() +
  ylim(-1, 30) +
  xlab("measurement location") +
  ylab("wear (mm)") +
  labs(col = "cumulative tonnage")

p_belt_ut_data

dev.off()
