library(ggplot2)

hweibull <- function(x, shape, scale) {
  h <- (shape / scale) * ((x / scale) ^ (shape - 1))
  return(h)
}

pdf(
  file.path("..", "..", "figures", "hazard_func_demo.pdf"),
  height = 3.5, width = 7 
)
ggplot() +
  xlim(0.1, 5) +
  geom_function(
    fun = hweibull,
    args = list(shape = 0.9, scale = 2),
    aes(colour = "shape = 0.9")
  ) +
  geom_function(
    fun = hweibull,
    args = list(shape = 1, scale = 2),
    aes(colour = "shape = 1")
  ) +
  geom_function(
    fun = hweibull,
    args = list(shape = 1.1, scale = 2),
    aes(colour = "shape = 1.1")
  ) +
  scale_colour_manual(
    values = c("shape = 0.9" = "green", "shape = 1" = "blue", "shape = 1.1" = "red"),
    name = "Shape Parameter",
    guide = guide_legend(override.aes = list(linetype = "solid"))
  ) +
  theme_minimal() +
  labs(y = "Hazard", x = "Exposure")
dev.off()

eta <- 12
ggplot() +
  xlim(0.1, 5) +
  geom_function(
    fun = hweibull,
    args = list(shape = 0.75, scale = eta),
    aes(colour = "shape = 0.9")
  ) +
  geom_function(
    fun = hweibull,
    args = list(shape = 1, scale = eta),
    aes(colour = "shape = 1")
  ) +
  geom_function(
    fun = hweibull,
    args = list(shape = 1.5, scale = eta),
    aes(colour = "shape = 1.1")
  ) +
  scale_colour_manual(
    values = c("shape = 0.9" = "green", "shape = 1" = "blue", "shape = 1.1" = "red"),
    name = "Shape Parameter",
    guide = guide_legend(override.aes = list(linetype = "solid"))
  ) +
  theme_minimal() +
  labs(y = "Hazard", x = "Exposure")
