library(dplyr)
library(ggplot2)
library(cowplot)

set.seed(222)
N_itr <- 4
N_units <- 10
N_obs <- 9

mu <- abs(rnorm(N_itr, 0.5, 0.2))
nu <- abs(rt(N_itr, df = 3) * 0.5)

z_itrs <- lapply(
  1:N_itr,
  function(itr) {
    z_itr <- lapply(
      1:N_units,
      function(j) {
        jumps <- rgamma(
          N_obs,
          shape = (0.1 / nu[itr]^2),
          rate = (1 / (nu[itr]^2 * mu[itr]))
        )
        z <- c(0, cumsum(jumps))
        return(z)
      }
    ) %>%
      unlist()
  }
) %>%
  unlist()

ppd_df <- data.frame(
  itr = factor(rep(1:N_itr, each = ((N_obs + 1) * N_units))),
  unit = factor(rep(rep(1:N_units, each = N_obs + 1), N_itr)),
  t = rep(seq(0, 0.9  , length.out = (N_obs + 1)), (N_units * N_itr)),
  z = z_itrs
)

p_ppc_units <- ppd_df %>%
  ggplot() +
  geom_line(
    aes(x = t, y = z, colour = itr, group = interaction(itr, unit))
  ) +
  xlim(0, 1) +
  xlab("Hundreds of thousands of cycles") +
  ylab("cumulative degradation") +
  facet_wrap(vars(itr)) +
  theme_minimal()
  
pdf(
  file.path("..", "..", "figures", "ch-5", "PPC_multi_unit.pdf"),
  height = 8, width = 8
)
p_ppc_units
dev.off()