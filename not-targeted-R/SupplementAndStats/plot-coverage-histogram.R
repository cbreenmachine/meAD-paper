library(tidyverse)

df <- read_csv("../../DataDerived/2023-02-01-AverageCoverage.csv", show_col_types = F)

df$weight <- df$N.positions / sum(df$N.positions)

df$weight %*% df$`100`

cov <- df %>%
  dplyr::select(-N.positions) %>%
  pivot_longer(-weight) %>%
  group_by(name) %>%
  summarize(Cov.bar = weight %*% value)

names(cov)[2] <- 'Mean effective coverage'

p <- cov %>%
  ggplot(aes(x = `Mean effective coverage`, y = after_stat(density))) +
  geom_histogram(bins = 25, fill = "black") +
  theme_bw() +
  ylab("Density") +
  xlim(c(20, 55))


cowplot::save_plot(filename = "../../Figs/EffectiveCoverage.png", p)
