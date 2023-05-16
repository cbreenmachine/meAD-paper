library(tidyverse)

df <- read_csv("../../DataDerived/DMRsParamSweep.csv", show=F)

p <- df %>%
  dplyr::mutate(MinCG = factor(MinCG)) %>%
  ggplot(aes(x = PercentSig, y = N.DMRs, color = MinCG)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  ylim(c(0, 30)) +
  ylab("Number of DMRs found") +
  xlab("Percent significant parameter") +
  labs(color = "Minimum CpGs\nparameter") +
  theme(plot.background = element_rect(fill = "white", color = "white"))


cowplot::save_plot(filename = "../../Figs/ParamSweep-DMRs.png", p)
