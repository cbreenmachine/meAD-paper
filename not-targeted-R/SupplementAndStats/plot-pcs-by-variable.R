# Plot four (or more) PCs
# Need to add variance explained
library(tidyverse)
library(cowplot)

load("../../DataRaw/2023-01-16-ImputationCheck/input-chr22.RData")

vv <- "sex"

df$diagnostic_group <- ifelse(df$diagnostic_group == "CONTROL", "Control", "LOAD")

plot_pcs <- function(vv, label.name){
  df %>%
    ggplot(aes(x = PC1, y = PC2, color = .data[[vv]])) +
    geom_point(alpha = 0.8) +
    theme_minimal() +
    labs(color = label.name) +
    theme(legend.position = c(0.75, 0.825),
          legend.background = element_rect(fill = alpha("white", 0.8),
                                           color = NA)) +
    xlim(c(-25, 75)) +
    ylim(c(-25, 35)) %>%
    return()
}


p1 <- plot_pcs("source", "Source") +
  scale_color_manual(values = c("#E1BE61", "#40B0A6"))

p2 <- plot_pcs("diagnostic_group", "Diagnostic group") +
  scale_color_manual(values = c("#994F00", "#006CD1"))

p3 <- plot_pcs("sex", "Sex") +
  scale_color_manual(values = c("#E66100", "#5D3A9B"))

p4 <- plot_pcs("race_primary", "Race") +
  scale_color_manual(values = c("#000000", "#e69f00", "#56b4e9"))


p <- cowplot::plot_grid(p1, p2, p3, p4, nrow = 2) +
  theme(panel.background = element_rect(fill = "white", color = "white")) +
  panel_border(remove = TRUE)


cowplot::save_plot(filename = "../../Figs/PCs/PanelPCs.png", p,
                   base_width = 10, base_height = 6)
