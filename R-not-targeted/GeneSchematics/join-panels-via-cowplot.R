
library(cowplot)

path1 <- "../../Figs/GeneSchematics/2023-03-08-panel-canonical.png"
path2 <- "../../Figs/GeneSchematics/2023-03-08-panel-highest_effect.png"

canonical <- ggdraw() + draw_image(path1)
high.effect <- ggdraw() + draw_image(path2)

stitched <- plot_grid(canonical, NULL, high.effect,
                      nrow = 1,
                      labels = c("A","", "B"),
                      rel_widths = c(1, -0.05, 1),
                      label_fontfamily = "Times",
                      label_colour = "black") +
  theme(plot.background = element_rect(fill = "white", color = "white"))

cowplot::save_plot(stitched, filename = "../../Figs/GeneSchematics/gene-schematics-stitched.png")
