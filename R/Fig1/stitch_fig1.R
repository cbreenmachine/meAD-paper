library(cowplot)
source("config.R")

volcano.file <- file.path(ODIR, "volcano.png")
sankey.file <- file.path(ODIR, "sankey.png")
pi.file <- file.path(ODIR, "pi.png")
cpg.file <- file.path(ODIR, "cpg_barchart.png")


volcano <- ggdraw() + draw_image(volcano.file)
sankey <- ggdraw() + draw_image(sankey.file)
pi <- ggdraw() + draw_image(pi.file)
cpg <- ggdraw() + draw_image(cpg.file)

topright <- plot_grid(pi, cpg, rel_widths = c(0.3, 0.7),
                      labels = c("B", "C"))

right <- plot_grid(topright, sankey,
                   ncol = 1,
                   rel_heights = c(0.5, 0.5),
                   labels = c(NA, "D"))
right


all <- plot_grid(volcano, right,
                 ncol = 2,
                 rel_widths = c(0.5, 0.5),
                 labels = c("A", NA)) +
  theme(plot.background = element_rect(fill = "white"))

save_plot(plot = all,
          file.path(ODIR, "fig1.png"),
          base_width = 8, base_height = 4)
