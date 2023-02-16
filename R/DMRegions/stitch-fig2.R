library(cowplot)

ODIR <- "../../Figs/Fig2/"

box.file <- file.path(ODIR, "boxplots-one-dmr.png")
go.file <- file.path(ODIR, "dmr-gene-ontology.png")
bar.file <- file.path(ODIR, "dmrs-hyper-hypo-bars.png")

box <- ggdraw() + draw_image(box.file)
go <- ggdraw() + draw_image(go.file)
bar <- ggdraw() + draw_image(bar.file) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_blank())

right <- cowplot::plot_grid(box, go,
                  ncol = 1,
                  rel_heights = c(0.5, 0.5),
                  labels = c("B", "C"),
                  label_fontfamily = "Times", label_colour = "black")

all <- plot_grid(bar, right,
                 ncol = 2,
                 rel_widths = c(0.5, 0.5),
                 labels = c("A", NA),
                 label_fontfamily = "Times", label_colour = "black") +
  panel_border(remove = TRUE)
  # theme(plot.background = element_rect(fill = "white", color = "white"),
  #       panel.background = element_blank())

save_plot(plot = all,
          file.path(ODIR, paste0(Sys.Date(), "-fig2.png")),
          base_width = 8, base_height = 4)

