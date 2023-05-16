library(showtext)
font_add(family = "Arial", regular = "Arial.ttf") ## here is the path to the font to add.
showtext.auto()
library(cowplot)

ifile.1 <- "../../Figs/PCHi-C/gene-ontology-pchic.png"
ifile.2 <- "../../Figs/PCHi-C/histogram-simulation-by-ct.png"

ifile.3 <- "../../Figs/UCSCBrowser/2023-03-21-B4GALT1-representative-interaction.pdf"

ofile <- "../../Figs/panel-pchic.png"


go <- ggdraw() + draw_image(ifile.1)
hist <- ggdraw() + draw_image(ifile.2)

interactions <- ggdraw() +
  draw_image(magick::image_read_pdf(ifile.3, density = 600))


# Start stitiching --------------------------------------------------------

bottom <- cowplot::plot_grid(go, hist, rel_widths = c(0.6, 0.4),
                             labels = c("B", "C"))

full <- cowplot::plot_grid(interactions, bottom, 
                           ncol = 1, 
                           rel_heights = c(0.5, 0.5),
                           labels = c("A", NULL)) +
  theme(plot.background = element_rect(fill = "white", color = "white"))
full


cowplot::save_plot(ofile, full)
