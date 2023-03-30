
library(showtext)
font_add(family = "Arial", regular = "Arial.ttf") ## here is the path to the font to add.
showtext.auto()
library(cowplot)

z <- ggdraw() +
  draw_image(magick::image_read_pdf("../Figs/UCSCBrowser/2023-03-14-B4GALT1-R01-report.pdf", density = 600))

print(z)
cowplot::save_plot("output.png", z, base_width = 1, base_height = .7, dpi = 1200)
