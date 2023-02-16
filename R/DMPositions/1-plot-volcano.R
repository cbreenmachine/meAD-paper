
suppressPackageStartupMessages({
  library(viridis)
  library(tidyverse)
  library(cowplot)
})

source("config.R"); gc()

# Base width and height to save as (inches)
BW <- 6
BH <- 6

# Directories
data <- data.table::fread(pvals.file, verbose = F)
data$lfdr <- data$lfdr.from.ss

######################################
######### PLOT VOLCANO ###############
######################################
plot_volcano <- function(data){
  # For thinning non-sig points
  set.seed(919)

  ix <- sample(which(data$y < -log10(0.2)))
  ix.to.ignore <- ix[1:floor(length(ix) * 0.9)]

  # Subset
  sub.df <- data[-ix.to.ignore, ]

  # Handle colors
  sub.df$color <- "grey"
  sub.df$color[sub.df$lfdr < LFDR.CUT & sub.df$pi.diff > 0] <- "hyper"
  sub.df$color[sub.df$lfdr < LFDR.CUT & sub.df$pi.diff < 0] <- "hypo"

  # Order of color and order of pallete
  sub.df$color <- factor(sub.df$color, levels = c("grey", "hyper", "hypo"))
  my.pal <- c("grey", C.HYPER, C.HYPO)

  #Plot volcano
  volcano <- sub.df %>%
    ggplot(aes(x = pi.diff, y = y, color = color)) +
    geom_point(size = 0.3) +
    xlab("LOAD effect size\n(difference in predicted methylation)") +
    ylab(expression(-log[10](lFDR))) +
    scale_x_continuous(breaks = round(seq(-0.15, 0.15, 0.05),2), lim = c(-0.15, 0.15)) +
    scale_color_manual(values = my.pal) +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05)) +
    theme(panel.grid.major.y = element_line(color = "grey", linewidth = 0.25),
          panel.grid.major.x = element_line(color = "grey", linewidth = 0.25),
          axis.text.x = element_text(angle = 0),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background = element_rect(fill = "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white"))
  volcano
}

# lFDR is what we'll use in pub
data$y <- -log10(data$lfdr)
lfdr.volcano <- plot_volcano(data)
  # scale_y_continuous(breaks = c(0, 5, 10, 15, 20), lim = c(0, 15)) +


ofile <- file.path(ODIR,  "volcano-lfdr.png")
cowplot::save_plot(plot = lfdr.volcano, filename=ofile,
                   base_width = BW, base_height = BH)


# Check p-values to make sure no numerical issues
data$y <- -log10(data$p.from.ss)
pp.volcano <- plot_volcano(data)


ofile <- file.path(ODIR,  "volcano-pvals.png")
cowplot::save_plot(plot = pp.volcano, filename=ofile,
                   base_width = BW, base_height = BH)



# DSS P-values ------------------------------------------------------------

# Check p-values to make sure no numerical issues
data$y <- -log10(data$p.from.DSS)
pp.volcano <- plot_volcano(data)

ofile <- file.path(ODIR,  "volcano-DSS-pvals.png")
cowplot::save_plot(plot = pp.volcano, filename=ofile,
                   base_width = BW, base_height = BH)
#END
