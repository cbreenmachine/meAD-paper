
suppressPackageStartupMessages({
  library(viridis)
  library(tidyverse)
  library(cowplot)
})

source("config.R")

# Base width and height to save as (inches)
BW <- 6
BH <- 6

# Directories
data <- data.table::fread(pvals.file)

######################################
######### PLOT VOLCANO ###############
######################################

# Thin non-significant points
data$y <- -log10(data$lfdr)
ix <- sample(which(data$y < -log10(0.01)))
ix.to.ignore <- ix[1:floor(length(ix) * 0.999)]

# Subset
sub.df <- data[-ix.to.ignore, ]

# Handle colors
sub.df$color <- "grey"
sub.df$color[sub.df$lfdr < LFDR.CUT & sub.df$diagnostic_group_coded > 0] <- "hyper"
sub.df$color[sub.df$lfdr < LFDR.CUT & sub.df$diagnostic_group_coded < 0] <- "hypo"

# Order of color and order of pallete
sub.df$color <- factor(sub.df$color, levels = c("grey", "hyper", "hypo"))
my.pal <- c("grey", C.HYPER, C.HYPO)

#Plot volcano
p.volcano <- sub.df %>%
  ggplot(aes(x = pi.diff, y = y, color = color)) +
  geom_jitter(size = 0.2, alpha = 1, height=1, width=0) +
  xlab("LOAD effect size\n(difference in adjusted methylation)") +
  ylab(expression(-log[10](lFDR))) +
  scale_x_continuous(breaks = seq(-1, 1, 0.25), lim = c(-1, 1)) +
  scale_y_continuous(breaks = c(0, 100, 200, 300), lim = c(0, 350)) +
  scale_color_manual(values = my.pal) +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(color = "grey", size=0.25),
        panel.grid.major.x = element_line(color = "grey", size=0.25),
        axis.text.x = element_text(angle = 0),
        axis.ticks = element_blank(),
        legend.position = "none")

ofile <- file.path(ODIR, "volcano.png")
cowplot::save_plot(plot = p.volcano, filename=ofile,
                   base_width = BW, base_height = BH)

