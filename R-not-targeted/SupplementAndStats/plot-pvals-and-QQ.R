# plot-pvals-and-QQ.R
# Three plots
# (A) Histogram of raw p-values from DSS
# (B) Histogram of 'adjusted' p-values after accounting for over-dispersion
# (C) QQ-plot with lambda

library(data.table)
library(fastqq)
library(cowplot)
library(fdrtool)

df <- fread("../../DataRaw/2023-02-14-Summaries-v6/pvals.bed")

# Extract the vector we need
pp.before <- df$p.from.DSS
pp.after <- df$p.from.ss

# Lambda_GC (genomic control)
# http://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html
lambda.before <- round(median(qchisq(1 - pp.before, 1)) / qchisq(0.5, 1), 1)
lambda.after <- round(median(qchisq(1 - pp.after, 1)) / qchisq(0.5, 1), 1)

hist(pp.before)
hist(pp.after)

plot_histogram <- function(xx, ti){
  data.frame(xx) %>%
    ggplot(aes(x = xx, y = after_stat(..density..))) +
    geom_histogram(bins = 50, fill = "darkgrey", color = "black", linewidth=0.1) +
    xlim(c(0,1)) +
    ylim(c(0, 1.7)) +
    xlab("P-values") +
    ylab("Density") +
    theme_minimal() +
    labs(subtitle = ti) +
    theme(plot.background = element_rect(fill = "white", color= "white"))
}

lambda.before
lambda.after

ti.before <- expression(paste("Raw P-values (", lambda[GC], " = 0.8)"))
ti.after <- expression(paste("Modeled P-values (", lambda[GC], " = 1.2)"))

plot.before <- plot_histogram(pp.before, ti.before)
plot.before
plot.after <- plot_histogram(pp.after, ti.after)

z <- cowplot::plot_grid(plot.before, plot.after,
                  labels = c("A", "B"))

cowplot::save_plot('../../Figs/histogram-pvals-before-after.png', z)


#END
