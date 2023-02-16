
library(data.table)
library(fastqq)
library(cowplot)

df <- fread("../../DataRaw/2023-01-24-Summaries-v4/pvals.bed")

# Plot plot plot
qq.plot <- "../../Figs/QQ-plot.png"

# Extract the vector we need
pp <- df$p.raw
ss <- df$stat
zz <- (ss - mean(ss)) / sd(ss)

zz.out <- fdrtool::fdrtool(zz, statistic = "normal", plot = F)
pp.out <- fdrtool::fdrtool(pp, statistic = "pvalue", plot = F)

hist(pp.out$pval)
hist(zz.out$pval)


chisq.stat <- (df$stat) ^ 2

# Lambda_GC (genomic control)
# http://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html
lambda <- median(chisq.stat) / qchisq(0.5, 1)

chisq.stat.crct <- chisq.stat / lambda

# From
stat.crct <- sign(df$stat) * chisq.stat.crct
pvals.crct = 2*pnorm(-abs(stat.crct))

hist(pvals.crct)


# FDR Tool --------------------------------------------------------------

out <- fdrtool::fdrtool(x = df$stat, statistic = "normal")
hist(out$pval)

# What is the genomic inflation after FDRtool?
out$pval

#


# QQ Plot -----------------------------------------------------------------


png(qq.plot,
    width = 900, height = 600,
    res = 120)
fastqq::qq(pp)
dev.off()



# Histogram ---------------------------------------------------------------

z <- df %>%
  ggplot(aes(x = p.raw)) +
  geom_histogram(bins = 100) +
  theme_minimal() +
  xlab("P-value")

p <- ggdraw() + draw_image(qq.plot)
out <- cowplot::plot_grid(p, z)

cowplot::save_plot(filename = "../../Figs/qq-with-pval-histogram.png", out)
