
compute_lambda_gc <- function(pp){
  round(median(qchisq(1 - pp, 1)) / qchisq(0.5, 1), 1)
}


make_title_for_pvals_plots <- function(desc, lambda.gc){
  bquote(.(desc)~~"("*lambda[GC]==.(lambda.gc)*")")
}

plot_hist_from_pvals <- function(pp){

  lambda.gc <- compute_lambda_gc(pp)
  ti <- make_title_for_pvals_plots("Modeled p-values", lambda.gc)

  hist(pp,
       probability = T,
       breaks = 50,
       main = as.expression(ti),
       xlab = expression(p))

  # data.frame(pp) %>%
  #   ggplot(aes(x = pp, y = after_stat(..density..))) +
  #   geom_histogram(bins = 50, fill = "darkgrey", color = "black", linewidth=0.1) +
  #   xlim(c(0,1)) +
  #   ylim(c(0, 1.7)) +
  #   xlab("P-values") +
  #   ylab("Density") +
  #   theme_minimal() +
  #   labs(subtitle = ti) +
  #   theme(plot.background = element_rect(fill = "white", color= "white"))
}


plot_qq_from_pvals <- function(pp){

  # Genomic in/de-flation
  lambda.gc <- compute_lambda_gc(pp)

  # Title
  ti <- make_title_for_pvals_plots("Raw p-values", lambda.gc)

  fastqq::qq(pp, main = as.expression(ti))
}


plot_pvals_routine <- function(pvals.data){
  pp.1 <- pvals.data$pval.Wald
  pp.2 <- pvals.data$pval

  # Run plotting (both coalculate lambda internally)

  par(mfrow = c(1,2))
  plot_qq_from_pvals(pp.1)
  plot_hist_from_pvals(pp.2)
  #
  #
  # plot_hist_from_pvals(pp.2)
  # hist.recorded <- recordPlot()
  #
  # plot_qq_from_pvals(pp.1)
  # qq.recorded <- recordPlot()
  # #
  # cowplot::plot_grid(ggdraw(qq.recorded),
  #                    ggdraw(hist.recorded),
  #                    labels = c("A", "B"),
  #                    rel_widths = c(0.5, 0.5))
}

run_and_save_pvals_plot <- function(pvals.data, file){
  pdf(file, width = 9, height = 5)
  plot_pvals_routine(pvals.data )
  dev.off()

  # Need to return a charcter vector
  file
}


