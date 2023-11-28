
plot_enhancer_enrichment_for_dmps <- function(enrichment.result, file) {
  # Histogram
  xx <- enrichment.result$simulated.df$N.dmps
  test.value <- enrichment.result$test.stat$N.dmps
  pval <- compute_empirical_p(xx, test.value)

  # Clean up
  if (pval == 0){
    pval.str <- "< 0.001"
  } else {
    pval.str <- round(pval, digits = 2)
  }

  z <- enrichment.result$simulated.df %>%
    ggplot(aes(x = N.dmps, y = after_stat(density))) +
    geom_histogram(bins = 25) +
    theme_minimal() +
    geom_vline(xintercept = enrichment.result$test.stat$N.dmps, color = "red") +
    xlab("Number of unique DMPs residing in feature") +
    ylab("Density") +
    labs(caption = paste0("Two-sided p-value: ", pval.str)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

  cowplot::save_plot(filename = file, z)
}
