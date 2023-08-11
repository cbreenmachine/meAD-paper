

plot_tech_comp_scatter <-  function(merged){
  # Small helper to plot a scatter plot
  merged %>%
    sample_n(10000) %>%
    ggplot(aes(x = EPIC, y = WGMS)) +
    geom_point(alpha = 0.1) +
    theme_minimal() %>%
    return()
}

