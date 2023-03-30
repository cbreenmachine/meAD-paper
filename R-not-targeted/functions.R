# General purpose functions
# Should be analysis agnostic


# General reading functions -----------------------------------------------

get_pvals_data <- function(ifile){
  # READ pvals from
  data.table::fread(ifile, verbose = F) %>%
    dplyr::mutate(lfdr = lfdr.from.ss,
                  y = -log10(lfdr)) %>%
    dplyr::select(-one_of(lfdr.from.zz))
}



# Essentially constants ---------------------------------------------------


get_hyper_hypo_colors <- function(){
  # Hyper and hypo colors, may be used outside volcano
  # so we'll keep it in the R/functions.R file
  return(list(hyper = "#0073C2FF", hypo = "#EFC000FF"))
}




# Volcano plotting --------------------------------------------------------


thin_large_pvals <- function(data){
  # Drop most of the non-significant points
  # For thinning non-sig points
  set.seed(919)

  ix <- sample(which(data$y < -log10(0.2)))
  ix.to.ignore <- ix[1:floor(length(ix) * 0.9)]

  # Subset
  data[-ix.to.ignore, ]
}

add_hyper_hypo_designation <- function(data, lfdr.cut=0.05){
  # Handle colors
  data$color <- "grey"
  data$color[sub.df$lfdr < lfdr.cut & sub.df$pi.diff > 0] <- "hyper"
  data$color[sub.df$lfdr < lfdr.cut & sub.df$pi.diff < 0] <- "hypo"

  # Order of color and order of pallete
  data$color <- factor(data$color, levels = c("grey", "hyper", "hypo"))
}

plot_volcano <- function(data){
  my.colors <- get_hyper_hypo_colors()
  my.pal <- c("grey", my.colors$hyper, my.colors$hypo)

  data %>%
      ggplot(aes(x = pi.diff, y = y, color = color)) +
      geom_point(size = 0.3) +
      xlab("LOAD effect size\n(difference in predicted methylation)") +
      ylab(expression(-log[10](lFDR))) +
      scale_x_continuous(breaks = round(seq(-0.15, 0.15, 0.05),2), lim = c(-0.15, 0.15)) +
      scale_color_manual(values = my.pal) +
      theme_minimal() +
      geom_hline(yintercept = -log10(0.05), alpha = 0.8) +
      theme(panel.grid.major.y = element_line(color = "grey", linewidth = 0.25),
            panel.grid.major.x = element_line(color = "grey", linewidth = 0.25),
            axis.text.x = element_text(angle = 0),
            axis.ticks = element_blank(),
            legend.position = "none",
            panel.background = element_rect(fill = "white", color = "white"),
            plot.background = element_rect(fill = "white", color = "white"))
}


volcano_routine <- function(data){
  data %>%
    thin_large_pvals() %>%
    add_hyper_hypo_designation() %>%
    plot_volcano()
}
