

my_color <- 'd3.scaleOrdinal() .domain(["base", "not.genic", "genic", "other"])
.range(["#444444", "#800000", "#F2BA49", "#FFFF9F"])'


plot_sankey <- function(data, N.dmps){
  # Need rownames and column called "value"

  # Add a name for "Starting DMPs" (sink)
  sink.name <- paste0("N = ", format(N.dmps, big.mark=","), " DMPs")
  nodes <- data.frame(name = c(sink.name, rownames(data)))
  nodes$group <- as.factor(c("base", "not.genic", "not.genic", "not.genic",
                          "genic", "genic", "not.genic", "other"))


  # Munging for Sankey
  data$source <- 0
  data$target <- 1:nrow(data)


  # Return the object
  sankey.plot <- sankeyNetwork(Links = data,
                            Nodes = nodes,
                            Source = "source", Target = "target",
                            Value = "value", NodeID = "name",
                            fontSize = 12, nodeWidth = 20,
                            colourScale=my_color, NodeGroup="group")

  sankey.plot
}

screenshot_sankey <- function(sankey.plot, filename){

  tmp.file <- ("tmp.sankey.html")
  saveNetwork(sankey.plot, file = tmp.file, selfcontained = F)

  # zoom increases the resolution it seems
  webshot::webshot(tmp.file, filename, vwidth = 500, vheight = 250, zoom=5)
  file.remove(tmp.file)

  return(filename)
}




# Volcano plotting --------------------------------------------------------


thin_large_pvals <- function(data){
  # Drop most of the non-significant points
  # For thinning non-sig points
  set.seed(919)

  ix <- sample(which(data$y < -log10(0.2)))
  ix.to.ignore <- ix[1:floor(length(ix) * 0.98)]

  # Subset
  data[-ix.to.ignore, ]
}

add_hyper_hypo_designation <- function(data, lfdr.cut=0.05){
  # Handle colors
  data$color <- "nonsig"
  data$color[data$lfdr < lfdr.cut & data$pi.diff > 0] <- "hyper"
  data$color[data$lfdr < lfdr.cut & data$pi.diff < 0] <- "hypo"

  # Order of color and order of pallete
  data$color <- factor(data$color, levels = c("nonsig", "hyper", "hypo"))
  data
}

plot_volcano <- function(data){
  my.colors <- get_hyper_hypo_colors()
  my.pal <- c("grey", my.colors$hyper, my.colors$hypo)

  data %>%
    drop_na() %>%
    ggplot(aes(x = pi.diff, y = y, color = color)) +
    geom_point(size = 0.3) +
    xlab("LOAD effect size\n(difference in predicted methylation)") +
    ylab(expression(-log[10](lFDR))) +
    scale_x_continuous(breaks = round(seq(-0.15, 0.15, 0.05), 2), lim = c(-0.15, 0.15)) +
    scale_color_manual(values = my.pal) +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05), alpha = 0.8) +
    theme(panel.grid.major.y = element_line(color = "grey", linewidth = 0.25),
          panel.grid.major.x = element_line(color = "grey", linewidth = 0.25),
          panel.grid.minor = element_blank(),
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



# Pi chart ----------------------------------------------------------------

dmp_pi_chart_routine <- function(data){
  my.colors <- get_hyper_hypo_colors()

  # How many DMPs
  N <- nrow(data)

  # Munge
  tally.df <- data %>%
    dplyr::mutate(direction = ifelse(pi.diff < 0, "Hypo", "Hyper")) %>%
    group_by(direction) %>%
    summarize(prop = round(100 * n() / N)) %>%
    arrange(-prop) %>%
    dplyr::mutate(ypos = cumsum(prop) - 0.5*prop) %>%
    dplyr::mutate(label = paste0(direction, "\n(", prop, "%)"))

  p.pi <- ggplot(tally.df, aes(x="", y=prop, fill=label)) +
    geom_bar(stat="identity", width=1.5, color="white") +
    coord_polar("y", start = 0, direction = -1) +
    theme_void() +
    theme(legend.position="none") +
    geom_text(aes(y = ypos, label = label), color = "white", size = 9) +
    scale_fill_manual(values = c(my.colors$hyper, my.colors$hypo))

  p.pi
}


# Plot CpG island, shore, shelf designation -------------------------------

plot_cpg_barchart <- function(data){
  data %>%
    filter(value >0.3) %>%
    rownames_to_column("Location") %>%
    dplyr::mutate(Percent = round(value)) %>%
    ggplot(aes(fill = Location, x = 0, y = Percent)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = NULL, values = pal_locuszoom()(4)) +
    coord_flip() +
    theme_void() +
    theme(plot.background = element_rect(fill = "transparent", colour = "transparent"),
          panel.background = element_rect(fill = "transparent", colour = "transparent"),
          plot.title = element_text(hjust = 0.5, size=20),
          legend.text=element_text(size=16),
          legend.position = "bottom",
          plot.margin = margin(t=5, r=5,b=5,l=5)) +
    guides(fill=guide_legend(ncol=2))

}




