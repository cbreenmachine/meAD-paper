

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



