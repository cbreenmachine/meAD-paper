



# Constants ---------------------------------------------------------------


GENIC.MAPPING <-
  data.frame(
    Annotation = c("1to5kb", "3UTRs", "5UTRs", "exons", "introns", "promoters", "Other"),
    Annotation.clean = c("1-5 kb", "3' UTR", "5' UTR", "Exon", "Intron", "Promoter", "Other")
  )


CPG.MAPPING <-
  data.frame(
    Annotation = c("hg38_cpg_inter", "hg38_cpg_islands", "hg38_cpg_shelves", "hg38_cpg_shores", "Other"),
    Annotation.clean = c("Open sea", "CpG island", "CpG shelf", "CpG shore", "Other")
  )

# Put factors in corret order, refelcting biology
GENIC.LEVELS <- c("1-5 kb", "Promoter", "5' UTR", "Exon", "Intron",  "3' UTR", "Other")
CPG.LEVELS <- c("Open sea", "CpG island", "CpG shelf", "CpG shore", "Other")



# Helper functions --------------------------------------------------------


annotate_dmps_to_intron_exon <- function(data.gr){
  my.anno <- annotatr::build_annotations(genome = 'hg38', annotations = 'hg38_basicgenes')
}


annotate_loci_to_genic_parts <- function(data.gr, annotate.to = 'hg38_basicgenes'){

  if ("hg38_basicgenes" %in% annotate.to){
    MAPPING <- GENIC.MAPPING
    LEVELS <- GENIC.LEVELS
  } else {
    MAPPING <- CPG.MAPPING
    LEVELS <- CPG.LEVELS
  }

  print(MAPPING)

  # hg38_basicgenes or hg38_cpg
  my.anno <- annotatr::build_annotations(genome = 'hg38', annotations = annotate.to)

  # Run the annotation function
  data.annotated <- annotatr::annotate_regions(data.gr, my.anno, quiet = T)

  # Post-process after linking a CpG to a genomic feature
  anno.df <- data.frame(data.annotated) %>%
    dplyr::transmute(cpg.locus = paste0(seqnames, "-", start), annot.type) %>% # unique id for each CpG
    distinct() %>%  # get rid of instances where one CpG resides in multiple introns
    dplyr::mutate(Annotation = str_remove(annot.type, "hg38_genes_"))

  # Count how many annotate to none of these
  none.of.the.above <- subsetByOverlaps(data.gr,
                                        data.annotated,
                                        invert = T)

  # How many went to none?
  N.other <- length(none.of.the.above)

  # Add an extra row to anno.df; annotatr does not (as best as I know)
  # natively support counting "nones"
  other.df <- data.frame(cpg.locus = rep("none", N.other),
                         annot.type = rep("none", N.other),
                         Annotation = rep("Other", N.other))


  # Get the counts
  anno.counts.df <- rbind(anno.df, other.df) %>%
    group_by(Annotation) %>%
    dplyr::summarize(value = round(100 * n() / length(data.gr), 1))

  # Clean some strings and put in correct order
  anno.counts.cleaned.df <- anno.counts.df %>%
    left_join(MAPPING, by = "Annotation") %>%
    dplyr::mutate(Annotation.clean = factor(Annotation.clean, levels = LEVELS)) %>%
    arrange(Annotation.clean)

  # Add a string that looks like 'Promoter (10 %)'
  anno.counts.cleaned.df %>%
    dplyr::mutate(label = paste0(Annotation.clean,  " (", value, "%)")) %>%
    column_to_rownames("label")
}



cpg_annotation_routine <- function(data.gr){
  annotate_loci_to_genic_parts(
    data.gr,
    annotate.to = c("hg38_cpg_inter", "hg38_cpg_islands", "hg38_cpg_shelves", "hg38_cpg_shores")
  )
}




# Plotting ----------------------------------------------------------------


my_color <- 'd3.scaleOrdinal() .domain(["base", "not.genic", "genic", "other"])
.range(["#444444", "#800000", "#F2BA49", "#FFFF9F"])'


plot_sankey <- function(data, N.dmps){
  # Need rownames and column called "value"

  # Add a name for "Starting DMPs" (sink)
  ti <- htmltools::HTML(paste0("<strong><center>N = ", format(N.dmps, big.mark=","), " DMPs</center></strong>"))
  sink.name <- ""
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

  sankey.plot <- htmlwidgets::prependContent(sankey.plot, ti)


  sankey.plot
}

screenshot_sankey <- function(sankey.plot, filename){

  tmp.file <- ("tmp.sankey.html")
  saveNetwork(sankey.plot, file = tmp.file, selfcontained = F)

  # zoom increases the resolution it seems
  webshot::webshot(tmp.file, filename, vwidth = 500, vheight = 300, zoom=5)
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
    xlab("AD methylation % minus no-AD methylation %") +
    ylab(expression(-log[10](lFDR))) +
    scale_x_continuous(breaks = round(seq(-0.15, 0.15, 0.05), 2), lim = c(-0.15, 0.15)) +
    scale_color_manual(values = my.pal) +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05), alpha = 0.8) +
    theme(panel.grid.major.y = element_line(color = "grey", linewidth = 0.1),
          panel.grid.major.x = element_line(color = "grey", linewidth = 0.1),
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
    dplyr::mutate(direction = ifelse(pi.diff < 0, "Hypomethylation", "Hypermethylation")) %>%
    group_by(direction) %>%
    summarize(prop = round(100 * n() / N)) %>%
    arrange(-prop) %>%
    dplyr::mutate(ypos = cumsum(prop) - 0.5*prop) %>%
    dplyr::mutate(label = paste0(direction, " (", prop, "%)"))

  p.pi <- ggplot(tally.df, aes(x="", y=prop, fill=label)) +
    geom_bar(stat="identity", width=1.5, color="white") +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = "white"),
          plot.title = element_text(hjust = 0.5, size=20),
          legend.text=element_text(size=16),
          legend.position = "bottom",
          plot.margin = margin(t=5, r=5,b=5,l=5)) +
    scale_fill_manual(name=NULL, values = c(my.colors$hyper, my.colors$hypo)) +
    guides(fill=guide_legend(ncol=1))

  p.pi
}


# Plot CpG island, shore, shelf designation -------------------------------

plot_cpg_pi_chart <- function(data){
  data %>%
    dplyr::filter(value > 0.1) %>%
    rownames_to_column("Location") %>%
    ggplot(aes(fill = Location, x = "", y = value)) +
    geom_bar(width=1, stat="identity", color="white") +
    scale_fill_manual(name = NULL, values = pal_locuszoom()(4)) +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = "white"),
          plot.title = element_text(hjust = 0.5, size=20),
          legend.text=element_text(size=16),
          legend.position = "bottom",
          plot.margin = margin(t=5, r=5,b=5,l=5)) +
    guides(fill=guide_legend(ncol=2, reverse = F))

}



# Stitch ------------------------------------------------------------------

stitch_birdseye_fig <- function(volcano.file, sankey.file, pi.file, cpg.file, out.file){

  volcano <- ggdraw() + draw_image(volcano.file)
  sankey <- ggdraw() + draw_image(sankey.file)
  pi <- ggdraw() + draw_image(pi.file)
  cpg <- ggdraw() + draw_image(cpg.file)

  topright <- plot_grid(pi, NULL, cpg, rel_widths = c(0.5, 0, 0.5),
                        labels = c("B","", "C"), nrow = 1)

  right <- plot_grid(topright, NULL, sankey,
                     ncol = 1,
                     rel_heights = c(0.5, 0, 0.5),
                     labels = c(NA, NA, "D"))

  z <- plot_grid(volcano, right,
            ncol = 2,
            rel_widths = c(0.5, 0.5),
            labels = c("A", NA)) +
    theme(plot.background = element_rect(fill = "white", color = "white")) +
    panel_border(remove = TRUE)

  cowplot::save_plot(z, filename = out.file)
  out.file

}


