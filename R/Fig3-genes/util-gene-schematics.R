suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    options(repos = BiocManager::repositories())
    library(EnsDb.Hsapiens.v86)
    library(ggbio)
    library(ggsci)
    library(GenomicRanges)
    library(scales)
    library(cowplot)
    library(argparse)
})



####################################
############ CONSTANTS #############
####################################
SZ <- 1.1
UP <- 5000 # upstream of TSS
DOWN <- 300 # downstream of Gene end

POINT.SZ <- 1.4
POINT.ALPHA <- 0.5

ALPHA <- 0.05 # significante level

C.GENE <- "#A73030FF" #maroonish
C.HYPER <- "#0073C2FF" #blue
C.HYPO <- "#EFC000FF" # Yellow

YLIM.EFFECT <- c(-1, 1)
YLIM.LFDR <- c(0, 150)

# Set plot style...
th <- theme(
  legend.position = "none",
  panel.border = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  plot.caption = element_text(hjust = 0),
  panel.background = element_rect(colour = NA),
  plot.background = element_rect(colour = NA)
)

theme_set(th)




####################################
######### Helper functions #########
####################################
add_tss_arrow <- function(p, xs, st) {
  sz <- 1.25
  ys <- c(0.725, 1.275)

  # How wide should arrow be
  delta.x <- xs[2] - xs[1]
  delta.y <- ys[2] - ys[1]

  # If top strand
  if (st == "+"){
    tmp.df <- data.frame(
      x1 = xs[1] + UP,
      x2 = xs[1] + UP + (0.1) * delta.x,
      y1 = 2,
      y2 = 2.35
    )
  } else if (st == "-"){
    tmp.df <- data.frame(
      x1 = xs[2] - UP,
      x2 = xs[2] - UP - (0.1) * delta.x,
      y1 = 2,
      y2 = 2.35
    )
  }

  p.new <- p +
    geom_segment(
      aes(x = x1, y = y1, xend = x1, yend = y2, color = NULL),
      data=tmp.df, size = sz,
      lineend = "square", linejoin = "round"
    ) +
    geom_segment(
      aes(x = x1, y = y2, xend = x2, yend = y2, color = NULL),
      data = tmp.df, size = sz,
      lineend = "butt", linejoin = "mitre",
      arrow = arrow(type = "closed",
                    length = unit(0.025, "npc"))
    )
  return(p.new)
}


add_promoter <- function(p, xs, st) {
    # adds green promoter box near TSS
    # p is an existing ggplot
    # xs is a length 2 vector telling start and stop of gene
    # st is strand ("+" or "-")
    sz <- 4
    ys <- c(1.725, 2.275)

    # How wide should arrow be
    delta.x <- xs[2] - xs[1]
    # delta.y <- ys[2] - ys[1]

    # Generate a small data frame representing
    # the x and y coordinates of the promoter
    if (st == "+"){
        tmp.df <- data.frame(
            x1 = xs[1], x2 = xs[1] + UP,
            y1 = ys[2], y2 = ys[1]
        )
    } else if (st == "-"){
        tmp.df <- data.frame(
            x1 = xs[2] - UP, x2 = xs[2],
            y1 = 1.75, y2 = 2.25
        )
    }

    p.new <- p +
        geom_rect(
            aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            fill = "springgreen4", alpha = 1, data=tmp.df
        )
    return(p.new)
}


get_gene_boundaries <- function(gene){
    # Tell me where the "gene" (including some UP nt upstream
    # and some DOWN bp downstream) starts and ends,
    # for plotting purposes

    # Top of plot, plots the gene introns/exons
    symbol.filt <- SymbolFilter(gene)

    # Derived values for marking TSS
    ix <- which(ensdb.genes$gene_name == gene)
    left <- start(ensdb.genes[ix])
    right <- end(ensdb.genes[ix])
    st <- as.character(strand(ensdb.genes[ix])[1])

    if (st == "+"){
        return(c(left - UP, right + DOWN))
    } else if (st == "-") {
        return(c(left - DOWN, right + UP))
    }
}



####################################
######### Plot functions ###########
####################################
plot_gene_model <- function(gene, ensdb.genes, ensdb){
    # Top of plot, plots the gene introns/exons
    symbol.filt <- SymbolFilter(gene)

    # Derived values for marking TSS
    ix <- which(ensdb.genes$gene_name == gene)

    # Another filter to handle wonky genes
    chr.filt <- SeqNameFilter(seqnames(ensdb.genes[ix]))
    left <- start(ensdb.genes[ix])
    right <- end(ensdb.genes[ix])
    st <- as.character(strand(ensdb.genes[ix])[1])

    ti <- bquote('Differential methylation in'~italic(.(gene)))

    tmp <- ensdb %>%
      autoplot(AnnotationFilterList(symbol.filt, chr.filt),
               stat="identity", label=F,
               color=  C.GENE, fill = C.GENE) +
      ggtitle(ti) +
      ylab("") +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_line(color = "grey", size=0.25),
        plot.caption = element_text(hjust = 0),
        plot.title = element_text(vjust = 0),
        axis.title.y=element_blank(),
        axis.ticks = element_blank()
      ) +
      ylim(c(0, 3.5))

    p <- tmp@ggplot

    p.2 <- add_promoter(p, c(left, right), st)
    p.3 <- add_tss_arrow(p.2, c(left, right), st)

  return(p.3)
}


plot_diffs <- function(gene, data.gr, ensdb.genes){
    #Index by gene name
    gene.ix <- which(ensdb.genes$gene_name == gene)

    #Find "rows" in GR that overlap the gene
    cpgs.ix <- queryHits(findOverlaps(data.gr, ensdb.genes[gene.ix, ], ignore.strand=T))
    data <- as.data.frame(data.gr[cpgs.ix, ])

    data %>%
        ggplot(aes(x = start, y = pi.diff, color = color)) +
        geom_point(size = POINT.SZ, alpha = POINT.ALPHA) +
        xlab("Genomic position") +
        ylab("Effect size") +
        scale_x_continuous(
          labels = unit_format(unit = "mb", scale=1e-6)
        ) +
        scale_y_continuous(
            limits = c(-1, 1),
            breaks = c(-1, 0, 1),
            minor_breaks = c(-0.5, 0.5),
        ) +
        scale_color_manual(values = my_pal) +
        theme(
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y.right = element_text(hjust= 0.5, angle = 90),
            axis.text.x = element_text(angle = 0),
            panel.grid.major.y = element_line(color = "grey", size=0.25),
            panel.grid.minor.y = element_line(color = "grey", size=0.2),
            panel.grid.major.x = element_line(color = "grey", size=0.25),
            plot.caption = element_text(hjust = 0),
            legend.position="none",
            axis.line.y.right = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "both"))
        )
}


#--> Bottom of plot, scatter of
plot_points <- function(gene, data.gr, ensdb.genes){

    y.lab <- bquote('-log'[10]*'(lFDR)')

    #Index by gene name
    gene.ix <- which(ensdb.genes$gene_name == gene)

    #Find "rows" in GR that overlap the gene
    cpgs.ix <- queryHits(findOverlaps(data.gr, ensdb.genes[gene.ix, ], ignore.strand=T))

    data <- as.data.frame(data.gr[cpgs.ix, ]) %>% distinct()
    N.Sig <- sum(abs(data$y) > -log10(ALPHA))

    p <- data  %>%
      ggplot(aes(x = start, y = abs(y), color = color)) +
      geom_point(size = POINT.SZ, alpha = POINT.ALPHA) +
      scale_y_continuous(breaks = c(0, 50, 100, 150, 200), lim = c(0, 250)) +
      geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
      geom_hline(yintercept = -log10(ALPHA), color = "black", alpha = 0.8) +
      xlab("") +
      ylab(y.lab) +
      scale_color_manual(values = my_pal) +
      annotate("text", -Inf, Inf, label = paste0("N=", N.Sig, " DMPs"), hjust = -0.3, vjust = 1.5) +
      theme(
        axis.text.y.right = element_text(angle = 90, hjust= 0.5),
        panel.grid.major.y = element_line(color = "grey", size=0.25),
        panel.grid.major.x = element_line(color = "grey", size=0.25),
        plot.caption = element_text(hjust = 0),
        legend.position= "top",
        legend.title = element_blank(),
        legend.background = element_rect(colour="white", fill=alpha("white", 0.5)),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()
      )

    return(p)
}


add_DMR <- function(p, left, right, pad=5000){
  # pad DMR (some are 10s of bases) to make it stand out
    p.ann <- p + annotate("rect", fill = "orange",
            xmin = left - pad, xmax = right + pad,
            ymin = -Inf, ymax = Inf,
            alpha = 0.3
        )
    return(p.ann)
}



# Helper function
stack_all_three <- function(a, b, c){
  p <- cowplot::plot_grid(a, NULL, b, NULL, c, ncol = 1, axis = "b",
                          align = "v", rel_heights = c(2, -0.5, 3, -0.1, 2)) +
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  return(p)
}


#--> Driver function
make_figure <- function(gene, data.gr, dmrs.df, ensdb.genes, ensdb, xlab=T, ylab=T, legend=T){

  bounds <- get_gene_boundaries(gene)

  # padding for DMR
  my.pad <- floor(diff(bounds) * 0.005)

  print("Building gene model")
  p1 <- plot_gene_model(gene, ensdb.genes, ensdb) +
    xlim(bounds) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

  print("Plotting lFDRs")
  p2 <- plot_points(gene, data.gr, ensdb.genes) +
    xlim(bounds)  +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

  # Conditoinal labels
  if (!ylab){p2 <- p2 + ylab("")}
  if (!legend){p2 <- p2 + theme(legend.position = "none")}

  print("Plotting differences")
  p3 <- plot_diffs(gene, data.gr, ensdb.genes) +
    xlim(bounds) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

  if (!xlab){p3 <- p3 + xlab("")}
  if (!ylab){p3 <- p3 + ylab("")}

if (gene %in% dmrs.df$nearest_gene_name){
    print("Adding DMR")
    keepix <- dmrs.df$nearest_gene_name == gene
    dmr.left <- dmrs.df$start[keepix]
    dmr.right <- dmrs.df$end[keepix]

    p2 <- add_DMR(p2, dmr.left, dmr.right, my.pad)
    p3 <- add_DMR(p3, dmr.left, dmr.right, my.pad)
}

  p <- stack_all_three(p1, p2, p3)
  return(p)
}

get_common_legend <- function(){

  # PANEL
  # Grab a legend from some gene that has all three point types
  tmp.plot <- plot_points("RBFOX1", data.gr, ensdb.genes) +
      theme(legend.position = "bottom",
            legend.background = element_rect(fill = "white", color="white"),
            plot.background = element_rect(fill= "white", color="white"),
            legend.key = element_rect(fill="white", color="white"),
            panel.background = element_rect(fill= "white", color="white"))

  common.legend <- cowplot::get_legend(tmp.plot)
  common.legend
}

