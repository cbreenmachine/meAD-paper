suppressPackageStartupMessages({
  library(viridis)
  library(tidyverse)
  library(data.table)
  library(GenomicRanges)
  library(cowplot)
  library(ggsci)
  library(ggbio)
  library(clusterProfiler)
  library(annotatr)
  library(networkD3)
  library(rbokeh)
  library(waffle)
})

# Constants
LFDR.CUT <- 0.05
C.HYPER <- "#0073C2FF"
C.HYPO <- "#EFC000FF"

MAPPING <-
  data.frame(
    Annotation = c("1to5kb", "3UTRs", "5UTRs", "exons", "introns", "promoters",
                   "inter", "islands", "shelves", "shores"),
    clean = c("1-5 kb", "3' UTR", "5' UTR", "Exon", "Intron", "Promoter",
              "Open sea", "CpG island", "CpG shelf", "CpG shore")
  )

# Directories
pvals.file <- "../../dataDerived/analysis-controlLOAD/test-diagnostic-group-coded/experimentSummary/pvals.bed"
odir <- "../../figs/2022-paper/"

dir.create(odir, showWarn=F, recurs=T)

# Something like dmps.all
prefix <- "fig1-"

# Read data and convert to GRanges
data <- fread(pvals.file)
dmps <- data %>% dplyr::filter(lfdr < LFDR.CUT)
dmps.gr <- dmps %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T, starts.in.df.are.0based = T)

seqlevelsStyle(dmps.gr) <- "UCSC"


annotate_regions <- function(regions.gr, anno.list='hg38_basicgenes'){

  # hg38_basicgenes or hg38_cpg
  my.anno <- build_annotations(genome = 'hg38', annotations = anno.list)

  # Run the annotation function
  regions.annotated <- annotatr::annotate_regions(regions.gr, my.anno, min=2)

  anno.df <- data.frame(regions.annotated) %>%
    dplyr::transmute(cpg.locus = paste0(seqnames, "-", start), annot.type) %>% # unique id for each CpG
    distinct() %>%  # get rid of instances where one CpG resides in multiple introns
    dplyr::mutate(Annotation = str_remove(str_remove(annot.type, "hg38_genes_"), "hg38_cpg_"))

  anno.counts.df <- anno.df %>%
    group_by(Annotation) %>%
    dplyr::summarize(value = round(100 * n() / length(regions.gr), 2)) %>%
    merge(MAPPING, by = "Annotation") %>%
    dplyr::mutate(label = paste0(clean,  " (", value, "%)")) %>%
    column_to_rownames("label") %>%
    dplyr::select(c("value"))

  anno.counts.df
}

# Annotation counts
ac.df <- annotate_regions(dmps.gr)
cpg.loc.df <- annotate_regions(dmps.gr, anno.list=c("hg38_cpg_inter", "hg38_cpg_islands", "hg38_cpg_shelves", "hg38_cpg_shores"))

################################
######### PLOT SANKEY ##########
################################
node.names <- data.frame(name = c(
  paste0("N = ", format(nrow(dmps), big.mark=","), " DMPs"),
  rownames(ac.df)))

ac.df$source <- 0
ac.df$target <- 1:nrow(ac.df)

p.sankey <- sankeyNetwork(Links = ac.df,
                          Nodes = node.names,
                          Source = "source", Target = "target",
                          Value = "value", NodeID = "name",
                          fontSize = 12, nodeWidth = 20)

widget2png(p.sankey, file.path(odir, paste0(prefix, "sankey.png")))


################################
######### PLOT WAFFLE ##########
################################
p.waffle <- cpg.loc.df %>%
  rownames_to_column("Location") %>%
  dplyr::mutate(Percent = round(value)) %>%
  dplyr::mutate(Percent = ifelse(Percent > 90, Percent-1, Percent)) %>%
  ggplot(aes(fill = Location, values = Percent)) +
  geom_waffle(nrow = 5, ncol = 20, size = 0.25, color="white") +
  scale_fill_manual(name = NULL, values = pal_locuszoom()(4)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", colour = "transparent"),
        panel.background = element_rect(fill = "transparent", colour = "transparent"),
        plot.title = element_text(hjust = 0.5, size=20),
        legend.text=element_text(size=16),
        legend.position = "bottom",
        plot.margin = margin(t=5, r=5,b=5,l=5)) +
  guides(fill=guide_legend(ncol=2))

ofile <- file.path(odir, paste0(prefix, "waffleDMPsAnno.png"))
cowplot::save_plot(plot = p.waffle, filename=ofile,
                   base_width = 6, base_height = 3)

######################################
######### PLOT VOLCANO ###############
######################################

# Thin non-significant points
data$y <- -log10(data$lfdr)
ix <- sample(which(data$y < -log10(0.01)))
ix.to.ignore <- ix[1:floor(length(ix) * 0.999)]

# Subset
sub.df <- data[-ix.to.ignore, ]

# Handle colors
sub.df$color <- "grey"
sub.df$color[sub.df$lfdr < LFDR.CUT & sub.df$diagnostic_group_coded > 0] <- "hyper"
sub.df$color[sub.df$lfdr < LFDR.CUT & sub.df$diagnostic_group_coded < 0] <- "hypo"

# Order of color and order of pallete
sub.df$color <- factor(sub.df$color, levels = c("grey", "hyper", "hypo"))
my.pal <- c("grey", C.HYPER, C.HYPO)

#Plot volcano
p.volcano <- sub.df %>%
  ggplot(aes(x = pi.diff, y = y, color = color)) +
  geom_jitter(size = 0.2, alpha = 1, height=1, width=0) +
  scale_color_manual(values = my.pal) +
  xlab("LOAD effect size\n(difference in adjusted methylation)") +
  ylab(expression(-log[10](lFDR))) +
  scale_x_continuous(breaks = seq(-1, 1, 0.25), lim = c(-1, 1)) +
  scale_y_continuous(breaks = c(0, 100, 200, 300), lim = c(0, 350)) +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(color = "grey", size=0.25),
        panel.grid.major.x = element_line(color = "grey", size=0.25),
        axis.text.x = element_text(angle = 0),
        axis.ticks = element_blank(),
        legend.position = "none")

ofile <- file.path(odir, paste0(prefix, "volcano.png"))
cowplot::save_plot(plot = p.volcano, filename=ofile,
                   base_width = 6, base_height = 6)



###################################
###################################


tally.df <- dmps %>%
  dplyr::mutate(direction = ifelse(diagnostic_group_coded < 0, "Hypo", "Hyper")) %>%
  group_by(direction) %>%
  summarize(prop = round(100 * n() / nrow(.))) %>%
  arrange(-prop) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5*prop) %>%
  dplyr::mutate(label = paste0(direction, "\n(", prop, "%)"))


p.pi <- ggplot(tally.df, aes(x="", y=prop, fill=label)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start= 0, direction = -1) +
  theme_void() +
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = label), color = "white", size=9) +
  scale_fill_manual(values = c(C.HYPER, C.HYPO))

ofile <- file.path(odir, paste0(prefix, "pi.png"))
cowplot::save_plot(plot = p.pi, filename=ofile,
                   base_width = 3, base_height = 3)
