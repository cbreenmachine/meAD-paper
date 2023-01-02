
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
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

source("config.R")

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
genic.df <- annotate_regions(dmps.gr)
cpg.df <- annotate_regions(dmps.gr,
                               anno.list=c("hg38_cpg_inter", "hg38_cpg_islands", "hg38_cpg_shelves", "hg38_cpg_shores"))

################################
######### PLOT SANKEY ##########
################################
node.names <- data.frame(name = c(
  paste0("N = ", format(nrow(dmps), big.mark=","), " DMPs"),
  rownames(genic.df)))

genic.df$source <- 0
genic.df$target <- 1:nrow(genic.df)

p.sankey <- sankeyNetwork(Links = genic.df,
                          Nodes = node.names,
                          Source = "source", Target = "target",
                          Value = "value", NodeID = "name",
                          fontSize = 12, nodeWidth = 20)

f1 <- file.path(ODIR, "sankey.html")
f2 <- file.path(ODIR, "sankey.png")

saveNetwork(p.sankey, file = f1, selfcontained = F)

# zoom increases the resolution it seems
webshot::webshot(f1, f2, vwidth = 500, vheight = 250, zoom=5)


# CpG island/shelf.shore --------------------------------------------------

p.bar <- cpg.df %>%
  rownames_to_column("Location") %>%
  dplyr::mutate(Percent = round(value)) %>%
  dplyr::mutate(Percent = ifelse(Percent > 90, Percent-1, Percent)) %>%
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

cowplot::save_plot(plot = p.bar, filename=file.path(ODIR, "cpg_barchart.png"),
                   base_width = 6, base_height = 3)



# Pi chart ----------------------------------------------------------------

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
  geom_text(aes(y = ypos, label = label), color = "white", size=7) +
  scale_fill_manual(values = c(C.HYPER, C.HYPO))

cowplot::save_plot(plot = p.pi, filename=file.path(ODIR, "pi.png"),
                   base_width = 3, base_height = 3)

