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
})

source("util-gene-schematics.R")
source("util-ensdb.R")

DIR <- "../../DataRaw/2023-02-14-Summaries-v6//"

IFILE <- file.path(DIR, "pvals.bed")
DMRS.PATH <- file.path(DIR, "DMRegions.bed")
ODIR <- "../../Figs/Fig3/"

df <- fread(IFILE)

dmrs.df <- fread(DMRS.PATH) %>%
  dplyr::filter(dist_to_nearest_gene == 0)


####################################
######### Data munging #############
####################################

df$y <- -log10(df$lfdr.from.ss)

sig_levels <- c("Not significant", "Significant (hyper)", "Significant (hypo)")
my_pal <- c("grey", C.HYPER, C.HYPO)

df$color <- sig_levels[1]
df$color[df$y > -log10(ALPHA) & df$diagnostic_group_coded > 0] <- sig_levels[2]
df$color[df$y > -log10(ALPHA) & df$diagnostic_group_coded < 0] <- sig_levels[3]

df$color <- factor(df$color, levels = sig_levels)

# Coerce to GRanges
data.gr <- makeGRangesFromDataFrame(df, keep.extra.columns=T, starts.in.df=T)
ensdb <- get_raw_ensdb()

ensdb.genes <- get_expanded_ensdb(UP, DOWN) %>% trim()
ensdb.genes <- ensdb.genes[seqnames(ensdb.genes) %in% paste0("chr", 1:22)]

seqlevelsStyle(data.gr) <- "UCSC"


wrapper <- function(g){
  ofile <- file.path(odir, paste0("schematic-", g, ".png"))
  p <- make_figure(g)
  cowplot::save_plot(ofile, p, base_width = 8, base_height = 5)
}


make_panel <- function(gene.list, data.gr, dmrs.df, ensdb.genes, ensdb){
  p1 <- make_figure(gene.list[1], data.gr, dmrs.df, ensdb.genes, ensdb, xlab = F, ylab=T, legend = F)
  p2 <- make_figure(gene.list[2], data.gr, dmrs.df, ensdb.genes, ensdb, xlab = F, ylab=F, legend = F)
  p3 <- make_figure(gene.list[3], data.gr, dmrs.df, ensdb.genes, ensdb, xlab = T, ylab=T, legend = F)
  p4 <- make_figure(gene.list[4], data.gr, dmrs.df, ensdb.genes, ensdb, xlab = T, ylab=F, legend = F)

  common.legend <- get_common_legend()

  panel <- cowplot::plot_grid(p1, p2, p3, p4, nrow = 2)
  cowplot::plot_grid(panel, common.legend, ncol = 1, rel_heights = c(2, .1)) +
    theme(
      legend.background = element_rect(fill = "white"),
      plot.background = element_rect(fill= "white", color="white"),
      panel.background = element_rect(fill= "white", color = "white")
    )
}

PANEL.WIDTH <- 13.5
PANEL.HEIGHT <- 7

gene.lists <- list(
  "paper" = c("HLA-DQA1", "HLA-DRB1", "MAPT", "UMAD1"),
  "ryrs" = c("RYR1", "RYR2", "RYR3", "MAPT"),
  "canonical" = c("APOE", "APP", "PSEN1", "PSEN2")
)


for (ii in 1:length(gene.lists)){
  n <- names(gene.lists)[ii]
  ll <- gene.lists[[ii]]
  ofile <- file.path(ODIR, paste0(Sys.Date(), "-panel-", n, ".png"))

  panel.out <- make_panel(ll,  data.gr, dmrs.df, ensdb.genes, ensdb)
  cowplot::save_plot(
      ofile,
      panel.out,
      base_height = PANEL.HEIGHT,
      base_width = PANEL.WIDTH)
}
#END
