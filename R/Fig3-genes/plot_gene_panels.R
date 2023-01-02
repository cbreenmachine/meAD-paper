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

source("geneSchematics.functions.R")
source("../get_expanded_ensdb.R")

parser <- ArgumentParser()
parser$add_argument("--ifile", default= "../../dataDerived/analysis-controlLOAD/test-diagnostic-group-coded/experimentSummary/pvals.bed", help='Where are the models stored')
parser$add_argument("--odir", default= "../../figs/2022-paper/geneModels/natGenPaper", help='Where are the models stored')
parser$add_argument("--dmrs_file", default= "../../dataDerived/analysis-controlLOAD/test-diagnostic-group-coded/experimentSummary/dmrs.all.bed", help='Where are the models stored')
parser$add_argument("--genes_file", default= "../../dataReference/2022-natGenPaper.txt", help='Where are the genes')
args <- parser$parse_args()

odir <- args$odir
dir.create(odir, showWarn=F, recursive=T)

genes.list <- read.table(args$genes_file)[[1]]
df <- fread(args$ifile)

dmrs.df <- fread(args$dmrs_file) %>% 
  dplyr::filter(dist_to_nearest_gene == 0)


####################################
######### Data munging #############
####################################

df$y <- -log10(df$lfdr)

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



p <- make_figure("HLA-DRB1", data.gr, dmrs.df, ensdb.genes, ensdb)

cowplot::save_plot("tmp.png", p)



wrapper <- function(g){
  ofile <- file.path(odir, paste0("geneModel-", g, ".png"))
  p <- make_figure(g)
  cowplot::save_plot(ofile, p, base_width = 8, base_height = 5)
}

# keepix <- (genes.list %in% ensdb.genes$gene_name)
# genes.list.sub <- genes.list[keepix]
# print(genes.list[!keepix])

# out <- lapply(X=genes.list.sub, FUN=wrapper)

# common.legend <- get_common_legend()


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

PANEL.WIDTH <- 13
PANEL.HEIGHT <- 7

gene.lists <- list(
  "ng1" = c("ANK3", "FERMT2", "HLA-DQA1", "INPP5D"),
  "ng2" = c("PLCG2", "UMAD1", "JAZF1", "SLC24A4"),
  "me1" = c("CEP112", "HLA-DRB1", "UMAD1", "RBFOX1"),
  "me2" = c("NRG1", "PTPRD", "NOCT", "HNRNPM")
)


for (ii in 1:length(gene.lists)){
  n <- names(gene.lists)[ii]
  ll <- gene.lists[[ii]]

  panel.out <- make_panel(ll,  data.gr, dmrs.df, ensdb.genes, ensdb)
  cowplot::save_plot(file.path(odir, paste0("2022-11-21-panel-", n, ".png")), 
      panel.out, 
      base_height = PANEL.HEIGHT, 
      base_width = PANEL.WIDTH)
}







panel.out <- make_panel(,  data.gr, dmrs.df, ensdb.genes, ensdb)
cowplot::save_plot("2022-11-21-panel-2.png", panel.out, 
  base_height = PANEL.HEIGHT, base_width = PANEL.WIDTH)
