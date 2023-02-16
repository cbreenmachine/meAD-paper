# tabulate-annotate-to-genes.R
# Takes a list of DMPs, runs gene ontologies
# for
# 1. Introns
# 2. Exons
# 3. Gene bodies
# 4. Genes +3kb, -200
# Data are assumed to be in hg38

suppressPackageStartupMessages({
  library(viridis)
  library(tidyverse)
  library(data.table)
  library(GenomicRanges)
  library(cowplot)
  library(ggsci)
  library(ggbio)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(Homo.sapiens)
  library(OrganismDbi)
  library(goplot)
})

source("config.R")


# Read data ---------------------------------------------------------------

# Read data and convert to GRanges
data <- fread(pvals.file, verbose = F)
dmps <- data %>%
  dplyr::mutate(lfdr = lfdr.from.ss) %>%
  dplyr::filter(lfdr < LFDR.CUT)
dmps.gr <- dmps %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T, starts.in.df.are.0based = T)

seqlevelsStyle(dmps.gr) <- "UCSC"


# Intersect ---------------------------------------------------------------

# Pull genes
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene,
               single.strand.genes.only = FALSE,
               columns = "gene_id")

# Introns and exons too
introns <- unlist(intronsByTranscript(TxDb.Hsapiens.UCSC.hg38.knownGene))
exons <- unlist(exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene))


get_first_feature <- function(tp = "exon"){
  # Returns a list of GRanges, where each "row"
  # is an intron / exon
  if (tp == "intron"){
    warning("not supported yet")
  } else if (tp == "exon"){
    by_gene <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")
  }

  # Grab the associated strand so we know whether to
  # grab the leftmost or rightmost value
  gene_strand <- (unlist(unique(strand(by_gene))))

  # Need to ignore a <10 genes with * strand (messes up downstream stuff)
  by_gene_sub <- by_gene[names(gene_strand)]

  # Grab the index
  ix <- ifelse(gene_strand == "+", 1, which.max(end(by_gene_sub))) %>% as.list()
  names(ix) <- names(by_gene_sub)

  first_in_each_gene <- unlist(by_gene_sub[ix], use.names = T)
  first_in_each_gene
}

first.exons <- get_first_feature("exon")


# Get the DMPs in these genic components
first.exon.dmps <- subsetByOverlaps(dmps.gr, first.exons)
intron.dmps <- subsetByOverlaps(dmps.gr, introns)
exon.dmps <- subsetByOverlaps(dmps.gr, exons)

# Histogram of pi.diffs ---------------------------------------------------

pis.df <-
  rbind(
    as.data.frame(first.exon.dmps[ ,'pi.diff']) %>% dplyr::mutate(type = "first exon"),
    as.data.frame(exon.dmps[ ,'pi.diff']) %>% dplyr::mutate(type = "all exons"),
    as.data.frame(intron.dmps[ ,'pi.diff']) %>% dplyr::mutate(type = "all introns")
    )

p <- pis.df %>%
  ggplot(aes(x = pi.diff, y = after_stat(density), fill = type)) +
  # geom_histogram(bins = 75, alpha = 0.5, position = "identity") +
  geom_density(alpha = 0.5, position = "identity", color = NA) +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0)) +
  xlab("Difference in predicted methylation\n(LOAD - control)") +
  ylab("Density") +
  labs(caption = "DMPs are filtered to include those in any exon,\nany intron, or the first exon only. There are far fewer\nDMPs in the first exon.")

cowplot::save_plot(file.path(ODIR, "intron-exon-comp-pi-diffs.png"), p)

# Workhorse function ------------------------------------------------------

get_gene_list <- function(my.dmps, out.type){
  # my.dmps : a subset of  DMPs residing in introns, exons, etc.
  # out.type : 'SYMBOL', 'ENSEMBL', etc.

  genes.w.dmp <- subsetByOverlaps(genes, my.dmps)

  # Gene IDs are consistent across Genome builds, so
  # no need to worry
  mapped.genes <- OrganismDbi::select(
    Homo.sapiens,
    names(genes.w.dmp),
    out.type,
    "GENEID")
  sort(unique(mapped.genes[[out.type]]))
}


gene.body.ids <- get_gene_list(dmps.gr, "ENSEMBL")
intron.ids <- get_gene_list(intronic.dmps, "ENSEMBL")
first.exon.ids <- get_gene_list(first.exon.dmps, "ENSEMBL")
exon.ids <- get_gene_list(exonic.dmps, "ENSEMBL")



# Gene ontologies ---------------------------------------------------------



# WHOLE GENE
first.exon.go.out <- goplot::run_gene_ontology(first.exon.ids)
p <- goplot::go_output_to_df(first.exon.go.out) %>%
  goplot::plot_go_barchart() +
  ggtitle("First exon-associated genes") +
  geom_hline(yintercept = -log10(0.05))

ofile <- file.path(ODIR, "GO-first-exons.png")
cowplot::save_plot(ofile, p, base_width = 14, base_height = 8)


# WHOLE GENE
gene.go.out <- goplot::run_gene_ontology(gene.body.ids)
p <- goplot::go_output_to_df(gene.go.out) %>%
  goplot::plot_go_barchart() +
  ggtitle("DMP-associated genes") +
  geom_hline(yintercept = -log10(0.05))

ofile <- file.path(ODIR, "GO-genes.png")
cowplot::save_plot(ofile, p, base_width = 14, base_height = 8)

#INTRONS
intron.go.out <- goplot::run_gene_ontology(intronic.ids)
p <- goplot::go_output_to_df(intron.go.out) %>%
  goplot::plot_go_barchart() +
  ggtitle("DMP-associated introns")

ofile <- file.path(ODIR, "GO-introns.png")
cowplot::save_plot(ofile, p, base_width = 14, base_height = 8)


#EXONS
exon.go.out <- goplot::run_gene_ontology(exonic.ids)
p <- goplot::go_output_to_df(exon.go.out) %>%
  goplot::plot_go_barchart() +
  ggtitle("DMP-associated exons")

ofile <- file.path(ODIR, "GO-exons.png")
cowplot::save_plot(ofile, p, base_width = 14, base_height = 8)

