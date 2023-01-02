library(tidyverse)
library(data.table)
library(poolr)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(goplot)

source("../functions_overlap_regions_compute_p.R")

# Reference data
genes <-
  EnsDb.Hsapiens.v86 |>
  filter(~ tx_biotype == "protein_coding") |>
  genes()

length(genes)
seqlevelsStyle(genes) <- "UCSC"


REF.DIR <- "../../DataReference/CTCF/"
ODIR <- "../../Figs/Fig4/"

pvals.df <- fread("../../DataRaw/pvals.bed")
meta.df <- read_tsv(file.path(REF.DIR, "metadata-v2.tsv"))


read_and_cast_to_granges <- function(ff){
  # Used for reference ENCODE data
  DT <- fread(ff)
  names(DT)[1:3] <- c("chr", "start", "end")

  DT %>%
    dplyr::select(c(chr, start, end)) %>%
    makeGRangesFromDataFrame(
      ignore.strand = T,
      starts.in.df.are.0based = T
      )
}


# Pvals needs to be made into GRanges just once
pvals <- makeGRangesFromDataFrame(pvals.df,
                                  keep.extra.columns = T,
                                  ignore.strand = T,
                                  starts.in.df.are.0based = T
)


# Tabulate everything -----------------------------------------------------

for (i in 2:nrow(meta.df)){
  ff <- file.path(REF.DIR, meta.df$ff[i])
  desc <- meta.df$biosample_term_name[i]

  ref <- read_and_cast_to_granges(ff)

  odir <- file.path(ODIR, desc)
  dir.create(odir, showWarnings = F, recursive = T)

  # Overlaps and pull out indices
  overlaps <- findOverlaps(pvals, ref)
  pvals.ix <- queryHits(overlaps)
  ref.ix <- subjectHits(overlaps)

  # This part takes a few minutes
  ref.summary <- run_min_routine(pvals, pvals.ix, ref, ref.ix)
  alpha <- 0.01 / length(ref.summary)

  print("Ran p-value routine")

  # Pull out genes and gene ids
  genes.ix <-
    GenomicRanges::nearest(
      ref.summary[ref.summary$p < alpha],
      genes
    )

  # Summary stats to keep track of
  N_domains <- length(ref)
  N_domains_with_CpG <- length(unique(ref.ix))
  N_CpGs_in_a_domain <-length(pvals.ix)
  N_sig_domains <- sum(ref.summary$p < alpha)
  N_unique_genes <- length(unique(genes[genes.ix]$symbol))

  summary.df <-
    rbind(N_domains, N_domains_with_CpG,
          N_CpGs_in_a_domain, N_sig_domains,
          N_unique_genes
          ) %>%
    as.data.frame() %>%
    rename(count = V1) %>%
    rownames_to_column("statistics")

  ofile <- file.path(odir, "summary_stats.csv")
  write_csv(summary.df, ofile)
  print("Wrote out summary statistics")

  # Gene Ontology
  out <- goplot::run_gene_ontology(genes[genes.ix]$gene_id)
  go.df <- goplot::go_output_to_df(out)

  ofile <- file.path(odir, "gene_ontology.csv")
  write_csv(go.df, ofile)

  # GO Plot
  go.plot <- goplot::plot_go_barchart(go.df) +
    ggtitle(desc)

  ofile <- file.path(odir, "gene_ontology.png")
  cowplot::save_plot(filename = ofile, go.plot,
                     base_width = 14, base_height = 9)
  print("Ran gene ontology")

  # Write significant genes
  sig_genes <- data.frame(significant_genes = genes[genes.ix]$symbol)
  ofile <- file.path(odir, "sig_genes.csv")
  write_csv(x = sig_genes, ofile)
}



