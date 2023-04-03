# _targets.R
library(targets)
library(tarchetypes)

source("R/functions.R")
source("R/functions_birdseye_analysis.R")
source("R/functions_birdseye_plot.R")


# Some constants
DATA.REFERENCE.DIR <- "./DataReference/"


# Set target options:
tar_option_set(
  packages = c("tidyverse",
               "data.table",
               "viridis",
               "GenomicRanges",
               "cowplot",
               "annotatr",
               "waffle",
               "networkD3",
               "TxDb.Hsapiens.UCSC.hg38.knownGene",
               "EnsDb.Hsapiens.v86",
               "ggsci",
               "harmonicmeanp"),
  format = "rds" # default storage format
)


# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Replace the target list below with your own:
list(

  ########## GET PUBLIC DATA ##########
  tar_target(chain.19to38, download_chain_from_ucsc(DATA.REFERENCE.DIR, "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz")),
  tar_target(chain.38to19, download_chain_from_ucsc(DATA.REFERENCE.DIR, "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz")),

  ######## DEFINE CONSTANTS ###########
  tar_target(natgen.genes, get_natgen_genes("DataReference/NatureGenetics042022_AD_genes.txt")),
  tar_target(pvals.ifile, "DataRaw/2023-02-14-Summaries-v6/pvals.bed", format = "file"),

  # General data wrangling
  tar_target(pvals.data, get_pvals_data(pvals.ifile)),
  tar_target(pvals.gr, to_granges(pvals.data)),
  tar_target(dmps.data, subset_dmps(pvals.data)),
  tar_target(dmps.gr, to_granges(dmps.data)),
  # Bird's eye plotting
  tar_target(volcano.plot,
             cowplot::save_plot(volcano_routine(pvals.data),
                                filename = "_targets/figs/fig2-birds-eye/dmps-volcano-scatter.png",
                                base_height = 6, base_width = 6),
             format="file"),
  tar_target(pi.plot,
             cowplot::save_plot(dmp_pi_chart_routine(dmps.data),
                                filename = "_targets/figs/fig2-birds-eye/dmps-hyper-hypo-pi.png",
                                base_height = 6, base_width = 6),
             format = "file"),
  tar_target(dmps.genic.df, annotate_loci_to_genic_parts(dmps.gr)),
  tar_target(dmps.cpg.df, cpg_annotation_routine(dmps.gr)),
  tar_target(cpg.island.plot,
             cowplot::save_plot(plot_cpg_barchart(dmps.cpg.df),
                                filename = "_targets/figs/fig2-birds-eye/dmps-cpg-island-barchart.png",
                                base_height = 3, base_width = 6),
             format = "file"),
  tar_target(genic.sankey.plot,
             screenshot_sankey(plot_sankey(dmps.genic.df, nrow(dmps.data)),
                               "_targets/figs/fig2-birds-eye/dmps-genic-sankey.png"),
             format="file"),

  # Promoters
  tar_target(proco.gene.bodies, get_protein_coding_gene_bodies(upstream = 3000, downstream = 200)),
  tar_target(proco.promoters, get_protein_coding_promoters(upstream = 10000, downstream = 500)),

  tar_target(gene.body.enrichment, harmonic_pvalue_routine(pvals.gr, proco.gene.bodies)),
  tar_target(promoter.enrichment, harmonic_pvalue_routine(pvals.gr, proco.promoters)),

  # Write some subsets
  tar_target(gene.body.natgen.table,
             my_write_csv(filter_by_natgen(gene.body.enrichment, natgen.genes),
                          file = "_targets/tables/gene-bodies-dm-subset-by-natgen.csv"),
             format = "file"),
  tar_target(promoter.natgen.table,
             my_write_csv(
               filter_by_natgen(promoter.enrichment, natgen.genes),
               file = "_targets/tables/promoters-dm-subset-by-natgen.csv"),
             format = "file"
             ),
  tar_target(gene.body.table,
             my_write_csv(dplyr::arrange(gene.body.enrichment, lfdr), "_targets/tables/gene-bodies-dm.csv"),
             format = "file"),
  tar_target(promoter.table,
             my_write_csv(dplyr::arrange(promoter.enrichment, lfdr),  "_targets/tables/promoters-dm.csv"),
             format = "file"),

  # Promoter analysis
  tar_target(interactions.hg19, clean_and_filter_interactions("./DataReference/PCHi-C/PCHiC_peak_matrix_cutoff5.txt")),
  tar_target(interactions.hg38, lift_promoter_capture_data_to_hg38(interactions.hg19, chain.19to38, return.granges = T))

  # PCHi-C
  # Get remote data, etc.
)
