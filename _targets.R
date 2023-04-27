# _targets.R
library(targets)
library(tarchetypes)

source("R/theme_meAD.R")
source("R/functions.R")
source("R/functions_summarize_DMPs.R")
source("R/functions_PCHiC.R")
source("R/functions_for_supplement.R")

#TODO: plots and saving is a little messy

# Some constants
DATA.REFERENCE.DIR <- "./DataReference/"

# lFDR Cutoff
ALPHA.CUT <- 0.05


# Set target options:
tar_option_set(
  packages = c("tidyverse",
               "readxl",
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
               "harmonicmeanp",
               "fastqq",
               "devEMF",
               "htmlwidgets"
  ),
  format = "rds" # default storage format
)


# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Replace the target list below with your own:
list(

  ########## GET PUBLIC DATA ##########
  tar_target(chain.19to38, download_chain_from_ucsc(DATA.REFERENCE.DIR, "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz")),
  tar_target(chain.38to19, download_chain_from_ucsc(DATA.REFERENCE.DIR, "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz")),
  # TODO: tar_target(diff.exp.genes, download_DE_genes)

  # Gene stuff
  tar_target(autosomal.symbols, get_autosomoal_gene_universe()),
  tar_target(autosomal.protein_coding.symbols, get_autosomoal_gene_universe(protein_coding = T)),

  ######## DEFINE CONSTANTS ###########
  tar_target(natgen.genes, get_natgen_genes("DataReference/NatureGenetics042022_AD_genes.txt")),
  tar_target(diff_exp.data,
             clean_differentially_expressed_genes(
               "DataReference/DEGenes/2020-Shigemizu-AD-RNAseq-DEGenes.xlsx",
               autosomal.symbols)
             ),
  tar_target(pvals.ifile, "DataRaw/2023-02-14-Summaries-v6/pvals.bed", format = "file"),

  ######### GENE STUFF ################


  ######### Wrangle DMPs #############
  tar_target(pvals.data, get_pvals_data(pvals.ifile)),
  tar_target(pvals.gr, to_granges(pvals.data)),
  tar_target(dmps.data, subset_dmps(pvals.data)),
  tar_target(dmps.gr, to_granges(dmps.data)),

  # Bird's eye plotting
  tar_target(dmps.volcano.plot,
             cowplot::save_plot(volcano_routine(pvals.data),
                                filename = "_targets/figs/dmps-volcano-scatter.png",
                                base_height = 6, base_width = 6),
             format="file"),
  tar_target(dmps.hyper.pi.plot,
             cowplot::save_plot(dmp_pi_chart_routine(dmps.data),
                                filename = "_targets/figs/dmps-hyper-hypo-pi-chart.png",
                                base_height = 4, base_width = 5),
             format = "file"),
  tar_target(dmps.genic.df, annotate_loci_to_genic_parts(dmps.gr)),
  tar_target(dmps.cpg.df, cpg_annotation_routine(dmps.gr)),
  tar_target(dmps.cpg.pi.plot,
             cowplot::save_plot(plot_cpg_pi_chart(dmps.cpg.df),
                                filename = "_targets/figs/dmps-cpg-island-pi-chart.png",
                                base_height = 4, base_width = 5),
             format = "file"),
  tar_target(dmps.genic.sankey.plot,
             screenshot_sankey(plot_sankey(dmps.genic.df, nrow(dmps.data)),
                               "_targets/figs/dmps-genic-sankey.png"),
             format="file"),


  # Pvalues plot
  tar_target(pvals.before.after.plot,
             run_and_save_pvals_plot(pvals.data, "_targets/figs/PvalsBeforeAfter.pdf"),
             format="file"),


  tar_target(birdseye.plot,
             stitch_birdseye_fig(
               dmps.volcano.plot,
               dmps.genic.sankey.plot,
               dmps.cpg.pi.plot,
               dmps.hyper.pi.plot,
               "_targets/figs/dmps-summary-panel.png"
               )
             ),

  # Promoters
  tar_target(proco.gene.bodies, get_protein_coding_gene_bodies(upstream = 3000, downstream = 200)),
  tar_target(proco.promoters, get_protein_coding_promoters(upstream = 10000, downstream = 500)),

  tar_target(gene.body.enrichment, harmonic_pvalue_routine(pvals.gr, proco.gene.bodies)),
  tar_target(promoter.enrichment, harmonic_pvalue_routine(pvals.gr, proco.promoters)),

  # Gene ontologiy for DM GENES
  tar_target(DMGenes.table,
             my_write_csv(dplyr::filter(gene.body.enrichment, lfdr < 0.05),
                          file = "_targets/tables/gene-ontology-DMGenes.csv"
                          ),
             format = "file"),
  tar_target(DMGenes.go.df, symbols_df_to_go_df_routine(DMGenes.table)),

  tar_target(DMGenes.gene.ont,
             cowplot::save_plot(plot_go_barchart(DMGenes.go.df, n=25),
                                filename = "_targets/figs/gene-ontology-DMGenes.png",
                                base_height = 7, base_width = 14)),

  # Write some subsets
  tar_target(gene.body.natgen.table,
             my_write_csv(filter_by_natgen(gene.body.enrichment, natgen.genes),
                          file = "_targets/tables/gene-bodies-dm-subset-by-natgen.csv"
              ),
             format = "file"),
  tar_target(promoter.natgen.table,
             my_write_csv(
               filter_by_natgen(promoter.enrichment, natgen.genes),
               file = "_targets/tables/promoters-dm-subset-by-natgen.csv"
              ),
             format = "file"
             ),
  tar_target(gene.body.table,
             my_write_csv(dplyr::arrange(gene.body.enrichment, lfdr), "_targets/tables/gene-bodies-dm.csv"),
             format = "file"),
  tar_target(promoter.table,
             my_write_csv(dplyr::arrange(promoter.enrichment, lfdr),  "_targets/tables/promoters-dm.csv"),
             format = "file"),

  # Promoter Capture Hi-C analysis
  tar_target(interactions.hg19, clean_interactions_data("./DataReference/PCHi-C/PCHiC_peak_matrix_cutoff5.txt")),
  tar_target(interactions.hg38, lift_promoter_capture_data_to_hg38(interactions.hg19, chain.19to38, return.granges = T)),
  tar_target(interactions.to.test, interactions.hg38[interactions.hg38$med.chicago > 5]),
  tar_target(interactions.with.dmp, combine_dmps_with_interactions(dmps.gr, interactions.to.test)),
  tar_target(interactions.summary, summarize_interactions_with_dmp(interactions.with.dmp)),
  tar_target(enhancer_genes.df, summarize_dm_enhancer_genes(interactions.summary)),
  # tar_target(gene_ont.dm_enhancers.plot,
  #            symbols_to_go_plot(unique(enhancer_genes.df$gene.name),
  #                               "_targets/figs/gene-ontology-dm-enhancers.png"),
  #            format = "file"),

  # PCHi-C Enrichment Analysis
  tar_target(test.enhancer.enrichment,
             test_enhancer_enrichment_for_dmps(pvals.gr,
                                               dmps.gr,
                                               interactions.to.test,
                                               B=10000)),
  tar_target(plot.enhancer.enrichment,
             plot_enhancer_enrichment_for_dmps(test.enhancer.enrichment,
                                               "_targets/figs/test-enhancer-enrichment.png"),
             format = "file"),




  # Integrate with expression data


  # Export UCSC data
  tar_target(dmps.lolly,
             format_and_write_ucsc_lolly(pvals.gr, "_targets/ucsc/DMPs.lolly.bb", ALPHA.CUT),
             format = "file")

  # PCHi-C
  # Get remote data, etc.
)
