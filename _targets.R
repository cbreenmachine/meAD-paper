# _targets.R
library(targets)
library(tarchetypes)

tar_source(files = "R")

# Some constants
DATA.REFERENCE.DIR <- "./DataReference/"
DATA.IDIR <- "DataRaw/2023-02-14-Summaries-v6/"

# lFDR Cutoff
ALPHA <- 0.05

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
               "htmlwidgets",
               "rjson"
  ),
  format = "rds" # default storage format
)


# Replace the target list below with your own:
list(

  tar_target(sample.ids, get_sample_ids_from_dss_inputs(
    file.path(DATA.IDIR, "DSS-outputs-chr22.RData"))),

  tar_target(apoe.df, get_apoe_allele_frequencies("DataRaw/masterSamplesheet.csv", sample.ids)),
  tar_target(apoe.test, test_apoe_allele_frequencies(apoe.df)),

  ########## GET PUBLIC DATA ##########
  tar_target(chain.19to38, download_chain_from_ucsc(DATA.REFERENCE.DIR, "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz")),
  tar_target(chain.38to19, download_chain_from_ucsc(DATA.REFERENCE.DIR, "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz")),

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
  tar_target(pvals.ifile, file.path(DATA.IDIR, "pvals.bed"), format = "file"),

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
                                base_height = 4.5, base_width = 5.5),
             format = "file"),
  tar_target(dmps.genic.df, annotate_loci_to_genic_parts(dmps.gr)),
  tar_target(dmps.cpg.df, cpg_annotation_routine(dmps.gr)),
  tar_target(dmps.cpg.pi.plot,
             cowplot::save_plot(plot_cpg_pi_chart(dmps.cpg.df),
                                filename = "_targets/figs/dmps-cpg-island-pi-chart.png",
                                base_height = 4.5, base_width = 5.5),
             format = "file"),
  tar_target(dmps.genic.sankey.plot,
             screenshot_sankey(plot_sankey(dmps.genic.df, nrow(dmps.data)),
                               "_targets/figs/dmps-genic-sankey.png"),
             format="file"),


# Sequencing quality plot -------------------------------------------------

tar_target(depth.df,
           munge_sequencing_depth("./DataSummaries/2023-01-30-AllDepth.txt", sample.ids)),
tar_target(effective.cov.df,
           munge_effective_coverage("./DataSummaries/2023-02-01-AverageCoverage.csv",
                                    sample.ids)),
tar_target(coverage.plot, plot_coverage_hist_routine(
  depth.df, effective.cov.df,
  "_targets/figs/seq-depth-effective-coverage-histograms.png"
), format = "file"),
tar_target(pvals.before.after.plot,
           run_and_save_pvals_plot(pvals.data, "_targets/figs/PvalsBeforeAfter.pdf"),
           format="file"),
tar_target(birdseye.plot,
           stitch_birdseye_fig(
             dmps.volcano.plot,
             dmps.genic.sankey.plot,
             dmps.hyper.pi.plot,
             dmps.cpg.pi.plot,
             "_targets/figs/dmps-summary-panel.png"
             )
           ),

# DM Genes (pc = protein coding and promoters) --------------------------------
# Check if there are DMPs in any of the Nature Genetics list of genes
tar_target(genes.expanded.25kb,
           get_gene_bodies(upstream = 25000,
                           downstream = 25000,
                           autosomes_only = T,
                           protein_coding_only = T)),
tar_target(NG.genes.expanded.25kb, filter_by_gene_symbols(genes.expanded.25kb, natgen.genes)),
tar_target(NG.genes.with.dmp, make_df_from_two_overlapping_granges(dmps.gr, NG.genes.expanded.25kb)),
tar_target(NG.tally.df, tally_dmps_in_genes(NG.genes.with.dmp)),

# Gene Body (GB)
tar_target(gene.bodies,
           get_gene_bodies(upstream = 3000,
                           downstream = 200,
                           autosomes_only = T,
                           protein_coding_only = T)),
tar_target(GB.with.dmp, make_df_from_two_overlapping_granges(dmps.gr, gene.bodies)),
tar_target(GB.tally.df, tally_dmps_in_genes(GB.with.dmp)),

# Baby bear: promoters only
tar_target(promoters, get_protein_coding_promoters(upstream = 10000, downstream = 500)),

# Combine P-values and such
tar_target(gene.body.enrichment, harmonic_pvalue_routine(pvals.gr, gene.bodies, ALPHA)),
tar_target(promoter.enrichment, harmonic_pvalue_routine(pvals.gr, promoters, ALPHA)),


# Gene ontologiy for DM GENES
tar_target(DMGenes.table,
           my_write_csv(dplyr::filter(gene.body.enrichment, lfdr < 0.01),
                        file = "_targets/tables/gene-ontology-DMGenes.csv"
                        ),
           format = "file"),
tar_target(DMGenes.go.df, symbols_df_to_go_df_routine(DMGenes.table)),

tar_target(DMGenes.gene.ont,
           cowplot::save_plot(plot_go_barchart(DMGenes.go.df, n = 25),
                              filename = "_targets/figs/gene-ontology-DMGenes.png",
                              base_height = 7, base_width = 14)),


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


# PCHi-C Analysis ---------------------------------------------------------

tar_target(interactions.hg19, clean_interactions_data("./DataReference/PCHi-C/PCHiC_peak_matrix_cutoff5.txt")),
tar_target(interactions.hg38, lift_promoter_capture_data_to_hg38(interactions.hg19, chain.19to38, return.granges = T)),
tar_target(interactions.to.test, interactions.hg38[interactions.hg38$med.chicago > 5]),
tar_target(interactions.with.dmp, combine_dmps_with_interactions(dmps.gr, interactions.to.test)),
tar_target(interactions.summary, summarize_interactions_with_dmp(interactions.with.dmp)),
tar_target(enhancer_genes.df, summarize_dm_enhancer_genes(interactions.summary)),

tar_target(test.enhancer.enrichment,
           test_enhancer_enrichment_for_dmps(pvals.gr,
                                             dmps.gr,
                                             interactions.to.test,
                                             B=10000)),
tar_target(plot.enhancer.enrichment,
           plot_enhancer_enrichment_for_dmps(test.enhancer.enrichment,
                                             "_targets/figs/test-enhancer-enrichment.png"),
           format = "file"),

tar_target(interactions.for.ucsc,
           export_significant_interactions_to_UCSC(interactions.with.dmp)),
tar_target(DE.vs.DME.test, test_ranks_of_pchic_rnaseq(diff_exp.data, interactions.with.dmp)),

# Export UCSC -------------------------------------------------------------
tar_target(interactions.for.ucsc.bed,
           format_and_write_ucsc_interactions(
             interactions.for.ucsc,
             "_targets/ucsc/interactions-with-DMP.hg38.bb"),
           format = "file"),
tar_target(dmps.lolly,
           format_and_write_ucsc_lolly(pvals.gr, "_targets/ucsc/DMPs-lolly.hg38.bb", ALPHA),
           format = "file"),



# Key summary stats -------------------------------------------------------

tar_target(pchic.stats, get_summary_stats_pchic(interactions.to.test, interactions.with.dmp)),
tar_target(reads.stats, get_all_stats_from_dir("DataSummaries/QCReports/")),

# Write out the data sets -------------------------------------------------

# 1. DMPs (good)
tar_target(table.s1.DMPs,
           process_and_write_dmps(
             dmps.gr,
             "_targets/tables/supplemental/S1-DMPs.xlsx"),
           format = "file"),

# 2. Nature genetics AD Risk Loci
tar_target(table.s2.NatGenLoci25kb,
           process_and_write_nature_genetics_list(
             NG.tally.df,
             natgen.genes,
             "_targets/tables/supplemental/S2-ADRiskLociNumberOfDMPs.xlsx"
           ),
           format = "file"),

# 3. All genes + Nature genetics annotation
tar_target(table.s3.DMGenes,
           process_and_write_DM_genes(
             gene.body.enrichment,
             "_targets/tables/supplemental/S3-DMGenes.xlsx"
           ),
           format = "file"),


# 4. Gene ontologies
tar_target(table.s4.GeneOntologies,
           process_and_write_gene_ontology_terms(
             DMGenes.go.df,
             "_targets/tables/supplemental/S4-DMGenes-GeneOntologies.xlsx"),
           format = "file"),

tar_target(table.s5.DMEnhancers,
           format_and_write_interactions(
             interactions.summary,
             "_targets/tables/supplemental/S5-DMEnhancers.xlsx"),
           format = "file")

)
