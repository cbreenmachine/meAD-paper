# _targets.R
library(targets)
library(tarchetypes)

tar_source(files = "R")

# Some constants
DATA.REFERENCE.DIR <- "./DataReference/"
DATA.IDIR <- "DataRaw/2023-02-14-Summaries-v6/"

# lFDR Cutoff
ALPHA <- 0.05
DMGENE.ALPHA <- 0.01

# Set target options:
tar_option_set(
  packages = c("tidyverse",
               "tibble",
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
               "rjson",
               "openxlsx",
               "chunkR"
  ),
  format = "rds" # default storage format
)


# Replace the target list below with your own:
preliminaries <- list(

  tar_target(sample.ids, get_sample_ids_from_dss_inputs(
    file.path(DATA.IDIR, "DSS-outputs-chr22.RData"))),

  tar_target(master.full.df, read_csv("DataRaw/masterSamplesheet.csv", show_col_types = F)),
  tar_target(master.df, munge_master_df(master.full.df, sample.ids)),

  tar_target(load.sample.ids, get_group_ids(master.df, "LOAD")),
  tar_target(control.sample.ids, get_group_ids(master.df, "CONTROL")),
  tar_target(mci.sample.ids, get_group_ids(master.full.df, "MCI")),

  # tar_target(apoe.df, get_apoe_allele_frequencies("DataRaw/masterSamplesheet.csv", sample.ids)),

  tar_target(ancestry.test, tabulate_and_test(master.df, "race_primary", "fisher")),
  tar_target(source.test, tabulate_and_test(master.df, "source", "chi")),
  tar_target(sex.test, tabulate_and_test(master.df, "sex", "chi")),
  tar_target(APOE.test, tabulate_and_test(master.df, "APOE_risk_allele", "chi")),

  tar_target(age.test, test_continuous(master.df, "age_at_visit")),
  tar_target(bmi.test, test_continuous(master.df, "bmi")),
  tar_target(education.test, test_continuous(master.df, "education")),

  # tar_target(apoe.test, test_apoe_allele_frequencies(apoe.df)),
  tar_target(madrid.data, read_madrid_data("./DataReference/madrid_cpgs_list.csv")),

  # Integrate with EPIC array later
  tar_target(array.gr, read_850k_array("DataRaw/array.M.hg38.bed")),
  tar_target(common.ids, intersect(sort(colnames(mcols(array.gr))), sample.ids)),

  # Missingness
  # tar_target(missing.chr1.load.df, read_and_join_m_cov_routine("DataRaw/2023-04-28-rawMCov/", load.sample.ids, 20)),
  # tar_target(missing.chr1.control.df, read_and_join_m_cov_routine("DataRaw/2023-04-28-rawMCov/", control.sample.ids, 20)),
  #
  # tar_target(miss.data, combine_missingness_with_pvals(missing.chr1.load.df, missing.chr1.control.df, pvals.gr)),
  # tar_target(miss.hexbin.plot, plot_missingness_hexbin(miss.data, "_targets/figs/missingness-hexbin.chr1.png")),

  ########## GET PUBLIC DATA ##########
  tar_target(chain.19to38, download_chain_from_ucsc(DATA.REFERENCE.DIR, "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz")),
  tar_target(chain.38to19, download_chain_from_ucsc(DATA.REFERENCE.DIR, "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz")),

  # Gene stuff
  tar_target(autosomal.symbols, get_autosomoal_gene_universe()),
  tar_target(autosomal.protein_coding.symbols, get_autosomoal_gene_universe(protein_coding = T)),

  ######## DEFINE CONSTANTS ###########
  tar_target(natgen.symbols, get_natgen_genes("DataReference/NatureGenetics042022_AD_genes.txt")),
  tar_target(diff.exp.data,
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
             format="file")
)

# Array Comparison --------------------------------------------------------


# Keep mapping sub-pipeline outside of the list
raw_wgms_vs_array <- list(
  tar_target(raw.chr1.wgms, join_raw_M_Cov_files("DataRaw/2023-04-28-raw-MCov-chr1/", common.ids)),
  tar_target(raw.chr1.wgms.gr, makeGRangesFromDataFrame(raw.chr1.wgms,
                                                        keep.extra.columns = T,
                                                        starts.in.df.are.0based = T)),
  tar_target(raw.chr1.merged, combine_850k_array_with_wgms(array.gr, raw.chr1.wgms.gr))
)


wgms_vs_array <-
  tar_map(
    unlist = F,
    values = data.frame(chr = paste0("chr", 1:22)),
    tar_target(merged,
               merge_array_and_wgms_subroutine(
                 array.gr,
                 paste0("DataRaw/2023-05-17-MCov-imputed/", chr)
               )),
    tar_target(sampled, sample_frac(merged, 0.1)),
    tar_target(stats, compute_wgms_vs_array_summary_stats(merged))
)


wgms_vs_array_combined <- list(
  tar_combine(wgms_vs_array.df, wgms_vs_array[["sampled"]]),
  tar_combine(wgms_vs_array.stats, wgms_vs_array[["stats"]])
)






# Sequencing quality plot -------------------------------------------------

sequencing_quality <- list(
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
             run_and_save_pvals_plot(pvals.data, "_targets/figs/pvals-before-after.pdf"),
             format="file"),
  tar_target(birdseye.plot,
             stitch_birdseye_fig(
               dmps.volcano.plot,
               dmps.genic.sankey.plot,
               dmps.hyper.pi.plot,
               dmps.cpg.pi.plot,
               "_targets/figs/dmps-summary-panel.png"
               )
             )
)

# DM Genes (pc = protein coding and promoters) --------------------------------
# Check if there are DMPs in any of the Nature Genetics list of genes

diff_methylation <- list(
  tar_target(genes.expanded.25kb,
             get_gene_bodies(upstream = 25000,
                             downstream = 25000,
                             autosomes_only = T,
                             protein_coding_only = T)),
  tar_target(NG.genes.expanded.25kb, filter_by_gene_symbols(genes.expanded.25kb, natgen.symbols)),
  tar_target(NG.genes.with.dmp, make_df_from_two_overlapping_granges(dmps.gr, NG.genes.expanded.25kb)),
  tar_target(NG.tally.df, tally_dmps_in_genes(NG.genes.with.dmp)),

  # Gene Body (GB)
  tar_target(gene.bodies,
             get_gene_bodies(upstream = 3000,
                             downstream = 200,
                             autosomes_only = T,
                             protein_coding_only = T)),
  tar_target(gene.bodies.with.dmp, make_df_from_two_overlapping_granges(dmps.gr, gene.bodies)),
  tar_target(gene.bodies.tally.df, tally_dmps_in_genes(gene.bodies.with.dmp)),
  tar_target(N.DMPs.not.in.genes, get_N_not_in_set(dmps.gr, gene.bodies)),

  # Baby bear: promoters only
  tar_target(promoters, get_protein_coding_promoters(upstream = 5000, downstream = 200)),

  # Combine P-values and such
  tar_target(gene.body.enrichment, harmonic_pvalue_routine(pvals.gr, gene.bodies, ALPHA)),
  tar_target(promoter.enrichment, harmonic_pvalue_routine(pvals.gr, promoters, ALPHA)),

  tar_target(dm.genes.df, dplyr::filter(gene.body.enrichment, lfdr < DMGENE.ALPHA,
                                        N.CpGs > 0)),
  tar_target(madrid.comp,
             compare_with_madrid_paper(
               dm.genes.df,
               madrid.data
              )
  ),


  tar_target(dm.genes.go.df, symbols_to_gene_ontology_routine(dm.genes.df$gene_name)),

  tar_target(dm.genes.gene.ont.df,
             cowplot::save_plot(plot_go_barchart(dm.genes.go.df, n = 25),
                                filename = "_targets/figs/gene-ontology-DMGenes.png",
                                base_height = 7, base_width = 14))
)


# PCHi-C Analysis ---------------------------------------------------------

pchic <- list(
  tar_target(interactions.hg19, clean_interactions_data("./DataReference/PCHi-C/PCHiC_peak_matrix_cutoff5.txt")),
  tar_target(interactions.hg38, lift_promoter_capture_data_to_hg38(interactions.hg19, chain.19to38, return.granges = T)),

  tar_target(enhancers.to.test, interactions.hg38[interactions.hg38$med.chicago > 5]),
  tar_target(baits.to.test, extract_promoters_from_interactions(enhancers.to.test)),
  tar_target(promoters.to.test, subsetByOverlaps(baits.to.test, promoters)),

  # Curate DMPs
  # TODO: turn this into mapped target, lots of duplication
  tar_target(dmps.in.enhancer.df, combine_dmps_with_interactions(dmps.gr, enhancers.to.test)),
  tar_target(dmps.in.promoter.df, combine_dmps_with_interactions(dmps.gr, promoters.to.test)),

  # Mostly questions from RA May 22, 2023
  tar_target(dmp.counts.in.enhancer.df, summarize_dmp_counts_in_enhancers(dmps.in.enhancer.df)),
  tar_target(dmp.counts.in.promoter.df, summarize_dmp_counts_in_promoters(dmps.in.promoter.df)),

  tar_target(shared.gene.symbols,
             get_common_genes_from_DM_interactions(
               dmps.in.enhancer.df,
               dmps.in.promoter.df)),

  tar_target(dm.enhancer.promoter.summary, summarize_counts_dm_interactions(dmps.in.enhancer.df, dmps.in.promoter.df)),

  # Gene ontology for genes with DM promoters and DM enhancers
  tar_target(dm.enhancer.promoter.gene.ont.df, symbols_to_gene_ontology_routine(shared.gene.symbols)),

  tar_target(enhancers.counts, count_enhancers_with_dmp(dmps.in.enhancer.df)),
  tar_target(enhancers.natgen.counts,
             count_enhancers_with_dmp(
               subset_interactions_by_diff_exp_data(dmps.in.enhancer.df, diff.exp.data$gene.name))),

  tar_target(enhancers.summary, summarize_interactions_with_dmp(dmps.in.enhancer.df)),
  tar_target(dm.enhancers.genes.df, summarize_dm_interaction_genes(enhancers.summary)),

  # Promoter stuff
  tar_target(promoters.counts, count_promoters_with_dmp(dmps.in.promoter.df)),
  tar_target(promoters.summary, summarize_interactions_with_dmp(dmps.in.promoter.df)),
  tar_target(dm.promoters.genes.df, summarize_dm_interaction_genes(promoters.summary)),

  tar_target(test.enhancer.enrichment,
             test_enrichment_for_dmps(pvals.gr,
                                      dmps.gr,
                                      enhancers.to.test,
                                      B=10000)),
  # tar_target(test.promoter.enrichment,
  #            test_enrichment_for_dmps(pvals.gr,
  #                                     dmps.gr,
  #                                     promoters.to.test,
  #                                     B=10000)),

  tar_target(plot.enhancer.enrichment,
             plot_enhancer_enrichment_for_dmps(test.enhancer.enrichment,
                                               "_targets/figs/test-enhancer-enrichment.png"),
             format = "file"),

  # tar_target(plot.promoter.enrichment,
  #            plot_enhancer_enrichment_for_dmps(test.promoter.enrichment,
  #                                              "_targets/figs/test-promoter-enrichment.png"),
  #            format = "file"),
  tar_target(DE.vs.DME.test, test_ranks_of_pchic_rnaseq(diff.exp.data, dmps.in.enhancer.df)),
  # tar_target(pchic.stats, get_summary_stats_pchic(enhancers.to.test, dmps.in.enhancer.df)),
  tar_target(reads.stats, get_all_stats_from_dir("DataSummaries/QCReports/")),
  tar_target(gene.list.by.diff.meth.df,
             curate_genes_by_dm_status(
               dm.genes.df$gene_name,
               dm.promoters.genes.df$gene.name,
               dm.enhancers.genes.df$gene.name,
               diff.exp.data$gene.name)
             )
)

# Export UCSC -------------------------------------------------------------

ucsc_exports <- list(
  tar_target(interactions.for.ucsc,
             export_significant_interactions_to_UCSC(enhancers.with.dmp)),
  tar_target(interactions.for.ucsc.bed,
             format_and_write_ucsc_interactions(
               interactions.for.ucsc,
               "_targets/ucsc/interactions-with-DMP.hg38.bb"),
             format = "file"),
  tar_target(dmps.lolly,
             format_and_write_ucsc_lolly(pvals.gr, "_targets/ucsc/DMPs-lolly.hg38.bb", ALPHA),
             format = "file")
)


# Write out the data sets -------------------------------------------------

supplemental_tables <- list(

  tar_target(paper.stats, get_paper_stats(gene.body.enrichment, DMGENE.ALPHA)),

  # 1. DMPs (good)
  tar_target(table.s1.DMPs,
             process_and_write_dmps(
               dmps.gr,
               "Table S1: DMPs.List of Differentially Methylated Positions (DMPs) with coordinates (hg38), effect sizes, and local False-Discovery Rates (lFDRs)",
               "_targets/tables/supplemental/S1-DMPs.xlsx"),
             format = "file"),

  # 2. Nature genetics AD Risk Loci
  tar_target(table.s2.NatGenLoci25kb,
             process_and_write_nature_genetics_list(
               NG.tally.df,
               natgen.symbols,
               "Table S2: List of 75 previously identified genetics risk loci with number of DMPs (if any) within 25,000 nucleotides of gene start/stop",
               "_targets/tables/supplemental/S2-ADRiskLociNumberOfDMPs.xlsx"
             ),
             format = "file"),

  # 3. All genes + Nature genetics annotation
  tar_target(table.s3.DMGenes,
             process_and_write_DM_genes(
               gene.body.enrichment,
               "Table S3: Differentially methylated genes with coordinates and lFDRs",
               "_targets/tables/supplemental/S3-DMGenes.xlsx"
             ),
             format = "file"),


  # 4. Gene ontologies
  tar_target(table.s4.GeneOntologies,
             process_and_write_gene_ontology_terms(
               dm.genes.go.df,
               "Table S4: Gene Ontologies for DM Genes",
               "_targets/tables/supplemental/S4-DMGenes-GeneOntologies.xlsx"),
             format = "file"),

  tar_target(table.s5.DMPromoters,
             format_and_write_dm_promoters(
               promoters.summary,
               "Table S5: Promoter-enhancer interactions with at least one DMP in promoter",
               "_targets/tables/supplemental/S5-DMPromoters.xlsx"),
             format = "file"),
  tar_target(table.s6.DMEnhancers,
             format_and_write_dm_enhancers(
               enhancers.summary,
               "Table S6: Promoter-enhancer interactions with at least one DMP in enhancer",
               "_targets/tables/supplemental/S6-DMEnhancers.xlsx"),
             format = "file"),

  tar_target(table.GenesWithAnyDMFeatures,
             my_write_csv(gene.list.by.diff.meth.df, "_targets/tables/supplemental/GenesWithAnyDMFeatures.csv"),
             format = "file")
)



list(preliminaries,
     raw_wgms_vs_array,
     wgms_vs_array,
     wgms_vs_array_combined,
     diff_methylation,
     pchic,
     sequencing_quality,
     supplemental_tables)

