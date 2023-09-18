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
               "gridExtra",
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
               "chunkR",
               "minfi",
               "IlluminaHumanMethylation450kanno.ilmn12.hg19",
               "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
  ),
  format = "rds" # default storage format
)


# Replace the target list below with your own:


# Preliminaries -----------------------------------------------------------

# 1. Getting data lifted to hg38
# 2. Performing quality controls
# 3. Tests on demographic data
# 4.

preliminaries <- list(

  tar_target(sample.ids, get_sample_ids_from_dss_inputs(
    file.path(DATA.IDIR, "DSS-outputs-chr22.RData"))),

  tar_target(master.full.df, read_csv("DataRaw/masterSamplesheet.csv", show_col_types = F)),
  tar_target(master.df, munge_master_df(master.full.df, sample.ids)),

  tar_target(load.sample.ids, get_group_ids(master.df, "LOAD")),
  tar_target(control.sample.ids, get_group_ids(master.df, "CONTROL")),
  tar_target(mci.sample.ids, get_group_ids(master.full.df, "MCI")),

  tar_target(ancestry.test, tabulate_and_test(master.df, "race_primary", "fisher")),
  tar_target(source.test, tabulate_and_test(master.df, "source", "chi")),
  tar_target(sex.test, tabulate_and_test(master.df, "sex", "chi")),
  tar_target(APOE.test, tabulate_and_test(master.df, "APOE_risk_allele", "chi")),

  tar_target(age.test, test_continuous(master.df, "age_at_visit")),
  tar_target(bmi.test, test_continuous(master.df, "bmi")),
  tar_target(education.test, test_continuous(master.df, "education")),

  # How many DMPs are shared with previously discovered DMPs?
  tar_target(dmps.array.gr, read_and_cast_madrid_data("./DataReference/madrid_cpgs_list.csv")),
  tar_target(wgms.sig.in.array.gr, subsetByOverlaps(pvals.gr, dmps.array.gr, minoverlap = 2)),
  tar_target(array.gene.symbols, get_array_genes_with_nearby_dmp(dmps.array.gr)),
  tar_target(wgms.vs.array.genes, compare_with_madrid_paper(dmps.array.gr, dm.genes.df)),
  tar_target(wgms.vs.array.genes.csv,
             my_write_csv(data.frame(shared_genes = wgms.vs.array.genes), "DataSummaries/2023-06-26-SharedGenesForKirk.csv")),

  # Integrate with EPIC array later
  tar_target(array.gr, read_850k_array("DataRaw/array.M.hg38.bed")),
  tar_target(common.ids, intersect(sort(colnames(mcols(array.gr))), sample.ids)),

  # Missingness - takes too long to run
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
  tar_target(natgen.symbols, get_natgen_genes("DataReference/NatureGenetics042022_AD_genes.txt",
                                              exclude_APP = T,
                                              exclude_IGH = T)),
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
             format="file"),
  tar_target(dmps.in.gene.stats, tally_dmps_in_out_genes(dmps.gr, gene.bodies))
)

# Array Comparison --------------------------------------------------------


# Keep mapping sub-pipeline outside of the list
# raw_wgms_vs_array <- list(
#   tar_target(raw.chr1.wgms, join_raw_M_Cov_files("DataRaw/2023-04-28-raw-MCov-chr1/", common.ids)),
#   tar_target(raw.chr1.wgms.gr, makeGRangesFromDataFrame(raw.chr1.wgms,
#                                                         keep.extra.columns = T,
#                                                         starts.in.df.are.0based = T)),
#   tar_target(raw.chr1.merged, combine_850k_array_with_wgms(array.gr, raw.chr1.wgms.gr))
# )



# WGMS vs Array Comparisons -----------------------------------------------

wgms_vs_array <-
  tar_map(
    unlist = F,
    values = data.frame(chr = paste0("chr", 1:22)),
    tar_target(merged,
               merge_array_and_wgms_subroutine(
                 array.gr,
                 paste0("DataRaw/2023-05-17-MCov-imputed/", chr)
               )),

    # Derived from merged
    tar_target(mid,
               dplyr::filter(merged,
                             WGMS > 0.1, WGMS < 0.9,
                             EPIC > 0.1, EPIC < 0.9)),
    tar_target(sampled, sample_frac(merged, 0.1)),
    tar_target(merged_no_dmps, remove_dmps_from_merged(merged, dmps.gr)),

    # Sampled estimates identical
    tar_target(stats, compute_wgms_vs_array_summary_stats(merged)),
    tar_target(stats_mid, compute_wgms_vs_array_summary_stats(mid)),
    tar_target(stats_no_dmps, compute_wgms_vs_array_summary_stats(merged_no_dmps)),
    tar_target(scatter, plot_tech_comp_scatter(merged))
)


wgms_vs_array_combined <- list(
  tar_combine(wgms_vs_array.df, wgms_vs_array[["sampled"]]),
  tar_combine(wgms_vs_array.stats, wgms_vs_array[["stats"]]),
  tar_combine(wgms_vs_array_mid.stats, wgms_vs_array[["stats_mid"]]),

  tar_combine(wgms_vs_array_no_dmps.stats, wgms_vs_array[["stats_no_dmps"]]),
  tar_target(N_loci_dropped,
             count_number_of_loci_dropped(
               wgms_vs_array.stats,
               wgms_vs_array_no_dmps.stats)
             )
)


wgms_vs_array_correlations <- list(
  tar_target(mid.avg.cor, compute_average_correlation(wgms_vs_array_mid.stats)),
  tar_target(avg.cor, compute_average_correlation(wgms_vs_array.stats))
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
  tar_target(dm.genes.go.df, symbols_to_gene_ontology_routine(dm.genes.df$gene_name)),

  tar_target(dm.genes.go.plot,
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

  # Ensure that baits are nearby genes (SK advised against this)
  # tar_target(promoters.to.test, subsetByOverlaps(baits.to.test, promoters)),
  tar_target(promoters.to.test, baits.to.test),

  # Curate DMPs
  # TODO: turn this into mapped target, lots of duplication because enhancer and promoter
  # stuff is paralleled
  tar_target(dmps.in.enhancer.df, combine_dmps_with_intervals(dmps.gr, enhancers.to.test)),
  tar_target(dmps.in.promoter.df, combine_dmps_with_intervals(dmps.gr, promoters.to.test)),

  # Mostly questions from RA May 22, 2023
  # Since more than one DMP can reside in an enhancer, we take the median effect size
  # and tally the number of DMPs. This gives you one row per enhancer/promoter respectively.
  tar_target(dmp.counts.in.enhancer.df, summarize_dmp_counts_in_pchic(dmps.in.enhancer.df, "oe.id", filter.zeros = T)),
  tar_target(dmp.counts.in.promoter.df, summarize_dmp_counts_in_pchic(dmps.in.promoter.df, "bait.id", filter.zeros = T)),

  # Extract Genes with DM enhancers, Genes with DM Promoters.
  tar_target(genes.w.dm.enhancers, get_unique_genes_from_baitName_col(dmp.counts.in.enhancer.df)),
  tar_target(genes.w.dm.promoters, get_unique_genes_from_baitName_col(dmp.counts.in.promoter.df)),

  tar_target(pchic.summary.stats,
             tabulate_pchic_findings(
               enhancers.to.test,
               baits.to.test,
               dmp.counts.in.enhancer.df,
               dmp.counts.in.promoter.df,
               genes.w.dm.enhancers,
               genes.w.dm.promoters
             )),
  tar_target(pchic.summary.stats.csv,
             my_write_csv(pchic.summary.stats, "_targets/tables/supplemental/PCHiC-integration-counts.csv"),
             format = "file"),

  # Now re-combine by interaction, integrate diff expression and GWAS by adding columns indicating
  # whether baitNames overlap with those respective sets
  tar_target(dmp.enhancer.promoter.integrated, combine_enhancer_promoter_gwas_diff_exp(
    dmp.counts.in.enhancer.df,
    dmp.counts.in.promoter.df,
    natgen.symbols,
    diff.exp.data$gene.name
  )),





  tar_target(genes.with.dm.enhancer.and.dm.promoter,
             get_common_genes_from_DM_interactions(
               dmps.in.enhancer.df,
               dmps.in.promoter.df)),

  tar_target(dm.interactions.summary, summarize_counts_dm_interactions(dmps.in.enhancer.df, dmps.in.promoter.df)),

  # Gene ontology for genes with DM promoters and DM enhancers
  tar_target(dm.enhancer.promoter.gene.ont.df, symbols_to_gene_ontology_routine(genes.with.dm.enhancer.and.dm.promoter)),

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

  # Summary stats (updated June 20, 2023)
  # To remedy some confusion with DM enhancers vs interactions with DM enhancers
  tar_target(test.enhancer.enrichment,
             test_enrichment_for_dmps(pvals.gr,
                                      dmps.gr,
                                      enhancers.to.test,
                                      B=10000)),
  tar_target(test.promoter.enrichment,
             test_enrichment_for_dmps(pvals.gr,
                                      dmps.gr,
                                      promoters.to.test,
                                      B=10000)),

  tar_target(plot.enhancer.enrichment,
             plot_enhancer_enrichment_for_dmps(test.enhancer.enrichment,
                                               "_targets/figs/test-enhancer-enrichment.png"),
             format = "file"),

  tar_target(DE.vs.DME.test, test_ranks_of_pchic_rnaseq(diff.exp.data, dmps.in.enhancer.df)),
  # tar_target(pchic.stats, get_summary_stats_pchic(enhancers.to.test, dmps.in.enhancer.df)),
  tar_target(reads.stats, get_all_stats_from_dir("DataSummaries/QCReports/")),
  tar_target(gene.list.by.diff.meth.df,
             curate_genes_by_dm_status(
               dm.genes.df$gene_name,
               dm.promoters.genes.df$gene.name,
               dm.enhancers.genes.df$gene.name,
               diff.exp.data$gene.name,
               natgen.symbols,
               array.gene.symbols)
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
               genes.w.dm.enhancers,
               genes.w.dm.promoters,
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
               diff.exp.data$gene.name,
               "Table S5: Promoter-enhancer interactions with at least one DMP in promoter",
               "_targets/tables/supplemental/S5-InteractionsWithDMPromoters.xlsx"),
             format = "file"),
  tar_target(table.s6.DMEnhancers,
             format_and_write_dm_enhancers(
               enhancers.summary,
               diff.exp.data$gene.name,
               "Table S6: Promoter-enhancer interactions with at least one DMP in enhancer",
               "_targets/tables/supplemental/S6-InteractionWithDMEnhancers.xlsx"),
             format = "file"),

  tar_target(table.GenesWithAnyDMFeatures,
             my_write_csv(gene.list.by.diff.meth.df, "_targets/tables/supplemental/GenesWithAnyDMFeatures.csv"),
             format = "file")
)

# How many DMPs are on each of the commonly used Illumina arrays?
tallies_from_array_techs <- list(
  tar_target(array.450.hg38, load_and_lift_450k(chain.19to38)),
  tar_target(array.EPIC.hg38, load_and_lift_EPIC(chain.19to38)),

  #
  tar_target(array.450.cpgs, subsetByOverlaps(pvals.gr, array.450.hg38)),
  tar_target(array.EPIC.cpgs, subsetByOverlaps(pvals.gr, array.EPIC.hg38)),

  #
  tar_target(array.450.dmps, subsetByOverlaps(dmps.gr, array.450.hg38)),
  tar_target(array.EPIC.dmps, subsetByOverlaps(dmps.gr, array.EPIC.hg38))
)



revisions <- list(
  # meQTL Comparison
  tar_target(assoc.df, read_GODMC_meQTLs("DataReference/GODMC-mQTL/assoc_meta_all.csv")),
  tar_target(godmc.hg19.gr, filter_450k_annotation(assoc.df$cpg)),
  tar_target(godmc.hg38.gr, unlist(rtracklayer::liftOver(godmc.hg19.gr, chain.19to38))),
  tar_target(dmp.and.mqtl, subsetByOverlaps(dmps.gr, godmc.hg38.gr)),

  tar_target(dmp.godmc.table,
             make_godmc_dmp_table(
               dmp.and.mqtl,
               array.450.dmps,
               godmc.hg38.gr,
               array.450.hg38)
             ),

  # Other Study's DMPs comparison
  tar_target(rou.dmps.hg19, read_roubroeks_data("DataReference/2020-Roubroeks/anova-dmps.xlsx")),
  tar_target(rou.dmps.hg38, unlist(rtracklayer::liftOver(rou.dmps.hg19, chain.19to38))),
  tar_target(rou.common.dmps, subsetByOverlaps(dmps.gr, rou.dmps.hg38)),

  # How many DMPs go to more than one gene
  tar_target(N.genes.per.dmp, count_n_dmps_multi_map_to_gene(dmps.gr, gene.bodies))
)

reports <- list(
  tar_render(report_post_fdr_stats, "2023-08-24-fdrtool-genomic-control.Rmd")
)


list(
  preliminaries,
  # raw_wgms_vs_array,
  wgms_vs_array,
  wgms_vs_array_combined,
  wgms_vs_array_correlations,
  diff_methylation,
  pchic,
  sequencing_quality,
  supplemental_tables,
  revisions,
  tallies_from_array_techs,
  reports
  )

