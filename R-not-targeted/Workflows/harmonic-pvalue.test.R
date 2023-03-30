library(EnsDb.Hsapiens.v86)
library(goplot)
library(fdrtool)
library(harmonicmeanp)


# Setup -------------------------------------------------------------------
db <- EnsDb.Hsapiens.v86

pvals <- fread("../../DataRaw/2023-02-14-Summaries-v6/pvals.bed") %>%
  dplyr::rename(pval = p.from.ss) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T, starts.in.df.are.0based = T)


db.genes <- genes(db)
db.genes <- db.genes[db.genes$gene_biotype == "protein_coding"]


db.promoters <- promoters(db, upstream = 25000)
length(db.promoters)
db.promoters <- db.promoters[db.promoters$tx_biotype == "protein_coding"]
length(unique(db.promoters$gene_id))

seqlevelsStyle(pvals) <- "NCBI"

# How far -----------------------------------------------------------------

aggregate_pvals <- function(pvals, features){
  overlaps <- findOverlaps(pvals, features)
  pvals.ix <- queryHits(overlaps)
  features.ix <- subjectHits(overlaps)

  result <- data.frame(pp = pvals$pval[pvals.ix],
                       gene.id = features$gene_id[features.ix])

  # Combine the p-values, make sure they're in the unit interval, drop
  # features with no CpGs
  aggregated <- result %>%
    group_by(gene.id) %>%
    summarize(
      k = length(pp),
      p.hmp = p.hmp(pp, L = k)) %>%
    dplyr::mutate(p.hmp = pmax(pmin(p.hmp, 1), 0)) %>%
    dplyr::filter(k > 0)

  # Add a local FDR
  pp <- aggregated$p.hmp
  fdr.out <- fdrtool::fdrtool(pp, statistic = "pvalue", plot = F)

  aggregated$lfdr <- fdr.out$lfdr
  aggregated
}


# Gene ontology -----------------------------------------------------------

promoters.agg <- aggregate_pvals(pvals, db.promoters)
genes.agg <- aggregate_pvals(pvals, db.genes)


promoter.dm.ids <- promoters.agg$gene.id[promoters.agg$lfdr < 0.01] %>% unique()
gene.dm.ids <- genes.agg$gene.id[genes.agg$lfdr < 0.01] %>% unique()

length(unique(promoter.dm.ids))
length(unique(gene.dm.ids))

length(intersect(promoter.dm.ids, gene.dm.ids)) / length(union(promoter.dm.ids, gene.dm.ids))

length(intersect(promoter.dm.ids, gene.dm.ids))

# Plot --------------------------------------------------------------------

# Gene ontology on top 1000 genes?
go.out <- goplot::run_gene_ontology(gene.dm.ids)

go.out %>%
  go_output_to_df() %>%
  plot_go_barchart()
