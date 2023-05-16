# fun.call-dm-features.R
# Functions to call differentially methylated features
#

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(parallel)
  library(harmonicmeanp)
})


library(EnsDb.Hsapiens.v86)

db <- EnsDb.Hsapiens.v86
# genes(db, filter = GeneBiotypeFilter('protein_coding'))

# Example CpGs
pvals <- fread("../../DataRaw/2023-02-14-Summaries-v6/pvals.bed") %>%
  dplyr::select(c(chr, start, end, p.from.ss)) %>%
  makeGRangesFromDataFrame(starts.in.df.are.0based = T,
                           keep.extra.columns = T)

# Example set of features
features <- promoters(db, filter = GeneBiotypeFilter('protein_coding'),
                      upstream = 25000, columns = c("gene_name", "gene_id"))


seqlevelsStyle(pvals) <- "NCBI"

# Steps...
# First,

mcols(subsetByOverlaps(pvals, features[100]))[ ,1]


subset_by_features <- function(data, features, ix){
  # First column should be p-values...
  pp <- mcols(subsetByOverlaps(data, features[ix]))[ ,1]
  pp
}

compute_hmp <- function(pp){
  # k is the number of p-values (numnber of tests)
  k <- length(pp)

  if (k > 0){
    harmonic.mean.p <- p.hmp(pp, L = k)
    return(c(harmonic.mean.p, k))
  } else {
    return(c(-1, 0))
  }
}

routine <- function(data, features, ix){
  pp <- subset_by_features(data, features, ix)
  out <- compute_hmp(pp)
  out
}



# Demonstration -----------------------------------------------------------

pp <- subset_by_features(pvals, features, 100)

# Harmoinc mean
harmonic.mean.p <- compute_hmp(pp)
harmonic.mean.p

# Or all in one spote
routine(pvals, features, 2)

# Wrap to use in parallel
wrapper <- function(ix){
  routine(pvals, features, ix)
}



# Parallelize -------------------------------------------------------------

out <- mclapply(1:5000, wrapper, mc.cores = 4)

bound.df <- as.data.frame(do.call(rbind, out))
names(bound.df) <- c("hmp", "k")

bound.df <- bound.df %>%
  filter(hmp >= 0) %>%
  filter(hmp <= 1)


# Fdr correction ----------------------------------------------------------

fdr.out <- fdrtool::fdrtool(bound.df$hmp, statistic = "pvalue", plot = F)

hist(fdr.out$lfdr)
