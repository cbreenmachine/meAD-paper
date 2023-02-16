# Results
# List of genes
# And how the DMR behaves?
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)
library(org.Hs.eg.db)
library(tidyverse)



# Variables ---------------------------------------------------------------
IFILE <- "../../DataRaw/2023-01-24-Summaries-v4/pvals.bed"
PCHIC.FILE <- "../../DataReference/PCHi-C/PCHiC_peak_matrix_cutoff5.txt"
REF.CHAIN.PATH <- "../../DataReference/hg38ToHg19.over.chain"
OFILE <- "../../Figs/PCHi-C/test-noverlaps.png"

# Enhancers
enhancers <- fread(PCHIC.FILE)  %>%
  makeGRangesFromDataFrame(seqnames.field = "oeChr",
                           start.field = "oeStart",
                           end.field = "oeEnd",
                           keep.extra.columns = T)

seqlevelsStyle(enhancers) <- "UCSC"

# CpGs need to be cast to hg19
pvals.hg38 <- fread(IFILE) %>%
  makeGRangesFromDataFrame(starts.in.df.are.0based = T,
                           keep.extra.columns = T)

chain <- import.chain(REF.CHAIN.PATH)
pvals.hg19 <- unlist(liftOver(pvals.hg38, chain))


# Significnat points ------------------------------------------------------
sig.ix <- pvals.hg19$lfdr < 0.05
N.DMPs <- sum(sig.ix)

# Overlap DMRs with PCHi-C ------------------------------------------------
overlap_and_tally <- function(pvals.sub){
  # Ensure that pvals.sub is the same length every time
  # Should be ~ 1300 (i.e. number of DMPs)
  overlaps <- findOverlaps(pvals.sub, enhancers)

  # Return the number of unique DMPs
  # residing in an enhancer
  # (this allows a DMP to be in multiple enhancers)
  length(unique(queryHits(overlaps)))
}

N.obs <- overlap_and_tally(pvals.hg19[sig.ix])


# Test --------------------------------------------------------------------

# Number of CpGs to samples from
G <- length(pvals.hg19)

# Number of simulations
B <- 1000
N.empirical <- rep(-1, B)

# Raleigh aka City of Oaks
set.seed(919)

for (b in 1:B){

  ix <- sample(1:G, N.DMPs, replace = F)
  N.empirical[b] <- overlap_and_tally(pvals.hg19[ix])
}


p <- 1 - (sum(N.obs < N.empirical) / B)

png(OFILE)
hist(N.empirical, xlab = paste0("One-sided p-value: ", p))
abline(v=N.obs,col="blue",lwd=2)

