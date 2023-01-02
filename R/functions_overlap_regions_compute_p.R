# functions_overlap_regions_compute_p.R
# Given two GenomicRegions (one with pvals)
# and one with regions to overlap,
# we can compute the p-values with min(p in that region)
# or run Fisher's method

library(tidyverse)
library(data.table)
library(poolr)
library(GenomicRanges)



# Tabulate everything -----------------------------------------------------

get_overlapped_ix <- function(pvals, ref){

  # This does all the computational work
  overlaps <- findOverlaps(pvals, ref)

  # Usually pvals.ix will look like
  # 1,2,3,4,5,6,...
  # and ref.ix will look like
  # 1,1,1,1, 3,3, 4, ...
  # indicating the first four entries of `pvals`
  # are in entry 1 of `ref`
  pvals.ix <- queryHits(overlaps)
  ref.ix <- subjectHits(overlaps)


  return(list(pvals.ix, ref.ix, ref.ix.range))
}





run_min_routine <-
  function(pvals,
           pvals.ix,
           ref,
           ref.ix){

    # Indices of the reference CTCF regions
    # that have at least one overlapping CpG
    # (what to loop over)
    ref.ix.range <- sort(unique(ref.ix))

    # z indexes the reference entries
    # Which CpGs are in the zth ref entry
    wrapper <- function(z){
      # used in lapply functions
      keep.ix <- (ref.ix == z)
      pp <- pvals[pvals.ix[keep.ix]]$p.corrected
      min(pp)
    }


    # Run the Fisher routine on all the CTCF binding domains
    ref.p <-
      do.call(rbind,
              parallel::mclapply(
                X=ref.ix.range,
                FUN=wrapper,
                mc.cores = 6
              )
      )

    # Include only entries with >=1 CpG
    out <- ref[ref.ix.range]

    # Add the p-values
    mcols(out)$p <- ref.p
    out
  }









