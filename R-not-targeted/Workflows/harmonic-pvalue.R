
# call-DM-associated-genes.R
#

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(parallel)
  library(harmonicmeanp)
  library(tidyverse)
})


fishers_method <- function(p.vec){
  # Runs fisher's method to aggregate p-values
  # It is well-known that Fisher's method is poorly calibrated when
  # tests are dependent (which is the case in methylation). So Fisher's method
  # provides a "poor baseline" to compare against

  # Args
  # p.vec is a vector of p-values

  # Returns a p-value
  k <- length(p.vec) # degrees of freedom
  chi.stat <- -2 * sum(log(p.vec))
  pchisq(chi.stat, df = 2*k)
}

dummy_function <- function(data){
  output_function <- function(){
    print(data)
  }
  return(output_function)
}

generate_wrapper <- function(loci.gr, features.gr){

  # Do some checks
  stopifnot("pval" %in% names(elementMetadata(loci.gr)))

  # Subset and get number of CpGs
  output_function <- function(i){
    cpgs.in.gene <- subsetByOverlaps(loci.gr, features.gr[i])
    k <- length(cpgs.in.gene)

    if (length(cpgs.in.gene) > 0){
      # Get vector of p-values
      pp <- cpgs.in.gene$pval
      return(c(p.hmp(pp, L = length(pp)), k = k))
    } else {
      return(c(p.hmp = -1, k = 0))
    }
  }

  return(output_function)
}

#END
