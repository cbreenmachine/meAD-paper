#

read_GODMC_meQTLs <- function(ifile){
  fread(ifile, select = c("cpg","snp","beta_a1","pval","cistrans","pval_are","pval_mre")) %>%
    dplyr::group_by(cpg) %>%
    dplyr::slice(1)
}





filter_450k_annotation <- function(valid.cpgs){

  # Load the manifest file and get it in the right format
  anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

  # Filter and add start/end
  anno.sig <- anno %>%
    dplyr::filter(Name %in% valid.cpgs) %>%
    dplyr::mutate(start = pos-1, end = pos+1)

  anno.sig.gr <- makeGRangesFromDataFrame(anno.sig, starts.in.df.are.0based = T)
  return(anno.sig.gr)
}


make_godmc_dmp_table <- function(
    dmp.and.mqtl,
    array.450.dmps,
    godmc.hg38.gr,
    array.450.hg38
    ){

  # Assume all 450k sites are included in universe
  N.universe <- length(array.450.hg38)

  x <- matrix(c(-1, -1, -1, -1), nrow = 2, ncol = 2)
  rownames(x) <- c("godmc_yes", "godmc_no")
  colnames(x) <- c("dmp_yes", "dmp_no")

  # Yes to both--DMP and meQTL
  x[1, 1] <- length(dmp.and.mqtl)

  # Yes to DMP, no to meQTL
  x[1, 2] <- length(array.450.dmps) - x[1, 1]

  # No to DMP, yes to meQTL
  x[2, 1] <- length(godmc.hg38.gr) - x[1, 1]

  # Inclusion-exclusion
  x[2, 2] <- length(array.450.hg38) - x[1, 2] - x[2, 1] + x[1, 1]

  tb <- as.table(x)
  tb
}

test_significance_of_meQTL_overlap <- function(N.meqtl.dmp, n.dmp, n.meqtl){
  # n.meqtl.dmp is the number of sites that are both DMPs and

  # See https://www.sciencedirect.com/science/article/pii/S0888754311001807?via%3Dihub
  # Bibikova 2014
  N.array <- 482421
  N.meqtl.dmp <- 10

}

compute_average_correlation <- function(df){
  (df$parameter %*% df$estimate) / sum(df$parameter)
}


count_n_dmps_multi_map_to_gene <- function(dmps, genes){
  seqlevelsStyle(dmps) <- "NCBI"
  seqlevelsStyle(genes) <- "NCBI"

  ov <- findOverlaps(dmps, genes)

  # Query is dmps in this case
  #check with max(queryHits(ov))
  data.frame(dmp.ix = queryHits(ov)) %>%
    dplyr::group_by(dmp.ix) %>%
    dplyr::summarize(N = n()) %>%
    return()
}
