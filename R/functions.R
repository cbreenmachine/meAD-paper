# General purpose functions
# Should be analysis agnostic


# General reading/casting functions -----------------------------------------------

get_pvals_data <- function(ifile){
  # READ pvals from
  data.table::fread(ifile, verbose = F) %>%
    dplyr::mutate(lfdr = lfdr.from.ss,
                  pval = p.from.ss,
                  y = -log10(lfdr)) %>%
    dplyr::select(-c(p.from.ss, p.from.zz, p.from.DSS,
                     lfdr.from.ss, lfdr.from.zz))
}


subset_dmps <- function(data, lfdr.cut=0.05){
  data %>%
    dplyr::filter(lfdr <= lfdr.cut)
}


to_granges <- function(data){
  out <- data %>%
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns=T,
      starts.in.df.are.0based = T
    )

  seqlevelsStyle(out) <- "UCSC"
  out
}


# Essentially constants ---------------------------------------------------


get_hyper_hypo_colors <- function(){
  # Hyper and hypo colors, may be used outside volcano
  # so we'll keep it in the R/functions.R file
  return(list(hyper = "#0073C2FF", hypo = "#EFC000FF"))
}

