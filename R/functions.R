# General purpose functions
# Should be analysis agnostic


# General reading/casting functions -----------------------------------------------

get_natgen_genes <- function(ifile){
  read_table(ifile, col_names = FALSE)$X1
}


my_write_csv <- function(data, file){
  write_csv(data, file = file)
  return(file)
}


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



# Manipulate GRanges ------------------------------------------------------

get_ensdb <- function(){
  # Make it super standard
  return(EnsDb.Hsapiens.v86)
}

expand_genes <- function(gr, upstream, downstream) {
  # https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
  strand_is_minus = strand(gr) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  start(gr)[on_plus] = start(gr)[on_plus] - upstream
  start(gr)[on_minus] = start(gr)[on_minus] - downstream
  end(gr)[on_plus] = end(gr)[on_plus] + downstream
  end(gr)[on_minus] = end(gr)[on_minus] + upstream
  return(gr)
}

get_protein_coding_gene_bodies <- function(upstream, downstream){
  db <- get_ensdb()
  genes <- genes(db)

  expand_genes(genes[genes$gene_biotype == "protein_coding"], upstream, downstream)
}

ids_to_symbols <- function(vv){
  genes.ref <- genes(get_ensdb())
  genes.ref[match(vv, genes.ref$gene_id)]$symbol
}

get_protein_coding_promoters <- function(upstream, downstream){
  db <- get_ensdb()
  promoters <- promoters(db, upstream, downstream)
  promoters <- promoters[promoters$tx_biotype == "protein_coding"]
  promoters$gene_name <- ids_to_symbols(promoters$gene_id)
  promoters
}


make_df_from_overlapping_granges <- function(range.1, range.2){
  # Common seqlevels
  seqlevelsStyle(range.1) <- "NCBI"
  seqlevelsStyle(range.2) <- "NCBI"

  # Overlap and get indices
  overlaps <- findOverlaps(range.1, range.2)
  ix.1 <- queryHits(overlaps)
  ix.2 <- subjectHits(overlaps)

  # We'll keep the seqnames for range.1, not range.2 (otherwise duplicate column names)
  left <- data.frame(range.1[ix.1, ]) %>%
    dplyr::select(-c(width, strand)) %>%
    dplyr::rename(chr = seqnames)

  right <- data.frame(range.2)[ix.2, -1:-3] %>%
    dplyr::select(-width)

  cbind(left, right)
}


# Harmonic P-value routine ------------------------------------------------

harmonic_pvalue_routine <- function(loci.gr, features.gr){
  df <- make_df_from_overlapping_granges(loci.gr, features.gr)

  hmp.df <- df %>%
    group_by(gene_name) %>%
    dplyr::summarize(
      k = n(),
      hm.pval = harmonicmeanp::p.hmp(pval, L = k)
    ) %>%
    dplyr::mutate(hm.pval = pmin(hm.pval, 1)) %>%
    dplyr::mutate(hm.pval = pmax(hm.pval, 0))

  # FDR calculation
  fdr.out <- fdrtool::fdrtool(x = hmp.df$hm.pval, statistic = "pvalue", plot = F)

  # Pack together
  hmp.df$lfdr <- fdr.out$lfdr
  hmp.df$qval <- fdr.out$qval

  # Output
  hmp.df
}



# -------------------------------------------------------------------------

filter_by_natgen<- function(data, genes){
  data %>%
    dplyr::filter(gene_name %in% genes) %>%
    arrange(lfdr)
}



# PCHi-C Data Wrangling ---------------------------------------------------

clean_and_filter_interactions <- function(file){

}
