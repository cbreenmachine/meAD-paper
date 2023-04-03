# General purpose functions
# Should be analysis agnostic


download_chain_from_ucsc <- function(dir, file){
  # Assign output name
  ofile.name <- file.path(dir, basename(file))

  # Download
  download.file(file, destfile = ofile.name)

  # Decompress
  R.utils::gunzip(ofile.name, remove = F, overwrite = T)

  # Return a chain object
  rtracklayer::import.chain(stringr::str_remove(ofile.name, "\\.gz"))
}


to_ucsc_format <- function(a, b, c){
  # Return chr1:1111-2222
  if (!str_detect(a, "chr")){a <- paste0("chr", a)}
  paste0(a, ":", b, "-", c)
}

to_ucsc_format_v <- Vectorize(to_ucsc_format)


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


make_df_from_two_overlapping_granges <- function(range.1, range.2){
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
  df <- make_df_from_two_overlapping_granges(loci.gr, features.gr)

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
  # Read as text file
  promoter_capture.data <- fread(file)

  # Filter out the following:
  # 1) oeName == "." indicates a bait-to-bait interaction
  # 2) baitName == NA indicates unknown promoter
  # 3) dist == NA indicates trans-chromosomal interaction
  # 4) Sex chromosomes (don't include for now)
  promoter_capture.filtered <- promoter_capture.data %>%
    dplyr::filter(oeName == ".",
                  !is.na(baitName),
                  !is.na(dist),
                  !(baitChr %in% c("X", "Y"))) %>%
    # Add identifiers for each bait and other end
    dplyr::mutate(bait.id = to_ucsc_format_v(baitChr, baitStart, baitEnd),
                  oe.id = to_ucsc_format_v(oeChr, oeStart, oeEnd)) %>%
    # Drop unneeded columns
    dplyr::select(-c(clusterID, clusterPostProb, oeName, oeID))

    # Get some summary statistics
    promoter_capture.filtered$min.chicago <- apply(X=promoter_capture.filtered[, 10:26], MARGIN=1, FUN=min)
    promoter_capture.filtered$med.chicago <- apply(X=promoter_capture.filtered[, 10:26], MARGIN=1, FUN=median)
    promoter_capture.filtered$mean.chicago <- apply(X=promoter_capture.filtered[, 10:26], MARGIN=1, FUN=mean)
    promoter_capture.filtered$max.chicago <- apply(X=promoter_capture.filtered[, 10:26], MARGIN=1, FUN=max)

    # Return this big fella
    promoter_capture.filtered
}

make_granges_with_common_field_prefix <- function(gr, prefix=""){
  makeGRangesFromDataFrame(gr,
                           seqnames.field =  paste0(prefix, "Chr"),
                           start.field = paste0(prefix, "Start"),
                           end.field = paste0(prefix, "End"),
                           keep.extra.columns = T)
}

liftover_wrapper <- function(gr, chain){
  seqlevelsStyle(gr) <- "UCSC"

  unlist(rtracklayer::liftOver(gr, chain)) %>%
    as.data.frame()
}

drop_granges_columns <- function(data){
  dplyr::select(data, -c(seqnames, start, end, width, strand))
}


lift_promoter_capture_data_to_hg38 <- function(interactions.hg19, chain, return.granges=F){
  # Need an id to re-combine later
  interactions.hg19$interaction.id <- 1:nrow(interactions.hg19)

  #TODO: how to reconcile differences here??
  baits.hg19 <- make_granges_with_common_field_prefix(interactions.hg19, prefix="bait")

  # Convert to UCSC style and liftOver
  baits.hg38 <-liftover_wrapper(baits.hg19, chain) %>%
    drop_granges_columns() %>%
    distinct()

  # Other ends
  other_ends.hg19 <- make_granges_with_common_field_prefix(interactions.hg19, prefix="oe")

  # Convert to UCSC style and liftOver
  other_ends.hg38 <- liftover_wrapper(other_ends.hg19, chain) %>%
    drop_granges_columns() %>%
    dplyr::select(c("interaction.id", "baitChr", "baitStart", "baitEnd")) %>%
    distinct()

  # Merge
  output <- dplyr::inner_join(baits.hg38, other_ends.hg38, by = "interaction.id")

  if (nrow(output) == length(unique(output$interaction.id))){
    if (return.granges){
      return(make_granges_with_common_field_prefix(output, prefix = "oe"))
    } else {
      return(output)
    }
  } else {
    warning("N row of interactions is not the same as unique interactions")
  }
}
