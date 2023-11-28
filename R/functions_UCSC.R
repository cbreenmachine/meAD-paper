# functions_UCSC.R
# Functions for exporintg to UCSC browser


ucsc.order.cols <- c("chrom", "chromStart", "chromEnd",
                     "name", "score", "value", "exp", "color",
                     "sourceChrom", "sourceStart", "sourceEnd", "sourceName", "sourceStrand",
                     "targetChrom", "targetStart", "targetEnd", "targetName", "targetStrand")




get_ucsc_seqlengths <- function(){
  # Need chromosome lengths for some GRanges functions
  db <- get_ensdb()
  lengths <- seqlengths(db)

  # Coinvert from 1 --> chr1
  names(lengths) <- paste0("chr", names(lengths))

  return(lengths[paste0("chr", 1:22)])
}



#TODO: refactor
format_and_write_ucsc_lolly <- function(data.gr, file, lfdr.cut, keep.nth=25){
  # data.gr is pvals.gr or something similar
  # alpha.cut is the level of significance
  # keep.nth is how much to thin by (keep every 25th non-significant point)

  colors <- get_hyper_hypo_colors(rgb = T)

  # Thin the points
  keepix <- c(
    which(data.gr$lfdr <= lfdr.cut),
    which(data.gr$lfdr > lfdr.cut)[c( rep(FALSE, keep.nth-1), TRUE)]
  )


  # Thin data
  out <- data.gr[keepix, ]

  # Manipulate the metadata columns to get in UCSC lolly format
  mdata.cleaned <- mcols(out) %>%
    as.data.frame() %>%
    dplyr::transmute(name = ".",
                     score = round(y),
                     thickStart = start(out),
                     thickEnd = end(out),
                     color = ifelse(pi.diff < 0, colors$hypo, colors$hyper),
                     lollySize = ifelse(lfdr < lfdr.cut, 4, 1))

  mdata.cleaned$color[mdata.cleaned$lollySize == 1] <- "220,220,220"

  mcols(out) <- mdata.cleaned
  genome(out) <- "hg38"
  seqlengths(out) <- get_ucsc_seqlengths()

  # Write and return
  rtracklayer::export.bb(out, file)
  return(file)
}




export_significant_interactions_to_UCSC <- function(interactions.with.dmp){

  interactions.with.dmp %>%
    separate(oe.id, into = c("oeChr", "oeStart", "oeEnd")) %>%
    dplyr::transmute(sourceChrom = oeChr,
                     sourceStart = as.numeric(oeStart),
                     sourceEnd = as.numeric(oeEnd),
                     targetChrom = paste0("chr", baitChr),
                     targetStart = as.numeric(baitStart),
                     targetEnd = as.numeric(baitEnd)) %>%
    dplyr::mutate(chrom = sourceChrom,
                  chromStart = floor((sourceStart + targetStart) / 2),
                  chromEnd = floor((sourceEnd + targetEnd) / 2)) %>%
    dplyr::mutate(value = 5, exp = ".", color = "255,0,0",
                  name = ".", score = 500,
                  sourceName = "OtherEnd:enhancer", sourceStrand = ".",
                  targetName = "Bait:promoter", targetStrand = ".") %>%
    dplyr::select(all_of(ucsc.order.cols)) %>%
    dplyr::arrange(chrom, chromStart)
}


format_and_write_ucsc_interactions <- function(interactions.for.ucsc, file){

  out <-
    interactions.for.ucsc %>%
    drop_na() %>% # need distance to be non-missing to visualize...
    makeGRangesFromDataFrame(seqnames.field = "chrom",
                             start.field = "chromStart",
                             end.field = "chromEnd",
                             keep.extra.columns = T,
                             ignore.strand = T)

  # Extra processing
  genome(out) <- "hg38"

  lengths <- get_ucsc_seqlengths()

  # Allows us to get the indices to "duplicate" the
  # lengths the appropriate number of times. E.g., chr1 is 1mb,
  # then if there are five itneractions, we have 1mb, 1mb, 1mb, 1mb, 1mb
  ix <- as.vector(match(seqnames(out), names(lengths)))

  # Which interactions are in correct range
  keep.ix <- ((end(out) < lengths[ix]) & (start(out) < lengths[ix])) &
    ((out$sourceEnd < lengths[ix]) & (out$sourceStart < lengths[ix]))

  # Subset
  out.trimmed <- out[keep.ix, ]
  seqlengths(out.trimmed) <- lengths

  seqlevelsStyle(out.trimmed) <- "UCSC"

  # Write write write
  rtracklayer::export.bb(out.trimmed, file)
  file
}



format_subset_of_interactions <- function(interactions, filename){
  # Just the area around b4galt1
  interactions %>%
    dplyr::filter(chrom == "chr9", chromStart > 32e6, chromEnd < 34e6) %>%
    write_tsv(filename)

  return(filename)
}

