library(EnsDb.Hsapiens.v86)

#TODO: fix out of range warning
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

get_raw_ensdb <- function(){
  # Ensembl database and filtered for coding genes only
  ensdb <- EnsDb.Hsapiens.v86
  seqlevelsStyle(ensdb) <- "UCSC"
  ensdb
}

get_expanded_ensdb <- function(up, down){
  #--> Ensembl database and filtered for coding genes only
  ensdb <- get_raw_ensdb()

  # Filter and expand ensdb
  ensdb.subset <- genes(ensdb, filter = GeneBiotypeFilter('protein_coding'))
  ensdb.expanded <- expand_genes(ensdb.subset, upstream=up, downstream=down)
  ensdb.expanded
}
