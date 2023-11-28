

# October 26, 2023 - DMPs with Gene Annotation ----------------------------

add_gene_annotations_to_dmps <- function(dmps.gr){
  genes <- genes(EnsDb.Hsapiens.v86)
  seqlevelsStyle(genes) <- "NCBI"
  seqlevelsStyle(dmps.gr) <- "NCBI"

  overlaps <- findOverlaps(dmps.gr, genes)

  dmps <- as.data.frame(dmps.gr)[queryHits(overlaps), ]
  genes <- as.data.frame(genes)[subjectHits(overlaps), ]

  return(cbind(dmps, genes))
}
