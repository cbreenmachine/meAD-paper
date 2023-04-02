


# Constants ---------------------------------------------------------------


GENIC.MAPPING <-
  data.frame(
    Annotation = c("1to5kb", "3UTRs", "5UTRs", "exons", "introns", "promoters", "Other"),
    Annotation.clean = c("1-5 kb", "3' UTR", "5' UTR", "Exon", "Intron", "Promoter", "Other")
  )


CPG.MAPPING <-
  data.frame(
    Annotation = c("hg38_cpg_inter", "hg38_cpg_islands", "hg38_cpg_shelves", "hg38_cpg_shores", "Other"),
    Annotation.clean = c("Open sea", "CpG island", "CpG shelf", "CpG shore", "Other")
  )

# Put factors in corret order, refelcting biology
GENIC.LEVELS <- c("1-5 kb", "Promoter", "5' UTR", "Exon", "Intron",  "3' UTR", "Other")
CPG.LEVELS <- c("Open sea", "CpG island", "CpG shelf", "CpG shore", "Other")



# Helper functions --------------------------------------------------------


annotate_loci_to_genic_parts <- function(data.gr, annotate.to = 'hg38_basicgenes'){

  if ("hg38_basicgenes" %in% annotate.to){
    MAPPING <- GENIC.MAPPING
    LEVELS <- GENIC.LEVELS
  } else {
    MAPPING <- CPG.MAPPING
    LEVELS <- CPG.LEVELS
  }

  print(MAPPING)

  # hg38_basicgenes or hg38_cpg
  my.anno <- annotatr::build_annotations(genome = 'hg38', annotations = annotate.to)

  # Run the annotation function
  data.annotated <- annotatr::annotate_regions(data.gr, my.anno, quiet = T)

  # Post-process after linking a CpG to a genomic feature
  anno.df <- data.frame(data.annotated) %>%
    dplyr::transmute(cpg.locus = paste0(seqnames, "-", start), annot.type) %>% # unique id for each CpG
    distinct() %>%  # get rid of instances where one CpG resides in multiple introns
    dplyr::mutate(Annotation = str_remove(annot.type, "hg38_genes_"))

  # Count how many annotate to none of these
  none.of.the.above <- subsetByOverlaps(data.gr,
                                        data.annotated,
                                        invert = T)

  # How many went to none?
  N.other <- length(none.of.the.above)

  # Add an extra row to anno.df; annotatr does not (as best as I know)
  # natively support counting "nones"
  other.df <- data.frame(cpg.locus = rep("none", N.other),
                         annot.type = rep("none", N.other),
                         Annotation = rep("Other", N.other))


  # Get the counts
  anno.counts.df <- rbind(anno.df, other.df) %>%
    group_by(Annotation) %>%
    dplyr::summarize(value = round(100 * n() / length(data.gr), 1))

  # Clean some strings and put in correct order
  anno.counts.cleaned.df <- anno.counts.df %>%
    left_join(MAPPING, by = "Annotation") %>%
    dplyr::mutate(Annotation.clean = factor(Annotation.clean, levels = LEVELS)) %>%
    arrange(Annotation.clean)

  # Add a string that looks like 'Promoter (10 %)'
  anno.counts.cleaned.df %>%
    dplyr::mutate(label = paste0(Annotation.clean,  " (", value, "%)")) %>%
    column_to_rownames("label")
}



cpg_annotation_routine <- function(data.gr){
  annotate_loci_to_genic_parts(
    data.gr,
    annotate.to = c("hg38_cpg_inter", "hg38_cpg_islands", "hg38_cpg_shelves", "hg38_cpg_shores")
  )
}

