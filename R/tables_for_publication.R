# tables_for_publication.R
# We include 5-6 tables in manuscript



process_and_write_dmps <- function(dmps.gr, desc, file){

  # Columns look like
  # Chrom, Start, End, AD - No AD, lFDR
  # chr1, 10, 12, 0.01, 0.02
  # ...

  data <- dmps.gr %>%
    as.data.frame() %>%
    dplyr::transmute(
      Chromosome = seqnames,
      "Start Coordinate (hg38)" = start,
      `End Coordinate (hg38)` = end,
      `Estimated Methylation Difference (AD minus no-AD)` = pi.diff,
      `Local False-Discovery Rate (lFDR)` = lfdr)

  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)

  return(file)
}


process_and_write_DM_genes <- function(gene.body.enrichment, desc, file){

  # Columns look like
  # Gene symbol, # CpGs, # DMPs, Harmonic Mean P-value, lFDR
  # BRCA1, 20, 1, 0.7, 0.99

  data <- gene.body.enrichment %>%
    transmute(
      `Gene Symbol` = gene_name,
      `Total Number of CpGs` = N.CpGs,
      `Number of Differentially Methylated Positions (DMPs)` = N.DMPs,
      `Harmonic Mean p-value` = HarmonicMeanPval,
      `Local False-Discovery Rate (lFDR)` = lfdr,
    )

  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)

  return(file)
}



process_and_write_gene_ontology_terms <- function(dm.genes.go.df, desc, file){

  # First convert ENSEMBL IDs to genes...
  data <- dm.genes.go.df %>%
    transmute(`Gene Ontology (GO) Domain` = Ontology,
              `GO Term ID` = ID,
              `GO Description` = Description,
              `Gene Symbols` = GeneSymbols,
              `Gene Ratio` = GeneRatio,
              `BG Ratio` = BgRatio,
              `Local False-Discovery Rate (lFDR)` = p.adjust)

  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)

  return(file)
}


format_and_write_dm_enhancers <- function(interactions.summary, de.genes, desc, file){

  data <- unnest_interactions_by_gene(interactions.summary) %>%
    dplyr::arrange(desc(k.dmps), gene.name) %>%
    transmute(
      `Number of Differentially Methylated Positions (DMPs)` = k.dmps,
      `Gene Symbol` = gene.name,
      `Promoter Locus (hg38)` = paste0(baitChr, ":", baitStart, "-", baitEnd),
      `Enhancer Locus (hg38)` = paste0("chr", oe.id),
      `Differentially Expressed` = ifelse(`Gene Symbol` %in% de.genes, "Differentially expressed", "Not differentially expressed")) %>%
    distinct() %>%
    arrange()


  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)
  return(file)
}


format_and_write_dm_promoters <- function(promoters.summary, de.genes, desc, file){

  data <- unnest_interactions_by_gene(promoters.summary) %>%
    dplyr::arrange(desc(k.dmps), gene.name) %>%
    transmute(
      `Number of Differentially Methylated Positions (DMPs)` = k.dmps,
      `Gene Symbol` = gene.name,
      `Promoter Locus (hg38)` = bait.id,
      `Enhancer Locus (hg38)` = paste0("chr", oeChr, ":", oeStart, "-", oeEnd),
      `Differentially Expressed` = ifelse(`Gene Symbol` %in% de.genes, "Differentially expressed", "Not differentially expressed")) %>%
    distinct()


  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)
  return(file)
}



process_and_write_nature_genetics_list <- function(NG.tally.df,
                                                   natgen.genes,
                                                   genes.w.dm.enhancers,
                                                   genes.w.dm.promoters, desc, file){


  # Need to expand current dataset
  genes.zero.dmps <- setdiff(natgen.genes, NG.tally.df$symbol)
  tmp <- data.frame(genes.zero.dmps, 0)
  colnames(tmp) <- colnames(NG.tally.df)

  # COmbine non-zero with zero counts
  data <- rbind(NG.tally.df, tmp) %>%
    dplyr::filter(!str_detect(symbol, "IGH")) %>%
    arrange(-N.DMPs) %>%
    dplyr::mutate(
      dm.enhancer = ifelse(symbol %in% genes.w.dm.enhancers, "1+ DMPs in 1+ enhancers", "None"),
      dm.promoter = ifelse(symbol %in% genes.w.dm.promoters, "1+ DMPs in 1+ promoters", "None"))
  head(data)

  colnames(data) <- c("AD Genomic Risk Locus Symbol",
                      "Number Of Differentially Methylated Positions (DMPs) Within 25kb",
                      "Differential Methylation in Enhancer",
                      "Differential Methylation in Promoter")

  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)

  return(file)

}

