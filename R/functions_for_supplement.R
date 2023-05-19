
compute_lambda_gc <- function(pp){
  round(median(qchisq(1 - pp, 1)) / qchisq(0.5, 1), 1)
}


make_title_for_pvals_plots <- function(desc, lambda.gc){
  bquote(.(desc)~~"("*lambda[GC]==.(lambda.gc)*")")
}

plot_hist_from_pvals <- function(pp){

  lambda.gc <- compute_lambda_gc(pp)
  ti <- make_title_for_pvals_plots("Modeled p-values", lambda.gc)

  hist(pp,
       probability = T,
       breaks = 50,
       main = as.expression(ti),
       xlab = expression(p))


}


plot_qq_from_pvals <- function(pp){

  # Genomic in/de-flation
  lambda.gc <- compute_lambda_gc(pp)

  # Title
  ti <- make_title_for_pvals_plots("Raw p-values", lambda.gc)

  fastqq::qq(pp, main = as.expression(ti))
}


plot_pvals_routine <- function(pvals.data){
  pp.1 <- pvals.data$pval.Wald
  pp.2 <- pvals.data$pval

  # Run plotting (both coalculate lambda internally)

  par(mfrow = c(1,2))
  plot_qq_from_pvals(pp.1)
  plot_hist_from_pvals(pp.2)

}

run_and_save_pvals_plot <- function(pvals.data, file){
  pdf(file, width = 9, height = 5)
  plot_pvals_routine(pvals.data )
  dev.off()

  # Need to return a charcter vector
  file
}



# Coverage histograms -----------------------------------------------------

munge_sequencing_depth <- function(file, sample.ids){
  read_csv(file,
           show_col_types = F,
           col_names = c("sample_id", "Mean sequencing depth")) %>%
    mutate(sample_id = as.character(sample_id)) %>%
    dplyr::filter(sample_id %in% sample.ids)
}


munge_effective_coverage <- function(file, sample.ids){

  df <- read_csv(file, show_col_types = F)

  # We'll need weighted averages based on number of positions per subject
  df$weight <- df$N.positions / sum(df$N.positions)

  # True sequencing coverage
  cov <- df %>%
    dplyr::select(-N.positions) %>%
    pivot_longer(-weight) %>%
    group_by(name) %>%
    summarize(Cov.bar = as.numeric(weight %*% value))

  names(cov) <-c("sample_id", "Mean effective coverage")
  cov %>%
    dplyr::filter(sample_id %in% sample.ids)
}



plot_coverage_hist_routine <- function(depth.df, cov.df, ofile){
  df <- full_join(cov.df, depth.df, by = "sample_id") %>%
    pivot_longer(-sample_id) %>%
    mutate(name = factor(name, levels = c("Mean sequencing depth", "Mean effective coverage")))

  # Plot this
  p <- df %>%
    ggplot(aes(x = value, y = after_stat(density))) +
    geom_histogram(bins = 20, color = "black", fill = "grey90") +
    theme_minimal() +
    facet_wrap(.~name, scales = "free") +
    ylab("Density") +
    xlab("") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill = "white", color = "white"))


  cowplot::save_plot(ofile, p)
  ofile

}




# QC Stuff from JSONS -----------------------------------------------------

get_stats_from_json <- function(file){
  z <- fromJSON(file = file)

  sample_id = tools::file_path_sans_ext(basename(file))

  # Get the four values (total (L, R), unmapped (L, R))
  vv <- as.numeric(unlist(z$Reads))

  if(length(z$Reads) > 0){
    data.frame(sample_id = sample_id,
               TotalForward = vv[1],
               TotalReverse = vv[2],
               UnmappedForward = vv[3],
               UnmappedReverse = vv[4])
  } else {
    NA
  }
}

get_all_stats_from_dir <- function(dir){
  all.files <- list.files(dir, full.names = T, pattern = "*json")

  result <- lapply(X=all.files, FUN=get_stats_from_json)
  result <- result[!is.na(result)]

  do.call(rbind, result) %>%
    mutate(TotalReads = TotalForward + TotalReverse,
           TotalMappedReads = TotalReads - UnmappedForward - UnmappedReverse,
           MappingPercentage = TotalMappedReads / TotalReads)
}




# Demographics table ------------------------------------------------------




# Write supplementary tables ----------------------------------------------


process_and_write_dmps <- function(dmps.gr, desc, file){

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

#
process_and_write_gene_ontology_terms <- function(DMGenes.go.df, desc, file){

  # First convert ENSEMBL IDs to genes...
  data <- convert_gene_ontology_ids_to_symbols(DMGenes.go.df) %>%
    transmute(`Gene Ontology (GO) Domain` = ONTOLOGY,
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


format_and_write_dm_enhancers <- function(interactions.summary, desc, file){

  data <- unnest_interactions_by_gene(interactions.summary) %>%
    dplyr::arrange(desc(k.dmps), gene.name) %>%
    transmute(
      `Number of Differentially Methylated Positions (DMPs)` = k.dmps,
      `Gene Symbol` = gene.name,
      `Promoter Locus (hg38)` = paste0(baitChr, ":", baitStart, "-", baitEnd),
      `Enhancer Locus (hg38)` = paste0("chr", oe.id)) %>%
    distinct() %>%
    arrange()


  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)
  return(file)
}


format_and_write_dm_promoters <- function(promoters.summary, desc, file){

  data <- unnest_interactions_by_gene(promoters.summary) %>%
    dplyr::arrange(desc(k.dmps), gene.name) %>%
    transmute(
      `Number of Differentially Methylated Positions (DMPs)` = k.dmps,
      `Gene Symbol` = gene.name,
      `Promoter Locus (hg38)` = bait.id,
      `Enhancer Locus (hg38)` = paste0("chr", oeChr, ":", oeStart, "-", oeEnd)) %>%
    distinct()


  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)
  return(file)
}



process_and_write_nature_genetics_list <- function(NG.tally.df, natgen.genes, desc,  file){


  # Need to expand current dataset
  genes.zero.dmps <- setdiff(natgen.genes, NG.tally.df$symbol)
  tmp <- data.frame(genes.zero.dmps, 0)
  colnames(tmp) <- colnames(NG.tally.df)

  # COmbine non-zero with zero counts
  data <- rbind(NG.tally.df, tmp) %>%
    dplyr::filter(!str_detect(symbol, "IGH")) %>%
    arrange(-N.DMPs)

  colnames(data) <- c("AD Genomic Risk Locus Symbol", "Number Of Differentially Methylated Positions (DMPs) Within 25kb")

  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)

  return(file)

}




# Array comparisons -------------------------------------------------------

read_one_raw_M_Cov <- function(file){
  out <- fread(file) %>%
    dplyr::filter(coverage > 5) %>%
    dplyr::transmute(chrom,
                     start = chromStart,
                     P = methylated / coverage,
                     sample)
  print(file)
  return(out)
}

join_raw_M_Cov_files <- function(dir, common.ids){
  all.files <- file.path(dir, paste0("chr1.", common.ids, ".bed"))
  data.table::rbindlist(lapply(X = all.files, FUN = read_one_raw_M_Cov)) %>%
    pivot_wider(id_cols = start, names_from = sample, values_from = P) %>%
    dplyr::mutate(chr = "chr1", end = start + 2)
}

chunker_wrapper <- function(file){
  chunker(file,
          sep = "\t",
          chunksize = 500000,
          has_colnames = T,
          has_rownames = F)
}


make_P_from_M_Cov <- function(M, Cov, chr){
  P <- M[, -1] / Cov[ ,-1]

  # Make this explicit
  P[Cov[ ,-1] == 0] <- NA

  P$chr <- chr
  P$start <- M$chromStart
  P$end <- M$chromStart + 2

  makeGRangesFromDataFrame(P,
                           keep.extra.columns = T,
                           starts.in.df.are.0based = T)
}

read_and_merge_M_Cov_files <- function(file.prefix){
  chr <- basename(file.prefix)

  M.file <- paste0(file.prefix, ".M.bed")
  Cov.file <- paste0(file.prefix, ".Cov.bed")

  M_chunker <- chunker_wrapper(M.file)
  Cov_chunker <- chunker_wrapper(Cov.file)


  P.all.grs <- list()
  ix <- 1

  while(next_chunk(M_chunker) & next_chunk(Cov_chunker)){

    M <- get_table(M_chunker)
    Cov <- get_table(Cov_chunker)

    P.all.grs[[ix]] <- make_P_from_M_Cov(M, Cov, chr)
    ix <- ix + 1
  }

  return(do.call(c, P.all.grs))

}


read_850k_array <- function(file){
  array <- fread(file)

  array.gr <- makeGRangesFromDataFrame(array,
                                       keep.extra.columns = T,
                                       starts.in.df.are.0based = F)

  array.gr
}


combine_850k_array_with_wgms <- function(array.gr, P.gr){

  # Get common columns
  cc <- intersect(names(mcols(P.gr)), names(mcols(array.gr)))

  df.1 <- subsetByOverlaps(P.gr, array.gr)[ ,cc] %>%
    as.data.frame() %>%
    dplyr::mutate(tech = "WGMS")
  df.2 <- subsetByOverlaps(array.gr, P.gr)[ ,cc] %>%
    as.data.frame() %>%
    dplyr::mutate(tech = "EPIC")

  out.df <- rbind(df.1, df.2) %>%
    dplyr::select(-c(end, width, strand)) %>%
    pivot_longer(starts_with("X"), names_to = "sample", values_to = "methylation", ) %>%
    pivot_wider(values_from = "methylation", names_from = "tech", values_fn = list) %>%
    unnest(cols = c(WGMS, EPIC)) %>%
    distinct() %>%
    dplyr::mutate(sample = as.numeric(str_remove(sample, "X")))

  return(out.df)
}


merge_array_and_wgms_subroutine <- function(array.gr, file.prefix){
  P.gr <- read_and_merge_M_Cov_files(file.prefix)

  # We sample because otherwise it's impossible to
  # compute a correlation on 850k * 84 entries
  combine_850k_array_with_wgms(array.gr, P.gr)
}


merge_array_and_wgms_routine <- function(array.gr, dir){
  all.file.prefixes <- file.path(dir, paste0("chr", 1:22))

  set.seed(919)

  out <- lapply(X = all.file.prefixes,
                FUN = merge_array_and_wgms_subroutine,
                array.gr = array.gr) %>% do.call(rbind)

}

compute_wgms_vs_array_summary_stats <- function(chr_merged){

  out <- data.frame(chr = chr_merged$seqnames[1],
                    rho = cor(chr_merged$WGMS, chr_merged$EPIC, use = "pairwise"),
                    N.loci = length(unique(chr_merged$start)))
  return(out)
}


plot_wgms_vs_array_hexbin <- function(wgms_vs_array.df){
  wgms_vs_array.df %>%
    drop_na() %>%
    sample_n(100000) %>%
    ggplot(aes(x = EPIC, y = WGMS)) +
    geom_hex(binwidth = c(0.01, 0.01), aes(fill = log(..density..))) +
    scale_fill_viridis(option = "plasma") +
    theme_classic()
}

