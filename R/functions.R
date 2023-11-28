# General purpose functions
# Should be analysis agnostic


# Constants ---------------------------------------------------------------

BLOOD.CELL.TYPES <- c("Mon", "Mac0", "Mac1", "Mac2", "Neu", "MK", "EP", "Ery",
                      "FoeT", "nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8",
                      "nB", "tB")

get_sample_ids_from_dss_inputs <- function(file){
  e <- new.env()
  load(file, envir = e)
  return(rownames(e$design.df))
}


run_chisq_test_from_counts <- function(N.in.a, N.in.b, N.in.both, N.total, labs){
  #
  # Args
  # N.in.a : number of counts in group A
  cell.11 <- N.total - N.in.a - N.in.b + N.in.both
  cell.12 <-  N.in.a - N.in.both
  cell.21 <- N.in.b - N.in.both
  cell.22 <- N.in.both

  X <- matrix(c(cell.11, cell.12, cell.21, cell.22), nrow = 2)
  colnames(X) <- paste0(c("No", "Yes"), labs[1])
  rownames(X) <- paste0(c("No", "Yes"), labs[2])

  chisq.test(X)

  prop.table(X, margin = 1)
}


get_apoe_allele_frequencies <- function(file, sample.ids){
  df <- read_csv(file, show_col_types = F) %>%
    dplyr::filter(sample_id %in% sample.ids) %>%
    dplyr::mutate(APOE.risk.allele = (apoe_e1 == 4) + ( apoe_e2 == 4)) %>%
    dplyr::select(diagnostic_group, sample_id, APOE.risk.allele, apoe_e1, apoe_e2) %>%
    drop_na()
}


test_apoe_allele_frequencies <- function(apoe.df){
  tb <- apoe.df %>%
    group_by(diagnostic_group, APOE.risk.allele) %>%
    summarize(count = n()) %>%
    pivot_wider(names_from = APOE.risk.allele, values_from = count) %>%
    column_to_rownames("diagnostic_group")

  return(list("Table" = tb, "Test" = chisq.test(tb)))
}


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


get_group_ids <- function(master.df, vv){
  sort(unique(master.df$sample_id[master.df$diagnostic_group == vv]))
}


munge_master_df <- function(master.full.df, sample.ids){
  master.full.df %>%
    dplyr::filter(sample_id %in% sample.ids) %>%
    group_by(sample_id) %>%  # one sample shows up multiple times
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::mutate(APOE_risk_allele = (apoe_e1 == 4) + ( apoe_e2 == 4))
}


tabulate_and_test <- function(master.df, var, test.type = "chi"){
  tb <- master.df %>%
    dplyr::count(diagnostic_group, !!rlang::sym(var)) %>%
    pivot_wider(names_from = "diagnostic_group", values_from = "n", values_fill = 0) %>%
    drop_na() %>%
    column_to_rownames(var)

  if (test.type == "chi"){
    test.result <- chisq.test(tb)
  } else if (test.type == "fisher"){
    test.result <- fisher.test(tb)
  }

  return(list("Table" = tb,
              "Test" = test.result))
}


stringify_mean_sd <- function(xx){
  paste0(round(mean(xx), 2), " (", round(sd(xx), 2), ")")
}


test_continuous <- function(master.df, var){

  if (!(var %in% colnames(master.df))){
    warning("Var not in master.df")
  }

  xx.control <- master.df %>%
    dplyr::filter(diagnostic_group == "CONTROL") %>%
    dplyr::pull(var) %>% as.numeric()

  xx.load <- master.df %>%
    dplyr::filter(diagnostic_group == "LOAD") %>%
    dplyr::pull(var) %>% as.numeric()

  list("Control: " = stringify_mean_sd(xx.control),
       "LOAD: " = stringify_mean_sd(xx.load),
       "Test: " = t.test(xx.control, xx.load))
}


# Missingness quanitifcation ----------------------------------------------


read_and_tally_subroutine <- function(files){
  all <- do.call(rbind,
                 lapply(X = files,
                        FUN = fread,
                        select = c("chromStart", "coverage", "sample")))

  wide.df <- pivot_wider(all,
                         "chromStart",
                         names_from = "sample",
                         values_from = "coverage")

  out <- data.frame(chr = "chr1",
             start = wide.df$chromStart,
             end = wide.df$chromStart + 2,
             N.missing = rowSums(is.na(wide.df)))

  print("Read and pivoted one chunk...")
  return(out)
}

read_and_join_m_cov_routine <- function(dir, sample.ids, chunk.size){
  all.files <- file.path(dir, paste0("chr1.", sample.ids, ".bed"))

  all.files.split <- split(all.files, ceiling(seq_along(all.files) / chunk.size))

  counts.df <- do.call(rbind, lapply(X = all.files.split, read_and_tally_subroutine)) %>%
    group_by(chr, start, end) %>%
    summarize(PercentMissing = sum(N.missing) / length(sample.ids))

  return(counts.df)
}


# Statistical functions ---------------------------------------------------

compute_empirical_p <- function(xx, test.value){
  positive.test.value <- abs(test.value)

  p1 <- (xx >= test.value) / length(xx)
  p2 <- (xx <= test.value) / length(xx)

  p <- 2 * min(p1, p2)
  return(p)
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





# Get public data ---------------------------------------------------------

clean_differentially_expressed_genes <- function(file, autosomal.symbols){
  readxl::read_xlsx(file, skip=1) %>%
    dplyr::rename(gene.name = "gene name",
                  gene.type = "gene type",
                  entrez.id = "Entrez Gene ID",
                  gene.id = "Ensemble ID",
                  expression.pval = "P-value") %>%
    dplyr::filter(gene.name %in% autosomal.symbols)
}

# Essentially constants ---------------------------------------------------


get_hyper_hypo_colors <- function(rgb = F){
  # Hyper and hypo colors, may be used outside volcano
  # so we'll keep it in the R/functions.R file

  if (rgb){
    return(list(hyper = "239,192,0", hypo = "0,119,154"))
  } else {
    return(list(hyper = "#0073C2FF", hypo = "#EFC000FF"))
  }
}



# Manipulate GRanges ------------------------------------------------------

get_ensdb <- function(){
  # Make it super standard
  return(EnsDb.Hsapiens.v86)
}

get_autosomoal_gene_universe <- function(protein_coding=F){
  # if protein_coding = T, return only protein_coding genes

  edb <- get_ensdb()

  all.genes <- genes(edb)
  autosomes <- all.genes[seqnames(all.genes) %in% 1:22]

  if (protein_coding){
    autosomes$symbol[autosomes$gene_biotype == "protein_coding"]
  } else {
    autosomes$symbol
  }
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

get_gene_bodies <- function(upstream, downstream, autosomes_only = T, protein_coding_only = T){
  db <- get_ensdb()
  genes <- genes(db)

  # Respective indices
  protein.coding.ix <- genes$gene_biotype == "protein_coding"
  autosome.ix <- as.vector(seqnames(genes) %in% 1:22)

  if (autosomes_only & protein_coding_only){
    out <- expand_genes(genes[protein.coding.ix & autosome.ix], upstream, downstream)
  } else if (autosomes_only & !protein_coding_only){
    out <- expand_genes(genes[autosome.ix], upstream, downstream)
  } else if (!autosomes_only & protein_coding_only){
    out <- expand_genes(genes[protein.coding.ix], upstream, downstream)
  } else {
    out <- expand_genes(genes, upstream, downstream )
  }

  seqlevelsStyle(out) <- "NCBI"
  out
}

ids_to_symbols <- function(vv){
  genes.ref <- ensembldb::genes(get_ensdb())
  genes.ref[match(vv, genes.ref$gene_id )]$symbol
}

symbols_to_ids <- function(vv){
  genes.ref <- ensembldb::genes(get_ensdb())
  ix <- match(vv, genes.ref$symbol)
  genes.ref[ix[!is.na(ix)]]$gene_id
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

combine_missingness_with_pvals <- function(missing.chr1.load.df,
                                           missing.chr1.control.df,
                                           pvals.gr){

  colnames(missing.chr1.load.df)[4] <- "AD"
  colnames(missing.chr1.control.df)[4] <- "NoAD"

  miss.df <- full_join(missing.chr1.load.df,
                               missing.chr1.control.df,
                               by = c("chr", "start", "end"))

  miss.gr <- makeGRangesFromDataFrame(miss.df, keep.extra.columns = T)
  make_df_from_two_overlapping_granges(miss.gr, pvals.gr)
}

plot_missingness_hexbin <- function(miss.data, file){
  z <- miss.data %>%
    pivot_longer(contains("AD"), values_to = "Percent missing", names_to = "Group") %>%
    ggplot(aes(x = `Percent missing`, y = -log10(lfdr))) +
      facet_wrap(.~Group) +
      stat_binhex() +
      theme_minimal() +
      scale_fill_viridis(option = "B", trans = "log", breaks = 10^(1:5)) +
      ylab(expression(-log[10](lFDR))) +
      geom_hline(yintercept = -log10(0.05), color = "grey40") +
    ggtitle("Percent missing values before imputation") +
    theme(plot.background = element_rect(fill = "white", color = "white"))

  cowplot::save_plot(file, z)
  return(file)
}


# Harmonic P-value routine ------------------------------------------------

harmonic_pvalue_routine <- function(loci.gr, features.gr, alpha){
  df <- make_df_from_two_overlapping_granges(loci.gr, features.gr)

  hmp.df <- df %>%
    group_by(gene_name) %>%
    dplyr::summarize(
      N.CpGs = n(),
      N.DMPs = sum(lfdr < alpha),
      HarmonicMeanPval = harmonicmeanp::p.hmp(pval, L = N.CpGs)
    ) %>%
    dplyr::mutate(HarmonicMeanPval = pmin(HarmonicMeanPval, 1)) %>%
    dplyr::mutate(HarmonicMeanPval = pmax(HarmonicMeanPval, 0))

  # FDR calculation
  fdr.out <- fdrtool::fdrtool(x = hmp.df$HarmonicMeanPval, statistic = "pvalue", plot = F)

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



# Gene ontology pipeline --------------------------------------------------

run_gene_ontology <- function(genes){
  clusterProfiler::enrichGO(
    genes,
    OrgDb = "org.Hs.eg.db",
    keyType = "ENSEMBL",
    ont = "ALL",
    pAdjustMethod = "fdr"
  )
}

go_output_to_df <- function(go.out){
  # Constant to define mapping from abbreviated terms to full ones used in plots
  MAPPING <- data.frame(
    ONTOLOGY = c("MF", "CC", "BP"),
    Ontology = c("Molecular\nfunction", "Cellular\ncomponent", "Biological\nprocess")
  )

  go.df <- data.frame(go.out)

  # Get the names right
  dplyr::left_join(go.df, MAPPING, by = "ONTOLOGY")
}


plot_go_barchart <- function(go.df, n=25){

  # Get top 25 by p-value, then arrang by gene set size
  subdata <- head(dplyr::arrange(go.df, p.adjust), n) %>%
    dplyr::arrange(-p.adjust) %>%
    dplyr::mutate(Count = as.numeric(Count))

  subdata$Description <- factor(subdata$Description, levels = subdata$Description)

  ggplot(data = subdata,
         aes(x = Description,
             y = Count,
             fill = Ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_meAD() +
    xlab("") +
    ylab("Number of genes") +
    labs(legend = "") +
    theme(legend.position = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_fill_manual(values = c("#B24745FF", "#DF8F44FF", "#00A1D5FF"))
}




symbols_df_to_go_df_routine <- function(file){

  # Assuming we read from file
  df <- read_csv(file, show_col_types = F)
  gene.ids <- symbols_to_ids(df$gene_name)

  # Gene ontology and data munging
  go.out <- run_gene_ontology(gene.ids)
  go_output_to_df(go.out)

}


convert_gene_ontology_ids_to_symbols <- function(DMGenes.go.df){
  # Return the same data frame but with gene symbols instead of gene ontologies
  DMGenes.go.df %>%
    mutate(gene.ids = str_split(geneID, "/")) %>%
    unnest(gene.ids) %>%
    mutate(gene.symbols = ids_to_symbols(gene.ids)) %>%
    group_by(ONTOLOGY, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust) %>%
    summarize(GeneSymbols = paste(gene.symbols, collapse = ";"))
}


symbols_to_gene_ontology_routine <- function(symbols){
  symbols_to_ids(symbols) %>%
    run_gene_ontology %>%
    go_output_to_df %>%
    convert_gene_ontology_ids_to_symbols %>%
    tidyr::separate(GeneRatio, c("Count", "BGCount"), remove = F) %>%
    dplyr::rename(Ontology = "ONTOLOGY")
}


# Export UCSC Genome Browser ----------------------------------------------



liftover_wrapper <- function(gr, chain){
  seqlevelsStyle(gr) <- "UCSC"

  unlist(rtracklayer::liftOver(gr, chain)) %>%
    as.data.frame()
}

filter_by_gene_symbols <- function(genes, inclusion.genes){
  out <- genes[genes$symbol %in% inclusion.genes]
  out
}


tally_dmps_in_genes <- function(genes.with.dmp){
  as.data.frame(genes.with.dmp) %>%
    group_by(symbol) %>%
    summarize(N.DMPs = n())
}



get_array_genes_with_nearby_dmp <- function(dmps.array.gr){
  sort(unique(unlist(str_split(dmps.array.gr$gene.symbols, ";"))))
}

compare_with_madrid_paper <- function(dmps.array.gr, dm.genes.df){
  array.genes <- get_array_genes_with_nearby_dmp(dmps.array.gr)
  wgms.genes <- dm.genes.df$gene_name

  zz <- intersect(array.genes, wgms.genes)

  return(zz)
}


create_wb_with_description <- function(desc){
  cc <- 1:6

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, 1)
  openxlsx::writeData(wb, sheet = 1,
                      startCol = 1,
                      startRow = 1,
                      desc,
                      headerStyle = )


  openxlsx::addStyle(wb, sheet = 1, cols = cc, rows = 1, createStyle(textDecoration = "Bold"))

  return(wb)
}

append_data_to_workbook <- function(wb, data){



  setColWidths(wb, sheet = 1, cols = 1:ncol(data), widths = pmax(nchar(colnames(data)) * 1.2, 8) )
  mergeCells(wb, sheet=1, cols = 1:ncol(data), rows=1)


  openxlsx::writeData(wb = wb, sheet=1, x = data, startRow = 2,
                      headerStyle = createStyle(textDecoration = "Bold",
                                                halign = "center"))
  return(wb)
}

tally_and_write_DMGenes_with_DMP_counts <- function(
    gene.body.enrichment,
    lfdr.cut,
    desc,
    file){

  data <- gene.body.enrichment %>%
    dplyr::filter(lfdr < lfdr.cut)

  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)

  return(file)
}


run_length_embed_dmps <- function(dmps.data, chr, dist){
  tmp <- dmps.data %>%
    dplyr::filter(chr == chr) %>%
    arrange(start)

  rle(diff(tmp$start) < dist)
}


get_N_not_in_set <- function(dmps.gr, features.gr){
  seqlevelsStyle(dmps.gr) <- "NCBI"
  seqlevelsStyle(features.gr) <- "NCBI"

  out <- subsetByOverlaps(x = dmps.gr, ranges = features.gr, invert = TRUE)
  return(length(out))
}

get_N_in_set <- function(dmps.gr, features.gr){
  seqlevelsStyle(dmps.gr) <- "NCBI"
  seqlevelsStyle(features.gr) <- "NCBI"

  out <- subsetByOverlaps(x = dmps.gr, ranges = features.gr, invert = FALSE)
  return(length(out))
}


tally_dmps_in_out_genes <- function(dmps.gr, features.gr){
  return(list(
    "N DMPs  in set:" = get_N_in_set(dmps.gr, features.gr),
    "N DMPs NOT in set:" = get_N_not_in_set(dmps.gr, features.gr)
  ))
}

#END
