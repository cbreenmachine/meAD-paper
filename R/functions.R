# General purpose functions
# Should be analysis agnostic


# Constants ---------------------------------------------------------------

BLOOD.CELL.TYPES <- c("Mon", "Mac0", "Mac1", "Mac2", "Neu", "MK", "EP", "Ery",
                      "FoeT", "nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8",
                      "nB", "tB")


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



# Statistical functions ---------------------------------------------------

compute_empirical_p <- function(xx, test.value){
  positive.test.value <- abs(test.value)

  mass <- sum(xx >= positive.test.value) + sum(xx <= -positive.test.value)
  mass / length(xx)
}


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

get_protein_coding_gene_bodies <- function(upstream, downstream){
  db <- get_ensdb()
  genes <- genes(db)

  expand_genes(genes[genes$gene_biotype == "protein_coding"], upstream, downstream)
}

ids_to_symbols <- function(vv){
  genes.ref <- ensembldb::genes(get_ensdb())
  genes.ref[match(vv, genes.ref$gene_id)]$symbol
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




# PCHi-C Analysis ---------------------------------------------------------


unnest_interactions <- function(interactions.gr){
  interactions.gr %>%
    dplyr::mutate(gene.name = str_split(baitName, ";")) %>%
    tidyr::unnest(gene.name)
}


combine_dmps_with_interactions <- function(dmps.gr, interactions.gr){
  # Overlap DMPs and interactions, and return a data frame.
  # There will be one row for each dmp by interaction (so some duplication if an)
  make_df_from_two_overlapping_granges(dmps.gr, interactions.gr) %>%
    group_by(interaction.id) %>%
    dplyr::mutate(
      median.pi.diff = median(pi.diff),
      k.dmps = n()) %>%
    ungroup()
}

summarize_interactions_with_dmp <- function(interactions.with.dmps) {
  interactions.with.dmps %>%
    dplyr::select(-c(chr, start, end, stat, pval, y, strand,
                     diagnostic_group_coded, pi.diff, lfdr)) %>%
    group_by(interaction.id) %>%
    dplyr::slice(1)
}


summarize_dm_enhancer_genes <- function(interactions.summary){
  interactions.summary %>%
    dplyr::mutate(gene.name = str_split(baitName, ";")) %>%
    unnest(gene.name) %>%
    distinct() %>%
    dplyr::select(gene.name, dist, bait.id, oe.id, median.pi.diff, k.dmps)
}


get_unique_genes_from_interactions <- function(interactions, protein_coding = F) {
  tmp <- unnest_interactions(interactions)

  out <- sort(unique(tmp$baitNameSplit))

  if (protein_coding) {

    genes <- genes(get_ensdb())
    filter.symbols <- genes$gene_name[genes$gene_biotype == "protein_coding"]

    out[out %in% filter.symbols]
  } else {
    out
  }
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


plot_go_barchart <- function(go.df, n=20){

  subdata <- head(arrange(go.df, -p.adjust), n)
  subdata$Description <- factor(subdata$Description, levels = subdata$Description)

  ggplot(data = subdata,
         aes(x = Description,
             y = -log10(p.adjust),
             fill = Ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal(base_size = 20) +
    xlab("") +
    ylab(expression(-log[10](p.adjusted))) +
    labs(legend = "") +
    theme(legend.position = "top",
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_fill_manual(values = c("red", "orange", "blue"))
}

symbols_to_go_plot <- function(symbols, file){

  ids <- symbols_to_ids(symbols)

  go.out <- run_gene_ontology(ids)
  go.df <- go_output_to_df(go.out)

  z <- plot_go_barchart(go.df)

  cowplot::save_plot(file, z, base_width = 12, base_height = 9)
}


# Export UCSC Genome Browser ----------------------------------------------

get_ucsc_seqlengths <- function(){
  good.chr <- paste0("chr", 1:22)

  db <- get_ensdb()
  seqlevelsStyle(db) <- "UCSC"

  # Usually pings to say "[some weird scaffold] didn't map"
  suppressWarnings(seqlengths(db)[good.chr])
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


liftover_wrapper <- function(gr, chain){
  seqlevelsStyle(gr) <- "UCSC"

  unlist(rtracklayer::liftOver(gr, chain)) %>%
    as.data.frame()
}


#END
