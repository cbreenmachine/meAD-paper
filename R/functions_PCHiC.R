


clean_interactions_data <- function(file){
  # Read as text file
  promoter_capture.data <- fread(file)
  
  # Filter out the following:
  # 1) oeName == "." indicates a bait-to-bait interaction
  # 2) baitName == NA indicates unknown promoter
  # 3) dist == NA indicates trans-chromosomal interaction
  # 4) Sex chromosomes (don't include for now)
  promoter_capture.cleaned <- promoter_capture.data %>%
    dplyr::filter(oeName == ".",
                  !is.na(baitName),
                  # !is.na(dist),
                  !(baitChr %in% c("X", "Y"))) %>%
    # Add identifiers for each bait and other end
    dplyr::mutate(bait.id = to_ucsc_format_v(baitChr, baitStart, baitEnd),
                  oe.id = to_ucsc_format_v(oeChr, oeStart, oeEnd)) %>%
    # Drop unneeded columns
    dplyr::select(-c(clusterID, clusterPostProb, oeName, oeID, baitID))
  
  # Get some summary statistics
  promoter_capture.cleaned$min.chicago <- apply(X=promoter_capture.cleaned[, ..BLOOD.CELL.TYPES], MARGIN=1, FUN=min)
  promoter_capture.cleaned$med.chicago <- apply(X=promoter_capture.cleaned[, ..BLOOD.CELL.TYPES], MARGIN=1, FUN=median)
  promoter_capture.cleaned$mean.chicago <- apply(X=promoter_capture.cleaned[, ..BLOOD.CELL.TYPES], MARGIN=1, FUN=mean)
  promoter_capture.cleaned$max.chicago <- apply(X=promoter_capture.cleaned[, ..BLOOD.CELL.TYPES], MARGIN=1, FUN=max)
  
  # Return this big fella
  promoter_capture.cleaned
}



make_granges_with_common_field_prefix <- function(gr, prefix=""){
  # Helper for interactions data. Basically allows you to specify the predix
  # so you can say "bait" and this function will figure out that
  # the output GRanges nees baitChr, baitStart, baitEnd
  
  makeGRangesFromDataFrame(gr,
                           seqnames.field =  paste0(prefix, "Chr"),
                           start.field = paste0(prefix, "Start"),
                           end.field = paste0(prefix, "End"),
                           keep.extra.columns = T)
}



drop_granges_columns <- function(data){
  # Remove these columns
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


subset_overlaps_return_unique <- function(features.gr, dmps.gr){
  seqlevelsStyle(dmps.gr) <- "UCSC"
  seqlevelsStyle(features.gr) <- "UCSC"
  
  subsetByOverlaps(features.gr, dmps.gr) %>%
    unique()
}




count_unique_overlaps <- function(dmps.gr, features.gr){
  # Helper function for null distribution (enrichment analysis)
  seqlevelsStyle(dmps.gr) <- "UCSC"
  seqlevelsStyle(features.gr) <- "UCSC"
  
  overlaps <- findOverlaps(dmps.gr, features.gr)
  
  # Tally
  return(list(
    "N.features" = length(unique(subjectHits(overlaps))),
    "N.dmps" = length(unique(queryHits(overlaps)))
  ))
}



# Enrichment Analysis -----------------------------------------------------

test_enhancer_enrichment_for_dmps <- function(null.space.gr, dmps.gr, features.gr, B=10000){
  
  N.simulated.dmps <- length(dmps.gr)
  
  # All possible DMPs in simulation
  null.ix <- 1:length(null.space.gr)
  
  N.enhancers <- rep(NA, B)
  N.dmps <- rep(NA, B)
  
  # Run the simulation
  for (b in 1:B){
    sim.ix <- sample(null.ix, size = N.simulated.dmps, replace = F)
    my.counts <- count_unique_overlaps(null.space.gr[sim.ix, ], features.gr)
    
    # Tally
    N.enhancers[b] <- my.counts$N.features
    N.dmps[b] <- my.counts$N.dmps
  }
  return(list("simulated.df" = data.frame(N.enhancers, N.dmps),
              "test.stat" = count_unique_overlaps(dmps.gr, features.gr)))
}

plot_enhancer_enrichment_for_dmps <- function(enrichment.result, file) {
  xx <- enrichment.result$simulated.df$N.dmps
  test.value <- enrichment.result$test.stat$N.dmps
  pval <- compute_empirical_p(xx, test.value)
  
  # Clean up
  if (pval == 0){
    pval.str <- "< 0.001"
  } else {
    pval.str <- round(pval, digits = 2)
  }
  
  z <- enrichment.result$simulated.df %>%
    ggplot(aes(x = N.dmps, y = after_stat(density))) +
    geom_histogram(bins = 25) +
    theme_minimal() +
    geom_vline(xintercept = enrichment.result$test.stat$N.dmps, color = "red") +
    xlab("Number of unique DMPs residing in enhancer") +
    ylab("Density") +
    labs(caption = paste0("Two-sided p-value: ", pval.str)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))
  
  cowplot::save_plot(filename = file, z)
}



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

