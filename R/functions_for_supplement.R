
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

  # data.frame(pp) %>%
  #   ggplot(aes(x = pp, y = after_stat(..density..))) +
  #   geom_histogram(bins = 50, fill = "darkgrey", color = "black", linewidth=0.1) +
  #   xlim(c(0,1)) +
  #   ylim(c(0, 1.7)) +
  #   xlab("P-values") +
  #   ylab("Density") +
  #   theme_minimal() +
  #   labs(subtitle = ti) +
  #   theme(plot.background = element_rect(fill = "white", color= "white"))
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
  #
  #
  # plot_hist_from_pvals(pp.2)
  # hist.recorded <- recordPlot()
  #
  # plot_qq_from_pvals(pp.1)
  # qq.recorded <- recordPlot()
  # #
  # cowplot::plot_grid(ggdraw(qq.recorded),
  #                    ggdraw(hist.recorded),
  #                    labels = c("A", "B"),
  #                    rel_widths = c(0.5, 0.5))
}

run_and_save_pvals_plot <- function(pvals.data, file){
  pdf(file, width = 9, height = 5)
  plot_pvals_routine(pvals.data )
  dev.off()

  # Need to return a charcter vector
  file
}



# Coverage histograms -----------------------------------------------------

munge_sequencing_depth <- function(file){
  read_csv(file,
           show_col_types = F,
           col_names = c("sample_id", "Mean sequencing depth")) %>%
    mutate(sample_id = as.character(sample_id))
}


munge_effective_coverage <- function(file){

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
  cov
}


plot_coverage_hist_routine <- function(depth.file, cov.file, ofile){

  depth.df <- munge_sequencing_depth(depth.file)
  cov.df <- munge_effective_coverage(cov.file)

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


process_and_write_dmps <- function(){
  return()
}

#
process_and_write_gene_ontology_terms <- function(DMGenes.go.df, file){
  # First convert ENSEMBL IDs to genes...
  tmp <- convert_gene_ontology_ids_to_symbols(DMGenes.go.df)
  writexl::write_xlsx(tmp, path = file)
  file
}
