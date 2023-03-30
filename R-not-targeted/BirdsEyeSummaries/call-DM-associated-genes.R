
# call-DM-associated-genes.R
#

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(parallel)
  library(harmonicmeanp)
  library(tidyverse)
})


# Get the function that allows us to expand up and downstream
source("../GeneSchematics/util-ensdb.R")
ensdb.genes <- get_expanded_ensdb(3000, 200)

IFILE <- "../../DataRaw/2023-02-14-Summaries-v6/pvals.bed"
OFILE <- "../../DataDerived/DMGenes-FishersP.RData"
ODIR <- "../../Figs/DM-associated-genes/"
NATGEN.FILE <- "../../DataReference/NatureGenetics042022_AD_genes.txt"


# Load data ---------------------------------------------------------------
pvals <- fread(IFILE)
pvals.gr <- makeGRangesFromDataFrame(pvals,
                                     keep.extra.columns = T,
                                     starts.in.df.are.0based = T)

# Just a vector
gwas.genes <- read.table(NATGEN.FILE)$V1


fishers_method <- function(p.vec){
  #p.vec is a vector of p-values
  k <- length(p.vec) # degrees of freedom
  chi.stat <- -2 * sum(log(p.vec))
  pchisq(chi.stat, df = 2*k)
}


wrapper <- function(i){
  if (i %% 1000 == 0){print(paste0("Iteration: ", i))}

  fishers.p <- NA
  min.p <- NA
  harmonic.mean.p <- NA

  # Subset and get number of CpGs
  cpgs.in.gene <- subsetByOverlaps(pvals.gr, ensdb.genes[i])
  k <- length(cpgs.in.gene)

  # Symbol (BRCA1)
  gn <- ensdb.genes$gene_name[i]

  # ENSEMBL ID (ENS0000wdeofwodifj)
  gn.id <- ensdb.genes$gene_id[i]

  if (length(cpgs.in.gene) > 0){
    # Get vector of p-values
    pp <- cpgs.in.gene$p.from.ss

    # Three methods of combining p-values
    min.p <- min(pp, na.rm = T)

    # Fisher's method
    fishers.p <- fishers_method(pp)

    # Harmoinc mean
    harmonic.mean.p <- p.hmp(pp, L = length(pp))

  }
  data.frame(fishers.p = fishers.p,
             min.p = min.p,
             harmonic.mean.p = harmonic.mean.p,
             k = k,
             gene.name = gn,
             gene.id = gn.id)
}


# Setup for parallelization -----------------------------------------------

# N.chunks <- 5
#
# # Split the 1,2,3,... into N.chunks (almost equal) sized vectors
# chunk.indices <- split(1:10, 1:N.chunks)
#
# run_thru_one_chunk <- function(chunk.ix){
#   out <- do.call(rbind, base::lapply(X=chunk.ix, FUN=wrapper))
#   row.names(out) <- NULL
#   out
# }
#
# # Demo on one "core" job
# run_thru_one_chunk(chunk.indices[[1]])

# Full fledge version -----------------------------------------------------

# chunk.indices <- split(1:1000, 1:N.chunks)
#
# # And demo on parallel version
# start <- Sys.time()
# result <- do.call(rbind, mclapply(X = chunk.indices,
#                                   FUN = run_thru_one_chunk,
#                                   mc.cores = N.chunks))
# Sys.time() - start



# Intensive part
B <- length(ensdb.genes)

# 1000 genes takes 5.5 mins
# 100 genes takes 40 secs with 6 cores; 1.25 mins with lapply; also 40 secs with 2 cores
# 20,000 should take ~2 hours

combp.result <- mclapply(X=1:B, FUN=wrapper, mc.cores = 2)

save(combp.result, file = OFILE)
# END

# Harmonic P-value distribution -------------------------------------------
#
# HMPs <- all.p$harmonic.mean.p[!is.na(all.p$harmonic.mean.p)]
#
# ix <- (HMPs > 0) & (HMPs < 1)
#
# hist(HMPs, breaks = 25)
#
# hmp.fdr <- fdrtool::fdrtool(HMPs[ix], statistic =  "pvalue")
#
# # Distribution of P-values
# png(file.path(ODIR, "2023-02-28-hist-fishersP-forSK.png"))
# hist(all.p$fishers.p,
#      main = "Distribution of combined P-values",
#      xlab = "Fisher's P-value",
#      ylab = "Density",
#      probability = T)
# dev.off()
#
# # How many sig at Bonferroni?
# p.vec <- P.df$fishers.p[!is.na(P.df$fishers.p)]
#
# ix <- p.vec < (0.05 / length(p.vec))
#
# N.sig <- sum(ix)
# N.sig
#
#
# # Too bad of a bias with gene length? -------------------------------------
#
# all.p$gene.width <- width(ensdb.genes)[1:B]
#
# p <- all.p %>%
#   drop_na() %>%
#   dplyr::mutate(Significance = ifelse(fishers.p < 0.01 / B,
#                                      "Significant\n(Bonferroni 0.01)",
#                                      "Not significant")) %>%
#   ggplot(aes(x = log10(gene.width), -log10(fishers.p), color = Significance)) +
#   geom_point() +
#   xlab("log10(Gene width)") +
#   ylab("-log10(Fisher's P)") +
#   theme_minimal() +
#   scale_color_manual(values = c("grey", "red")) +
#   ggtitle("Fisher's p vs gene width")
#
#
# cowplot::save_plot(
#   filename = file.path(ODIR, "2023-02-28-scatter-gene-width-vs-fishersP-forSK.png"),
#   p)
#
# # Compare with GWAS -------------------------------------------------------
#
# intersect(all.p$gene.name[ix], gwas.genes)
