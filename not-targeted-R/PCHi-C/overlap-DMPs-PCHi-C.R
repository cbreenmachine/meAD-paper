# Results
# List of genes
# And how the DMR behaves?

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(liftOver)
  library(ggbump)
  library(org.Hs.eg.db)
  library(goplot)
  library(tibble)
  library(tidygraph)
  library(kableExtra)
  library(pheatmap)
  library(ggsci)
  library(ggraph)
  library(igraph)
})


# Notes from 2/21/23
# When min.chicago > 5; we get the result in meAD-paper
# But this is very stringent.

ALPHA <- 0.05


# Get valid gene symbols --------------------------------------------------
gene.ids <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene,
                  single.strand.genes.only = FALSE,
                  columns = "gene_id")

VALID.SYMBOLS <- OrganismDbi::select(Homo.sapiens, names(gene.ids), "SYMBOL", "GENEID")$SYMBOL



# All cell types available
cell.types <- c("Mon", "Mac0", "Mac1", "Mac2",
                "Neu", "MK", "EP", "Ery", "FoeT", "nCD4", "tCD4",
                "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB")

# Ones we're interested in
sub.cell.types <- c("Neu", "Mon",
                    "nCD4", "tCD4", "aCD4", "naCD4",
                    "nCD8", "tCD8",
                    "nB", "tB")


# baits -- parts of genome mostly in promoters of protein coding genes
# other ends -- "discoveries" where they're in close spatial proximity; likely enhancers
# Map DMRs to other ends

# Lots of baths to read in
REF.CHAIN.PATH <- "../../DataReference/hg38ToHg19.over.chain"
DE.PATH <- "../../DataReference/DEGenes/2020-Shigemizu-AD-RNAseq-DEGenes.xlsx"
IFILE <- "../../DataRaw/2023-02-14-Summaries-v6/pvals.bed"

FIG.DIR <- "../../Figs/PCHi-C/"
DATA.ODIR <- "../../DataDerived/"


# DE Genes ----------------------------------------------------------------
# differentially expressed
expression.df <- readxl::read_xlsx(DE.PATH, skip = 1) %>%
  dplyr::rename(DE.FDR = FDR,
                DE.P = `P-value`,
                gene.name = `gene name`)


# PCHi-C Data Read-in and Munging -----------------------------------------

# Now the PCHI-C data (don't lift; keep in hg19)
# 1. Drop the cell types we don't care about (see `sub.cell.types` for types of interest)
# 2. Compute a max chicago score for remaining cells and subset if
# it's large (> 10)
promcap.raw <- fread("../../DataReference/PCHi-C/PCHiC_peak_matrix_cutoff5.txt")
promcap <- promcap.raw %>%
  dplyr::select(-base::setdiff(cell.types, sub.cell.types)) %>%
  dplyr::mutate(bait.id = paste0("chr", baitChr,":", baitStart, "-", baitEnd),
                oe.id = paste0("chr", oeChr,":", oeStart, "-", oeEnd),
                ) %>%
  rowwise() %>%
  dplyr::mutate(min.chicago = pmin(!!!rlang::syms(sub.cell.types)),
                med.chicago = median(!!!rlang::syms(sub.cell.types)))

# For a given link (row), what percentage are signficant?
X <- dplyr::select(promcap, all_of(sub.cell.types))
promcap$percent.sig <- rowSums(X >= 5) / ncol(X)

promcap.filtered <- promcap %>%
  dplyr::filter(percent.sig >= 1)

# Don't map baits (genic loci),
# instead map the "other ends", which are likely enhancers
promcap.gr <- makeGRangesFromDataFrame(promcap.filtered,
                                       seqnames.field = "oeChr",
                                       start.field = "oeStart",
                                       end.field = "oeEnd",
                                       keep.extra.columns = T)

# Rename to chr1 format
seqlevelsStyle(promcap.gr) <- "UCSC"


# Plot significant cell types ---------------------------------------------
tmp <- data.frame(p.sig = colSums(promcap[ ,sub.cell.types] >= 5) / nrow(promcap)) %>%
  rownames_to_column("cell.type")

tmp$cell.type <- factor(tmp$cell.type, levels = tmp$cell.type[order(tmp$p.sig, decreasing = T)])

z <- tmp %>%
  ggplot(aes(x = cell.type, y = p.sig, fill = p.sig)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = 'none',
        plot.background = element_rect(fill = "white", color = "white")) +
  xlab("Cell type") +
  ylab("Percent significant") +
  ggtitle("Percent of interactions considered significant by cell type")

cowplot::save_plot(file.path(FIG.DIR, "barchart-percent-sig-by-celltype.png"), z)


# Read and clean DMPs -----------------------------------------------------
pvals.df <- fread(IFILE) %>%
  dplyr::select(-c(p.from.zz, p.from.DSS, lfdr.from.zz))

pvals.gr <- pvals.df %>%
  makeGRangesFromDataFrame(starts.in.df.are.0based = T,
                           keep.extra.columns = T)

chain <- import.chain(REF.CHAIN.PATH)
pvals.lifted.gr <- unlist(liftOver(pvals.gr, chain))


# Functions ------------------------------------------------------------

find_and_curate_overlaps <- function(dmps, pchic){
  overlaps <- findOverlaps(dmps, pchic)

  dmp.ix <- queryHits(overlaps)
  pchich.ix <- subjectHits(overlaps)

  # This is the data we want from the dmps
  a <- as.data.frame(dmps[dmp.ix])

  # Data from PCHi-C with some munging
  b <- as.data.frame(pchic[pchich.ix]) %>%
    dplyr::transmute(oe.chr = seqnames,
                     oe.start = start,
                     oe.end = end,
                     oe.id,
                     bait.chr = paste0("chr", baitChr),
                     bait.start = baitStart,
                     bait.end = baitEnd,
                     bait.id,
                     linked.gene = baitName,
                     dist,
                     !!!rlang::syms(sub.cell.types))

  result <- cbind(a, b) %>%
    dplyr::mutate(Chr = as.numeric(stringr::str_remove(seqnames, "chr"))) %>%
    arrange(Chr, start)

  return(result)
}


get_unique_genes_from_overlaps <- function(result, filter = F){
  # filter subsets by VALID.IDS
  symbols <- unlist(stringr::str_split(result$linked.gene, ";"))

  if (filter){
    sort(unique(symbols[symbols %in% VALID.SYMBOLS]))
  } else {
    sort(unique(symbols))
  }
}


# Non-parameteric test ----------------------------------------------------

# Two tests going on here
# 1) Are DMPs overlapping enhancers more than randomly selected CpGs?
# 2) Are DMPs linking to more/fewr genes than randomly selecte CpGs?

# We'll sample N.dmps number of CpGs
N.dmps <- sum(pvals.df$lfdr.from.ss < ALPHA)

# How many possible genes?
dummy <- find_and_curate_overlaps(promcap.gr[, -1:-25], promcap.gr)
all.possible.genes <- get_unique_genes_from_overlaps(dummy)

G <- nrow(pvals.df) # number of CpGs (used for index)
B <- 100 # number of simulations

# Where we'll store the randomly generated values
P.uniq.genes <- rep(NA, B)
P.enhancer.dmps <- rep(NA, B)

set.seed(919)
for (b in 1:B){

  if (b %% 1000 == 0){print(paste0("Iteration: ", b))}

  # Sample G CpGs from all of the 25+ million possible
  ix <- sample(1:G, size = N.dmps)
  dmps.sim <- pvals.gr[ix]

  # Overlap with PCHi-C
  rr <- find_and_curate_overlaps(dmps.sim, promcap.gr)

  overlapped.genes <- get_unique_genes_from_overlaps(rr, filter = T)
  P.uniq.genes[b] <- length(overlapped.genes) / length(all.possible.genes)
  P.enhancer.dmps[b] <- nrow(rr) / length(promcap.gr)
}

# Test statistic
result <- find_and_curate_overlaps(pvals.lifted.gr[pvals.lifted.gr$lfdr.from.ss < ALPHA],
                                   promcap.gr)

T.uniq.genes <- length(get_unique_genes_from_overlaps(result, filter = T))
T.enhancer.dmps <- nrow(result)

p.uniq.genes <- 1 - (sum(T.uniq.genes / length(all.possible.genes) <= P.uniq.genes) / B)
p.enhancer.dmps <- 1 - (sum(T.enhancer.dmps / length(promcap.gr) >= P.enhancer.dmps) / B)


# Plot the result with p-values -------------------------------------------
sim.df <- data.frame(
  Number = c(P.uniq.genes ,
             P.enhancer.dmps),
  Feature = c(rep("Unique genes", B), rep("DMPs in enhancer", B))
  )

teststat.df <- data.frame(Observed = c(T.uniq.genes / length(all.possible.genes),
                                       T.enhancer.dmps/ length(promcap.gr)),
                          Feature = c("Unique genes", "DMPs in enhancer"))


cap <- "Distribution of simulated counts of (i) DMPs overlapping enhancers and
(ii) unique protein coding genes linked by PCHi-C dataset. Red lines indicate
observed values. One sided tests reject the null for both count distributions
at level 0.01. This indicates that DMPs are enriched in enhancers, but the
number of linked genes is under (i.e. multiple DMPs are linking to the same
genes, either because they reside in the same enhancer or because multiple
enhancers are in close proximity to promoter.)"

z <- sim.df %>%
  ggplot(aes(x = Number, y = after_stat(..density..))) +
  geom_histogram(bins = 50, color = "black", fill="darkgrey") +
  facet_wrap(.~Feature, scales = "free") +
  geom_vline(data = teststat.df,
             aes(xintercept = Observed),
             color = "red") +
  theme_minimal() +
  xlab("Count") +
  ylab("Density") +
  # labs(caption = cap) +
  theme(legend.position = "none",
        plot.caption = element_text(hjust = 0),
        plot.background = element_rect(fill = "white", color = "white"))
z

# cowplot::save_plot(filename = file.path(FIG.DIR, "enrichment-test-histograms.png"),
#                    z, base_width = 7)


# Integrate with Differential Expression ----------------------------------

result.2 <-  result %>%
  mutate(gene.name = strsplit(linked.gene, ";")) %>%
  unnest(gene.name)

ofile <- file.path(DATA.ODIR, "genes-with-DMEnhancers.csv")
write_csv(result.2, ofile)

ex.vs.meth.df <- result.2 %>%
  inner_join(expression.df, by = "gene.name") %>%
  group_by(gene.name) %>%
  summarize(
    DM.effect = median(pi.diff),
    DE.effect = logFC,
    DM.direction = ifelse(DM.effect > 0, "Hyper", "Hypo"),
    DE.direction = Regulated,
    DM.abs.effect = (abs(DM.effect)),
    DE.abs.effect = (abs(DE.effect))) %>%
  distinct()



# Add ranking columns
set.seed(1234)
ex.vs.meth.df$`Differential expression rank` <- rank(ex.vs.meth.df$DE.abs.effect, ties.method = "random")
ex.vs.meth.df$`Differential methylation rank` <- rank(ex.vs.meth.df$DM.abs.effect,  ties.method = "random")

# Test direction
tb <- table(ex.vs.meth.df$DM.direction, ex.vs.meth.df$DE.direction)
tb
fisher.test(tb)

# Test magnitude
kp.test <- cor.test(ex.vs.meth.df$DM.abs.effect,
                    ex.vs.meth.df$DE.abs.effect,
                    method="kendall")

kp.p <- signif(kp.test$p.value, 2)

# There's ranking here!!!
p <- ex.vs.meth.df %>%
  pivot_longer(cols = contains("Diff"),
               values_to = "rank",
               names_to = "type")  %>%
  ggplot(aes(x = type, y = rank,
             group = gene.name,
             color = gene.name,
             label = gene.name)) +
  geom_bump() +
  geom_label(hjust = 1.05, size = 1.5) +
  theme_minimal() +
  scale_y_reverse() +
  theme(legend.position = "none",
        plot.caption = element_text(hjust = 0),
        plot.background = element_rect(fill = "white", color="white")) +
  xlab("") +
  ylab("Rank")
  cap <-    "
      Forty-three unique genes have both (i) DMP(s) in enhancer
      and (ii) differential expression. Shown are the rankings of the
      absolute value of logFC (DE rank) and ranking of absolute
      value of LOAD methylation effect (DMP rank) for corresponding
      genes. A Kendall's rank test gives a p-value < 0.001, indicating
      that rank correlation is non-zero."


p

cowplot::save_plot(filename = file.path(FIG.DIR, "rank-bump-plot.png"),
                   p, base_width = 4.5, base_height = 7)


ofile <- file.path(FIG.DIR, "genes-DMEnhancers-and-DiffExpression.csv")
write_csv(ex.vs.meth.df, ofile)

# WashU -------------------------------------------------------------------


format_for_washu <- function(rr){
  N <- nrow(rr) * 2

  output <- data.frame(matrix(NA, nrow = N, ncol = 4))

  for (i in seq(1, N, by=2)){
    output[i, ] <- rr[i, c("bait.chr", "start", "end", "bait.id")]
  }

  for (i in seq(2, N, by=2)){
    bb <- paste0(rr[i-1, "seqnames"], ":", rr[i-1, "start"], "-", rr[i-1, "end"] )
    output[i, ] <- c(rr[i-1, c("bait.chr", "oe.start", "oe.end")], bb)
  }
  output
}

washu.output <- format_for_washu(result) %>%
  drop_na() %>%
  dplyr::arrange(X1, X2) %>%
  dplyr::mutate(X4 = paste0(X4, ";55"),
                X5 = 1:nrow(.),
                V6 = ".")
washu.output %>%
  fwrite(sep = "\t", file = file.path(DATA.ODIR, "interactions.forWashU.txt"), col.names = F)

washu.output %>%
  dplyr::filter(X1 == "chr17") %>%
  pull(X2) %>% summary()


# Format for UCSC ---------------------------------------------------------

ucsc.order.cols <- c("chrom", "chromStart", "chromEnd",
                     "name", "score", "value", "exp", "color",
                     "sourceChrom", "sourceStart", "sourceEnd", "sourceName", "sourceStrand",
                     "targetChrom", "targetStart", "targetEnd", "targetName", "targetStrand")

ucsc.df <- result %>%
  dplyr::transmute(sourceChrom = oe.chr, sourceStart = oe.start, sourceEnd = oe.end,
                targetChrom = bait.chr, targetStart = bait.start, targetEnd = bait.end) %>%
  dplyr::mutate(chrom = sourceChrom,
                chromStart = floor((sourceStart + targetStart) / 2),
                chromEnd = floor((sourceEnd + targetEnd) / 2)) %>%
  dplyr::mutate(value = 5, exp = ".", color = 0,
                name = ".", score = 500,
                sourceName = "OtherEnd:enhancer", sourceStrand = ".",
                targetName = "Bait:promoter", targetStrand = ".") %>%
  dplyr::select((ucsc.order.cols)) %>%
  dplyr::filter(chrom != "chr11") %>% dplyr::filter(chrom != "chr21")



ofile <- file.path(DATA.ODIR, "interactions.ucsc.bed")

cat(
  file = ofile,
  'track type=interact name="PCHi-C Integration" description="Interactions" useScore=on maxHeightPixels=200:100:50 visibility=full\n')
write.table(ucsc.df, file = ofile, append = T,
            row.names = F, col.names = F, quote = F)

setdiff(result$linked.gene, ex.vs.meth.df$gene.name)


# Plot boxplots by gene ---------------------------------------------------

result.2 %>%
  dplyr::filter(linked.gene == "SORT1") %>%
  dplyr::select(lfdr.from.ss)


# Hive plot ---------------------------------------------------------------

result.2$node.DM <- paste0("(DM) ", result.2$gene.name)
expression.df$node.DE <- paste0("(DE) ", expression.df$gene.name)

xx <- sort(unique(result.2$node.DM))
yy <- sort(unique(expression.df$node.DE))

nodes <- data.frame(name = c(xx, yy),
                    type = c(rep("Methy", length(xx)), rep("Expr", length(yy))))

connections <- result.2 %>%
  inner_join(expression.df, by = "gene.name") %>%
  dplyr::select(starts_with("node"))
names(connections) <- c("from", "to")

gene.graph <- graph_from_data_frame(connections, vertices = nodes)


gene.graph %>%
  ggraph('hive', axis = type, sort.by = degree) +
  geom_edge_hive() +
  geom_axis_hive()

# Gene ontology for DMP linked genes --------------------------------------

ensembl.ids <-
  OrganismDbi::select(Homo.sapiens, result.2$gene.name,
                      "ENSEMBL", "SYMBOL")$ENSEMBL
ensembl.ids <- ensembl.ids[!is.na(ensembl.ids)]

go.out <- goplot::run_gene_ontology(ensembl.ids)
z <- go.out %>%
  goplot::go_output_to_df() %>%
  goplot::plot_go_barchart() +
  geom_hline(yintercept = -log10(0.05))

cowplot::save_plot(filename = file.path(FIG.DIR, "2023-02-15-GO-PCHiC-genes.png"),
                   z, base_width = 14, base_height = 8)


# Heatmap prep ------------------------------------------------------------
X.df <- de.pc.df %>%
  tidyr::drop_na() %>%
  group_by(`gene name`, Regulated) %>%
  arrange(-abs(pi.diff)) %>%
  slice(1)

# What to make of this?
# Certain amount of hyper methylated DMPs must be annotating to
# multiple genes

# PHEATMAP ----------------------------------------------------------------

X <- X.df %>%
  arrange(`gene name`) %>%
  column_to_rownames("gene name") %>%
  dplyr::select(all_of(sub.cell.types)) %>%
  base::as.data.frame()

Y <- ifelse(X < 5, 0, 1)


anno <- base::data.frame(DMP = ifelse(X.df$pi.diff< 0, "Hypo", "Hyper"))
anno$DE <- ifelse(X.df$Regulated == "UP", "Up", "Down")

rownames(anno) <- rownames(anno)
rownames(anno) <- rownames(Y)

# CHICAGO score too high
plot_and_save <- function(a, b){
  ofile <- file.path(FIG.DIR, paste0("pchic-heatmap-", a, "-thru-", b, ".png"))

  png(ofile, width = 1080, height = 960, res = 120)
  pheatmap(Y[a:b, ],
           cluster_cols = T,
           color = c("white", "black"),
           annotation_row = anno[a:b, ],
           cluster_rows = F,
           treeheight_row = 0, treeheight_col = 0)
  dev.off()
}

plot_and_save(1, 50)
plot_and_save(51, 100)


plot_and_save(201, 300)
plot_and_save(301, nrow(X))

ofile <- file.path(FIG.DIR, "pchic-heatmap-all-genes-clustered-xy.png" )
png(ofile, width = 1080, height = 960, res = 120)
pheatmap(Y,
         cluster_cols = T,
         cluster_rows = T,
         color = c("white", "black"),
         annotation_row = anno,
         show_rownames = F,
         treeheight_row = 0, treeheight_col = 0)
dev.off()

