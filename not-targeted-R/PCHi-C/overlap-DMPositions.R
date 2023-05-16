# 1-overlap-DMPositions.R
# Three ouputs (stored in one RData object)
# 1) Interactions w/ any number of DMPs
# 2) Simulation results (is the above significant)
# 3)

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



# Parameters --------------------------------------------------------------


MIN.CHICAGO <- 5
ALPHA <- 0.05


# Paths -------------------------------------------------------------------

# Data
IFILE <- "../../DataRaw/2023-02-14-Summaries-v6/pvals.bed"
PCHIC.FILE <- "../../DataReference/PCHi-C/PCHiC_peak_matrix_cutoff5.txt"
EXPRESSION.FILE <- "../../DataReference/DEGenes/2020-Shigemizu-AD-RNAseq-DEGenes.xlsx"

# Utility
REF.CHAIN.PATH <- "../../DataReference/hg38ToHg19.over.chain"

# Outputs
ODIR <- "../../DataDerived/Outputs-PCHi-C/"

# Get valid gene symbols --------------------------------------------------
gene.ids <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene,
                  single.strand.genes.only = FALSE,
                  columns = "gene_id")



# Constants ---------------------------------------------------------------

VALID.SYMBOLS <- OrganismDbi::select(Homo.sapiens, names(gene.ids), "SYMBOL", "GENEID")$SYMBOL

# All cell types available
ALL.CELL.TYPES <- c("Mon", "Mac0", "Mac1", "Mac2",
                "Neu", "MK", "EP", "Ery", "FoeT", "nCD4", "tCD4",
                "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB")

# Ones we're interested in
# SUB.CELL.TYPES <- c("Neu", "Mon",
#                     "nCD4", "tCD4", "aCD4", "naCD4",
#                     "nCD8", "tCD8",
#                     "nB", "tB")
SUB.CELL.TYPES <- ALL.CELL.TYPES

dir.create(ODIR, showW = F)

# Setup -------------------------------------------------------------------
# Load data and perform light processing
# Differentialyl Expressed genes
expression.df <- readxl::read_xlsx(EXPRESSION.FILE, skip = 1) %>%
  dplyr::rename(DE.FDR = FDR,
                DE.P = `P-value`,
                gene.name = `gene name`)


# PCHi-C Data Read-in and Munging -----------------------------------------

# Now the PCHI-C data (don't lift; keep in hg19)
# 1. Drop the cell types we don't care about (see `sub.cell.types` for types of interest)
# 2. Compute a max chicago score for remaining cells and subset
promcap.raw <- fread(PCHIC.FILE)

promcap <- promcap.raw %>%
  dplyr::select(-base::setdiff(ALL.CELL.TYPES, SUB.CELL.TYPES)) %>%
  dplyr::mutate(bait.id = paste0("chr", baitChr,":", baitStart, "-", baitEnd),
                oe.id = paste0("chr", oeChr,":", oeStart, "-", oeEnd)) %>%
  rowwise() %>%
  dplyr::mutate(min.chicago = pmin(!!!rlang::syms(SUB.CELL.TYPES)),
                med.chicago = median(!!!rlang::syms(SUB.CELL.TYPES)))

# For a given link (row), what percentage are signficant?
X <- dplyr::select(promcap, all_of(SUB.CELL.TYPES))
promcap$percent.sig <- rowSums(X >= 5) / ncol(X)

# Keep only ubiquitous enhancers
promcap.filtered <- promcap %>%
  dplyr::filter(percent.sig >= 0.5)

# Don't map baits (genic loci),
# instead map the "other ends", which are likely enhancers
promcap.filt.gr <- makeGRangesFromDataFrame(promcap.filtered,
                                       seqnames.field = "oeChr",
                                       start.field = "oeStart",
                                       end.field = "oeEnd",
                                       keep.extra.columns = T)

# Rename to chr1 format
seqlevelsStyle(promcap.filt.gr) <- "UCSC"


# Read and clean DMPs -----------------------------------------------------

chain <- import.chain(REF.CHAIN.PATH)

# 1. Read only the necessary columns for pvals
# 2. Conver to GRanges
# 3. Lift from hg38 to hg19
pvals.lifted.gr <- fread(IFILE) %>%
  dplyr::select(-c(p.from.zz, p.from.DSS, lfdr.from.zz)) %>%
  makeGRangesFromDataFrame(starts.in.df.are.0based = T,
                           keep.extra.columns = T) %>%
  liftOver(chain) %>%
  unlist()


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
                     !!!rlang::syms(SUB.CELL.TYPES))

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

# Test 1
# Null: the number of DMPs residing in enhancers is more than CpGs being drawn randomly
dmps.gr <- pvals.lifted.gr[pvals.lifted.gr$lfdr.from.ss < ALPHA]

# We'll sample N.dmps number of CpGs
N.dmps <- length(dmps.gr)

dmps.in.enhancer <- find_and_curate_overlaps(dmps.gr, promcap.filt.gr)
N.dmps.in.enhancer <- nrow(dmps.in.enhancer)

G <- length(pvals.lifted.gr) # number of CpGs (used for index)
B <- 1000 # number of simulations

test1.draws <- rep(NA, B)

set.seed(919)
for (b in 1:B){

  if (b %% 1000 == 0){print(paste0("Iteration: ", b))}

  # Sample G CpGs from all of the 25+ million possible
  ix <- sample(1:G, size = N.dmps)
  random.dmps <- pvals.lifted.gr[ix]

  # Overlap with PCHi-C
  test1.draws[b] <- nrow(find_and_curate_overlaps(random.dmps, promcap.filt.gr))
}

N.dmps.in.enhancer
hist(test1.draws)



# Integrate with Differential Expression ----------------------------------

# When we split
dmps.in.enhancer$dmp.id <- 1:nrow(dmps.in.enhancer)


dmps.in.gene.enhancer <-
  dmps.in.enhancer %>%
  mutate(gene.name = strsplit(linked.gene, ";")) %>%
  unnest(gene.name)

# Join with expression data
expression.methy.df <- dmps.in.gene.enhancer %>%
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
expression.methy.df$`Differential expression rank` <- rank(ex.vs.meth.df$DE.abs.effect, ties.method = "random")
expression.methy.df$`Differential methylation rank` <- rank(ex.vs.meth.df$DM.abs.effect,  ties.method = "random")

# Test direction
tb <- table(ex.vs.meth.df$DM.direction, ex.vs.meth.df$DE.direction)
tb
fisher.test(tb)

# Test magnitude
kp.test <- cor.test(expression.methy.df$DM.abs.effect,
                    expression.methy.df$DE.abs.effect,
                    method="kendall")

kp.p <- signif(kp.test$p.value, 2)
kp.p


# Save --------------------------------------------------------------------

ofile <- file.path(ODIR, "meth-w-expression-integrated.RData")
save(expression.methy.df, dmps.in.gene.enhancer,
     kp.p, ALPHA, file = ofile)

#END

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
