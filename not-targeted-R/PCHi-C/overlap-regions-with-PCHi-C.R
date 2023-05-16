# Results
# List of genes
# And how the DMR behaves?
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)
library(UpSetR)
library(org.Hs.eg.db)
library(goplot)
library(BioCircos)
library(networkD3)
library(rbokeh)
library(kableExtra)

# baits -- parts of genome mostly in promoters of protein coding genes
# other ends -- "discoveries" where they're in close spatial proximity; likely enhancers
# Map DMRs to other ends

REF.CHAIN.PATH <- "../../DataReference/hg38ToHg19.over.chain"
NATGEN.PATH <- "../../DataReference/NatureGenetics042022_AD_genes.txt"
DMR.PATH <- "../../DataRaw/2023-01-24-Summaries-v4/DMRegions-pis.bed"
OFILE <- "../../DataDerived/DMRegions-PCHi-C-result.csv"
FIG.DIR <- "../../Figs/PCHi-C/"

# Now the PCHI-C data (don't lift; keep in hg19)
data <- fread("../../DataReference/PCHi-C/PCHiC_peak_matrix_cutoff5.txt")  %>%
  dplyr::mutate(bait.id = paste0("chr", baitChr,":", baitStart, "-", baitEnd),
                oe.id = paste0("chr", oeChr,":", oeStart, "-", oeEnd))

# genes <- sort(as.vector(read.table(NATGEN.PATH)$V1))
dmrs <- fread(DMR.PATH)

keep.cols <- c("chr.dmr", "start.dmr", "end.dmr",
               "nCG", "areaStat", "persistent.effect",
               "one.large.effect", "n.sig", "n.sig.hyper",
               "nearest_gene_name", "nearest_gene_id", "dist_to_nearest_gene")

dmrs.gr <- dmrs %>%
  dplyr::select(all_of(keep.cols)) %>%
  dplyr::mutate(DMR.id = paste0(chr.dmr,":", start.dmr, "-", end.dmr)) %>%
  distinct() %>% # important; just unique DMRs
  makeGRangesFromDataFrame(starts.in.df.are.0based = T,
                           seqnames.field = "chr.dmr",
                           start.field = "start.dmr",
                           end.field = "end.dmr",
                           keep.extra.columns = T)

chain <- import.chain(REF.CHAIN.PATH)
dmrs.lifted.gr <- unlist(liftOver(dmrs.gr, chain))

# Don't map baits (genic loci),
# instead map the "other ends", which are likely enhancers
data.gr <- makeGRangesFromDataFrame(data,
                                    seqnames.field = "oeChr",
                                    start.field = "oeStart",
                                    end.field = "oeEnd",
                                    keep.extra.columns = T)

# Rename to chr1 format
seqlevelsStyle(data.gr) <- "UCSC"

# Overlap DMRs with PCHi-C ------------------------------------------------
overlaps <- findOverlaps(dmrs.lifted.gr, data.gr)

print(paste0("N DMRs considered: " , length(dmrs.gr)))
N <- length(unique(queryHits(overlaps)))
print(paste0("N unique DMRs overlapping an 'other end': ", N))
N <- length(unique(subjectHits(overlaps)))

# Sankey of DMR - Overlapped - Gene ---------------------------------------
#
# links <-
#   data.frame(
#     source = c("DMRs", "DMRs", "Linked"),
#     target = c("Not linked", "Linked", "Genes"),
#     value = c(20, 7, 41))
#
# # From these flows we need to create a node data frame: it lists every entities involved in the flow
# nodes <- data.frame(
#   name=c(as.character(links$source),
#          as.character(links$target)) %>% unique()
# )
#
# links$IDsource <- match(links$source, nodes$name)-1
# links$IDtarget <- match(links$target, nodes$name)-1
#
#
# p.sankey <- sankeyNetwork(Links = links, Nodes = nodes,
#               Source = "IDsource", Target = "IDtarget",
#               Value = "value", NodeID = "name",
#               sinksRight = F)
#
#
# p.sankey
# f1 <- file.path( "sankey.html")
# f2 <- file.path( "sankey.png")
#
# saveNetwork(p.sankey, file = f1, selfcontained = F)

# Data Munging ------------------------------------------------------------
dmr.ix <- queryHits(overlaps); dmr.ix
data.ix <- subjectHits(overlaps); data.ix
# bait.genes <- sort(unique(data.gr[data.ix]$"baitName"))

# Join
a <- as.data.frame(dmrs.lifted.gr[dmr.ix]) %>%
  dplyr::transmute(dmr.chr = seqnames, dmr.start = start, dmr.end = end, DMR.id)
b <- as.data.frame(data.gr[data.ix]) %>%
  dplyr::transmute(oe.chr = seqnames, oe.start = start, oe.end = end, oe.id,
                   bait.chr = paste0("chr", baitChr), bait.start = baitStart,
                   bait.end = baitEnd, bait.id,
                   linked.gene = baitName, dist,
                   Mon, Mac0, Mac1, Mac2, Neu, MK, EP, Ery, FoeT, nCD4, tCD4,
                   aCD4, naCD4, nCD8, tCD8, nB, tB)

result <- cbind(a, b)
write_csv(result, OFILE)

# Use the interactions later
# interactions <- as.data.frame(data.gr[data.ix])

# Table of DMRs versus linked genes ---------------------------------------
# result %>%
#   group_by(DMR.id) %>%
#   summarize(Genes = str_replace_all(toString(sort(unique(baitName))), ";", ", ")) %>%
#   head(5) %>%
#   kable()%>%
#   kable_styling("striped", font_size = 14) %>%
#   save_kable(file = "../../Figs/PCHi-C/dmr-gene-links-table.png")


# Genes identified as significant by PCHi-C analysis
# symbols <- sort(unique(unlist(str_split(result$baitName, ";"))))
# symbols
#
# # Now map them to ensembl
# ensembls <- AnnotationDbi::mapIds(org.Hs.eg.db, symbols, 'ENSEMBL', 'SYMBOL')
# ensembls
#
# go.out <- goplot::run_gene_ontology(ensembls)
# go.df <- goplot::go_output_to_df(go.out)
# # go.df$Description <- sapply(str_split(go.df$Description, ","), '[[', 1)
# go.df$Description <- sapply(str_split(go.df$Description, " and "), '[[', 1)
#
# p <- goplot::plot_go_barchart(go.df) +
#   geom_hline(yintercept = -log10(0.05))
#
# cowplot::save_plot(filename = "../../Figs/Fig2/pchic-gene-ontology.png", p,
#                    base_width = 12, base_height = 6)

# Post-processing ---------------------------------------------------------



# What do we want out of this analysis?
# celltypes <- c("Mon", "Mac0", "Mac1", "Mac2", "Neu",
#                "MK", "EP", "Ery", "FoeT", "nCD4", "tCD4",
#                "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB")
#
# result.2 <-
#   result %>%
#   pivot_longer(cols = all_of(celltypes),
#                names_to = "celltype",
#                values_to = "CHICAGO")
#
# out <-
#   result.2 %>%
#   group_by(seqnames, start, end) %>%
#   summarize(median.CHICAGO = median(CHICAGO),
#             mean.CHICAGO = mean(CHICAGO),
#             nCG = first(nCG),
#             n.sig = first(n.sig),
#             n.hyper = first(n.sig.hyper),
#             baitName = first(baitName)) %>%
#   arrange(-mean.CHICAGO)
#
# write_tsv(out, file = OFILE)


# Top genes ---------------------------------------------------------------

# UpSet Plot Of Genes -----------------------------------------------------
# upset.df <- as.data.frame(ifelse(result[ , celltypes] < 5, 0, 1))
# upset.df$Name <- make.names(result$baitName, unique = T)
# # upset.df[ ,"Name"] <- as.vector(make.names(result$baitName, unique = F))
# # rownames(upset.df) <- make.names(result$baitName, unique = F)
#
# UpSetR::upset(data = upset.df,
#               # sets = celltypes,
#               nsets = 5,
#               order.by = "freq",
#               )
#
# tmp <- upset.df %>%
#   pivot_longer(all_of(celltypes)) %>%
#   group_by(Name) %>%
#   summarize(n.sig = sum(value)) %>%
#   arrange(-n.sig) %>%
#   head(20) %>%
#   arrange(n.sig)
#
# tmp$Name <- factor(tmp$Name, levels = tmp$Name)
#
# p <- tmp %>%
#   dplyr::mutate(Name = as.factor(Name)) %>%
#   ggplot(aes(x = Name, y = n.sig)) +
#   geom_col() +
#   coord_flip() +
#   theme_minimal() +
#   ylab("Number of significant cell types") +
#   xlab("Promoter gene name") +
#   ggtitle("Top 20 promoter-enhancer interactions") +
#   theme(plot.background = element_rect(fill = "white"))
#
# cowplot::save_plot(
#   filename = file.path(FIG.DIR, paste0(Sys.Date(), "-numsig-by-ct.png")),
#   base_width = 9,
#   p)
#
#
#
#
# # Format interactions -----------------------------------------------------
#
# score <- 55
#
# left.bed <- interactions[ ,c("seqnames", "start", "end")]
# left.bed$inter <- paste0(interactions$seqnames, ":",
#                          interactions$baitStart, "-",
#                          interactions$baitEnd,
#                          ",", score)
#
#
# right.bed <- interactions[ ,c("seqnames", "baitStart", "baitEnd")]
# names(right.bed) <- c("seqnames", "start", "end")
# right.bed$inter <- paste0(interactions$seqnames, ":",
#                           interactions$start, "-",
#                           interactions$end,
#                           ",", score)
#
# N <- nrow(left.bed)
# orderix <- rep(1:N, each = 2) + rep(c(0, N), times = N)
# pchic <- rbind(left.bed, right.bed)[orderix, ]
#
#
# write_tsv(pchic, col_names = F,
#           file = "../../DataDerived/interactions.txt")
#
# #END
#
#
