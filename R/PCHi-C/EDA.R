# Results
# List of genes
# And how the DMR behaves?
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)

REF.CHAIN.PATH <- "../../DataReference/hg19ToHg38.over.chain"
NATGEN.PATH <- "../../DataReference/NatureGenetics042022_AD_genes.txt"
DMR.PATH <- "../../DataRaw/2022-12-28-ExperimentSummary-v1/DMR.pis.bed"
OFILE <- "../../DataDerived/DMRs-overlapped-PCHi-C.bed"

genes <- sort(as.vector(read.table(NATGEN.PATH)$V1))

dmrs <- fread(DMR.PATH)

keep.cols <- c("chr.dmr", "start.dmr", "end.dmr",
               "nCG", "areaStat", "persistent.effect",
               "one.large.effect", "n.sig", "n.sig.hyper",
               "nearest_gene_name", "nearest_gene_id", "dist_to_nearest_gene")

dmrs.gr <- dmrs %>%
  dplyr::select(all_of(keep.cols)) %>%
  distinct() %>% # important; just unique DMRs
  makeGRangesFromDataFrame(starts.in.df.are.0based = T,
                           seqnames.field = "chr.dmr",
                           start.field = "start.dmr",
                           end.field = "end.dmr",
                           keep.extra.columns = T)

data <- fread("../../DataReference/PCHi-C/PCHiC_peak_matrix_cutoff5.txt")
names(data)

# Baits are essentially genes
length(unique(data.gr$baitName))

# Plug in different genes and you'll find them here
sum(str_detect(unique(data.gr$baitName), "RBFOX"), na.rm = T)

chain <- import.chain(REF.CHAIN.PATH)

# Don't map baits (genic loci),
# instead map the "other ends", which are regulatory
data.gr <- makeGRangesFromDataFrame(data,
                                    seqnames.field = "oeChr",
                                    start.field = "oeStart",
                                    end.field = "oeEnd",
                                    keep.extra.columns = T)
# Rename to chr1 format
seqlevelsStyle(data.gr) <- "UCSC"
data.lifted.gr <- unlist(liftOver(data.gr, chain))


# Overlap DMRs with PCHi-C ------------------------------------------------
overlaps <- findOverlaps(dmrs.gr, data.lifted.gr)

print(paste0("N DMRs considered: " , length(dmrs.gr)))
N <- length(unique(queryHits(overlaps)))
print(paste0("N unique DMRs overlapping an 'other end': ", N))

dmr.ix <- queryHits(overlaps)
data.ix <- subjectHits(overlaps)
bait.genes <- sort(unique(data.lifted.gr[data.ix]$"baitName"))

# Join
a <- as.data.frame(dmrs.gr[dmr.ix])
b <- as.data.frame(data.lifted.gr[data.ix])[, -1:-5]

result <- cbind(a, b)
dplyr::select(result, c(nearest_gene_name, baitName)) %>%
  arrange(nearest_gene_name) %>%
  distinct()

result %>%
  dplyr::filter(baitName == "BNIP3")

# What to plot?
# We could plot box plots like in the Bayes project. These would be at the location of the DMR
# Certainly interested in HLA-DQA1
# BNIP3 (dysregulation implicated in AD in rats)

# What do we want out of this analysis?

celltypes <- c("Mon", "Mac0", "Mac1", "Mac2", "Neu",
               "MK", "EP", "Ery", "FoeT", "nCD4", "tCD4",
               "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB")

result.2 <-
  result %>%
  pivot_longer(cols = all_of(celltypes),
               names_to = "celltype",
               values_to = "CHICAGO")

out <-
  result.2 %>%
  group_by(seqnames, start, end) %>%
  summarize(median.CHICAGO = median(CHICAGO),
            mean.CHICAGO = mean(CHICAGO),
            nCG = first(nCG),
            n.sig = first(n.sig),
            n.hyper = first(n.sig.hyper),
            baitName = first(baitName)) %>%
  arrange(-mean.CHICAGO)

write_tsv(out, file = OFILE)

#END
