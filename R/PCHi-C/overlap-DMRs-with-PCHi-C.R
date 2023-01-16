# Results
# List of genes
# And how the DMR behaves?
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)
library(UpSetR)

REF.CHAIN.PATH <- "../../DataReference/hg38ToHg19.over.chain"

NATGEN.PATH <- "../../DataReference/NatureGenetics042022_AD_genes.txt"
DMR.PATH <- "../../DataRaw/2022-12-28-ExperimentSummary-v1/DMR.pis.bed"
OFILE <- "../../DataDerived/DMRs-overlapped-PCHi-C.bed"
FIG.DIR <- "../../Figs/PCHi-C/"

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

chain <- import.chain(REF.CHAIN.PATH)

dmrs.lifted.gr <- unlist(liftOver(dmrs.gr, chain))


# Now the PCHI-C data (don't lift; keep in hg19)
data <- fread("../../DataReference/PCHi-C/PCHiC_peak_matrix_cutoff5.txt")

# Don't map baits (genic loci),
# instead map the "other ends", which are regulatory
data.gr <- makeGRangesFromDataFrame(data,
                                    seqnames.field = "oeChr",
                                    start.field = "oeStart",
                                    end.field = "oeEnd",
                                    keep.extra.columns = T)

# Rename to chr1 format
seqlevelsStyle(data.gr) <- "UCSC"

# Plug in different genes and you'll find them here
sum(str_detect(unique(data.gr$baitName), "RBFOX"), na.rm = T)




# Overlap DMRs with PCHi-C ------------------------------------------------
overlaps <- findOverlaps(dmrs.lifted.gr, data.gr)

print(paste0("N DMRs considered: " , length(dmrs.gr)))
N <- length(unique(queryHits(overlaps)))
print(paste0("N unique DMRs overlapping an 'other end': ", N))

dmr.ix <- queryHits(overlaps)
data.ix <- subjectHits(overlaps)
bait.genes <- sort(unique(data.gr[data.ix]$"baitName"))

# Join
a <- as.data.frame(dmrs.lifted.gr[dmr.ix])
b <- as.data.frame(data.gr[data.ix])[, -1:-5]

# Use the interactions later
interactions <- as.data.frame(data.lifted.gr[data.ix])

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



# UpSet Plot Of Genes -----------------------------------------------------
upset.df <- as.data.frame(ifelse(result[ , celltypes] < 5, 0, 1))
upset.df$Name <- make.names(result$baitName, unique = T)
# upset.df[ ,"Name"] <- as.vector(make.names(result$baitName, unique = F))
# rownames(upset.df) <- make.names(result$baitName, unique = F)

UpSetR::upset(data = upset.df,
              # sets = celltypes,
              nsets = 5,
              order.by = "freq",
              )

tmp <- upset.df %>%
  pivot_longer(all_of(celltypes)) %>%
  group_by(Name) %>%
  summarize(n.sig = sum(value)) %>%
  arrange(-n.sig) %>%
  head(20) %>%
  arrange(n.sig)

tmp$Name <- factor(tmp$Name, levels = tmp$Name)

p <- tmp %>%
  dplyr::mutate(Name = as.factor(Name)) %>%
  ggplot(aes(x = Name, y = n.sig)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  ylab("Number of significant cell types") +
  xlab("Promoter gene name") +
  ggtitle("Top 20 promoter-enhancer interactions") +
  theme(plot.background = element_rect(fill = "white"))

cowplot::save_plot(
  filename = file.path(FIG.DIR, paste0(Sys.Date(), "-numsig-by-ct.png")),
  base_width = 9,
  p)




# Format interactions -----------------------------------------------------

score <- 55

left.bed <- interactions[ ,c("seqnames", "start", "end")]
left.bed$inter <- paste0(interactions$seqnames, ":",
                         interactions$baitStart, "-",
                         interactions$baitEnd,
                         ",", score)


right.bed <- interactions[ ,c("seqnames", "baitStart", "baitEnd")]
names(right.bed) <- c("seqnames", "start", "end")
right.bed$inter <- paste0(interactions$seqnames, ":",
                          interactions$start, "-",
                          interactions$end,
                          ",", score)

N <- nrow(left.bed)
orderix <- rep(1:N, each = 2) + rep(c(0, N), times = N)
pchic <- rbind(left.bed, right.bed)[orderix, ]


write_tsv(pchic, col_names = F,
          file = "../../DataDerived/interactions.txt")

#END


