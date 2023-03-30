
# Output directory
pvals.file <- "../../DataRaw/2023-02-14-Summaries-v6/pvals.bed"
ODIR <- "../../Figs/DMPositions/"

dir.create(ODIR, showWarn = F, recurs = T)

print(paste0("Reading in ", pvals.file))

# Constants
LFDR.CUT <- 0.05
C.HYPER <- "#0073C2FF"
C.HYPO <- "#EFC000FF"

# When annotating, the text is lower case, abbreviated
MAPPING <-
  data.frame(
    Annotation = c("1to5kb", "3UTRs", "5UTRs", "exons", "introns", "promoters",
                   "inter", "islands", "shelves", "shores", "Other"),
    clean = c("1-5 kb", "3' UTR", "5' UTR", "Exon", "Intron", "Promoter",
              "Open sea", "CpG island", "CpG shelf", "CpG shore", "Other")
  )
