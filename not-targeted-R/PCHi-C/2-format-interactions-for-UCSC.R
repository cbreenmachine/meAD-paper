# 2-format-interactions-for-UCSC.R
# Take the PCHi-C data integrated with DMPs and subset
# pvals to have parallel tracks (-log10(lFDR))
# and

# Output is hg19!!!

# After running, go to DataExport/UCSC and run `bash process.all.sh`


suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(GenomicRanges)
  library(liftOver)
})

pvals.file <- "../../DataRaw/2023-02-14-Summaries-v6/pvals.bed"
REF.CHAIN.PATH <- "../../DataReference/hg38ToHg19.over.chain"

DIR <- "../../DataDerived/Outputs-PCHi-C/"
ODIR <- "../../DataExport/UCSC-Lolly-w-PCHi-C/"

ifile <- file.path(DIR, "meth-w-expression-integrated.RData")

pad <- 1e5


# Load and process --------------------------------------------------------

load(ifile)

chain <- import.chain(REF.CHAIN.PATH)

pvals.lifted.gr <- fread(pvals.file) %>%
  dplyr::select(-c(p.from.zz, p.from.DSS, lfdr.from.zz)) %>%
  makeGRangesFromDataFrame(starts.in.df.are.0based = T,
                           keep.extra.columns = T) %>%
  liftOver(chain) %>%
  unlist()



# PCHi-C Result -----------------------------------------------------------

ucsc.order.cols <- c("chrom", "chromStart", "chromEnd",
                     "name", "score", "value", "exp", "color",
                     "sourceChrom", "sourceStart", "sourceEnd", "sourceName", "sourceStrand",
                     "targetChrom", "targetStart", "targetEnd", "targetName", "targetStrand")

ucsc.df <- dmps.in.gene.enhancer %>%
  dplyr::transmute(sourceChrom = oe.chr, sourceStart = oe.start, sourceEnd = oe.end,
                   targetChrom = bait.chr, targetStart = bait.start, targetEnd = bait.end) %>%
  dplyr::mutate(chrom = sourceChrom,
                chromStart = floor((sourceStart + targetStart) / 2),
                chromEnd = floor((sourceEnd + targetEnd) / 2)) %>%
  dplyr::mutate(value = 5, exp = ".", color = "255,0,0",
                name = ".", score = 500,
                sourceName = "OtherEnd:enhancer", sourceStrand = ".",
                targetName = "Bait:promoter", targetStrand = ".") %>%
  dplyr::select((ucsc.order.cols)) %>%
  dplyr::filter(!(chrom == "chr21" & chromStart == 69336611)) # one bad record

# Extend the interval of the interaction by [pad] in both directions
intervals.gr <- makeGRangesFromDataFrame(ucsc.df,
                                         seqnames.field = "chrom",
                                         start.field = "chromStart",
                                         end.field = "chromEnd",
                                         keep.extra.columns = T) %>%
  IRanges::resize(width = GenomicRanges::width(.) + pad*2, fix = "center")


# Now Subset pvals
pvals.sub.df <- subsetByOverlaps(pvals.lifted.gr, intervals.gr) %>%
  as.data.frame() %>%
  dplyr::rename(chrom = seqnames)

# # Write files chromosome by chromosome ------------------------------------
#
# write_interactions <- function(data, ofile){
#   # First line is hashed
#   cat(
#     file = ofile,
#     'track type=interact name="PCHi-C Integration" description="Interactions" useScore=on maxHeightPixels=200:100:50 visibility=full\n')
#
#   # Then append
#   write.table(data,
#               file = ofile,
#               append = T,
#               row.names = F,
#               col.names = F,
#               quote = F)
# }


# Big Lolly ---------------------------------------------------------------
format_for_lolly <- function(data){
  out <- data %>%
    dplyr::transmute(chrom,
                     start,
                     end,
                     name = ".",
                     score = round(-log10(lfdr.from.ss)),
                     strand = ".",
                     thickStart = start,
                     thickEnd = end,
                     color = ifelse(pi.diff < 0, "0,119,154", "239,192,0"),
                     lollySize = ifelse(lfdr.from.ss < 0.05, 5, 1))

  # Make non-sig points grey
  out[out$pValueLog < -log10(0.05), "color"] <- "220,220,220"
  out
}


# Subset of P-values
out <- format_for_lolly(pvals.sub.df)
ofile <- file.path(ODIR,  "lolly.bed")
write_tsv(out, ofile, col_names = F)

# All p-values (thinned)
out <- pvals.lifted.gr %>%
  as.data.frame() %>%
  dplyr::rename(chrom = seqnames) %>%
  format_for_lolly()

# Keep all significant loci,
# Every 25th non-significant
keepix <- c(
  which(pvals.lifted.gr$lfdr.from.ss <= 0.05),
  which(pvals.lifted.gr$lfdr.from.ss > 0.05)[c( rep(FALSE, 24), TRUE)]
)

length(keepix)

ofile <- file.path(ODIR,  "all.cpgs.lolly.bed")
write_tsv(out[keepix, ], ofile, col_names = F)

# Interactions ---------------------------------------------------------------

ofile <- file.path(ODIR, "interactions.bed")
write_tsv(ucsc.df, ofile, col_names = F)

# END

