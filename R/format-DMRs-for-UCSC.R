library(data.table)

IFILE <- "../../DataRaw/2022-12-28-ExperimentSummary-v1/DMR.pis.bed"
OFILE <- "../../DataDerived/2023-01-02-DMRs.bed"

data <- fread(IFILE) %>%
  dplyr::select(1:3) %>%
  dplyr::distinct()

colnames(data) <- NULL

sink(OFILE)
cat('browser position chr22:20100000-20140000\n')
cat('track name=DMRs description="DMRs colored blue" color=0,0,255,\n')
cat('#chrom chromStart chromEnd')
print(data, row.names = FALSE)
sink()

