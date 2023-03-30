
array.gr <-
  fread("../../DataReference/madrid_cpgs_list.csv", skip = 1) %>%
  dplyr::select(-c(Chromosome, `Relation to Island`)) %>%
  dplyr::filter(`aLIS P-value` < 0.05) %>%
  separate(col = hg38_coordinates, into = c("chr", "start"), sep = "\\:") %>%
  dplyr::mutate(
    start = as.numeric(start) - 2,
    end = as.numeric(start) +4) %>%
  makeGRangesFromDataFrame()


overlaps <- findOverlaps(array.gr, pvals.gr[pvals.gr$lfdr.from.ss <0.05])


pvals.df %>%
  dplyr::filter(lfdr.from.ss < 0.05) %>%
  dplyr::select(chr, start, end, diagnostic_group_coded) %>%
  fwrite("2023-02-16-DMPsForAndy.csv")
