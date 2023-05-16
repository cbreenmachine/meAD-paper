library(tidyverse)
library(data.table)

df <- fread("../../DataRaw/2023-01-24-Summaries-v4/pvals.bed")

out <- data.frame(
  N.Sig = sum(df$lfdr <= 0.05),
  N.Total = nrow(df)
)
