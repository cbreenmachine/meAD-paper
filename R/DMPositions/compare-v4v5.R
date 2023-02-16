library(data.table)

xx <- fread("../../DataRaw/2023-01-24-Summaries-v4/pvals.bed")$stat
yy <- fread("../../DataRaw/2023-02-09-Summaries-v5/pvals.bed")$stat

all(xx == yy)
