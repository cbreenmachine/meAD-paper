
library(tidyverse)

ifile <- "../../DataReference/CTCF/metadata.tsv"
ofile <- "../../DataReference/CTCF/metadata-v2.tsv"

meta.df <- read_tsv(ifile, show_col_types = F)


sub.df <- meta.df %>%
  dplyr::select(c(
    `File format type`,
    `Experiment accession`,
    `File download URL`,
    `Biosample term name`
  )) %>%
  dplyr::filter(`File format type` == "idr_ranked_peak") %>%
  dplyr::mutate(file_name =
                  str_split_fixed(`File download URL`, "download\\/", 2)[ ,2] %>%
                  str_remove("\\.gz")
                ) %>%
  janitor::clean_names() %>%
  dplyr::select(-c("file_download_url", "file_format_type"))

one_only.df <- sub.df %>%
  group_by(biosample_term_name) %>%
  summarize(ff = dplyr::first(file_name))

write_tsv(one_only.df, file = ofile)

# all_files <- file.path("../../DataReference/CTCF/", one_only.df$ff)

