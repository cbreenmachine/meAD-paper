# make-demographics-table.R
# Spits out latex to use in a demographics table

library(tidyverse)
library(kableExtra)

# Does not matter which one
OFILE <- "../../DataDerived/DemographicsTable.tex"
dmrs <- read_table("../../DataRaw/2023-01-24-Summaries-v4/DMRegions-CpGs.bed", show_col_types = F)

valid.ids <- as.numeric(names(dmrs))
valid.ids <- valid.ids[!is.na(valid.ids)]

df <- read_csv("../../DataRaw/masterSamplesheet.csv", show_col = F) %>%
  dplyr::filter(sample_id %in% valid.ids) %>%
  group_by(sample_id) %>%
  slice(1)

# Helper function to get values from numeric
get_var <- function(vv, group){
  x <- df[df$diagnostic_group == group, ][[vv]]
  x[!is.na(x)]
}


vec_to_string <- function(x){
  dd <- 1

  # Controls summaries
  x.bar <- round(mean(x, na.rm = T), dd)
  x.sd <- round(sd(x), dd)

  paste0(x.bar, " (", x.sd, ")")
}

p_to_string <- function(p){
  if (p < 0.01) {
    p.string <- paste0("$p<", 10 ^ ceiling(log10(p)), "$")
  } else {
    p.string <- paste0("$p=", round(p, 2), "$")
  }
  p.string
}

summarize_numeric <- function(vv, clean.vv){

  x <- get_var(vv, "CONTROL")
  y <- get_var(vv, "LOAD")

  test.out <- t.test(x, y)

  stat.str <- paste0("$t=", round(test.out$statistic, 2), "$")
  p.str <- p_to_string(test.out$p.value)
  last.str <- paste0(stat.str, "; ", p.str)

  data.frame(Phenotype = clean.vv,
             Control = vec_to_string(x),
             LOAD = vec_to_string(y),
             test = last.str)
}


# Proof of concept --------------------------------------------------------
# Chi-squared test

tt <- df %>%
  group_by(diagnostic_group, .["sex"]) %>%
  summarize(Count = n()) %>%
  pivot_wider(names_from = "sex", values_from = "Count") %>%
  column_to_rownames("diagnostic_group")


summarize_categorical <- function(vv){

  tt <- df %>%
    group_by(diagnostic_group, .[vv]) %>%
    summarize(Count = n()) %>%
    pivot_wider(names_from = vv, values_from = "Count", values_fill = 0) %>%
    column_to_rownames("diagnostic_group") %>%
    t() %>%
    as.data.frame()

  tt

  # Run Chi-squared test and clean the result
  test.out <- fisher.test(tt)

  # stat.str <- paste0("X=", round(test.out$statistic, 2))
  p.str <- p_to_string(test.out$p.value)
  # last.str <- paste0(stat.str, "; ", p.str)
  last.str <- p.str

  out <- tt %>%
    rownames_to_column(vv)
  out$test <- c(last.str, rep("", nrow(tt) - 1))

  names(out)[1] <- "Phenotype"
  names(out)[2] <- "Control"
  rbind(c("", NA, NA, ""), out)
}

opts <- options(knitr.kable.NA = "")

output <-
  rbind(
    summarize_numeric("age_at_visit", "Age (yrs)"),
    summarize_numeric("education", "Education (yrs)"),
    summarize_numeric("bmi", "Body Mass Index"),
    summarize_categorical("source"),
    summarize_categorical("race_primary"),
    summarize_categorical("sex")
  ) %>%
  kableExtra::kable(format = 'latex', booktabs = T, escape = F)

writeLines(output, OFILE)
