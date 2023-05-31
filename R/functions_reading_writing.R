

# General reading/casting functions -----------------------------------------------

get_natgen_genes <- function(ifile, exclude_APP = T, exclude_IGH = T){
  zz <- read_table(ifile, col_names = FALSE)$X1

  if (exclude_APP){
    zz <- zz[!(zz == "APP")]
  }

  if (exclude_IGH){
    zz <- zz[!stringr::str_detect(zz, "IGH")]
  }
  zz

}


my_write_csv <- function(data, file){
  write_csv(data, file = file)
  return(file)
}


get_pvals_data <- function(ifile){
  # READ pvals from
  data.table::fread(ifile, verbose = F) %>%
    dplyr::mutate(lfdr = lfdr.from.ss,
                  pval = p.from.ss,
                  pval.Wald = p.from.DSS,
                  y = -log10(lfdr)) %>%
    dplyr::select(-c(p.from.ss, p.from.zz, p.from.DSS,
                     lfdr.from.ss, lfdr.from.zz))
}

read_and_cast_madrid_data <- function(file){
  read_csv(file, show_col_types = F, skip = 1) %>%
    separate(hg38_coordinates, into = c("chr", "start")) %>%
    dplyr::mutate(start = as.numeric(start)) %>%
    dplyr::mutate(start = start - 1,
                  end = start + 2,
                  strand = Strand) %>%
    drop_na() %>%
    dplyr::transmute(chr, start, end,
                     gene.symbols = `Gene Symbol(s)`,
                     effect.size = -1 * `Mean Beta Difference (Control-AD)`) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T,
                             ignore.strand = T,
                             seqnames.field = "chr")
}

