# BEGIN
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(ComplexHeatmap)
    library(viridis)
    library(circlize)
    library(ggsci)
})

ODIR <- "../../Figs/Fig3/"
SAMPLES.PATH <- "../../DataRaw/masterSamplesheet.csv"
PIS.PATH <- "../../DataRaw/2022-12-28-ExperimentSummary-v1/DMR.pis.bed"
N.COL <- 100 # number of breaks in color map


# Read data ---------------------------------------------------------------
samples.df <- read_csv(SAMPLES.PATH, show_col_types = F)

#Load and filter DMRs
data <- fread(PIS.PATH) %>%
  unite("dmr.id", all_of(c("chr.dmr", "start.dmr")), sep = "-")


# We'll pull out the sample pis,
# but need to know which columns are from samples
# versus what describes the DMR
sample.cols <- as.character(as.numeric(names(data)))
sample.cols <- sample.cols[!is.na(sample.cols)]

# Control vs LOAD status
ix <- match(sample.cols, samples.df$sample_id)
load.status <- substr(samples.df$diagnostic_group[ix], 1, 1)

prepare_data <- function(data, dmr.name){
  # Filter and grab the right columns
  X <- data %>%
    dplyr::filter(dmr.id == dmr.name) %>%
    dplyr::select(all_of(sample.cols)) %>%
    as.matrix() %>%
    t()
  rownames(X) <- NULL
  X
}



make_annotation <- function(load.status){
    HeatmapAnnotation(
        DG = load.status,
        which = "row",
        col = list(DG = c("C"="#868686FF", "L"="#8843F2"))
    )
}


make_title <- function(data, dmr.name){
  dplyr::filter(data, dmr.id == dmr.name)$nearest_gene_name[1]
}

# Colors need to be set manually
# colors <- circlize::colorRamp2(breaks = 0:N.COL / N.COL, colors = magma(N.COL+1))

# Interpolate red - white - blue
color.fun <- circlize::colorRamp2(c(0, 0.5, 1), colors = c("#377EB8", "white", "#E41A1C"))
colors <- color.fun( 0:N.COL / N.COL)

plot_hm_matrix <- function(X, load.status, ti){
    load.annotation <- make_annotation(load.status)
    n.cpgs <- ncol(X)

    # Space out CpGs if there aren't too many
    if (n.cpgs > 12){
        n.cpgs <- 1
    }


    ComplexHeatmap::Heatmap(
        X,
        col = color.fun,
        column_order = 1:ncol(X), #turn off clustering on CpGs
        heatmap_legend_param =
          list(title = "Pi",
               at = c(0, 0.25, 0.5, 0.75, 1)
               ),
        left_annotation = load.annotation,
        column_title = ti,
        row_km = 2, # KNN for rows (samples)
        column_split = 1:n.cpgs,
        border = F
      )

}


pipeline <- function(dmr.id){
  print(dmr.id)
  X <- prepare_data(data, dmr.id)
  ti <- make_title(data, dmr.id)
  plot_hm_matrix(X, load.status, ti)
}





wrapper <- function(x){
  # Pulls out gene name to name the output file appropriately
  # Wraps `pipeline()` and saves png
  dmr.name <- x
  gene.name <- unique(dplyr::filter(data, dmr.id == dmr.name)$nearest_gene_name)
  ofile <- file.path(ODIR, paste0(dmr.name, "-", gene.name, ".png"))
  print(ofile)

  # Old school save
  png(ofile, width = 8, height = 8, unit="in", res = 150)
  print(pipeline(dmr.name))
  dev.off()

}

dmr.range <- unique(data$dmr.id)
lapply(X = dmr.range, FUN = wrapper)

#END
