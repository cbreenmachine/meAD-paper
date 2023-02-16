
library(tidyverse)
library(goplot)
library(kableExtra)


IFILE <- "../../DataRaw/2023-01-24-Summaries-v4/DMRegions-pis.bed"
ODIR <- "../../Figs/DMRegions/"

# Read data ---------------------------------------------------------------

df <- read_csv(IFILE, show_col_types = F) %>%
  dplyr::mutate(
    DMR.id = paste0(chr.dmr, ":", start.dmr-1, "-", end.dmr)
    )

ss.df <- read_csv("../../DataRaw/masterSamplesheet.csv", show_col = F) %>%
  dplyr::filter(sample_id %in% names(df)) %>%
  mutate(sample_id = as.character(sample_id))

load.samples <- ss.df$sample_id[ss.df$diagnostic_group == "LOAD"]
control.samples <- ss.df$sample_id[ss.df$diagnostic_group == "CONTROL"]

length(unique(df$DMR.id))

df$pi.diff <-
  rowMeans(df[ , as.character(load.samples)]) -
  rowMeans(df[ , as.character(control.samples)])


data <- df %>%
  dplyr::mutate(direction = sign(pi.diff)) %>%
  group_by(DMR.id, direction) %>%
  summarize(N = sum(direction))

ll <- data %>%
  group_by(DMR.id) %>%
  summarize(Total = sum(abs(N))) %>%
  arrange(-Total) %>%
  pull(DMR.id)

data$DMR.id <- factor(data$DMR.id, levels = ll)
data$direction <- factor(data$direction, levels = c("-1", "1"))

p <- data %>%
  ggplot(aes(x = DMR.id, y = N, fill = direction)) +
    geom_bar(position="stack", stat="identity") +
  xlab("DMR location") +
  ylab("Number of CpGs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 75, vjust = 0.9, hjust=1),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = "white")) +
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF")) +
  ylim(c(-15, 15))

cowplot::save_plot(filename = "../../Figs/Fig2/dmrs-hyper-hypo-bars.png", p)


# Plot boxplots of DMRs ---------------------------------------------------

my.dmr <- unique(df$DMR.id)[[2]]

plot_boxes <- function(xx){

  # Filter and pivot to get in ggplot-friendly form
  data <- df %>%
    dplyr::filter(DMR.id == xx) %>%
    pivot_longer(cols = as.character(c(load.samples, control.samples)),
                 names_to = "sample_id", values_to = "predicted_methylation") %>%
    left_join(dplyr::select(ss.df, sample_id, diagnostic_group), by = "sample_id")

  # "x" coordinate needs to be factor
  data$predicted_methylation
  data$start <- factor(data$start)

  p <- data %>%
    ggplot(aes(x = start, y = predicted_methylation,
               fill = diagnostic_group)) +
    geom_boxplot(outlier.size = 0.5) +
    # ylim(c(0, 0.3)) +
    xlab("Position") +
    ylab("Predicted methylation") +
    labs(fill='Diagnostic group') +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 75, vjust = 0.95, hjust=1),
          # legend.position = c(0.88, 0.8),
          plot.background = element_rect(fill = "white", color = "white")) +
    ggtitle(xx) +
    scale_fill_npg()
  p
}

# Get a list of things to iterate over
uniq.dmrs <- unique(df$DMR.id)

for (xx in uniq.dmrs){
  ofile <- str_replace_all(file.path(ODIR, paste0("dmr-boxplot-", xx, ".png")), ":", "-")
  p <- plot_boxes(xx)
  cowplot::save_plot(filename = ofile, p)
}

# Summary Tables ----------------------------------------------------------

genic.df <- df %>%
  dplyr::filter(dist_to_nearest_gene == 0) %>%
  dplyr::transmute(
    "Gene symbol" = nearest_gene_name,
    "DMR coordinates" = paste0(chr.dmr, ":", start.dmr, "-", end.dmr),
    "Number of CpGs" = nCG) %>%
  distinct() %>%
  arrange(-`Number of CpGs`)


ofile <- file.path(ODIR, "dmrs-in-genes-table.png")
genic.df %>%
  kable() %>%
  kable_styling("striped", full_width = F, font_size = 18) %>%
  save_kable(ofile)


# Out of curiosity, are any really close to a gene?
ofile <- file.path("dmrs-not-in-genes.png")
df %>%
  dplyr::filter(dist_to_nearest_gene != 0) %>%
  dplyr::select(nearest_gene_name, dist_to_nearest_gene) %>%
  distinct() %>%
  kable() %>%
  kable_styling("striped", full_width = F, font_size = 18) %>%
  save_kable(ofile)

# Gene ontology of the (now 11) DMRs in genes -----------------------------
gene.ids <- genic.df$nearest_gene_id

go.out <- goplot::run_gene_ontology(gene.ids)
go.df <- goplot::go_output_to_df(go.out)

p <- goplot::plot_go_barchart(go.df)
cowplot::save_plot(filename = file.path(ODIR, "dmr-gene-ontology.png",
                   p, base_height = 7, base_width = 12)



