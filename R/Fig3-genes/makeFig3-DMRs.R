
suppressPackageStartupMessages({
    library(viridis)
    library(tidyverse)
    library(data.table)
    library(argparse)
    library(cowplot)
    library(ggsci)
    library(ggbio)
    library(EnsDb.Hsapiens.v86)
    library(clusterProfiler)
    library(annotatr)
    library(networkD3)
    library(rbokeh)
    library(waffle)
}) 

source("../functions.R")

parser <- ArgumentParser(description='')
parser$add_argument('--ifile', default = "../../dataDerived/analysis-controlLOAD/test-diagnostic-group-coded/experimentSummary/dmrs.all.bed", help = "input BED file")
parser$add_argument('--odir', default = "../../figs/2022-paper/dmrGenes/", help = "input BED file")
args <- parser$parse_args()

odir <- args$odir
dir.create(odir, showWarn=F, recurs=T)

# Read data and convert to GRanges
data <- fread(args$ifile)
data$class <- "mixed"

hyperix <- (data$n.sig.hyper / data$n.sig) > 0.7
hypoix <- (data$n.sig.hyper / data$n.sig) < 0.3
data$class[hyperix] <- "hyper"
data$class[hypoix] <- "hypo"

data %>% 
  group_by(class) %>% 
  summarize(count = n())

data.gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T, starts.in.df.are.0based = T)
hyper.gr <- makeGRangesFromDataFrame(data[data$class == "hyper", ], keep.extra.columns=T, starts.in.df.are.0based = T)
hypo.gr <- makeGRangesFromDataFrame(data[data$class == "hypo", ], keep.extra.columns=T, starts.in.df.are.0based = T)

seqlevelsStyle(data.gr) <- "UCSC"
seqlevelsStyle(hyper.gr) <- "UCSC"
seqlevelsStyle(hypo.gr) <- "UCSC"


annotate_regions <- function(regions.gr, anno.list='hg38_basicgenes'){

  # hg38_basicgenes or hg38_cpg
  my.anno <- build_annotations(genome = 'hg38', annotations = anno.list)

  # Run the annotation function
  regions.annotated <- annotatr::annotate_regions(regions.gr, my.anno, min=2)

  anno.df <- data.frame(regions.annotated) %>%
    dplyr::transmute(cpg.locus = paste0(seqnames, "-", start), annot.type) %>% # unique id for each CpG
    distinct() %>%  # get rid of instances where one CpG resides in multiple introns
    dplyr::mutate(Annotation = str_remove(str_remove(annot.type, "hg38_genes_"), "hg38_cpg_"))
    
  anno.counts.df <- anno.df %>% 
    group_by(Annotation) %>%
    dplyr::summarize(value = round(100 * n() / nrow(data), 2)) %>%
    dplyr::mutate(Annotation = paste0(str_to_title(Annotation),  " (", value, "%)")) %>%
    column_to_rownames("Annotation")

  anno.counts.df
}

# Annotation counts
ac.df <- annotate_regions(data.gr)
cpg.loc.df <- annotate_regions(data.gr, anno.list=c("hg38_cpg_inter", "hg38_cpg_islands", "hg38_cpg_shelves", "hg38_cpg_shores"))

################################
######### PLOT SANKEY ##########
################################
node.names <- data.frame(name = c(
  paste0("N = ", format(nrow(data), big.mark=","), " DMPs"), 
  rownames(ac.df)))

ac.df$source <- 0
ac.df$target <- 1:nrow(ac.df)

p <- sankeyNetwork(Links = ac.df, 
            Nodes = node.names, 
            Source = "source", Target = "target", 
            Value = "value", NodeID = "name",
            fontSize = 12, nodeWidth = 20,
            )

widget2png(p, file.path(odir, paste0("sankey.png")))



################################
######### PLOT WAFFLE ##########
################################
p <- cpg.loc.df %>%
  rownames_to_column("Location") %>%
  dplyr::mutate(Percent = round(value)) %>%
  dplyr::mutate(Percent = ifelse(Percent > 90, Percent-1, Percent)) %>%
  ggplot(aes(fill = Location, values = Percent)) +
    geom_waffle(nrow = 5, ncol = 20, size = 0.25, color="white") +
    scale_fill_manual(name = NULL, values = pal_locuszoom()(4)) +
  # coord_equal() +
  ggtitle("Annotation of DMPs") +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", colour = "transparent"),
        panel.background = element_rect(fill = "transparent", colour = "transparent"),
        plot.title = element_text(hjust = 0.5, size=20),
        legend.text=element_text(size=16),
        legend.position = "bottom",
        plot.margin = margin(t=5, r=5,b=5,l=5)) +
  guides(fill=guide_legend(ncol=2))


ofile <- file.path(odir, paste0("waffleAnno.png"))
cowplot::save_plot(plot = p, filename=ofile, 
                   base_width = 7, base_height = 5)





######################################
######### ANNOTATE TO GENES ##########
######################################

#TODO: fix out of range warning
expand_genes <- function(gr, upstream, downstream) {
# https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
  strand_is_minus = strand(gr) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  start(gr)[on_plus] = start(gr)[on_plus] - upstream
  start(gr)[on_minus] = start(gr)[on_minus] - downstream
  end(gr)[on_plus] = end(gr)[on_plus] + downstream
  end(gr)[on_minus] = end(gr)[on_minus] + upstream
  return(gr)
}

#--> Ensembl database and filtered for coding genes only
ensdb <- EnsDb.Hsapiens.v86
seqlevelsStyle(ensdb) <- "UCSC"

ensdb.subset <- genes(ensdb, filter = GeneBiotypeFilter('protein_coding'))
ensdb.expanded <- expand_genes(ensdb.subset, upstream=5000, downstream=300)

gr_to_gene_hits <- function(my.gr){

  # Overlap our loci and ENSDB
  overlaps <-findOverlaps(my.gr, ensdb.expanded)
  gr.ix <- queryHits(overlaps)
  ensdb.ix <- subjectHits(overlaps)

  # Pull out gene name, gene id from hits in ensdb
  gene.hits <- data.frame(mcols(ensdb.expanded[ensdb.ix]))
  list("gene_name" = sort(unique(gene.hits$gene_name)),
       "gene_ids" <- sort(unique(gene.hits$gene_id)))
}

all.hits <- gr_to_gene_hits(data.gr)
hyper.hits <- gr_to_gene_hits(hyper.gr)
hypo.hits <- gr_to_gene_hits(hypo.gr)

write.table(all.hits[[1]], file.path(odir, "allGeneHits.txt"), row.names=F, col.names=F, quote=F)
write.table(hyper.hits[[1]], file.path(odir, "hyperGeneHits.txt"), row.names=F, col.names=F, quote=F)
write.table(hypo.hits[[1]], file.path(odir, "hypoGeneHits.txt"), row.names=F, col.names=F, quote=F)


########################################
########## GENE ONTOLOGY ###############
########################################

gene_ids_to_go_df <- function(gene.ids){
  clusterProfiler::enrichGO(gene.ids, OrgDb = 'org.Hs.eg.db', 
    keyType = "ENSEMBL", ont="ALL", pAdjustMethod="fdr" ) %>%
    process_GO_output()
}

all.go.df <- gene_ids_to_go_df(all.hits[[2]])
hyper.go.df <- gene_ids_to_go_df(hyper.hits[[2]])
hypo.go.df <- gene_ids_to_go_df(hypo.hits[[2]])

plot_and_save <- function(go.df, oname){
  p <- plot_GO(go.df)
  ofile <- file.path(odir, paste0(oname,  "-geneOnt.png"))
  cowplot::save_plot(ofile, p, base_width=15.5, base_height=9)
}

plot_and_save(all.go.df, "all")
plot_and_save(hyper.go.df, "hyper")
plot_and_save(hypo.go.df, "hypo")


summary(data$n.sig)
