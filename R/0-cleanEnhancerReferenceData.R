library(liftOver)
library(data.table)
library(tidyverse)

ofile <- "../DataReference/enhancer_promoter_map.hg38.BED11"

# Read in data
enhancer.df <- fread("../DataReference/4-enhancerAnalysis/blood_enhancers_FANTOM.bed3",
                     col.names = c("chr", "start", "end"))

enhancer_promoter.df <- fread("../DataReference/4-enhancerAnalysis/hg19_enhancer_promoter_correlations_distances_cell_type.txt")
chain.filepath <- "../DataReference/hg19ToHg38.over.chain"

# Canonical code to download and use chain if needed. Make this into a gist or something!!!
if (file.exists(chain.filepath)){
  chain <- import.chain(chain.filepath)
} else {
  download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
                destfile = paste0(chain.filepath, ".gz"))
  R.utils::gunzip(paste0(chain.filepath, ".gz"))
}


# Clean enhancer_promoter
# Both are in hg19, make sure we don't look track of the mapping
# We'll split enhancer chr,start,stop and lift
# Do the same with promoter
# Then join them at the end
enhancer_promoter.df$ID <- 1:nrow(enhancer_promoter.df)
head(enhancer_promoter.df)

enhancer.map <- enhancer_promoter.df %>%
  tidyr::separate(col = enhancer, into = c("chr", "start", "end"), sep = "\\:|\\-") %>%
  dplyr::select(c("chr", "start", "end", "ID"))

head(enhancer.map)

# Strand imputation is only needed because of quirk in `separate()`
promoter.map <- enhancer_promoter.df %>%
  tidyr::separate(col = promoter, into = c("chr", "start", "end", "strand"), sep = "\\:|\\-|\\.\\.|\\,") %>%
  dplyr::mutate(strand = ifelse(strand == "+", "+", "-")) %>% # fill missing values (about half are empty strings at this stage)
  dplyr::select(c("chr", "start", "end","strand", "ID"))

head(promoter.map)

############# liftOver ################
enhancer.38 <- makeGRangesFromDataFrame(enhancer.map, keep.extra.columns = T) %>%
  liftOver(chain = chain) %>%
  unlist() %>%
  as.data.frame() %>%
  rename(e.chr = "seqnames",
         e.start = "start",
         e.end = "end",
         e.width = "width")

# Rename columns so we can distingush which start/end we're taling about
promoter.38 <- makeGRangesFromDataFrame(promoter.map, keep.extra.columns = T) %>%
  liftOver(chain = chain) %>%
  unlist() %>%
  as.data.frame() %>%
  rename(p.chr = "seqnames",
         p.start = "start",
         p.end = "end",
         p.width = "width",
         p.strand = "strand")

# Pull it all back together
map.df <- inner_join(enhancer.38, promoter.38, by = "ID") %>%
  dplyr::mutate(p.start = p.start - 1,
                e.start = e.start - 1)

head(map.df)

fwrite(x = map.df, file = ofile, sep= "\t")
