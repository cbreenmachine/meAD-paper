
library(tidyverse)x

# Email from March 3, 2022
# Hi Ligia,
#
# To follow up on your texts from yesterday, I don't have sex or control/LOAD/MCI split on any of the samplesheets you've given me.
# I'll email Lindsay today or tomorrow about getting the rest of the phenotypes. However, I noticed your numbers add to 379--
# I think you need to remove 203 and 208 from your tally!
#
# Best,
# Coleman

df <- read_csv("../../DataRaw/masterSamplesheet.csv", show = F) %>%
  dplyr::select(c(sample_id, study_id)) %>%
  dplyr::filter(!sample_id %in% c("203", "208")) %>%
  group_by(study_id) %>%
  slice(1) %>%
  arrange(sample_id)

df %>% write_csv("../../DataDerived/2023-02-21-sampleIDs-forLigia.csv")


