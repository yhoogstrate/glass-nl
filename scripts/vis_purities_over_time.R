#!/usr/bin/env R

# load libs ----

library(tidyverse)
source('scripts/R/job_gg_theme.R')
source('scripts/R/youri_gg_theme.R')


# load data ----


source("scripts/load_tumour_purities.R")


# fig ----


plt <- metadata.glass.per.resection %>% 
  dplyr::select(Sample_Name, resection, methylation.purity.absolute, `dna.shallow.ACE.purity` , `dna.shallow.ACE.purity.below.1`, dna.wes.VAF_IDH ,dna.purity.manual.Erik  ) %>% 
  dplyr::mutate(dna.wes.VAF_IDH = dna.wes.VAF_IDH *2) %>% 
  tidyr::pivot_longer(cols = -c(Sample_Name, resection)) %>% 
  dplyr::rename(Purity.method = name) %>% 
  dplyr::mutate(reliability = ifelse(Purity.method %in% c('methylation.purity.absolute','dna.wes.VAF_IDH'),"most reliable","not most reliable"))


ggplot(plt, aes(x=resection, y=value, fill = resection)) +
  facet_grid(cols = vars(Purity.method),  scales = "free", space="free") +
  geom_violin(alpha=0.1) + 
  ggbeeswarm::geom_beeswarm(pch=21,cex=1.2) +
  youri_gg_theme





