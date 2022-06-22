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



plt <- metadata.glass.per.resection %>% 
  dplyr::select(Sample_Name, resection, methylation.purity.absolute, `dna.shallow.ACE.purity` , `dna.shallow.ACE.purity.below.1`, dna.wes.VAF_IDH ,dna.purity.manual.Erik  ) %>% 
  dplyr::mutate(dna.wes.VAF_IDH = dna.wes.VAF_IDH *2) %>% 
  tidyr::pivot_longer(cols = -c(Sample_Name, resection)) %>% 
  dplyr::rename(Purity.method = name) %>% 
  dplyr::mutate(reliability = ifelse(Purity.method %in% c('methylation.purity.absolute','dna.wes.VAF_IDH'),"most reliable","not most reliable")) %>% 
  dplyr::filter(Purity.method == "dna.wes.VAF_IDH")


ggplot(plt, aes(x=resection, y=value, fill = resection)) +
  facet_grid(cols = vars(Purity.method),  scales = "free", space="free") +
  geom_violin(alpha=0.1,draw_quantiles = c(0.6), col="black") + 
  ggbeeswarm::geom_beeswarm(pch=21,cex=1.2) +
  youri_gg_theme +
  ylim(0,1.4) +
  ggsignif::geom_signif( comparisons = list(c("S1" , "S2"),c("S2" , "S3"),c("S1" , "S3")), y_position=c(1.125,1.2,1.325),    test="wilcox.test", col="black" )+
  labs(y="fraction tumor cells")






