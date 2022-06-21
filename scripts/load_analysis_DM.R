#!/usr/bin/env R

# These are the differential methylated regions and probes from the main analysis from Wies

# all these data are hg19(!)


#DMPs_IvR_FDR <- read.csv('data/DM_regions_probes_Wies/DMPs_IvR_FDR 0.05_10052022.csv')

DMRs_IvR_FDR <- read.csv('data/DM_regions_probes_Wies/DMRs_IvR_FDR 0.05_20052022.csv') %>% 
  dplyr::mutate(sign = ifelse(meandiff < 9, -1, 1)) %>% 
  dplyr::mutate(logp = sign * -log10(Fisher)) %>% 
  dplyr::arrange(start, end)
