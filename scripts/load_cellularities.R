#!/usr/bin/env R


cnv.cellularities <- read.delim('output/tables/dna-seq/CellularitiesManuallyCurated.xlsx.txt') %>% 
  dplyr::mutate(names = gsub("-","_",names,fixed=T)) %>% 
  dplyr::mutate(Sample_ID = gsub("^(.+_.+)_I.+$","\\1",names)) %>% 
  dplyr::mutate(names = NULL)

