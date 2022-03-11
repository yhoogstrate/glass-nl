#!/usr/bin/env R 

# load libs ----


# load data ----

if(!exists('expression.glass.vst')) {
  source('scripts/load_rna-counts.R')
}

if(!exists("metadata.glass.per.resection") | "dna.wes.VAF_IDH" %in% names(metadata.glass.per.resection) == F) {
  source('scripts/load_tumour_purities.R')
}


# calc corr for ACE method ----

tmp.metadata <- metadata.glass.per.resection %>%
  dplyr::filter(excluded == F) %>% 
  dplyr::filter(!is.na(dna.shallow.ACE.purity))

tmp.data <- expression.glass.vst %>% 
  dplyr::select(tmp.metadata$genomescan.sid)


stopifnot(tmp.metadata$genomescan.sid == colnames(tmp.data))


cor.test(tmp.metadata$dna.shallow.ACE.purity, as.numeric(tmp.data[1,]))


fun <- function(vec) {
  t <- cor.test(tmp.metadata$dna.shallow.ACE.purity, as.numeric(vec))
  
  return( t$statistic)
  
}

apply(tmp.data[3,],1,fun)


