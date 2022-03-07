#!/usr/bin/env R


if(!exists("metadata.glass.per.resection")) {
  warning('metadata was not loaded')
  
  source('scripts/load_metadata.R')
}


expression.glass <- read.delim('data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T) %>% 
  `colnames<-`(gsub("^X.+new.","",colnames(.))) %>% 
  `colnames<-`(gsub(".Aligned.+bam$","",colnames(.))) %>% 
  `colnames<-`(gsub(".","-",colnames(.),fixed=T)) 


stopifnot(metadata.glass.per.resection$genomescan.sid %in% colnames(expression.glass)) # all metadata included samples must exist expression data






