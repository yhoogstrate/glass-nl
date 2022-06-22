#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(recursiveCorPlot)


# load data ----

 
source('scripts/R/chrom_sizes.R')
 


# check DGE info
if("padj.partially.paired.exon" %in% colnames(expression.glass.exon.metadata) == F) {
  warning('DGE analysis results were not loaded')
  
  #source('scripts/analysis_DGE.R') 
  source('scripts/load_analysis_DGE.R')
}



if(!exists('cycling.cell.markers')) {
  warning('cycling cell marker genes were not loaded')
   %>% %>% 