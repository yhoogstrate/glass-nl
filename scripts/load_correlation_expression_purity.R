#!/usr/bin/env R



# if read counts and per-gene metadata are not loaded, load it
if(!exists('expression.glass.exon.metadata')) {
  warning('read count data was not loaded, but is now')
  
  source('scripts/load_rna-counts.R')
}


stopifnot(nrow(expression.glass.exon.metadata) == 22346)

expression.glass.exon.metadata <- expression.glass.exon.metadata |> 
  dplyr::left_join(readRDS("cache/correlation_expression_purity.Rds"), by=c('gene_uid'='gene_uid'), suffix = c("", "")) # keep=F

stopifnot(nrow(expression.glass.exon.metadata) == 22346)



