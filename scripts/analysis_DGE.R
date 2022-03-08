#!/usr/bin/env R


# load libs ----


source('scripts/R/youri_gg_theme.R')
library(DESeq2)



# load data ----


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if(!exists("expression.glass.vst")) {
  source('scripts/load_rna-counts.R')
}


# 1. unpaired ----



tmp.data <- expression.glass %>% 
  dplyr::select(all_of(
  c(metadata.glass.per.patient$genomescan.sid.I, metadata.glass.per.patient$genomescan.sid.R) %>% 
    purrr::keep(~ !is.na(.))))


tmp.metadata <- data.frame(genomescan.sid = colnames(tmp.data) ) %>% 
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid')) %>% 
  dplyr::select(genomescan.sid, Sample_Sex,Sample_Type)


stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)


dds <- DESeqDataSetFromMatrix(countData = tmp.data,
                              colData = tmp.metadata,
                              design= ~ Sample_Type)



dds <- DESeq(dds)
res <- results(dds) %>% 
  as.data.frame() %>% 
  dplyr::filter(!is.na(padj)) %>% 
  dplyr::arrange(pvalue,padj)


head(res)


dim(res %>%  dplyr::filter(padj < 0.01))


a = res %>% head(n=50) %>%  rownames() # geeft iets aan signaal, lijkt

View(res)


b <- a %>%
  as.data.frame %>%
  dplyr::rename('gene_uid' = '.') %>% 
  dplyr::left_join(expression.glass.metadata, by=c('gene_uid'='gene_uid'))

View(b)



