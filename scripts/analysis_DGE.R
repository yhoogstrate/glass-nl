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


# 


dds <- DESeqDataSetFromMatrix(countData = expression.glass,
                              colData = metadata.glass.per.resection %>%  dplyr::filter(excluded == F) %>%  dplyr::select(genomescan.sid, 
                                                                                                                          Sample_Sex,Sample_Type),
                              design= ~ Sample_Type)
#                              design= ~ GLASS_ID + Sample_Type)
#                            design= ~ Sample_Sex)


dds <- DESeq(dds)
res <- results(dds) %>% 
  as.data.frame() %>% 
  dplyr::filter(!is.na(padj)) %>% 
  dplyr::arrange(pvalue,padj)
dim(res %>%  dplyr::filter(padj < 0.01))

a = res %>% head(n=50) %>%  rownames() # geeft iets aan signaal, lijkt

View(res)


b <- a %>%
  as.data.frame %>%
  dplyr::rename('gene_uid' = '.') %>% 
  dplyr::left_join(expression.glass.metadata, by=c('gene_uid'='gene_uid'))

View(b)



