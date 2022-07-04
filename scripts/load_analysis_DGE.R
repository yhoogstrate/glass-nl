#!/usr/bin/env R

## loads data from analysis_DGE.R

stopifnot(file.exists(c("cache/res.paired.a.exon.Rds","cache/res.paired.a.exon.design.Rds")))

# load data ----

# if read counts and per-gene metadata are not loaded, load it
if(!exists('expression.glass.exon.metadata')) {
  warning('read count data was not loaded')
  
  source('scripts/load_rna-counts.R')
}




# load from cache
res.paired.a.exon <- readRDS("cache/res.paired.a.exon.Rds")
res.paired.a.exon.design <- readRDS("cache/res.paired.a.exon.design.Rds")
res.paired.a.covar.regression <- readRDS("cache/res.paired.a.covar.regression.Rds")


# append to gene metadata

# n sign
res.paired.a.exon %>%
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  dplyr::filter(padj < 0.01) %>% 
  dim


nrow(expression.glass.exon.metadata) == nrow(expression.glass.exon)
stopifnot(rownames(expression.glass.exon) == expression.glass.exon.metadata$gene_uid)

# append results
expression.glass.exon.metadata <- expression.glass.exon.metadata %>% 
  dplyr::left_join(
    res.paired.a.exon %>% 
      dplyr::select(gene_uid, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>% 
      tibble::column_to_rownames('gene_uid') %>% 
      `colnames<-`(paste0(colnames(.),".partially.paired.exon")) %>% 
      tibble::rownames_to_column('gene_uid'),
    by=c('gene_uid'='gene_uid'),suffix = c("", "")
  )  %>% 
  dplyr::left_join(
    res.paired.a.covar.regression %>% 
      dplyr::rename(coef.chemotherapy = status.chemo_chemo_vs_no.chemo) %>% 
      dplyr::rename(coef.radiotherapy = status.radio_radio_vs_no.radio) %>% 
      dplyr::rename(coef.grading = status.grading_Recurrent.High.Grade_vs_Recurrent.Low.Grade) %>% 
      dplyr::rename(coef.resection = Sample_Type_recurrent_vs_initial)
      
    , by=c('gene_uid'='gene_uid'),suffix = c("", "")
  )


nrow(expression.glass.exon.metadata) == nrow(expression.glass.exon)
stopifnot(rownames(expression.glass.exon) == expression.glass.exon.metadata$gene_uid)

# expression.glass.exon.metadata$log2FoldChange.partially.paired.exon
rm(res.paired.a.exon)
rm(res.paired.a.covar.regression)

