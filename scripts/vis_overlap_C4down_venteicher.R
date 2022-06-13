#!/usr/bin/env R

# load data ----

source('scripts/load_hclust.R')
source('scripts/R/lgg.transcriptional.programs.Venteicher.R')

# vis1 ----


intersect(dge.partially.paired.clusters %>%  dplyr::filter(down.1 | down.2) %>% dplyr::pull(gene_name),
          lgg.transcr.prog.astro)

intersect(dge.partially.paired.clusters %>%  dplyr::filter(down.1 | down.2) %>% dplyr::pull(gene_name),
          lgg.transcr.prog.oligo)

intersect(dge.partially.paired.clusters %>%  dplyr::filter(down.1 | down.2) %>% dplyr::pull(gene_name),
          lgg.transcr.prog.stemn)


# 



metadata <- readRDS("/tmp/glass-metadata-per-resection.Rds")
expr.data <- readRDS("/tmp/glass-readcounts-vst-transformed.Rds")

stopifnot(metadata$genomescan.sid == colnames(expr.data))


tmp <- dge.partially.paired.clusters %>%  dplyr::filter(down.1 | down.2) %>% dplyr::pull(gene_name)

plt <- expr.data %>%
  tibble::rownames_to_column('gene_symbol') %>%
  dplyr::filter(gsub("^[^_]+_","",gene_symbol) %in% c("ID3","DOCK7","ATP1A2","AGT",tmp)) %>% 
  dplyr::mutate(astr = gsub("^[^_]+_","",gene_symbol) %in% c("ID3","DOCK7","ATP1A2","AGT") ) %>%
  dplyr::arrange(-astr) %>%
  dplyr::mutate(astr=NULL) %>%
  head(n=135) %>%
  #tibble::column_to_rownames('gene_symbol') %>%
  dplyr::mutate(gene_symbol = NULL) %>%
  t() %>%
  cor


corrplot::corrplot(plt,tl.cex=0.75)
