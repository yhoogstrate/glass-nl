#!/usr/bin/env

# load libs ----

# load data ----


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if(!exists("expression.glass.exon.vst")) {
  source('scripts/load_rna-counts.R')
}

# check DGE info
if("padj.partially.paired.exon" %in% colnames(expression.glass.exon.metadata) == F) {
  warning('DGE analysis results were not loaded')
  
  #source('scripts/analysis_DGE.R') 
  source('scripts/load_analysis_DGE.R')
}


# load hclust ----

# full gene_uid labels
tmp <- expression.glass.exon.metadata %>%
  dplyr::filter(abs(log2FoldChange.partially.paired.exon) > 0.75 & padj.partially.paired.exon < 0.01) %>%
  dplyr::select('gene_uid') %>% 
  dplyr::mutate(gene_name = gsub("^.+_","",gene_uid))


dge.partially.paired.h <- readRDS('cache/h.Rds')
dge.partially.paired.h$labels <- gsub("ARHGAP11B.2", "ARHGAP11B",dge.partially.paired.h$labels)
dge.partially.paired.h$labels.uid <- data.frame(gene_name = dge.partially.paired.h$label) %>% 
  dplyr::left_join(tmp, by=c('gene_name' = 'gene_name')) %>% 
  dplyr::pull(gene_uid)

stopifnot(length(dge.partially.paired.h$labels) == length(dge.partially.paired.h$labels.uid))


dge.partially.paired.h.up <- dendextend::prune(dge.partially.paired.h,cutree(dge.partially.paired.h,2) %>% purrr::keep(function(x) x == "1") %>%  names )
dge.partially.paired.h.down <- dendextend::prune(dge.partially.paired.h,cutree(dge.partially.paired.h,2) %>% purrr::keep(function(x) x == "2") %>%  names )


stopifnot(length(dge.partially.paired.h$order) == 604)

#dge.partially.paired.clusters %>% dplyr::filter(grepl("ARHGA",gene_name) | is.na(gene_name) | is.na(gene_uid))


dge.partially.paired.clusters <- 
  rbind(
    data.frame(cluster = cutree(dge.partially.paired.h.down, 2)) %>% 
      dplyr::mutate(cluster = paste0('down.',cluster))
    ,
    data.frame(cluster = cutree(dge.partially.paired.h.up, 3)) %>% 
      dplyr::mutate(cluster = paste0('up.',cluster))
  ) %>% 
  dplyr::mutate(val=T) %>% 
  tibble::rownames_to_column('gene_name') %>% 
  dplyr::left_join(tmp, by=c('gene_name'='gene_name')) %>% 
  dplyr::mutate(marker=T) %>% 
  tidyr::pivot_wider(names_from = cluster,values_from = marker, values_fill=F) %>% 
  dplyr::mutate(val = NULL) %>% 
  dplyr::mutate(down = down.1 | down.2) %>% 
  dplyr::select(gene_name, gene_uid  , up.1 , up.2, up.3, down , down.1, down.2 )


stopifnot(nrow(dge.partially.paired.clusters) == 604)


rm(dge.partially.paired.h.up, dge.partially.paired.h.down)
rm(tmp)


tmp <- data.frame(
  gene_name = dge.partially.paired.h$labels,
  hclust_order = dge.partially.paired.h$order,
  hclust_rank = order(dge.partially.paired.h$order)
)



dge.partially.paired.clusters <- dge.partially.paired.clusters |> 
  dplyr::left_join(tmp, by=c('gene_name'='gene_name'),suffix=c('','')) |> 
  dplyr::arrange(hclust_rank)

rm(tmp)

stopifnot(nrow(dge.partially.paired.clusters) == 604)
stopifnot(is.na(dge.partially.paired.clusters$hclust_rank) == F)


stopifnot(dge.partially.paired.clusters |> head(n=1) |> dplyr::pull(gene_name) == "CHL1-AS2")
stopifnot(dge.partially.paired.clusters |> tail(n=1) |> dplyr::pull(gene_name) == "COL6A2")



