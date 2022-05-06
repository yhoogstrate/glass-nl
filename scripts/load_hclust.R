#!/usr/bin/env


if(!exists("expression.glass.exon.vst")) {
  source('scripts/load_rna-counts.R')
}



dge.partially.paired.h <- readRDS('cache/h.Rds')
dge.partially.paired.h.up <- dendextend::prune(dge.partially.paired.h,cutree(dge.partially.paired.h,2) %>% purrr::keep(function(x) x == "1") %>%  names )
dge.partially.paired.h.down <- dendextend::prune(dge.partially.paired.h,cutree(dge.partially.paired.h,2) %>% purrr::keep(function(x) x == "2") %>%  names )


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
  dplyr::mutate(marker=T) %>% 
  tidyr::pivot_wider(names_from = cluster,values_from = marker, values_fill=F) %>% 
  dplyr::mutate(val = NULL) %>% 
  dplyr::mutate(gene_name = gsub('^ARHGAP11B.2$','ARHGAP11B',gene_name)) %>% 
  dplyr::mutate(down = down.1 | down.2) %>% 
  dplyr::select(gene_name  , up.1 , up.2, up.3, down , down.1, down.2 ) %>% 
  dplyr::left_join(expression.glass.exon.metadata %>% dplyr::select('gene_uid','gene_name'), by=c('gene_name'='gene_name'))



rm(dge.partially.paired.h.up, dge.partially.paired.h.down)


