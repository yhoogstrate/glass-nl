#!/usr/bin/env R




if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if(!exists("expression.glass.exon.vst")) {
  source('scripts/load_rna-counts.R')
}





# 1. take same DGE input data and (re-)generate hclust object ----


clusters <- data.frame(cluster = cutree(readRDS('cache/h.Rds'), k=5)) %>% 
  dplyr::mutate(cluster = paste0('c',cluster)) %>% 
  dplyr::mutate(val=T) %>% 
  tibble::rownames_to_column('gene_name') %>% 
  dplyr::mutate(marker=T) %>% 
  tidyr::pivot_wider(names_from = cluster,values_from = marker, values_fill=F) %>% 
  dplyr::mutate(gene_name = gsub('^ARHGAP11B.2$','ARHGAP11B',gene_name)) %>% 
  dplyr::left_join(expression.glass.exon.metadata %>% dplyr::select('gene_uid','gene_name'), by=c('gene_name'='gene_name'))




# 2. take all expression data ----

# take all samples of sufficient quality, also the R2/R3's not taken into account in the DE analysis
metadata <- metadata.glass.per.resection %>%
  dplyr::filter(excluded == F)


cp1 <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid) %>%  # do for all samples w/ enough quality - also the R2/R3 that are not that relevant for DE
  dplyr::filter(rownames(.) %in% (clusters %>% dplyr::filter(c1 == T) %>%  dplyr::pull(gene_uid)))

cp2 <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid) %>%  # do for all samples w/ enough quality - also the R2/R3 that are not that relevant for DE
  dplyr::filter(rownames(.) %in% (clusters %>% dplyr::filter(c2 == T) %>%  dplyr::pull(gene_uid)))

cp3 <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid) %>%  # do for all samples w/ enough quality - also the R2/R3 that are not that relevant for DE
  dplyr::filter(rownames(.) %in% (clusters %>% dplyr::filter(c3 == T) %>%  dplyr::pull(gene_uid)))

cp4 <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid) %>%  # do for all samples w/ enough quality - also the R2/R3 that are not that relevant for DE
  dplyr::filter(rownames(.) %in% (clusters %>% dplyr::filter(c4 == T | c5 == T) %>%  dplyr::pull(gene_uid)))

cp4a <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid) %>%  # do for all samples w/ enough quality - also the R2/R3 that are not that relevant for DE
  dplyr::filter(rownames(.) %in% (clusters %>% dplyr::filter(c4 == T) %>%  dplyr::pull(gene_uid)))

cp4b <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid) %>%  # do for all samples w/ enough quality - also the R2/R3 that are not that relevant for DE
  dplyr::filter(rownames(.) %in% (clusters %>% dplyr::filter(c5 == T) %>%  dplyr::pull(gene_uid)))





pc1 <- prcomp(t(cp1))
#plot(pc1)
#plot(pc1$x[,1:2])


pc2 <- prcomp(t(cp2))
#plot(pc2)
#plot(pc2$x[,1:2])


pc3 <- prcomp(t(cp3))
#plot(pc3)
#plot(pc3$x[,1:2])


pc4 <- prcomp(t(cp4))
#plot(pc4)
#plot(pc4$x[,1:2])


pc4a <- prcomp(t(cp4a))
#plot(pc4a)
#plot(pc4a$x[,1:2])


pc4b <- prcomp(t(cp4b))
#plot(pc4b)
#plot(pc4b$x[,1:2])


rm(cp1, cp2, cp3, cp4, cp4a, cp4b)


transcriptional.signatures <- 
  pc1$x %>% as.data.frame(stringsAsFactor=F) %>% dplyr::select("PC1") %>% dplyr::rename(lts1 = PC1) %>% tibble::rownames_to_column('genomescan.sid') %>% 
  dplyr::left_join(pc2$x %>% as.data.frame(stringsAsFactor=F) %>% dplyr::select("PC1") %>% dplyr::rename(lts2 = PC1) %>% tibble::rownames_to_column('genomescan.sid'), by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(pc3$x %>% as.data.frame(stringsAsFactor=F) %>% dplyr::select("PC1") %>% dplyr::rename(lts3 = PC1) %>% tibble::rownames_to_column('genomescan.sid'), by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(pc4$x %>% as.data.frame(stringsAsFactor=F) %>% dplyr::select("PC1") %>% dplyr::rename(lts4 = PC1) %>% tibble::rownames_to_column('genomescan.sid'), by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(pc4a$x %>% as.data.frame(stringsAsFactor=F) %>% dplyr::select("PC1") %>% dplyr::rename(lts4a = PC1) %>% tibble::rownames_to_column('genomescan.sid'), by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(pc4b$x %>% as.data.frame(stringsAsFactor=F) %>% dplyr::select("PC1") %>% dplyr::rename(lts4b = PC1) %>% tibble::rownames_to_column('genomescan.sid'), by=c('genomescan.sid'='genomescan.sid')) 


rm(pc1, pc2, pc3, pc4, pc4a, pc4b)
rm(clusters,metadata)


saveRDS(transcriptional.signatures, file="cache/transcriptional.signatures.Rds")




  