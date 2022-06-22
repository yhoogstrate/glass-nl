#!/usr/bin/env R




if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if(!exists("expression.glass.exon.vst")) {
  source('scripts/load_rna-counts.R')
}


if(!exists("dge.partially.paired.clusters")) {
  source('scripts/load_hclust.R')
}





# 1. take all expression data ----

# take all samples of sufficient quality, also the R2/R3's not taken into account in the DE analysis
metadata <- metadata.glass.per.resection %>%
  dplyr::filter(excluded == F)


up1 <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid) %>%  # do for all samples w/ enough quality - also the R2/R3 that are not that relevant for DE
  dplyr::filter(rownames(.) %in% (dge.partially.paired.clusters %>% dplyr::filter(up.1 == T) %>%  dplyr::pull(gene_uid)))

up2 <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid) %>%  # do for all samples w/ enough quality - also the R2/R3 that are not that relevant for DE
  dplyr::filter(rownames(.) %in% (dge.partially.paired.clusters %>% dplyr::filter(up.2 == T) %>%  dplyr::pull(gene_uid)))

up3 <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid) %>%  # do for all samples w/ enough quality - also the R2/R3 that are not that relevant for DE
  dplyr::filter(rownames(.) %in% (dge.partially.paired.clusters %>% dplyr::filter(up.3 == T) %>%  dplyr::pull(gene_uid)))

down <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid) %>%  # do for all samples w/ enough quality - also the R2/R3 that are not that relevant for DE
  dplyr::filter(rownames(.) %in% (dge.partially.paired.clusters %>% dplyr::filter(down.1 == T | down.2 == T) %>%  dplyr::pull(gene_uid)))

down.1 <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid) %>%  # do for all samples w/ enough quality - also the R2/R3 that are not that relevant for DE
  dplyr::filter(rownames(.) %in% (dge.partially.paired.clusters %>% dplyr::filter(down.1 == T) %>%  dplyr::pull(gene_uid)))

down.2 <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid) %>%  # do for all samples w/ enough quality - also the R2/R3 that are not that relevant for DE
  dplyr::filter(rownames(.) %in% (dge.partially.paired.clusters %>% dplyr::filter(down.2 == T) %>%  dplyr::pull(gene_uid)))





pc.up1 <- prcomp(t(up1))
c <- cor(t(up1), as.data.frame(pc.up1$x)$PC1)
if(median(c) < 0) {
  pc.up1$x[,1] <- pc.up1$x[,1] * -1
}


pc.up2 <- prcomp(t(up2))
c <- cor(t(up2), as.data.frame(pc.up2$x)$PC1)
if(median(c) < 0) {
  pc.up2$x[,1] <- pc.up2$x[,1] * -1
}



pc.up3 <- prcomp(t(up3))
c <- cor(t(up3), as.data.frame(pc.up3$x)$PC1)
if(median(c) < 0) {
  pc.up3$x[,1] <- pc.up3$x[,1] * -1
}


pc.down <- prcomp(t(down))
c <- cor(t(down), as.data.frame(pc.down$x)$PC1)
if(median(c) < 0) {
  pc.down$x[,1] <- pc.down$x[,1] * -1
}


pc.down.1 <- prcomp(t(down.1))
c <- cor(t(down.1), as.data.frame(pc.down.1$x)$PC1)
if(median(c) < 0) {
  pc.down.1$x[,1] <- pc.down.1$x[,1] * -1
}


pc.down.2 <- prcomp(t(down.2))
c <- cor(t(down.2), as.data.frame(pc.down.2$x)$PC1)
if(median(c) < 0) {
  pc.down.2$x[,1] <- pc.down.2$x[,1] * -1
}

#cor(t(cp1)[,1], as.data.frame(pc1$x)$PC1)
rm(up1, up2, up3, down, down.1, down.2)


transcriptional.signatures <- 
                   pc.up1$x %>% as.data.frame(stringsAsFactor=F) %>% dplyr::select("PC1") %>% dplyr::rename(lts.up1 = PC1) %>% tibble::rownames_to_column('genomescan.sid') %>% 
  dplyr::left_join(pc.up2$x %>% as.data.frame(stringsAsFactor=F) %>% dplyr::select("PC1") %>% dplyr::rename(lts.up2 = PC1) %>% tibble::rownames_to_column('genomescan.sid'), by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(pc.up3$x %>% as.data.frame(stringsAsFactor=F) %>% dplyr::select("PC1") %>% dplyr::rename(lts.up3 = PC1) %>% tibble::rownames_to_column('genomescan.sid'), by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(pc.down$x %>% as.data.frame(stringsAsFactor=F) %>% dplyr::select("PC1") %>% dplyr::rename(lts.down = PC1) %>% tibble::rownames_to_column('genomescan.sid'), by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(pc.down.1$x %>% as.data.frame(stringsAsFactor=F) %>% dplyr::select("PC1") %>% dplyr::rename(lts.down.a = PC1) %>% tibble::rownames_to_column('genomescan.sid'), by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(pc.down.2$x %>% as.data.frame(stringsAsFactor=F) %>% dplyr::select("PC1") %>% dplyr::rename(lts.down.b = PC1) %>% tibble::rownames_to_column('genomescan.sid'), by=c('genomescan.sid'='genomescan.sid')) 


rm(pc.up1, pc.up2, pc.up3, pc.down, pc.down.1, pc.down.2)
rm(c, metadata)


saveRDS(transcriptional.signatures, file="cache/transcriptional.signatures.Rds")




  
