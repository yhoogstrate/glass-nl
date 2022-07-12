#!/usr/bin/env R


#' Clusters patients at RNA level, superviser /w DE genes
#' 


# load libs ----

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(viridis)


# load data  ----

if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if('mean.DNA.methylation.signature' %in% colnames(metadata.glass.per.resection) == F) {
  print("Loading methylation analysis")
  source('scripts/load_analysis_DM.R')
}


if('IDH.mutation.WES' %in% colnames(metadata.glass.per.resection) == F) {
  print("Loading genomic alterations and up2date IDH calls")
  source('scripts/load_genomic_alterations.R')
}


if(!exists('dge.partially.paired.clusters')) { # obtain DE genes
  source('scripts/load_hclust.R')
}


# 1. obtain vst counts ----

metadata <- metadata.glass.per.resection %>% 
  dplyr::filter(excluded == F) %>%
  dplyr::arrange(-lts.up1)


# 2.  ----


plt <- expression.glass.exon.vst %>% 
  dplyr::select(metadata$genomescan.sid) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::filter(gene_uid %in% dge.partially.paired.clusters$gene_uid) %>% 
  tibble::column_to_rownames('gene_uid')

colnames(plt) <- data.frame(genomescan.sid = colnames(plt)) %>% 
  dplyr::left_join(
    metadata.glass.per.resection %>% 
      dplyr::select(genomescan.sid, Sample_Name),
    by=c('genomescan.sid'='genomescan.sid')
  ) %>% 
  dplyr::pull(Sample_Name)


# metadata for patients
plt.x <- metadata %>%
  #tibble::column_to_rownames('genomescan.sid') %>% 
  tibble::column_to_rownames('Sample_Name') %>% 
  dplyr::select(lts.up1, mean.DNA.methylation.signature ,methylation.sub.diagnosis, Sample_Type, CDKN2AB, IDH.mutation) %>% 
  
  dplyr::rename(RNA.cell.cycling.signature = lts.up1) %>% 
  
  dplyr::mutate(Sample_Type = as.character(Sample_Type)) %>% 
  dplyr::mutate(Sample_Type = ifelse(Sample_Type == "X", NA, Sample_Type)) %>% 
  
  dplyr::mutate(methylation.sub.diagnosis = ifelse(methylation.sub.diagnosis %in% c("A_IDH", "A_IDH_HG", "na") == F, "other", methylation.sub.diagnosis)) %>% 
  dplyr::mutate(methylation.sub.diagnosis = as.factor(methylation.sub.diagnosis))


# metadata for genes / signatures
plt.y <- data.frame(gene_uid = rownames(plt)) %>% 
  dplyr::left_join(
    dge.partially.paired.clusters %>%
      dplyr::filter(!duplicated(gene_name)) %>%
      dplyr::mutate(gene_name = NULL) %>% 
      dplyr::mutate(cluster = case_when(
        up.1 == T ~ "up.1",
        up.2 == T ~ "up.2",
        up.3 == T ~ "up.3",
        down.1 == T ~ "down.1",
        down.2 == T ~ "down.2",
        T ~ "?"
      )) %>% 
      dplyr::select(gene_uid, cluster),
    by=c('gene_uid'='gene_uid')
  ) %>%  tibble::column_to_rownames('gene_uid')



colnames(plt) == rownames(plt.x)
rownames(plt) == rownames(plt.y)



plt.y.h <- dge.partially.paired.h
plt.y.h$labels <- plt.y.h$labels.uid


# non linear color palette: https://slowkow.com/notes/pheatmap-tutorial/

quantile_breaks <- function(xs, n = 10) {
 breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
 breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(as.matrix(t(scale(t(plt), center=T, scale=T))), n = 21)




pheatmap::pheatmap(t(scale(t(plt), center=T, scale=T)),
                   annotation_col = plt.x, 
                   annotation_row = plt.y,
                   cluster_cols = F,
                   cluster_rows = plt.y.h,
                   
                   color = inferno(length(mat_breaks) - 1),
                   breaks = mat_breaks,
                   labels_row = rep("",nrow(plt)) # gene names are to much detail
)



#'@todo add IDH1/2/R132H
#'@todo add ATRX mut
#'@todo add DAXX meth


ComplexHeatmap::pheatmap(t(scale(t(plt), center=T, scale=T)),
                          annotation_col = plt.x, 
                          annotation_row = plt.y,
                          cluster_cols = F,
                          cluster_rows = plt.y.h,
                          
                          color = inferno(length(mat_breaks) - 1),
                          breaks = mat_breaks,
                          labels_row = rep("",nrow(plt)) # gene names are to much detail
)




