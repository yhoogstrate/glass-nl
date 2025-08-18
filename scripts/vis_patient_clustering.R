#!/usr/bin/env R


#' Clusters patients at RNA level, superviser /w DE genes
#' 


# load libs ----

library(tidyverse)
#library(pheatmap)
#library(RColorBrewer)
#library(viridis)


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
  dplyr::arrange(-lts.up1) # %>%
  #dplyr::arrange(mean.DNA.methylation.signature)


## add cnv metadata ----

metadata.cnv <- cnv2 %>% 
  tibble::rownames_to_column('segment.id') %>% 
  dplyr::mutate(chr = gsub("^([^:]+).+$","\\1",segment.id)) %>% 
  dplyr::mutate(start = as.numeric(gsub("^[^:]+:([^\\-]+)\\-.+$","\\1",segment.id))) %>% 
  dplyr::mutate(end = as.numeric(gsub("^[^:]+:[^\\-]+\\-(.+)$","\\1",segment.id)))


# MSH4 // loss // chr == "chr1" &  start < 76320465 & end > 76320504 // chr1:76300001-76400000
metadata.cnv %>% 
  dplyr::filter(chr == "chr1" &  start < 76320465 & end > 76320504) %>% 
  dplyr::mutate(label = "MSH4 (loss)")


# PDGFRA // gain // chr == "chr4" &  start < 55161200 & end > 55161239
metadata.cnv %>% 
  dplyr::filter(chr == "chr4" &  start < 55161200 & end > 55161239) %>% 
  dplyr::mutate(label = "PDGFRA (gain)")


# MASTL/RAB18/MXK // gain // chr == "chr10" &  start < 27830903 & end > 27831192
metadata.cnv %>% 
  dplyr::filter(chr == "chr10" &  start < 27830903 & end > 27831192) %>% 
  dplyr::mutate(label = "")


# CDK4 // gain // chr == "chr12" &  start < 58142912 & end > 58143373
metadata.cnv %>% 
  dplyr::filter(chr == "chr12" &  start < 58142912 & end > 58143373) %>% 
  dplyr::mutate(label = "CDK4 (gain)")


# CCNE1 // gain // chr == "chr19" &  start < 39502161 & end > 39502733
metadata.cnv %>% 
  dplyr::filter(chr == "chr19" & start < 30308415 & end > 30308454) %>% 
  dplyr::mutate(label = "CCNE1 (gain)")

# RB1 // loss // chr == "chr13" &  start <= 48474808 & end >= 48476242
metadata.cnv %>% 
  dplyr::filter(chr == "chr13" &  start <= 48474808 & end >= 48476242) %>% 
  dplyr::mutate(label = "")

# PTEN // loss // chr == "chr10" &  start < 87863980 & end > 87865285
metadata.cnv %>% 
  dplyr::filter(chr == "chr10" &  start < 87863980 & end > 87865285) %>% 
  dplyr::mutate(label = "")

# TCFL2 // loss // chr == "chr10" &  start < 113033733 & end > 113076506
metadata.cnv %>% 
  dplyr::filter(chr == "chr10" &  start < 113033733 & end > 113076506) %>% 
  dplyr::mutate(label = "")

# MYC // gain // chr == "chr8" &  start < 128752748 & end > 128752787
metadata.cnv %>% 
  dplyr::filter(chr == "chr8" &  start < 128752748 & end > 128752787) %>% 
  dplyr::mutate(label = "")

# SMARCA2 // loss // chr == "chr9" &  start < 2137095 & end > 2138078
metadata.cnv %>% 
  dplyr::filter(chr == "chr9" &  start < 2137095 & end > 2138078) %>% 
  dplyr::mutate(label = "")

## add tp53 / ATRX mut data ----


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



stopifnot(colnames(plt) == rownames(plt.x))
stopifnot(rownames(plt) == rownames(plt.y))



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




