#!/usr/bin/env R

# Clusters patients at RNA level, superviser /w DE genes

# load data  ----

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


# metadata for patients
plt.x <- metadata %>%
  tibble::column_to_rownames('genomescan.sid') %>% 
  dplyr::select(methylation.sub.diagnosis,mean.DNA.methylation.signature , Sample_Type, CDKN2AB) %>% 
  
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
  )  %>%  tibble::column_to_rownames('gene_uid')  %>% 
  dplyr::mutate(srt = runif(n())) %>% 
  dplyr::arrange(srt) %>% 
  dplyr::mutate(srt = NULL)


plt.y.h <- dge.partially.paired.h

plt.y.h$labels %in% 



# quantile_breaks <- function(xs, n = 10) {
#   breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
#   breaks[!duplicated(breaks)]
# }
# 
# mat_breaks <- quantile_breaks(as.matrix(plt), n = 11)

# 
# pheatmap::pheatmap(plt, annotation_col = m, clustering_distance_rows = 'correlation'
#                    ,clustering_distance_cols = 'correlation'
#                    )

pheatmap::pheatmap(t(scale(t(plt), center=T, scale=T)),
                   annotation_col = plt.x, 
                   annotation_row = plt.y,
                   cluster_cols = F,
                   cluster_rows = 
                   
                   #color = inferno(length(mat_breaks) - 1)
                   )

#p$tree_col
