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

# 
# h2a = metadata %>%
#   dplyr::select(genomescan.sid , lts.up1) %>%
#   #dplyr::arrange(lts.up1) %>% 
#   #dplyr::mutate(lts.up1 = order(lts.up1)) %>% 
#   dplyr::mutate(lts.up1.rep = lts.up1 + runif(n()) * 0.0001) %>% 
#   tibble::column_to_rownames('genomescan.sid') %>% 
#   t() %>% 
#   as.data.frame
# 
# h2b = metadata %>%
#   dplyr::select(genomescan.sid , lts.up1) %>%
#   #dplyr::arrange(lts.up1) %>% 
#   dplyr::mutate(lts.up1 = -rank(lts.up1)) %>% 
#   dplyr::mutate(lts.up1.rep = lts.up1 + runif(n()) * 0.0001) %>% 
#   tibble::column_to_rownames('genomescan.sid') %>% 
#   t() %>% 
#   as.data.frame


#plt <- plt %>% dplyr::select(colnames(h2a))

# 
# #phe
# h1 = hclust(dist(t(h2a))^2, method="single")
# #h1$order <- 
# plot(h1)
# pheatmap::pheatmap(h2b,  cluster_cols = h1)



m = metadata %>%
  tibble::column_to_rownames('genomescan.sid') %>% 
  dplyr::select(methylation.sub.diagnosis,mean.DNA.methylation.signature , Sample_Type, CDKN2AB) %>% 
  
  dplyr::mutate(Sample_Type = as.character(Sample_Type)) %>% 
  dplyr::mutate(Sample_Type = ifelse(Sample_Type == "X", NA, Sample_Type)) %>% 
  
  dplyr::mutate(methylation.sub.diagnosis = ifelse(methylation.sub.diagnosis %in% c("A_IDH", "A_IDH_HG", "na") == F, "other", methylation.sub.diagnosis)) %>% 
  dplyr::mutate(methylation.sub.diagnosis = as.factor(methylation.sub.diagnosis))


quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(as.matrix(plt), n = 11)


n = data.frame(gene_uid = rownames(plt)) %>% 
  dplyr::left_join(
    dge.partially.paired.clusters %>% dplyr::filter(!duplicated(gene_name)) %>%  dplyr::mutate(gene_name = NULL),
    by=c('gene_uid'='gene_uid')
  )  %>%  tibble::column_to_rownames('gene_uid') %>% 
  dplyr::mutate(up.1 = as.character(up.1)) %>% 
  dplyr::select(up.1)
  
  

pheatmap::pheatmap(plt, annotation_col = m, clustering_distance_rows = 'correlation'
                   ,clustering_distance_cols = 'correlation'
                   )

pheatmap::pheatmap(t(scale(t(plt), center=T, scale=T)),
                   annotation_col = m, 
                   annotation_row = n,
                   cluster_cols = F,
                   cluster_rows = dge.partially.paired.h
                   
                   #color = inferno(length(mat_breaks) - 1)
                   )

#p$tree_col
