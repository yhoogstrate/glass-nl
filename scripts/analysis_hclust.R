#!/usr/bin/env R


# load libs ----

library(recursiveCorPlot)


# load data ----


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}


# if read counts and per-gene metadata are not loaded, load it
if(!exists('expression.glass.exon.metadata')) {
  warning('read count data was not loaded')
  
  source('scripts/load_rna-counts.R')
}


if("padj.partially.paired.exon" %in% colnames(expression.glass.exon.metadata) == F) {
  warning('paired exon-count DGE results were not loaded')
  
  source('scripts/load_analysis_DGE.R')
}



# do hclust ----


tmp.signi <- expression.glass.exon.metadata %>%
  dplyr::filter(abs(log2FoldChange.partially.paired.exon) > 0.75 & padj.partially.paired.exon < 0.01) %>%
  dplyr::pull(gene_uid)

# 
# # these contain interim resections
# tmp.patients <- metadata.glass.per.resection %>%
#  dplyr::filter(excluded == F) %>%
#  dplyr::pull(genomescan.sid)

# these contain those used for DE
tmp.patients <- metadata.glass.per.patient %>%
  dplyr::select(genomescan.sid.I, genomescan.sid.R) %>%
  tidyr::pivot_longer(cols=c('genomescan.sid.I', 'genomescan.sid.R')) %>%
  tidyr::drop_na(value) %>% 
  dplyr::pull('value')

# # these contain those used for DE
# tmp.patients <- metadata.glass.per.patient %>%
#   dplyr::pull(genomescan.sid)


plt <- expression.glass.exon.vst %>% 
  dplyr::select(all_of(tmp.patients)) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::filter(gene_uid %in% tmp.signi) %>% 
  tibble::column_to_rownames('gene_uid')


#rm(tmp.signi, tmp.patients)

labels <- t(plt) %>% as.data.frame %>%  dplyr::select(`ENSG00000233542_AL391845.2`) %>%  dplyr::mutate(`ENSG00000233542_AL391845.2` = "A")

a = recursiveCorPlot::recursiveCorPlot(plt, labels , 7 , 7)


#a[[1]]$data
write.table(a[[2]]$data, file="output/tables/recursiveCorr.txt", sep="\t", quote = F, row.names = F)



#saveRDS(plt, file="/tmp/glass_nl__hclust.Rds")


#h <- recursiveCorPlot::recursiveCorPlot(cp, cpm, 2 ,2)
#ggsave("/tmp/glass-supervised.png",height=20 * 1.3,width=30 * 1.3)
