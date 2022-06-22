#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(DESeq2)



# load data ----

# if metadata was not loaded, load it
if(!exists("metadata.glass.per.resection")) {
  warning('metadata was not loaded')
  
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


if(!exists('chrs_hg38_s')) {
  source('scripts/R/chrom_sizes.R')
}


source('scripts/R/youri_gg_theme.R')



# chromosome plot [all DGE] ----

plt <- expression.glass.exon.metadata %>% 
  dplyr::left_join(chrs_hg38_s, by=c('gene_chr'='chr')) %>% 
  dplyr::mutate(x = gene_chr_center_loc + pos) %>% 
  dplyr::mutate(gene_chr = factor(gene_chr, levels=gtools::mixedsort(unique(as.character(gene_chr))) )) %>% 
  dplyr::mutate(significant = padj.partially.paired.exon < 0.01 & abs(log2FoldChange.partially.paired.exon) > 0.75)



ggplot(plt, aes(x=gene_chr_center_loc / 1000000,y=stat.partially.paired.exon,col=gene_chr)) + 
  facet_grid(cols = vars(gene_chr), scales = "free", space="free") +
  geom_point(pch=19,cex=0.2) +
  geom_point(data = subset(plt, significant==T), pch=21,cex=0.8,col='black',fill=NA) +
  geom_smooth(se=F,col="black", lwd=0.7) +
  youri_gg_theme + labs(x=NULL)


## export for circos plot Wies

write.csv(
  plt %>%
    dplyr::select(gene_id, gene_name, gene_chr, V4, V5, gene_strand, gene_chr_center_loc, gene_loc, contains('partially.paired.exon'), significant) %>% 
    dplyr::rename(gene_id_gencode34 = gene_id) %>% 
    dplyr::rename(gene_name_gencode34 = gene_name) %>% 
    dplyr::rename(gene_chr_hg38 = gene_chr) %>% 
    dplyr::rename(gene_start_hg38 = V4) %>% 
    dplyr::rename(gene_end_hg38 = V5) %>% 
    dplyr::rename(gene_strand_hg38 = gene_strand) %>% 
    dplyr::rename(gene_center_hg38 = gene_chr_center_loc) %>% 
    dplyr::rename(gene_center_text_hg38 = gene_loc)
  , "tmp/2022-05-27_glass-nl_dge_results_paired.csv")


write.table(
  plt %>%
    dplyr::select(gene_chr, V4, V5, gene_id) %>% 
    dplyr::rename(gene_id_gencode34 = gene_id) %>% 
    dplyr::rename(gene_chr_hg38 = gene_chr) %>% 
    dplyr::rename(gene_start_hg38 = V4) %>% 
    dplyr::rename(gene_end_hg38 = V5)
  , "tmp/2022-05-27_glass-nl_dge_results_paired.bed",sep=" ",quote=F,col.names = F,
  row.names=F)

# ^ file kan zo in ucsc liftover



# chromosome plot [all DGE + hist & cycling] ----

plt <- expression.glass.exon.metadata %>% 
  dplyr::left_join(chrs_hg38_s, by=c('gene_chr'='chr')) %>% 
  dplyr::mutate(x = gene_chr_center_loc + pos) %>% 
  dplyr::mutate(gene_chr = factor(gene_chr, levels=gtools::mixedsort(unique(as.character(gene_chr))) )) %>% 
  dplyr::mutate(significant = padj.partially.paired.exon < 0.01 & abs(log2FoldChange.partially.paired.exon) > 0.75)



ggplot(plt, aes(x=gene_chr_center_loc / 1000000,y=stat.partially.paired.exon,col=gene_chr)) + 
  facet_grid(cols = vars(gene_chr), scales = "free", space="free") +
  geom_point(pch=19,cex=0.2) +
  geom_point(data = subset(plt, significant==T), pch=21,cex=0.8,col='black',fill=NA) +
  geom_smooth(se=F,col="black", lwd=0.7) +
  youri_gg_theme +
  labs(x=NULL)




