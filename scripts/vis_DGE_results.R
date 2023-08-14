#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(recursiveCorPlot)


# load data ----


source('scripts/R/chrom_sizes.R')
source('scripts/R/chrom_sizes.R')
 


# check DGE info
if("padj.partially.paired.exon" %in% colnames(expression.glass.exon.metadata) == F) {
  warning('DGE analysis results were not loaded')
  
  #source('scripts/analysis_DGE.R') 
  source('scripts/load_analysis_DGE.R')
}



#if(!exists('cycling.cell.markers')) {
#  warning('cycling cell marker genes were not loaded')
#   %>% %>% 



## plot cycling genes ----


source('scripts/load_hclust.R')



# plt <-  dge.partially.paired.clusters |> 
#   dplyr::mutate(x = rev(hclust_rank)) |> 
#   dplyr::arrange(x) |> 
#   dplyr::left_join(
#     expression.glass.exon.metadata,
#     by=c('gene_uid'='gene_uid'), suffix=c('','')
#   ) |> 
#   dplyr::mutate(`chr6 histone locus`= ifelse(
#     (chr.hg19 == "chr6" | gene_chr == "chr6") & grepl("^H[0-9]+[A-Z]",gene_name)
#     ,T,F)) |> 
#   dplyr::select(gene_name, gene_uid, x, 
#                 G1.S.tirosh,
#                 G2.M.tirosh,
#                 cycling.melanoma.tirosh,
#                 G1.S.neftel,
#                 G2.M.neftel,
#                 `chr6 histonel locus`
#                 )



plt <- readRDS(file="/tmp/vis_DGE_results__cycling_genes__plt.Rds")


plt.long <- plt |>
  tidyr::pivot_longer(cols = -c(gene_name, gene_uid, x), values_to = "col") |>
  dplyr::mutate(col = ifelse(!is.na(col) & col != FALSE, "y","n"))


ggplot(plt.long, aes(x=x, y=name, col = col, alpha=col)) +
  geom_point(pch="|",size=4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(c(0,604)) +
  scale_color_manual(values=c('y'='red','n'='gray80')) +
  scale_alpha_manual(values=c('y'=1,'n'=0.25))

ggsave("output/figures/vis_DGE_results__cycling_genes.pdf",width=8.5,height=2)
ggsave("output/figures/vis_DGE_results__cycling_genes.svg",width=8.5,height=2)






