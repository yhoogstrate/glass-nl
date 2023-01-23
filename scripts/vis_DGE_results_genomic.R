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
  dplyr::mutate(significant = padj.partially.paired.exon < 0.01 & abs(log2FoldChange.partially.paired.exon) > 0.75) |> 
  dplyr::filter(gene_chr %in% c('chrX', 'chrY') == F)



ggplot(plt, aes(x=gene_chr_center_loc / 1000000,y=stat.partially.paired.exon,col=gene_chr)) + 
  facet_grid(cols = vars(gene_chr), scales = "free", space="free") +
  geom_point(pch=19,cex=0.2) +
  geom_point(data = subset(plt, significant==T), pch=21,cex=0.8,col='black',fill=NA) +
  geom_smooth(se=F,col="black", lwd=0.7) +
  youri_gg_theme + 
  theme_bw + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  labs(x=NULL)


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


# volcano pim ----


plt <- expression.glass.exon.metadata |> 
  dplyr::mutate(log10.padj = 1 - log10(padj.partially.paired.exon)) |> 
  dplyr::left_join(
    dge.partially.paired.clusters |> dplyr::mutate(gene_name=NULL), by=c('gene_uid'='gene_uid'),suffix=c('','')
  ) |> 
  dplyr::mutate(col = case_when(
    !is.na(up.1) & up.1 == T ~ "up-1 (Cell cycling + histones)",
    !is.na(up.2) & up.2 == T ~ "up-2 (Col/ECM + VEGFA + CD248)",
    !is.na(up.3) & up.3 == T ~ "up-3 (Fuzzy)",
    !is.na(down.1) & down.1 == T ~ "down-1 (AC/AC-like?)",
    !is.na(down.2) & down.2 == T ~ "down-2 (AC/AC-like? + Neuron?)",
    T ~ "-"
  )) |> 
  dplyr::mutate(hist = grepl("^H[0-9]+[A-Z]",gene_name)) |> 
  dplyr::mutate(hox = grepl("^HOX",gene_name))


ggplot(plt, aes(y=log10.padj, x=log2FoldChange.partially.paired.exon,label=gene_name, fill=col)) + 
  geom_point(data = subset(plt, col == "-"),col="gray") +
  geom_hline(yintercept = 1 - log10(0.01), lty=2, lwd=0.5) + 
  geom_vline(xintercept = -0.75, lty=2, lwd=0.5) + 
  geom_vline(xintercept = 0.75, lty=2, lwd=0.5) + 
  geom_point(data = subset(plt, col != "-"),size=2,pch=21) +
  ggrepel::geom_text_repel(data = subset(plt, abs(log2FoldChange.partially.paired.exon) > 3)) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(subtitle="DGE primary -> recurrence GLASS-NL (patient paired analysis)") + 
  xlim(-3,3.5)
ggsave("Rplot1.pdf",width=6*1.3,height=5*1.3)


ggplot(plt, aes(y=log10.padj, x=log2FoldChange.partially.paired.exon,label=gene_name, fill=col)) + 
  geom_point(data = subset(plt, col == "-"),col="gray") +
  geom_hline(yintercept = 1 - log10(0.01), lty=2, lwd=0.5) + 
  geom_vline(xintercept = -0.75, lty=2, lwd=0.5) + 
  geom_vline(xintercept = 0.75, lty=2, lwd=0.5) + 
  geom_point(data = subset(plt, col != "-"),size=2,pch=21) +
  ggrepel::geom_text_repel(data = subset(plt, abs(log2FoldChange.partially.paired.exon) > 3)) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(subtitle="DGE primary -> recurrence GLASS-NL (patient paired analysis)")
ggsave("Rplot2.pdf",width=10*1.3,height=5*1.3)


ggplot(plt, aes(y=log10.padj, x=log2FoldChange.partially.paired.exon,label=gene_name, fill=hist)) + 
  geom_point(data = subset(plt, hist == F),col="gray") +
  geom_hline(yintercept = 1 - log10(0.01), lty=2, lwd=0.5) + 
  geom_vline(xintercept = -0.75, lty=2, lwd=0.5) + 
  geom_vline(xintercept = 0.75, lty=2, lwd=0.5) + 
  geom_point(data = subset(plt, hist == T),size=2,pch=21) +
  ggrepel::geom_text_repel(data = subset(plt, abs(log2FoldChange.partially.paired.exon) > 3)) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(subtitle="DGE primary -> recurrence GLASS-NL (patient paired analysis)",fill="Histon gene") + 
  xlim(-3,3.5)
ggsave("Rplot3.pdf",width=6*1.3,height=5*1.3)


ggplot(plt, aes(y=log10.padj, x=log2FoldChange.partially.paired.exon,label=gene_name, fill=hox)) + 
  geom_point(data = subset(plt, hox == F),col="gray") +
  geom_hline(yintercept = 1 - log10(0.01), lty=2, lwd=0.5) + 
  geom_vline(xintercept = -0.75, lty=2, lwd=0.5) + 
  geom_vline(xintercept = 0.75, lty=2, lwd=0.5) + 
  geom_point(data = subset(plt, hox == T),size=2,pch=21) +
  ggrepel::geom_text_repel(data = subset(plt, hox==T)) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(subtitle="DGE primary -> recurrence GLASS-NL (patient paired analysis)",fill="HOX gene")
  #xlim(-3,3.5)
ggsave("Rplot4.pdf",width=10*1.3,height=5*1.3)



saveRDS(plt, file="plt.glass.Rds")

