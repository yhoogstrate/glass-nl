#!/usr/bin/env R


# load libs ----


source('scripts/R/youri_gg_theme.R')
library(DESeq2)
library(EnhancedVolcano)


# load data ----


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if(!exists("expression.glass.vst")) {
  source('scripts/load_rna-counts.R')
}


# 1. unpaired ----


tmp.metadata <- metadata.glass.per.patient %>%
  dplyr::select(genomescan.sid.I, genomescan.sid.R) %>%
  tidyr::pivot_longer(cols=c('genomescan.sid.I', 'genomescan.sid.R')) %>%
  tidyr::drop_na(value) %>%
  dplyr::mutate(name=NULL) %>%
  dplyr::rename(genomescan.sid = value) %>%
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(metadata.glass.per.patient %>% dplyr::select(GLASS_ID, patient.correction.id), by=c('GLASS_ID'='GLASS_ID')) %>% 
  dplyr::arrange(Sample_Type, genomescan.sid)


tmp.data <- expression.glass %>%
  dplyr::select(all_of( tmp.metadata.2$genomescan.sid ))



stopifnot(colnames(tmp.data.2) == tmp.metadata.2$genomescan.sid)



dds.unpaired.a <- DESeqDataSetFromMatrix(countData = tmp.data,
                              colData = tmp.metadata,
                              design= ~ Sample_Type)



dds.unpaired.a <- DESeq(dds.unpaired.a)
res.unpaired.a <- results(dds.unpaired.a) %>% 
  as.data.frame(stringsAsFactors=F) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::filter(!is.na(padj)) %>% 
  dplyr::arrange(pvalue,padj) %>% 
  dplyr::left_join(expression.glass.metadata %>% dplyr::select(gene_uid, gene_name, gene_type, gene_strand, gene_loc),by=c('gene_uid'='gene_uid'))



res.unpaired.a %>%
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  dplyr::filter(padj < 0.01) %>% 
  dim


#saveRDS(res.unpaired.a, "cache/res.unpaired.a.Rds")
res.unpaired.a <- readRDS("cache/res.unpaired.a.Rds")




# 2a. paired [incomplete as separate group] ----


tmp.metadata <- metadata.glass.per.patient %>%
  dplyr::select(genomescan.sid.I, genomescan.sid.R) %>%
  tidyr::pivot_longer(cols=c('genomescan.sid.I', 'genomescan.sid.R')) %>%
  tidyr::drop_na(value) %>%
  dplyr::mutate(name=NULL) %>%
  dplyr::rename(genomescan.sid = value) %>%
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(metadata.glass.per.patient %>% dplyr::select(GLASS_ID, patient.correction.id), by=c('GLASS_ID'='GLASS_ID')) %>% 
  dplyr::arrange(Sample_Type, genomescan.sid)


tmp.data <- expression.glass %>%
  dplyr::select(all_of( tmp.metadata.2$genomescan.sid ))



stopifnot(colnames(tmp.data.2) == tmp.metadata.2$genomescan.sid)



dds.paired.a <- DESeqDataSetFromMatrix(countData = tmp.data,
                              colData = tmp.metadata,
                              design= ~ patient.correction.id + Sample_Type)



dds.paired.a <- DESeq(dds.paired.a)
res.paired.a <- results(dds.paired.a) %>% 
  as.data.frame() %>% 
  dplyr::filter(!is.na(padj)) %>% 
  dplyr::arrange(pvalue,padj) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::left_join(expression.glass.metadata %>% 
                     dplyr::select(gene_uid, gene_name, gene_type, gene_strand, gene_chr, gene_chr_center_loc, gene_loc),by=c('gene_uid'='gene_uid'))




res.paired.a %>%
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  dplyr::filter(padj < 0.01) %>% 
  dim


#saveRDS(res.paired.a, "cache/res.paired.a.Rds")
res.paired.a <- readRDS("cache/res.paired.a.Rds")



# 2b. paired [only paired] ----


tmp.metadata <- metadata.glass.per.patient %>% 
  dplyr::filter(pair.status == "complete") %>% 
  dplyr::select(genomescan.sid.I, genomescan.sid.R) %>% 
  tidyr::pivot_longer(cols=c('genomescan.sid.I', 'genomescan.sid.R')) %>% 
  dplyr::mutate(name=NULL) %>% 
  dplyr::rename(genomescan.sid = value) %>% 
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(metadata.glass.per.patient %>% dplyr::select(GLASS_ID, patient.correction.id), by=c('GLASS_ID'='GLASS_ID'))


tmp.data <- expression.glass %>% 
  dplyr::select(all_of( tmp.metadata$genomescan.sid ))


stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)



dds.paired.b <- DESeqDataSetFromMatrix(countData = tmp.data,
                              colData = tmp.metadata,
                              design= ~ patient.correction.id + Sample_Type)



dds.paired.b <- DESeq(dds.paired.b)
res.paired.b <- results(dds.paired.b) %>% 
  as.data.frame() %>% 
  dplyr::filter(!is.na(padj)) %>% 
  dplyr::arrange(pvalue,padj) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::left_join(expression.glass.metadata %>% 
                     dplyr::select(gene_uid, gene_name, gene_type, gene_strand, gene_chr, gene_chr_center_loc, gene_loc),by=c('gene_uid'='gene_uid'))




res.paired.b %>%
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  dplyr::filter(padj < 0.01) %>% 
  dim



# small visualisations for checking etc. ----


## PCA ----



tmp.plt.data <- expression.glass.vst %>%
  dplyr::select(
    c(metadata.glass.per.patient$genomescan.sid.I, metadata.glass.per.patient$genomescan.sid.R) %>% 
      purrr::keep(~ !is.na(.))
  ) %>% 
  dplyr::filter(rownames(.) %in% (res %>%  dplyr::slice_head(n=250) %>% dplyr::pull(gene_uid) )) %>% 
  t() %>% 
  prcomp() %>% 
  purrr::pluck('x') %>% 
  as.data.frame(stringsAsFactor=F) %>% 
  tibble::rownames_to_column('genomescan.sid') %>% 
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid'))



ggplot(tmp.plt.data,aes(x=PC1,y=PC2,col=Sample_Type, group=GLASS_ID, label=Sample_Name)) +
  geom_point() +
  youri_gg_theme +
  ggrepel::geom_text_repel(size=3, col="gray80")





## Another PCA of some sort ----


tmp.plt.data <- expression.glass.vst %>%
  dplyr::select(
    c(metadata.glass.per.patient$genomescan.sid.I, metadata.glass.per.patient$genomescan.sid.R) %>% 
      purrr::keep(~ !is.na(.))
  ) %>% 
  dplyr::filter(rownames(.) %in% (res.paired %>%  dplyr::slice_head(n=250) %>% dplyr::pull(gene_uid) )) %>% 
  t() %>% 
  prcomp() %>% 
  purrr::pluck('x') %>% 
  as.data.frame(stringsAsFactor=F) %>% 
  tibble::rownames_to_column('genomescan.sid') %>% 
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid'))



ggplot(tmp.plt.data,aes(x=PC1,y=PC2,col=Sample_Type, group=GLASS_ID, label=Sample_Name)) +
  geom_line(col="gray90") + 
  geom_point() +
  youri_gg_theme
#ggrepel::geom_text_repel(size=3, col="gray80")



## plot integration 1 ~ 2a ----


p1 <- res.unpaired.a %>%
  dplyr::mutate(stat.unpaired = stat) %>%
  dplyr::select(gene_uid, stat.unpaired,gene_name) %>% 
  dplyr::left_join(
    res.paired.a %>%
      dplyr::mutate(stat.paired = stat) %>%
      dplyr::select(gene_uid, stat.paired)
    , by=c('gene_uid'='gene_uid')
  ) %>% 
  dplyr::filter(!is.na(stat.unpaired) & !is.na(stat.paired)) %>% 
  dplyr::mutate(dist = abs(stat.unpaired - stat.paired))


ggplot(p1, aes(x= stat.unpaired, y=stat.paired,label=gene_name)) +
  geom_point() +
  youri_gg_theme + 
  ggpubr::stat_cor() +
  ggrepel::geom_text_repel(data = subset(p1, dist > 5))



p1 <- res.unpaired.a %>%
  dplyr::mutate(stat.unpaired = stat) %>%
  dplyr::select(gene_uid, stat.unpaired,gene_name) %>% 
  dplyr::left_join(
    res.paired.b %>%
      dplyr::mutate(stat.paired = stat) %>%
      dplyr::select(gene_uid, stat.paired)
    , by=c('gene_uid'='gene_uid')
  ) %>% 
  dplyr::filter(!is.na(stat.unpaired) & !is.na(stat.paired)) %>% 
  dplyr::mutate(dist = abs(stat.unpaired - stat.paired))


ggplot(p1, aes(x= stat.unpaired, y=stat.paired,label=gene_name)) +
  geom_point() +
  youri_gg_theme + 
  ggpubr::stat_cor() +
  ggrepel::geom_text_repel(data = subset(p1, dist > 5))




## chromosome plot

source('scripts/R/chrom_sizes.R')

plt <- res.paired %>% 
  dplyr::left_join(chrs_hg38_s, by=c('gene_chr'='chr')) %>% 
  dplyr::mutate(x = gene_chr_center_loc + pos) %>% 
  dplyr::mutate(gene_chr = factor(gene_chr, levels=gtools::mixedsort(unique(as.character(gene_chr))) ))


ggplot(plt , aes(x=x,y=stat,col=gene_chr)) + 
  geom_point(pch=19,cex=0.2) +
  geom_smooth() +
  youri_gg_theme

ggplot(plt, aes(x=gene_chr_center_loc / 1000000,y=stat,col=gene_chr)) + 
  facet_grid(cols = vars(gene_chr), scales = "free", space="free") +
  geom_point(pch=19,cex=0.2) +
  geom_smooth() +
  youri_gg_theme


plt$gene_chr = factor(plt$gene_chr, levels = unique(plt$gene_chr))

## enhancedvolcano's ----

EnhancedVolcano(res.unpaired.a,
                lab = res.unpaired.a$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.01)


EnhancedVolcano(res.paired.a,
                lab = res.paired.a$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.01)


EnhancedVolcano(res.paired.b,
                lab = res.paired.b$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.01)


# geiserplot x tpc ----


plt <- res.paired.a %>% 
  dplyr::left_join(
    expression.glass.metadata %>%
      dplyr::select(gene_uid, cor.t.dna.shallow.ACE.purity, cor.t.dna.wes.VAF_IDH, cor.t.methylation.purity.absolute),
    by=c('gene_uid'='gene_uid'))


ggplot(plt, aes(x=log2FoldChange, y=cor.t.dna.shallow.ACE.purity)) +
  geom_point(pch=19,cex=0.5) +
  xlim(-3,3)
ggplot(plt, aes(x=log2FoldChange, y=cor.t.dna.wes.VAF_IDH)) +
  geom_point(pch=19,cex=0.5) +
  xlim(-3,3)
ggplot(plt, aes(x=log2FoldChange, y=cor.t.methylation.purity.absolute)) +
  geom_point(pch=19,cex=0.5) +
  xlim(-3,3)




## comparison methylering ----

plt1 <- expression.glass.vst %>% 
  dplyr::filter(grepl("HOXD11",rownames(.))) %>% 
  t() %>% 
  as.data.frame %>% 
  tibble::rownames_to_column('genomescan.sid') %>% 
  dplyr::left_join(tmp.metadata, by=c('genomescan.sid'='genomescan.sid')) %>% 
  tidyr::drop_na(Sample_Type) 



ggplot(plt1, aes(x = Sample_Type, y= ENSG00000128713_HOXD11)) +
  ggbeeswarm::geom_quasirandom()


plt2 <- rowSums(t(expression.glass)) %>% 
  as.data.frame %>% 
  tibble::rownames_to_column('genomescan.sid') %>% 
  dplyr::left_join(tmp.metadata, by=c('genomescan.sid'='genomescan.sid')) %>% 
  tidyr::drop_na(Sample_Type) 


ggplot(plt2, aes(x = Sample_Type, y= `.`)) +
  ggbeeswarm::geom_quasirandom()



p3 <- plt1 %>% 
  dplyr::left_join(plt2, by=c('genomescan.sid'='genomescan.sid'))


ggplot(p3, aes(x = `.`, y=ENSG00000128713_HOXD11, col=Sample_Type.x)) +
  geom_point()



sign <- res.paired %>% 
  dplyr::filter(padj < 0.01) %>% 
  dplyr::filter(abs(log2FoldChange) >= 0.75)

cp <- expression.glass.vst %>%
  dplyr::select(tmp.metadata$genomescan.sid) %>% 
  dplyr::filter(rownames(.) %in% sign$gene_uid) %>% 
  `rownames<-`(gsub("ENSG00000284906_ARHGAP11B","ENSG00000284906_ARHGAP11B.2",rownames(.),fixed=T)) %>% 
  `rownames<-`(gsub("^ENS.+_","",rownames(.)))


cpm <- data.frame(gid=rownames(cp)) %>% 
  dplyr::mutate(HOX = grepl("^HOX", gid)) %>% 
  dplyr::mutate(COL = ifelse(grepl("^COL", gid),"red","green")) %>% 
  tibble::column_to_rownames('gid')


recursiveCorPlot::recursiveCorPlot(cp, cpm, 2 ,2)

#ggsave("/tmp/glass-supervised.png",height=20 * 1.3,width=30 * 1.3)



"CD248" %in% rownames(cp)




