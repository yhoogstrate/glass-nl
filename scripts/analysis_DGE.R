#!/usr/bin/env R



# load libs ----


source('scripts/R/youri_gg_theme.R')


library(DESeq2)
library(EnhancedVolcano)


# load data ----


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if(!exists("expression.glass.exon.vst")) {
  source('scripts/load_rna-counts.R')
}




# 1. unpaired exon ----



if(file.exists("cache/res.unpaired.a.exon.Rds")) {
  
  print("Loading 'cache/res.unpaired.a.exon.Rds' from cache")
  res.unpaired.a.exon <- readRDS("cache/res.unpaired.a.exon.Rds")
  
} else {

  tmp.metadata <- metadata.glass.per.patient %>%
    dplyr::select(genomescan.sid.I, genomescan.sid.R) %>%
    tidyr::pivot_longer(cols=c('genomescan.sid.I', 'genomescan.sid.R')) %>%
    tidyr::drop_na(value) %>%
    dplyr::mutate(name=NULL) %>%
    dplyr::rename(genomescan.sid = value) %>%
    dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid')) %>%
    dplyr::left_join(metadata.glass.per.patient %>% dplyr::select(GLASS_ID, patient.correction.id), by=c('GLASS_ID'='GLASS_ID')) %>% 
    dplyr::arrange(Sample_Type, genomescan.sid)
  
  
  tmp.data <- expression.glass.exon %>%
    dplyr::select(all_of( tmp.metadata$genomescan.sid ))
  
  
  
  stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)
  
  
  
  dds.unpaired.a.exon <- DESeqDataSetFromMatrix(countData = tmp.data,
                                colData = tmp.metadata,
                                design= ~ Sample_Type)
  
  
  
  dds.unpaired.a.exon <- DESeq(dds.unpaired.a.exon)
  res.unpaired.a.exon <- results(dds.unpaired.a.exon) %>% 
    as.data.frame(stringsAsFactors=F) %>% 
    tibble::rownames_to_column('gene_uid') %>% 
    dplyr::filter(!is.na(padj)) %>% 
    dplyr::arrange(pvalue,padj) %>% 
    dplyr::left_join(expression.glass.exon.metadata %>% dplyr::select(gene_uid, gene_name, gene_type, gene_strand, gene_loc),by=c('gene_uid'='gene_uid'))
  

  saveRDS(res.unpaired.a.exon, "cache/res.unpaired.a.exon.Rds")

}  


# n-sig
res.unpaired.a.exon %>%
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  dplyr::filter(padj < 0.01) %>% 
  dim


# append results
expression.glass.exon.metadata <- expression.glass.exon.metadata %>% 
  dplyr::left_join(
    res.unpaired.a.exon %>% 
      dplyr::select(gene_uid, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>% 
      tibble::column_to_rownames('gene_uid') %>% 
      `colnames<-`(paste0(colnames(.),".unpaired.exon")) %>% 
      tibble::rownames_to_column('gene_uid'),
    by=c('gene_uid'='gene_uid'),suffix = c("", "")
  ) 



rm(res.unpaired.a.exon) # cleanup, is joined to gene metadata anyway





# 2a. paired exon [incomplete as separate group] ----


if(file.exists("cache/res.paired.a.exon.Rds")) {
  
  # load from cache
  print("Loading 'cache/res.paired.a.exon.Rds' from cache")
  res.paired.a.exon <- readRDS("cache/res.paired.a.exon.Rds")
  
} else {
  
  tmp.metadata <- metadata.glass.per.patient %>%
    dplyr::select(genomescan.sid.I, genomescan.sid.R) %>%
    tidyr::pivot_longer(cols=c('genomescan.sid.I', 'genomescan.sid.R')) %>%
    tidyr::drop_na(value) %>%
    dplyr::mutate(name=NULL) %>%
    dplyr::rename(genomescan.sid = value) %>%
    dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid')) %>%
    dplyr::left_join(metadata.glass.per.patient %>% dplyr::select(GLASS_ID, patient.correction.id), by=c('GLASS_ID'='GLASS_ID')) %>% 
    dplyr::arrange(Sample_Type, genomescan.sid)
  
  
  tmp.data <- expression.glass.exon %>%
    dplyr::select(all_of( tmp.metadata$genomescan.sid ))
  
  
  stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)
  
  
  
  dds.paired.a.exon <- DESeqDataSetFromMatrix(countData = tmp.data,
                                colData = tmp.metadata,
                                design= ~ patient.correction.id + Sample_Type)
  
  
  
  dds.paired.a.exon <- DESeq(dds.paired.a.exon)
  res.paired.a.exon <- results(dds.paired.a.exon) %>% 
    as.data.frame() %>% 
    dplyr::filter(!is.na(padj)) %>% 
    dplyr::arrange(pvalue,padj) %>% 
    tibble::rownames_to_column('gene_uid') %>% 
    dplyr::left_join(expression.glass.exon.metadata %>% 
                       dplyr::select(gene_uid, gene_name, gene_type, gene_strand, gene_chr, gene_chr_center_loc, gene_loc),by=c('gene_uid'='gene_uid'))

  
  saveRDS(res.paired.a.exon, "cache/res.paired.a.exon.Rds")
}


# n sign
res.paired.a.exon %>%
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  dplyr::filter(padj < 0.01) %>% 
  dim



# append results
expression.glass.exon.metadata <- expression.glass.exon.metadata %>% 
  dplyr::left_join(
    res.paired.a.exon %>% 
      dplyr::select(gene_uid, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>% 
      tibble::column_to_rownames('gene_uid') %>% 
      `colnames<-`(paste0(colnames(.),".partially.paired.exon")) %>% 
      tibble::rownames_to_column('gene_uid'),
    by=c('gene_uid'='gene_uid'),suffix = c("", "")
  ) 




rm(res.paired.a.exon) # cleanup, is joined to gene metadata anyway




# 2b. paired exon [only paired] ----
# 
# 
# tmp.metadata <- metadata.glass.per.patient %>% 
#   dplyr::filter(pair.status == "complete") %>% 
#   dplyr::select(genomescan.sid.I, genomescan.sid.R) %>% 
#   tidyr::pivot_longer(cols=c('genomescan.sid.I', 'genomescan.sid.R')) %>% 
#   dplyr::mutate(name=NULL) %>% 
#   dplyr::rename(genomescan.sid = value) %>% 
#   dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid')) %>%
#   dplyr::left_join(metadata.glass.per.patient %>% dplyr::select(GLASS_ID, patient.correction.id), by=c('GLASS_ID'='GLASS_ID'))
# 
# 
# tmp.data <- expression.exon.glass %>% 
#   dplyr::select(all_of( tmp.metadata$genomescan.sid ))
# 
# 
# stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)
# 
# 
# 
# dds.paired.b <- DESeqDataSetFromMatrix(countData = tmp.data,
#                               colData = tmp.metadata,
#                               design= ~ patient.correction.id + Sample_Type)
# 
# 
# 
# dds.paired.b <- DESeq(dds.paired.b)
# res.paired.b <- results(dds.paired.b) %>% 
#   as.data.frame() %>% 
#   dplyr::filter(!is.na(padj)) %>% 
#   dplyr::arrange(pvalue,padj) %>% 
#   tibble::rownames_to_column('gene_uid') %>% 
#   dplyr::left_join(expression.glass.metadata %>% 
#                      dplyr::select(gene_uid, gene_name, gene_type, gene_strand, gene_chr, gene_chr_center_loc, gene_loc),by=c('gene_uid'='gene_uid'))
# 
# 
# 
# 
# res.paired.b %>%
#   dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
#   dplyr::filter(padj < 0.01) %>% 
#   dim
# 



# 3. unpaired gene ----


if(file.exists("cache/res.unpaired.a.gene.Rds")) {
  
  print("Loading 'cache/res.unpaired.a.gene.Rds' from cache")
  res.unpaired.a.gene <- readRDS("cache/res.unpaired.a.gene.Rds")
  
} else {

  tmp.metadata <- metadata.glass.per.patient %>%
    dplyr::select(genomescan.sid.I, genomescan.sid.R) %>%
    tidyr::pivot_longer(cols=c('genomescan.sid.I', 'genomescan.sid.R')) %>%
    tidyr::drop_na(value) %>%
    dplyr::mutate(name=NULL) %>%
    dplyr::rename(genomescan.sid = value) %>%
    dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid')) %>%
    dplyr::left_join(metadata.glass.per.patient %>% dplyr::select(GLASS_ID, patient.correction.id), by=c('GLASS_ID'='GLASS_ID')) %>% 
    dplyr::arrange(Sample_Type, genomescan.sid)
  
  
  tmp.data <- expression.glass.gene %>%
    dplyr::select(all_of( tmp.metadata$genomescan.sid ))
  
  
  
  stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)
  
  
  
  dds.unpaired.a.gene <- DESeqDataSetFromMatrix(countData = tmp.data,
                                                colData = tmp.metadata,
                                                design= ~ Sample_Type)
  
  
  
  dds.unpaired.a.gene <- DESeq(dds.unpaired.a.gene)
  res.unpaired.a.gene <- results(dds.unpaired.a.gene) %>% 
    as.data.frame(stringsAsFactors=F) %>% 
    tibble::rownames_to_column('gene_uid') %>% 
    dplyr::filter(!is.na(padj)) %>% 
    dplyr::arrange(pvalue,padj) %>% 
    dplyr::left_join(expression.glass.gene.metadata %>% dplyr::select(gene_uid, gene_name, gene_type, gene_strand, gene_loc),by=c('gene_uid'='gene_uid'))
  
  
  saveRDS(res.unpaired.a.gene, "cache/res.unpaired.a.gene.Rds")
  
  
}  


# n-sig
res.unpaired.a.gene %>%
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  dplyr::filter(padj < 0.01) %>% 
  dim


# append results
expression.glass.gene.metadata <- expression.glass.gene.metadata %>% 
  dplyr::left_join(
    res.unpaired.a.gene %>% 
      dplyr::select(gene_uid, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>% 
      tibble::column_to_rownames('gene_uid') %>% 
      `colnames<-`(paste0(colnames(.),".unpaired.gene")) %>% 
      tibble::rownames_to_column('gene_uid'),
    by=c('gene_uid'='gene_uid'),suffix = c("", "")
  ) 



rm(res.unpaired.a.gene) # cleanup, is joined to gene metadata anyway



# 4a. paired gene [incomplete as separate group] ----


if(file.exists("cache/res.paired.a.gene.Rds")) {
  
  # load from cache
  print("Loading 'cache/res.paired.a.gene.Rds' from cache")
  res.paired.a.gene <- readRDS("cache/res.paired.a.gene.Rds")
  

} else {
  
  tmp.metadata <- metadata.glass.per.patient %>%
    dplyr::select(genomescan.sid.I, genomescan.sid.R) %>%
    tidyr::pivot_longer(cols=c('genomescan.sid.I', 'genomescan.sid.R')) %>%
    tidyr::drop_na(value) %>%
    dplyr::mutate(name=NULL) %>%
    dplyr::rename(genomescan.sid = value) %>%
    dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid')) %>%
    dplyr::left_join(metadata.glass.per.patient %>% dplyr::select(GLASS_ID, patient.correction.id), by=c('GLASS_ID'='GLASS_ID')) %>% 
    dplyr::arrange(Sample_Type, genomescan.sid)
  
  
  tmp.data <- expression.glass.gene %>%
    dplyr::select(all_of( tmp.metadata$genomescan.sid ))
  
  
  stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)
  
  
  
  dds.paired.a.gene <- DESeqDataSetFromMatrix(countData = tmp.data,
                                              colData = tmp.metadata,
                                              design= ~ patient.correction.id + Sample_Type)
  
  
  
  dds.paired.a.gene <- DESeq(dds.paired.a.gene)
  res.paired.a.gene <- results(dds.paired.a.gene) %>% 
    as.data.frame() %>% 
    dplyr::filter(!is.na(padj)) %>% 
    dplyr::arrange(pvalue,padj) %>% 
    tibble::rownames_to_column('gene_uid') %>% 
    dplyr::left_join(expression.glass.gene.metadata %>% 
                       dplyr::select(gene_uid, gene_name, gene_type, gene_strand, gene_chr, gene_chr_center_loc, gene_loc),by=c('gene_uid'='gene_uid'))
  
  
  saveRDS(res.paired.a.gene, "cache/res.paired.a.gene.Rds")
  
  
  
}


# n sign
res.paired.a.gene %>%
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  dplyr::filter(padj < 0.01) %>% 
  dim



# append results
expression.glass.gene.metadata <- expression.glass.gene.metadata %>% 
  dplyr::left_join(
    res.paired.a.gene %>% 
      dplyr::select(gene_uid, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>% 
      tibble::column_to_rownames('gene_uid') %>% 
      `colnames<-`(paste0(colnames(.),".partially.paired.gene")) %>% 
      tibble::rownames_to_column('gene_uid'),
    by=c('gene_uid'='gene_uid'),suffix = c("", "")
  ) 




rm(res.paired.a.gene) # cleanup, is joined to gene metadata anyway


# small plot on concordange gene/exon ----


plt <- dplyr::inner_join(
  expression.glass.exon.metadata,
  expression.glass.gene.metadata,
  by=c('gene_uid'='gene_uid') ) %>%
  dplyr::mutate(signi.paired = padj.partially.paired.exon < 0.01 | padj.partially.paired.gene < 0.01) %>%
  dplyr::mutate(signi.unpaired = padj.unpaired.exon < 0.01 | padj.unpaired.gene < 0.01)



p1 <- ggplot(plt, aes(x=stat.partially.paired.exon,
                y=stat.partially.paired.gene,col=signi.paired)) +
  geom_point(pch=19,cex=0.5) +
  xlim(-8.5,8.5) +
  ylim(-8.5,8.5) +
  geom_abline(intercept = 0,col="gray60",lty=1,lwd=0.5) +
  geom_smooth(method="lm",col="black",lwd=0.5) +
  theme_bw()




p2 <- ggplot(plt, aes(x=stat.unpaired.exon,
                      y=stat.unpaired.gene,col=signi.unpaired)) +
  geom_point(pch=19,cex=0.5) +
  xlim(-8.5,8.5) +
  ylim(-8.5,8.5) +
  geom_abline(intercept = 0,col="gray60",lty=1,lwd=0.5) +
  geom_smooth(method="lm",col="black",lwd=0.5) +
  theme_bw()

p1+p2



p1 <- ggplot(plt, aes(x=-log10(padj.unpaired.exon),
                y=-log10(padj.unpaired.gene),
                col=signi.unpaired)) +
  geom_point(pch=19,cex=0.5) +
  xlim(-0.1,22) +
  ylim(-0.1,22) +
  geom_abline(intercept = 0,col="gray30",lty=1,lwd=0.5) +
  #geom_smooth(method="lm",col="black",lwd=0.5) +
  theme_bw()


p2 <- ggplot(plt, aes(x=-log10(padj.partially.paired.exon),
                y=-log10(padj.partially.paired.gene),
                col=signi.paired)) +
  geom_point(pch=19,cex=0.5) +
  xlim(-0.1,11) +
  ylim(-0.1,11) +
  geom_abline(intercept = 0,col="gray30",lty=1,lwd=0.5) +
  #geom_smooth(method="lm",col="black",lwd=0.5) +
  theme_bw()


p1 + p2




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
  dplyr::mutate(dist = abs(stat.unpaired - stat.paired)) %>% 
  dplyr::mutate(proteomics.rf.imp = gene_name %in%
                  c("GCLC","SORD","SNX18","AGRN","METTL7A","BAG6","FCGBP",
                    "NADK2","ACAA2","SEMA4B") |
                  grepl("ENSG00000001084",gene_uid) |
                  grepl("ENSG00000140263",gene_uid) |
                  grepl("ENSG00000152620",gene_uid) |
                  grepl("ENSG00000167315",gene_uid)
                )


ggplot(p1, aes(x= stat.unpaired, y=stat.paired,label=gene_name,col=proteomics.rf.imp)) +
  geom_point(data = p1 %>% dplyr::filter(proteomics.rf.imp == F)) +
  geom_point(data = p1 %>% dplyr::filter(proteomics.rf.imp == T)) +
  youri_gg_theme + 
  #ggpubr::stat_cor() +
  ggrepel::geom_text_repel(data = subset(p1, dist > 5 | proteomics.rf.imp)) +
  scale_color_manual(values=c('TRUE'='black','FALSE'='red'))







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






### chromomsome 9 ----

plt <- res.paired.a %>%
  dplyr::filter(gene_chr == "chr9") %>% 
  dplyr::left_join(chrs_hg38_s, by=c('gene_chr'='chr')) %>% 
  dplyr::mutate(x = gene_chr_center_loc + pos) %>% 
  dplyr::mutate(gene_chr = factor(gene_chr, levels=gtools::mixedsort(unique(as.character(gene_chr))) )) %>% 
  dplyr::mutate(m.l.padj = -log10(padj)) %>% 
  tidyr::pivot_longer(cols=c(m.l.padj, log2FoldChange, stat)) %>% 
  dplyr::mutate(name = factor(name, levels=c('m.l.padj', 'log2FoldChange', 'stat') ))


ggplot(plt, aes(x=gene_chr_center_loc / 1000000,y=value,col=gene_chr)) + 
  facet_grid(cols = vars(gene_chr), rows=vars(name), scales = "free", space="free_x") +
  geom_smooth(alpha=0.5, se = FALSE, col="gray60") +
  geom_point(pch=19,cex=0.2) +
  youri_gg_theme +
  labs(x=NULL,y="DESeq2 wald statistic") +
  geom_vline(xintercept = 36, alpha=0.4)
#geom_vline(xintercept = 33.67, alpha=0.4)
#geom_vline(xintercept = 36, alpha=0.4)



### chromomsome 11 ----


plt <- res.paired.a %>%
  dplyr::filter(gene_chr == "chr11") %>% 
  dplyr::left_join(chrs_hg38_s, by=c('gene_chr'='chr')) %>% 
  dplyr::mutate(x = gene_chr_center_loc + pos) %>% 
  dplyr::mutate(gene_chr = factor(gene_chr, levels=gtools::mixedsort(unique(as.character(gene_chr))) )) %>% 
  dplyr::mutate(m.l.padj = -log10(padj)) %>% 
  tidyr::pivot_longer(cols=c(m.l.padj, log2FoldChange, stat)) %>% 
  dplyr::mutate(name = factor(name, levels=c('m.l.padj', 'log2FoldChange', 'stat') ))


ggplot(plt, aes(x=gene_chr_center_loc / 1000000,y=value,col=gene_chr)) + 
  facet_grid(cols = vars(gene_chr), rows=vars(name), scales = "free", space="free_x") +
  #geom_point(data = plt %>% dplyr::filter(gene_name == "CCND1"), pch=19,cex=1.2,col="black") +
  geom_point(data = plt %>% dplyr::filter(gene_name %in% c("CD248","CDCA5")), pch=19,cex=1.2,col="black") +
  geom_point(pch=19,cex=0.2) +
  geom_smooth(alpha=0.5, se = FALSE, col="gray60") +
  youri_gg_theme + 
  labs(x=NULL)




res.paired.a %>%
  dplyr::filter(gene_chr == "chr11") %>%
  dplyr::filter(stat > 4) %>%
  dplyr::filter(gene_chr_center_loc > 50000000 & gene_chr_center_loc < 75000000) %>% 
  dplyr::arrange(gene_chr_center_loc)


### chromomsome 19 ----


plt <- res.paired.a %>%
  dplyr::filter(gene_chr == "chr19") %>% 
  dplyr::left_join(chrs_hg38_s, by=c('gene_chr'='chr')) %>% 
  dplyr::mutate(x = gene_chr_center_loc + pos) %>% 
  dplyr::mutate(gene_chr = factor(gene_chr, levels=gtools::mixedsort(unique(as.character(gene_chr))) )) %>% 
  dplyr::mutate(m.l.padj = -log10(padj)) %>% 
  tidyr::pivot_longer(cols=c(m.l.padj, log2FoldChange, stat)) %>% 
  dplyr::mutate(name = factor(name, levels=c('m.l.padj', 'log2FoldChange', 'stat') ))


ggplot(plt, aes(x=gene_chr_center_loc / 1000000,y=value,col=gene_chr)) + 
  facet_grid(cols = vars(gene_chr), rows=vars(name), scales = "free", space="free_x") +
  geom_vline(xintercept = 19.937823 , alpha=0.4) +
  geom_vline(xintercept = 24.219001, alpha=0.4) +
  geom_point(pch=19,cex=0.2) +
  geom_smooth(alpha=0.5, se = FALSE, col="gray60") +
  youri_gg_theme + 
  labs(x=NULL)




res.paired.a %>%
  dplyr::filter(gene_chr == "chr19") %>%
  dplyr::filter(stat > 4) %>%
  dplyr::filter(gene_chr_center_loc > 50000000 & gene_chr_center_loc < 75000000) %>% 
  dplyr::arrange(gene_chr_center_loc)



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



if("cor.t.methylation.purity.absolute" %in% names(expression.glass.metadata) == F) {
  source('scripts/load_correlation_expression_purity.R') 
}




plt <- res.paired.a %>% 
  dplyr::left_join(
    expression.glass.metadata %>%
      dplyr::select(gene_uid, cor.t.dna.shallow.ACE.purity, cor.t.dna.wes.VAF_IDH, cor.t.methylation.purity.absolute, cor.t.dna.purity.manual.Erik),
    by=c('gene_uid'='gene_uid')) %>% 
  dplyr::select(gene_name, log2FoldChange, cor.t.dna.purity.manual.Erik,  cor.t.dna.wes.VAF_IDH, cor.t.methylation.purity.absolute) %>% 
  tidyr::pivot_longer(cols = -c(gene_name, log2FoldChange)) %>% 
  dplyr::rename(method = name) %>% 
  dplyr::rename(purity = value) %>% 
  dplyr::arrange(abs(log2FoldChange))



ggplot(plt, aes(x=log2FoldChange, y=purity, fill = method)) +
  facet_grid(cols = vars(method),  scales = "free", space="free") +
  geom_point(pch=21, fill=alpha("white",0.6), col="gray10", cex=2) +
  xlim(-2.5,2.5)+
  youri_gg_theme








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





# DE invasive ~ expansive ----


tmp.metadata <- metadata.glass.per.resection %>% 
  dplyr::filter(excluded == F & !is.na(imaging.growth_pattern) & imaging.growth_pattern %in% c('Mostly expansive','Mostly invasive')) %>% 
  dplyr::mutate(imaging.growth_pattern = factor(gsub(" ",".",imaging.growth_pattern),levels=c('Mostly.expansive','Mostly.invasive')))


tmp.data <- expression.glass.exon %>%
  dplyr::select(all_of( tmp.metadata$genomescan.sid ))


stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)


dds <- DESeq2::DESeqDataSetFromMatrix(countData = tmp.data,
                                         colData = tmp.metadata,
                                         design= ~ imaging.growth_pattern) %>% 
  DESeq2::DESeq()

res <- DESeq2::results(dds) %>% 
  as.data.frame(stringsAsFactors=F) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::filter(!is.na(padj)) %>% 
  dplyr::arrange(pvalue,padj) %>% 
  dplyr::left_join(expression.glass.exon.metadata %>% dplyr::select(gene_uid, gene_name, gene_type, gene_strand, gene_loc),by=c('gene_uid'='gene_uid'))


res %>%
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  dplyr::filter(padj < 0.01) %>% 
  dim


EnhancedVolcano(res,
                lab = res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.01)



plt <- metadata.glass.per.resection %>% 
#plt <- tmp.metadata %>% 
  dplyr::filter(excluded == F) %>% 
  dplyr::left_join(
    expression.glass.vst %>%
      t() %>%
      as.data.frame %>%
      tibble::rownames_to_column('genomescan.sid') %>% 
      dplyr::select( genomescan.sid , ENSG00000167244_IGF2, ENSG00000128710_HOXD10), by=c('genomescan.sid'='genomescan.sid')
  )



ggplot(plt, aes(x=time.resection.until.last.event, y=ENSG00000167244_IGF2, col=imaging.growth_pattern)) +
  geom_point() +
  theme_bw()

ggplot(plt, aes(x=time.resection.until.last.event, y=ENSG00000128710_HOXD10, col=imaging.growth_pattern)) +
  geom_point() +
  theme_bw()



## tmp ----


plt.tmp <- res.paired.a %>% 
  dplyr::left_join(
    data.frame(pbc.obj$importance) %>% tibble::rownames_to_column('gene_uid'),
    by=c('gene_uid' = 'gene_uid')
  ) %>% 
  dplyr::mutate(m.l.padj = -log10(padj))  %>%
  dplyr::mutate(pbc.obj.importance = log((pbc.obj.importance * 1000) + 1 ))



ggplot(plt.tmp, aes(x = m.l.padj, y= pbc.obj.importance) ) +
  geom_point(pch=19, cex=0.5) +
  xlim(0,5.5)  +
  labs(x = "adjusted p-value time difference", y="importance predicting survival") +
  theme_bw()


## tmp2 ----


plt.tmp <- res %>% 
  dplyr::left_join(
    data.frame(pbc.obj$importance) %>% tibble::rownames_to_column('gene_uid'),
    by=c('gene_uid' = 'gene_uid')
  ) %>% 
  dplyr::mutate(m.l.padj = -log10(padj))  %>%
  dplyr::mutate(pbc.obj.importance = log((pbc.obj.importance * 1000) + 1 ))



ggplot(plt.tmp, aes(x = m.l.padj, y= pbc.obj.importance) ) +
  geom_point(pch=19, cex=0.5) +
  xlim(0,2)  +
  labs(x = "adjusted p-value time difference", y="importance predicting survival") +
  theme_bw()




