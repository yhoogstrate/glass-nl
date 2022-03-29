#!/usr/bin/env R


# load libs ----


source('scripts/R/youri_gg_theme.R')
source('scripts/R/chrom_sizes.R')


library(DESeq2)
library(EnhancedVolcano)


# load data ----


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if(!exists("expression.glass.vst")) {
  source('scripts/load_rna-counts.R')
}

if("cor.t.methylation.purity.absolute" %in% names(expression.glass.metadata) == F) {
  source('scripts/load_correlation_expression_purity.R') 
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
  dplyr::select(all_of( tmp.metadata$genomescan.sid ))



stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)



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
  dplyr::select(all_of( tmp.metadata$genomescan.sid ))



stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)



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




## chromosome plot ----


plt <- res.paired.a %>% 
  dplyr::left_join(chrs_hg38_s, by=c('gene_chr'='chr')) %>% 
  dplyr::mutate(x = gene_chr_center_loc + pos) %>% 
  dplyr::mutate(gene_chr = factor(gene_chr, levels=gtools::mixedsort(unique(as.character(gene_chr))) ))


# ggplot(plt , aes(x=x,y=stat,col=gene_chr)) + 
#   geom_point(pch=19,cex=0.2) +
#   geom_smooth() +
#   youri_gg_theme

ggplot(plt, aes(x=gene_chr_center_loc / 1000000,y=stat,col=gene_chr)) + 
  facet_grid(cols = vars(gene_chr), scales = "free", space="free") +
  geom_point(pch=19,cex=0.2) +
  geom_smooth(se=F,col="black", lwd=0.7) +
  youri_gg_theme + labs(x=NULL)


#plt$gene_chr = factor(plt$gene_chr, levels = unique(plt$gene_chr))


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


# recursiveCorPlot ----


sign <- res.paired.a %>% 
  dplyr::filter(padj < 0.01) %>% 
  dplyr::filter(abs(log2FoldChange) >= 0.75)

cp <- expression.glass.vst %>%
  dplyr::select(tmp.metadata$genomescan.sid) %>% 
  dplyr::filter(rownames(.) %in% sign$gene_uid) %>% 
  `rownames<-`(gsub("ENSG00000284906_ARHGAP11B","ENSG00000284906_ARHGAP11B.2",rownames(.),fixed=T)) %>% 
  `rownames<-`(gsub("^ENS.+_","",rownames(.)))


dim(cp)
# find n PCA - 5 features?!
pca <- prcomp(t(cp))
plot(pca)
dev.off()


pcs <- pca$rotation[,1:5] %>%
  as.data.frame() %>% 
  dplyr::mutate(main.PC = {names(.)[max.col(.)]}) %>% 
  dplyr::select(main.PC) %>% 
  tibble::rownames_to_column('gene_name') %>% 
  dplyr::mutate(PC1 = ifelse(main.PC == "PC1", NA , F )) %>% 
  dplyr::mutate(PC2 = ifelse(main.PC == "PC2", NA , F )) %>% 
  dplyr::mutate(PC3 = ifelse(main.PC == "PC3", NA , F )) %>% 
  dplyr::mutate(PC4 = ifelse(main.PC == "PC4", NA , F )) %>% 
  dplyr::mutate(PC5 = ifelse(main.PC == "PC5", NA , F )) %>% 
  dplyr::mutate(main.PC = NULL)




cpm <- data.frame(gid=rownames(cp)) %>% 
  dplyr::mutate(HOX = grepl("^HOX", gid)) %>% 
  dplyr::mutate(COL = ifelse(grepl("^COL", gid),"red","green")) %>%
  dplyr::left_join(pcs, by=c('gid'='gene_name')) %>% 
  tibble::column_to_rownames('gid')


h <- recursiveCorPlot::recursiveCorPlot(cp, cpm, 2 ,2)

ggsave("/tmp/glass-supervised.png",height=20 * 1.3,width=30 * 1.3)



"CD248" %in% rownames(cp)



# DE invasive ~ expansive ----


tmp.metadata <- metadata.glass.per.resection %>% 
  dplyr::filter(excluded == F & !is.na(imaging.growth_pattern) & imaging.growth_pattern %in% c('Mostly expansive','Mostly invasive')) %>% 
  dplyr::mutate(imaging.growth_pattern = factor(gsub(" ",".",imaging.growth_pattern),levels=c('Mostly.expansive','Mostly.invasive')))


tmp.data <- expression.glass %>%
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
  dplyr::left_join(expression.glass.metadata %>% dplyr::select(gene_uid, gene_name, gene_type, gene_strand, gene_loc),by=c('gene_uid'='gene_uid'))


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




