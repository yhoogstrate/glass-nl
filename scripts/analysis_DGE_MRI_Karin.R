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




