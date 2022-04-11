#!/usr/bin/env R

# load libs ----


source('scripts/R/youri_gg_theme.R')


# load data ----


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if(!exists("expression.glass.vst")) {
  source('scripts/load_rna-counts.R')
}


# simple PCA ----


plt <- metadata.glass.per.resection %>% 
  dplyr::filter(excluded == F)


stopifnot(plt$genomescan.sid == colnames(expression.glass.vst))


plt.pca <- expression.glass.vst %>%
  dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>%  # use MAD instead of SD to find most variable genes
  dplyr::arrange(-mad) %>% # arrange in asc? order
  dplyr::slice_head(n=1000) %>%  # pick top 1000
  #dplyr::filter(rownames(.) %in% a) %>% 
  dplyr::mutate(mad=NULL) %>% # remove the sd to obtain original vst matrix
  t() %>% # transpose, to PCA the genes rather than the patients
  prcomp %>% # PCA
  purrr::pluck('x') %>%  # take coordinates
  as.data.frame(stringsAsFactor=F) %>% # transform back from matrix to data.frame 
  dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
  tibble::rownames_to_column('genomescan.sid') %>% 
  dplyr::left_join(metadata.glass.per.resection , by=c('genomescan.sid'='genomescan.sid')) %>% 
  dplyr::mutate(low.q.suspects = genomescan.sid %in% c('104059-002-054','104059-003-035','104059-002-103','104059-002-176'))




# ggplot(plt.pca, aes(x = PC1, y=PC2, col=low.q.suspects, group=GLASS_ID, label=Sample_Name)) +
#   geom_point() +
#   youri_gg_theme +
#   ggrepel::geom_text_repel(size=3, col="gray80")
# 
# 
# plt.pca %>% dplyr::filter(PC2 > 25) %>%  dplyr::pull(genomescan.sid)
# 
# 

# 
# ggplot(plt.pca, aes(x = PC1, y=PC2, col=low.q.suspects, group=GLASS_ID, label=Sample_Name)) +
#   geom_point() +
#   youri_gg_theme +
#   ggrepel::geom_text_repel(size=3, col="gray80")
# 


ggplot(plt.pca, aes(x = PC1, y=PC2, col=resection, group=GLASS_ID, label=Sample_Name)) +
  geom_point() +
  youri_gg_theme +
  ggrepel::geom_text_repel(size=3, col="gray80")


ggsave("output/figures/vis_unsupervised_expression_01.png",width=8,height=6)




ggplot(plt.pca, aes(x = PC3, y=PC4, col=Sample_Sex, group=GLASS_ID, label=Sample_Name)) +
  geom_line() +
  geom_point() +
  youri_gg_theme +
  ggrepel::geom_text_repel(size=3, col="gray80")

ggsave("output/figures/vis_unsupervised_expression_02.png",width=8,height=6)


# gender PCA exon + gene-body ----


stopifnot(plt$genomescan.sid == colnames(expression.glass.exon.vst))
stopifnot(plt$genomescan.sid == colnames(expression.glass.gene.vst))


## exon part ----


plt.pca.exon <- expression.glass.exon.vst %>%
  dplyr::select(metadata.glass.per.resection %>% 
                  dplyr::filter(excluded == F) %>%
                  dplyr::pull(genomescan.sid)) %>% 
  dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>%  # use MAD instead of SD to find most variable genes
  dplyr::arrange(-mad) %>% # arrange in asc? order
  dplyr::slice_head(n=1000) %>%  # pick top 1000
  #dplyr::filter(rownames(.) %in% a) %>% 
  dplyr::mutate(mad=NULL) %>% # remove the sd to obtain original vst matrix
  t() %>% # transpose, to PCA the genes rather than the patients
  prcomp %>% # PCA
  purrr::pluck('x') %>%  # take coordinates
  as.data.frame(stringsAsFactor=F) %>% # transform back from matrix to data.frame 
  dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
  `colnames<-`(paste0(colnames(.),".exon")) %>% 
  tibble::rownames_to_column('genomescan.sid') 



## gene-body part ----


plt.pca.gene <- expression.glass.gene.vst %>%
  dplyr::select(metadata.glass.per.resection %>% 
                  dplyr::filter(excluded == F) %>%
                  dplyr::pull(genomescan.sid)) %>% 
  dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>%  # use MAD instead of SD to find most variable genes
  dplyr::arrange(-mad) %>% # arrange in asc? order
  dplyr::slice_head(n=1000) %>%  # pick top 1000
  #dplyr::filter(rownames(.) %in% a) %>% 
  dplyr::mutate(mad=NULL) %>% # remove the sd to obtain original vst matrix
  t() %>% # transpose, to PCA the genes rather than the patients
  prcomp %>% # PCA
  purrr::pluck('x') %>%  # take coordinates
  as.data.frame(stringsAsFactor=F) %>% # transform back from matrix to data.frame 
  dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
  `colnames<-`(paste0(colnames(.),".gene")) %>% 
  tibble::rownames_to_column('genomescan.sid')


## integrate ----


stopifnot(plt.pca.exon$genomescan.sid == plt.pca.gene$genomescan.sid)


plt.pca <- dplyr::left_join(plt.pca.exon, plt.pca.gene, by=c('genomescan.sid'='genomescan.sid')) %>% 
  dplyr::left_join(metadata.glass.per.resection , by=c('genomescan.sid'='genomescan.sid'))


p1 <- ggplot(plt.pca, aes(x=PC4.exon, y=PC3.exon, col=Sample_Sex, group=GLASS_ID, label=Sample_Name)) +
  geom_line() +
  geom_point() +
  youri_gg_theme +
  ggrepel::geom_text_repel(size=3, col="gray80")

p2 <- ggplot(plt.pca, aes(x=PC4.gene, y=PC3.gene, col=Sample_Sex, group=GLASS_ID, label=Sample_Name)) +
  geom_line() +
  geom_point() +
  youri_gg_theme +
  ggrepel::geom_text_repel(size=3, col="gray80")


p1 + p2 




##  concorande exon & gene-body ----

plt <- dplyr::left_join(
  data.frame(exon = rowSums(expression.glass.exon)) %>% 
    tibble::rownames_to_column('gene_uid'),
  data.frame(gene = rowSums(expression.glass.gene)) %>% 
    tibble::rownames_to_column('gene_uid'),
  by=c('gene_uid'='gene_uid')
)


ggplot(plt, aes(x=exon, y=gene)) +
  geom_point()



