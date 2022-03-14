#!/usr/bin/env R 

# load libs ----


# load data ----

if(!exists('expression.glass.vst')) {
  source('scripts/load_rna-counts.R')
}

if(!exists("metadata.glass.per.resection") | "dna.wes.VAF_IDH" %in% names(metadata.glass.per.resection) == F) {
  source('scripts/load_tumour_purities.R')
}


# calc corr for ACE method ----


tmp.metadata <- metadata.glass.per.resection %>%
  dplyr::filter(excluded == F) %>% 
  dplyr::filter(!is.na(dna.shallow.ACE.purity)) %>% 
  dplyr::select(genomescan.sid, dna.shallow.ACE.purity)

tmp.data <- expression.glass.vst %>% 
  dplyr::select(tmp.metadata$genomescan.sid)


stopifnot(tmp.metadata$genomescan.sid == colnames(tmp.data))


tmp.out <- data.frame(cor.t.dna.shallow.ACE.purity = apply(tmp.data,1, function(vec) {return( cor.test(tmp.metadata %>% dplyr::pull(dna.shallow.ACE.purity), as.numeric(vec))$statistic) })) %>% 
  tibble::rownames_to_column('gene_uid')


tmp.out %>% dplyr::arrange(cor.t.dna.shallow.ACE.purity) %>% head(n=10)
tmp.out %>% dplyr::arrange(-cor.t.dna.shallow.ACE.purity) %>% head(n=10)


expression.glass.metadata <- expression.glass.metadata %>% 
  dplyr::left_join(tmp.out , by=c('gene_uid'='gene_uid'), keep=F,suffix = c("", "")) # force overwrite


rm(tmp.metadata, tmp.data, tmp.out)



# calc corr for ACE method (< 1.0) ----

tmp.metadata <- metadata.glass.per.resection %>%
  dplyr::filter(excluded == F) %>% 
  dplyr::filter(!is.na(dna.shallow.ACE.purity.below.1)) %>% 
  dplyr::select(genomescan.sid, dna.shallow.ACE.purity.below.1)

tmp.data <- expression.glass.vst %>% 
  dplyr::select(tmp.metadata$genomescan.sid)


stopifnot(tmp.metadata$genomescan.sid == colnames(tmp.data))


tmp.out <- data.frame(cor.t.dna.shallow.ACE.purity.below.1 = apply(tmp.data,1, function(vec) {return( cor.test(tmp.metadata %>% dplyr::pull(dna.shallow.ACE.purity.below.1), as.numeric(vec))$statistic) })) %>% 
  tibble::rownames_to_column('gene_uid')


tmp.out %>% dplyr::arrange(cor.t.dna.shallow.ACE.purity.below.1) %>% head(n=10)
tmp.out %>% dplyr::arrange(-cor.t.dna.shallow.ACE.purity.below.1) %>% head(n=10)


expression.glass.metadata <- expression.glass.metadata %>% 
  dplyr::left_join(tmp.out , by=c('gene_uid'='gene_uid'), keep=F,suffix = c("", "")) # force overwrite


rm(tmp.metadata, tmp.data, tmp.out)


# calc corr for Erik/Manual method ----

tmp.metadata <- metadata.glass.per.resection %>%
  dplyr::filter(excluded == F) %>% 
  dplyr::filter(!is.na(dna.purity.manual.Erik)) %>% 
  dplyr::select(genomescan.sid, dna.purity.manual.Erik)


tmp.data <- expression.glass.vst %>% 
  dplyr::select(tmp.metadata$genomescan.sid)


stopifnot(tmp.metadata$genomescan.sid == colnames(tmp.data))


tmp.out <- data.frame(cor.t.dna.purity.manual.Erik = apply(tmp.data,1, function(vec) {return( cor.test(tmp.metadata %>% dplyr::pull(dna.purity.manual.Erik), as.numeric(vec))$statistic) })) %>% 
  tibble::rownames_to_column('gene_uid')


tmp.out %>% dplyr::arrange(cor.t.dna.purity.manual.Erik) %>% head(n=10)
tmp.out %>% dplyr::arrange(-cor.t.dna.purity.manual.Erik) %>% head(n=10)


expression.glass.metadata <- expression.glass.metadata %>% 
  dplyr::left_join(tmp.out , by=c('gene_uid'='gene_uid'), keep=F,suffix = c("", "")) # force overwrite


rm(tmp.metadata, tmp.data, tmp.out)



# calc corr for VAF method ----

tmp.metadata <- metadata.glass.per.resection %>%
  dplyr::filter(excluded == F) %>% 
  dplyr::filter(!is.na(dna.wes.VAF_IDH))


tmp.data <- expression.glass.vst %>% 
  dplyr::select(tmp.metadata$genomescan.sid)


stopifnot(tmp.metadata$genomescan.sid == colnames(tmp.data))


tmp.out <- data.frame(cor.t.dna.wes.VAF_IDH = apply(tmp.data,1, function(vec) {return( cor.test(tmp.metadata$dna.wes.VAF_IDH, as.numeric(vec))$statistic) })) %>% 
  tibble::rownames_to_column('gene_uid')



tmp.out %>% dplyr::arrange(cor.t.dna.wes.VAF_IDH) %>% head(n=10)
tmp.out %>% dplyr::arrange(-cor.t.dna.wes.VAF_IDH) %>% head(n=20)


expression.glass.metadata <- expression.glass.metadata %>% 
  dplyr::left_join(tmp.out , by=c('gene_uid'='gene_uid'), keep=F,suffix = c("", "")) # force overwrite


rm(tmp.metadata, tmp.data, tmp.out)




# calc corr for Meth/RFpurity method ----


tmp.metadata <- metadata.glass.per.resection %>%
  dplyr::filter(excluded == F) %>% 
  dplyr::filter(!is.na(methylation.purity.absolute))


tmp.data <- expression.glass.vst %>% 
  dplyr::select(tmp.metadata$genomescan.sid)


stopifnot(tmp.metadata$genomescan.sid == colnames(tmp.data))


tmp.out <- data.frame(cor.t.methylation.purity.absolute = apply(tmp.data,1, function(vec) {return( cor.test(tmp.metadata$methylation.purity.absolute, as.numeric(vec))$statistic) })) %>% 
  tibble::rownames_to_column('gene_uid')



tmp.out %>% dplyr::arrange(cor.t.methylation.purity.absolute) %>% head(n=10)
tmp.out %>% dplyr::arrange(-cor.t.methylation.purity.absolute) %>% head(n=20)


expression.glass.metadata <- expression.glass.metadata %>% 
  dplyr::left_join(tmp.out , by=c('gene_uid'='gene_uid'), keep=F,suffix = c("", "")) # force overwrite


rm(tmp.metadata, tmp.data, tmp.out)




