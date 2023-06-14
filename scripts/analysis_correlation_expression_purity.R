#!/usr/bin/env R 

# load libs ----


# load data ----

if(!exists('expression.glass.exon.vst')) {
  source('scripts/load_rna-counts.R')
}

if(!exists("metadata.glass.per.resection") | "dna.wes.VAF_IDH" %in% names(metadata.glass.per.resection) == F) {
  source('scripts/load_tumour_purities.R')
}


# calc correlations ----


# generic output object
out <- data.frame(gene_uid = rownames(expression.glass.exon.vst))


## for ACE method ----


# tmp.metadata <- metadata.glass.per.resection %>%
#   dplyr::filter(excluded == F) %>% 
#   dplyr::filter(!is.na(dna.shallow.ACE.purity)) %>% 
#   dplyr::select(genomescan.sid, dna.shallow.ACE.purity)
# 
# tmp.data <- expression.glass.exon.vst %>% 
#   dplyr::select(tmp.metadata$genomescan.sid)
# 
# 
# stopifnot(tmp.metadata$genomescan.sid == colnames(tmp.data))
# 
# 
# tmp.out <- data.frame(cor.t.dna.shallow.ACE.purity = apply(tmp.data,1, function(vec) {return( cor.test(tmp.metadata %>% dplyr::pull(dna.shallow.ACE.purity), as.numeric(vec))$statistic) })) %>% 
#   tibble::rownames_to_column('gene_uid')
# 
# 
# # tmp.out %>% dplyr::arrange(cor.t.dna.shallow.ACE.purity) %>% head(n=10)
# # tmp.out %>% dplyr::arrange(-cor.t.dna.shallow.ACE.purity) %>% head(n=10)
# 
# 
# expression.glass.exon.metadata <- expression.glass.exon.metadata %>% 
#   dplyr::left_join(tmp.out , by=c('gene_uid'='gene_uid'), keep=F,suffix = c("", "")) # force overwrite
# 
# 
# rm(tmp.metadata, tmp.data, tmp.out)



## for ACE method (< 1.0) ----

# tmp.metadata <- metadata.glass.per.resection %>%
#   dplyr::filter(excluded == F) %>% 
#   dplyr::filter(!is.na(dna.shallow.ACE.purity.below.1)) %>% 
#   dplyr::select(genomescan.sid, dna.shallow.ACE.purity.below.1)
# 
# tmp.data <- expression.glass.exon.vst %>% 
#   dplyr::select(tmp.metadata$genomescan.sid)
# 
# 
# stopifnot(tmp.metadata$genomescan.sid == colnames(tmp.data))
# 
# 
# tmp.out <- data.frame(cor.t.dna.shallow.ACE.purity.below.1 = apply(tmp.data,1, function(vec) {return( cor.test(tmp.metadata %>% dplyr::pull(dna.shallow.ACE.purity.below.1), as.numeric(vec))$statistic) })) %>% 
#   tibble::rownames_to_column('gene_uid')
# 
# 
# 
# nrow(expression.glass.exon.metadata) == nrow(expression.glass.exon)
# stopifnot(rownames(expression.glass.exon) == expression.glass.exon.metadata$gene_uid)
# 
# expression.glass.exon.metadata <- expression.glass.exon.metadata %>% 
#   dplyr::left_join(tmp.out , by=c('gene_uid'='gene_uid'), keep=F,suffix = c("", "")) # force overwrite
# 
# nrow(expression.glass.exon.metadata) == nrow(expression.glass.exon)
# stopifnot(rownames(expression.glass.exon) == expression.glass.exon.metadata$gene_uid)
# 
# # expression.glass.exon.metadata %>% dplyr::arrange(cor.t.dna.shallow.ACE.purity.below.1) %>% head(n=10)
# # expression.glass.exon.metadata %>% dplyr::arrange(-cor.t.dna.shallow.ACE.purity.below.1) %>% head(n=10)
# # 
# # 
# # expression.glass.exon.metadata %>% dplyr::arrange(cor.t.dna.shallow.ACE.purity.below.1) %>% head(n=10) %>% dplyr::pull(gene_name)
# # expression.glass.exon.metadata %>% dplyr::arrange(-cor.t.dna.shallow.ACE.purity.below.1) %>% head(n=10) %>% dplyr::pull(gene_name)
# # 


#rm(tmp.metadata, tmp.data, tmp.out)


## Erik/Manual method ----

# tmp.metadata <- metadata.glass.per.resection %>%
#   dplyr::filter(excluded == F) %>% 
#   dplyr::filter(!is.na(dna.purity.manual.Erik)) %>% 
#   dplyr::select(genomescan.sid, dna.purity.manual.Erik)
# 
# 
# tmp.data <- expression.glass.exon.vst %>% 
#   dplyr::select(tmp.metadata$genomescan.sid)
# 
# 
# stopifnot(tmp.metadata$genomescan.sid == colnames(tmp.data))
# 
# 
# tmp.out <- data.frame(cor.t.dna.purity.manual.Erik = apply(tmp.data,1, function(vec) {return( cor.test(tmp.metadata %>% dplyr::pull(dna.purity.manual.Erik), as.numeric(vec))$statistic) })) %>% 
#   tibble::rownames_to_column('gene_uid')
# 
# 
# nrow(expression.glass.exon.metadata) == nrow(expression.glass.exon)
# stopifnot(rownames(expression.glass.exon) == expression.glass.exon.metadata$gene_uid)
# 
# expression.glass.exon.metadata <- expression.glass.exon.metadata %>% 
#   dplyr::left_join(tmp.out , by=c('gene_uid'='gene_uid'), keep=F,suffix = c("", "")) # force overwrite
# 
# 
# nrow(expression.glass.exon.metadata) == nrow(expression.glass.exon)
# stopifnot(rownames(expression.glass.exon) == expression.glass.exon.metadata$gene_uid)
# 
# # expression.glass.exon.metadata %>% dplyr::arrange(cor.t.dna.purity.manual.Erik) %>% head(n=10)
# # expression.glass.exon.metadata %>% dplyr::arrange(-cor.t.dna.purity.manual.Erik) %>% head(n=10)
# # 
# # 
# # expression.glass.exon.metadata %>% dplyr::arrange(cor.t.dna.purity.manual.Erik) %>% head(n=10) %>% dplyr::pull(gene_name)
# # expression.glass.exon.metadata %>% dplyr::arrange(-cor.t.dna.purity.manual.Erik) %>% head(n=10) %>% dplyr::pull(gene_name)
# 
# 
# rm(tmp.metadata, tmp.data, tmp.out)



## IDH VAF method ----


tmp.metadata <- metadata.glass.per.resection |> 
  dplyr::filter(excluded == F) |> 
  dplyr::filter(!is.na(genomescan.sid)) |> # proteomics-only sample
  dplyr::filter(!is.na(dna.wes.VAF_IDH))


tmp.data <- expression.glass.exon.vst |> 
  dplyr::select(tmp.metadata$genomescan.sid)


stopifnot(tmp.metadata$genomescan.sid == colnames(tmp.data))



tmp.out <- do.call(
  rbind,
  pbapply::pbapply(tmp.data, 1, function(vec) {
    c <- cor.test(tmp.metadata$dna.wes.VAF_IDH, as.numeric(vec))

    return(list(
      "cor.dna.wes.VAF_IDH.t.statistic" = as.numeric(c$statistic),
      "cor.dna.wes.VAF_IDH.p.value" = as.numeric(c$p.value),
      "cor.dna.wes.VAF_IDH.cor.estimate" = as.numeric(c$estimate)
    ))
  })
) |>
  as.data.frame() |> 
  tibble::rownames_to_column("gene_uid")


nrow(expression.glass.exon.metadata) == nrow(expression.glass.exon)
stopifnot(rownames(expression.glass.exon) == expression.glass.exon.metadata$gene_uid)


out <- out |> # keep=F
   dplyr::left_join(tmp.out , by=c('gene_uid'='gene_uid'),suffix = c("", "")) # force overwrite



rm(tmp.metadata, tmp.data, tmp.out)




# calc corr for Meth/RFpurity method ----


tmp.metadata <- metadata.glass.per.resection |> 
  dplyr::filter(excluded == F) |> 
  dplyr::filter(!is.na(genomescan.sid)) |> # proteomics-only sample
  dplyr::filter(!is.na(methylation.purity.absolute))


tmp.data <- expression.glass.exon.vst |> 
  dplyr::select(tmp.metadata$genomescan.sid)


stopifnot(tmp.metadata$genomescan.sid == colnames(tmp.data))



tmp.out <- do.call(
  rbind,
  pbapply::pbapply(tmp.data, 1, function(vec) {
    c <- cor.test(tmp.metadata$methylation.purity.absolute, as.numeric(vec))
    
    return(list(
      "cor.methylation.purity.absolute.t.statistic" = as.numeric(c$statistic),
      "cor.methylation.purity.absolute.p.value" = as.numeric(c$p.value),
      "cor.methylation.purity.absolute.cor.estimate" = as.numeric(c$estimate)
    ))
  })
) |>
  as.data.frame() |> 
  tibble::rownames_to_column("gene_uid")


nrow(expression.glass.exon.metadata) == nrow(expression.glass.exon)
stopifnot(rownames(expression.glass.exon) == expression.glass.exon.metadata$gene_uid)


out <- out |> # keep=F
  dplyr::left_join(tmp.out , by=c('gene_uid'='gene_uid'),suffix = c("", "")) # force overwrite


rm(tmp.metadata, tmp.data, tmp.out)




# export -----


out <- out |> 
  tibble::column_to_rownames('gene_uid') |> 
  dplyr::mutate_all(as.numeric) |> 
  tibble::rownames_to_column('gene_uid')


saveRDS(out, file="cache/correlation_expression_purity.Rds")





