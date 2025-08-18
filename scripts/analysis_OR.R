#!/usr/bin/env R

if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if(!exists("expression.glass.exon.vst")) {
  source('scripts/load_rna-counts.R')
}



# get meth data ----



# get RNA data ----


a <- expression.glass.exon |>
  tibble::rownames_to_column('gene_uid') |> 
  dplyr::filter(grepl("_OR[0-9]", gene_uid)) |> 
  tibble::column_to_rownames('gene_uid') |> 
  as.matrix() |> 
  rowMeans() |> 
  data.frame() |> 
  dplyr::rename(mean_count = 1) |> 
  tibble::rownames_to_column('gene_uid')


b <- readRDS("cache/res.paired.a.exon.Rds")


plt <- a |> 
  dplyr::left_join(b, by=c('gene_uid'='gene_uid'))


ggplot(plt, aes(y=mean_count, x=log2FoldChange)) +
  geom_point() +
  ylim(0,35)


out = metadata.glass.per.resection |> dplyr::select(genomescan.sid, methylation.sid) |> dplyr::filter(!is.na(genomescan.sid) | !is.na(methylation.sid))
write.csv(out, file="/tmp/GLASS-NL_identifiers.csv")



