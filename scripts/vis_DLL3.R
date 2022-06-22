#!/usr/bin/env R

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


# make plot ----


metadata <- metadata.glass.per.resection %>%
  dplyr::filter(excluded != T )

expr.data <- expression.glass.exon.vst %>%
  dplyr::select(metadata$genomescan.sid)

stopifnot(metadata$genomescan.sid == colnames(expr.data))


saveRDS(metadata,"/tmp/glass-metadata-per-resection.Rds")
saveRDS(expr.data,"/tmp/glass-readcounts-vst-transformed.Rds")





metadata <- readRDS("/tmp/glass-metadata-per-resection.Rds")
expr.data <- readRDS("/tmp/glass-readcounts-vst-transformed.Rds")

stopifnot(metadata$genomescan.sid == colnames(expr.data))


plt <- expr.data %>%
  tibble::rownames_to_column('gene_symbol') %>%
  dplyr::filter(grepl("DLL3", gene_symbol)) %>%
  tibble::column_to_rownames('gene_symbol') %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column('genomescan.sid') %>%
  dplyr::left_join(metadata,by=c('genomescan.sid'='genomescan.sid'))


ggplot(plt, aes(x = resection, y=ENSG00000090932_DLL3, group = GLASS_ID)) + 
  geom_line(alpha=0.6) +
  geom_point() + 
  theme_bw()


ggsave("/tmp/GLASS-DLL3-expression.pdf",width=10,height=6)


