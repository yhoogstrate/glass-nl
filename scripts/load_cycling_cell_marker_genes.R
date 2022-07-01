#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(readODS)


# load data ----


# to match identifiers


# if read counts and per-gene metadata are not loaded, load it
if(!exists('expression.glass.exon.metadata')) {
  warning('read count data was not loaded')
  
  source('scripts/load_rna-counts.R')
}


# cycling reading ----

cycling.cell.markers <- readODS::read_ods('data/scRNA-cycling-genes.ods',col_names=F) %>% 
  dplyr::mutate(A = ifelse(duplicated(A),NA,A)) %>% 
  dplyr::mutate(B = ifelse(duplicated(B),NA,B)) %>% 
  dplyr::mutate(C = ifelse(duplicated(C),NA,C)) %>% 
  dplyr::mutate(E = ifelse(duplicated(E),NA,E)) %>% 
  dplyr::mutate(F = ifelse(duplicated(F),NA,F)) %>% 
  dplyr::rename(G1.S.tirosh = A) %>% 
  dplyr::rename(G2.M.tirosh = B) %>% 
  dplyr::rename(cycling.melanoma.tirosh = C) %>% 
  dplyr::rename(G1.S.neftel = E) %>% 
  dplyr::rename(G2.M.neftel = F) %>% 
  dplyr::mutate(D = NULL) %>% 
  dplyr::mutate(G = NULL) %>% 
  dplyr::mutate(H = NULL) %>%
  tidyr::pivot_longer(cols=colnames(.)) %>% 
  tidyr::drop_na(value) %>% 
  dplyr::filter(value %in% c('10.1126/science.aad0501','10.1016/j.cell.2019.06.024','G1/S phase-specific','G2/M phase-specific','melanoma cell cycle genes','G1/S phase-specific','G2/M phase-specific') == F) %>% 
  dplyr::mutate(marker=T) %>% 
  dplyr::mutate(name = as.factor(name)) %>% 
  dplyr::mutate(value = as.factor(value)) %>% 
  dplyr::rename(source = name) %>% 
  dplyr::rename(gene_name = value) %>% 
  pivot_wider(names_from = source, values_from = marker, values_fill = F) %>% 
  dplyr::mutate(gene_name = as.character(gene_name))




# cycling.cell.markers$gene_name[cycling.cell.markers$gene_name %in% expression.glass.exon.metadata$gene_name == F]

# KIAA0101 -> ENSG00000166803 -> PCLAF
#expression.glass.exon.metadata %>%
#  dplyr::filter(grepl("ENSG00000166803",gene_id)) %>%
#  dplyr::select(gene_id, gene_name, gene_uid)
cycling.cell.markers <- cycling.cell.markers %>% 
  dplyr::mutate(gene_name = gsub("^KIAA0101$","PCLAF",gene_name))


# HIST1H4C -> H4C3
cycling.cell.markers <- cycling.cell.markers %>% 
  dplyr::mutate(gene_name = gsub("^HIST1H4C$","H4C3",gene_name))


# MLF1IP -> ENSG00000151725 -> CENPU
#expression.glass.exon.metadata %>%
#  dplyr::filter(grepl("ENSG00000151725",gene_id)) %>%
#  dplyr::select(gene_id, gene_name, gene_uid)
cycling.cell.markers <- cycling.cell.markers %>% 
  dplyr::mutate(gene_name = gsub("^MLF1IP$","CENPU",gene_name))


# FAM64A -> ENSG00000129195 -> PIMREG
# expression.glass.exon.metadata %>%
#  dplyr::filter(grepl("ENSG00000129195",gene_id)) %>%
#  dplyr::select(gene_id, gene_name, gene_uid)
cycling.cell.markers <- cycling.cell.markers %>% 
  dplyr::mutate(gene_name = gsub("^FAM64A$","PIMREG",gene_name))


# HN1 -> ENSG00000189159 -> JPT1
# expression.glass.exon.metadata %>%
#  dplyr::filter(grepl("ENSG00000189159",gene_id)) %>%
#  dplyr::select(gene_id, gene_name, gene_uid)
cycling.cell.markers <- cycling.cell.markers %>% 
  dplyr::mutate(gene_name = gsub("^HN1$","JPT1",gene_name))


# H2AFZ -> ENSG00000164032 -> H2AZ1
# expression.glass.exon.metadata %>%
#  dplyr::filter(grepl("ENSG00000164032",gene_id)) %>%
#  dplyr::select(gene_id, gene_name, gene_uid)
cycling.cell.markers <- cycling.cell.markers %>% 
  dplyr::mutate(gene_name = gsub("^H2AFZ$","H2AZ1",gene_name))


# TMEM194A -> ENSG00000166881 -> NEMP1
# expression.glass.exon.metadata %>%
#  dplyr::filter(grepl("ENSG00000166881",gene_id)) %>%
#  dplyr::select(gene_id, gene_name, gene_uid)
cycling.cell.markers <- cycling.cell.markers %>% 
  dplyr::mutate(gene_name = gsub("^TMEM194A$","NEMP1",gene_name))


stopifnot(duplicated(cycling.cell.markers$gene_name) == F)
stopifnot(cycling.cell.markers$gene_name %in% expression.glass.exon.metadata$gene_name)
#cycling.cell.markers$gene_name[cycling.cell.markers$gene_name %in% expression.glass.exon.metadata$gene_name == F]



cycling.cell.markers <- cycling.cell.markers %>% 
  dplyr::left_join(
    expression.glass.exon.metadata %>% dplyr::select(gene_name, gene_uid), by=c('gene_name'='gene_name')
  ) %>% 
  tibble::column_to_rownames('gene_name')


stopifnot(duplicated(cycling.cell.markers$gene_uid) == F)


nrow(expression.glass.exon.metadata) == nrow(expression.glass.exon)
stopifnot(rownames(expression.glass.exon) == expression.glass.exon.metadata$gene_uid)

expression.glass.exon.metadata <- expression.glass.exon.metadata %>% 
  dplyr::left_join(cycling.cell.markers , by=c('gene_uid'='gene_uid'), keep=F,suffix = c("", "")) # force overwrite

nrow(expression.glass.exon.metadata) == nrow(expression.glass.exon)
stopifnot(rownames(expression.glass.exon) == expression.glass.exon.metadata$gene_uid)



rm(cycling.cell.markers)




