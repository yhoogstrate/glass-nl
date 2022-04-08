#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(readODS)


# load data ----


cycling.cell.markers <- readODS::read_ods('data/scRNA-cycling-genes.ods',col_names=F) %>% 
  dplyr::mutate(A = ifelse(duplicated(A),NA,A)) %>% 
  dplyr::mutate(B = ifelse(duplicated(B),NA,B)) %>% 
  dplyr::mutate(C = ifelse(duplicated(C),NA,C)) %>% 
  dplyr::mutate(E = ifelse(duplicated(E),NA,E)) %>% 
  dplyr::mutate(F = ifelse(duplicated(F),NA,F)) %>% 
  dplyr::rename(G1.S.tirosh = A) %>% 
  dplyr::rename(G2.M.tirosh = B) %>% 
  dplyr::rename(melanoma.tirosh = C) %>% 
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
  pivot_wider(names_from = source, values_from = marker, values_fill = F)



