#!usr/bin/env R

# load libs ----

# obtain metadata ---


if(!exists('dge.partially.paired.clusters')) {
  source('scripts/load_hclust.R')
}


if(!exists('expression.proteomics.metadata')) {
  source('scripts/load_proteomics_expression.R')
}


if(!exists('theme_cellpress')) {
  source('scripts/R/youri_gg_theme.R')
}


# DPE: primary - recurrence ----
## limma - imputed [cycling] ----


tmp.metadata <- readRDS('cache/meta.proteomics.proteomics.primary__recurrence.Rds')
plt <- readRDS('cache/res.proteomics.primary__recurrence.Rds') |> 
  tibble::rownames_to_column('hugo_symbol') |> 
  dplyr::mutate(cell.cycling = hugo_symbol %in% c("ANLN","ANP32E","ARHGAP11A","ARL6IP1","ASF1B","ASPM","ATAD2","AURKA","AURKB","BIRC5","BLM","BRIP1","BUB1","BUB1B","CASP8AP2","CBX5","CCNA2","CCNB1","CCNB2","CCNE2","CDC20","CDC25B","CDC25C","CDC45","CDC6","CDCA2","CDCA3","CDCA5","CDCA7","CDCA8","CDK1","CDKN3","CENPA","CENPE","CENPF","CENPK","CENPM","CENPW","CHAF1B","CKAP2","CKAP2L","CKAP5","CKS1B","CKS2","CLSPN","CTCF","DEK","DHFR","DLGAP5","DNMT1","DSCC1","DSN1","DTL","DTYMK","DUT","E2F8","ECT2","EXO1","EZH2","FABP5","FAM64A","FANCI","FEN1","FOXM1","G2E3","GAS2L3","GINS2","GMNN","GPSM2","GTSE1","H2AFZ","HAT1","HELLS","HIST1H4C","HJURP","HMGB1","HMGB2","HMGB3","HMMR","HN1","KIAA0101","KIF11","KIF20A","KIF20B","KIF22","KIF23","KIF2C","KIF4A","KIFC1","KNSTRN","KPNA2","LBR","LMNB1","MAD2L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MELK","MKI67","MLF1IP","MND1","MSH2","MXD3","MZT1","NASP","NCAPD2","NDC80","NEK2","NUDT1","NUF2","NUSAP1","OIP5","ORC6","PBK","PCNA","PHF19","PKMYT1","PLK1","POLA1","POLA2","POLD3","PRIM1","PSRC1","PTTG1","RACGAP1","RAD51","RAD51AP1","RANBP1","RANGAP1","REEP4","RFC2","RFC3","RFC4","RFC5","RNASEH2A","RPA2","RPL39L","RRM1","RRM2","SAE1","SDF2L1","SHCBP1","SLBP","SMC4","SNRNP25","SPAG5","TACC3","TCF19","TIMELESS","TIPIN","TK1","TMEM106C","TMEM194A","TMPO","TOP2A","TPX2","TROAP","TTK","TUBA1B","TUBA1C","TUBB4B","TUBB6","TUBG1","TYMS","UBE2C","UBE2T","UBR7","UHRF1","UNG","USP1","VRK1","WDR34","WDR76","ZWILCH","ZWINT"))

ggplot(plt, aes(x=logFC, y=-log10(adj.P.Val), col=cell.cycling, label=hugo_symbol)) +
  geom_point(data = plt |>  dplyr::filter(cell.cycling==F),cex=0.5) +
  geom_point(data = plt |>  dplyr::filter(cell.cycling==T)) +
  geom_hline(yintercept = -log10(0.05),col="black", lty=2, size=0.5/2.14) +
  ggrepel::geom_text_repel(data = plt |>  dplyr::filter(cell.cycling==T), col="black",alpha=0.55,
                           nudge_y=-0.25,
                           segment.size=0.5 / 2.14,
                           size=7 * (3.88/11)
  ) +
  annotate(geom="text", x=-0.75, y=-log10(0.065), label="Padj = 0.05", color="black",
           size=7 * (3.88/11)
  ) +
  labs(x = "logFC proteomics",
       y = "-log10(adjusted P-value limma)",
       caption = paste0("n=",length(tmp.metadata$Sample_Type),
                        "  (initial: ",sum(tmp.metadata$Sample_Type == "initial"),
                        ", recurrent: ",sum(tmp.metadata$Sample_Type == "recurrent")
                        ,")"),
       title = "Differential proteomics: primary  -  recurrence"
  ) +
  theme_cellpress + # in compliance with most journals
  theme(legend.position = 'bottom', 
        legend.key.height = unit(0, 'cm'))



ggsave("output/figures/vis_DPE_primary_recurrence.pdf", width=8.5/2, height=11/3)


rm(plt, tmp.metadata)

## limma - imputed [collagen/c1?] ----


tmp.metadata <- readRDS('cache/meta.proteomics.proteomics.primary__recurrence.Rds')
plt <- readRDS('cache/res.proteomics.primary__recurrence.Rds') |> 
  tibble::rownames_to_column('hugo_symbol__proteomics') |> 
  dplyr::left_join(
    expression.proteomics.metadata |> 
      tibble::rownames_to_column('hugo_symbol__proteomics') |> 
      dplyr::select(hugo_symbol__proteomics, hugo_symbol_rna_prot_shared),
    by=c('hugo_symbol__proteomics'='hugo_symbol__proteomics')
  ) |> 
  dplyr::left_join(
    dge.partially.paired.clusters, by=c('hugo_symbol_rna_prot_shared'='gene_name')
  ) |> 
  dplyr::mutate(up.2 = ifelse(is.na(up.2), F, up.2))
stopifnot(nrow(plt) == 3216)


ggplot(plt, aes(x=logFC, y=-log10(adj.P.Val), col=up.2, label=hugo_symbol__proteomics)) +
  geom_point(data = plt |>  dplyr::filter(up.2==F),cex=0.5) +
  geom_point(data = plt |>  dplyr::filter(up.2==T)) +
  geom_hline(yintercept = -log10(0.05),col="black", lty=2, size=0.5/2.14) +
  ggrepel::geom_text_repel(data = plt |>  dplyr::filter(up.2==T), col="black",alpha=0.55,
                           nudge_y=-0.25,
                           segment.size=0.5 / 2.14,
                           size=7 * (3.88/11)
  ) +
  annotate(geom="text", x=-0.75, y=-log10(0.065), label="Padj = 0.05", color="black",
           size=7 * (3.88/11)
  ) +
  labs(x = "logFC proteomics",
       y = "-log10(adjusted P-value limma)",
       caption = paste0("n=",length(tmp.metadata$Sample_Type),
                        "  (initial: ",sum(tmp.metadata$Sample_Type == "initial"),
                        ", recurrent: ",sum(tmp.metadata$Sample_Type == "recurrent")
                        ,")"),
       title = "Differential proteomics: primary  -  recurrence"
  ) +
  theme_cellpress + # in compliance with most journals
  theme(legend.position = 'bottom', 
        legend.key.height = unit(0, 'cm'))



ggsave("output/figures/vis_DPE_primary_recurrence__c2_ECM.pdf", width=8.5/2, height=11/3)


rm(plt, tmp.metadata)





# DPE: A_IDH - A_IDH_HG ----
## limma - imputed ----


tmp.metadata <- readRDS('cache/meta.proteomics.meth_a_idh__a_idh_hg.Rds.Rds')
plt <- readRDS('cache/res.proteomics.meth_a_idh__a_idh_hg.Rds') |> 
  tibble::rownames_to_column('hugo_symbol') |> 
  dplyr::mutate(cell.cycling = hugo_symbol %in% c("ANLN","ANP32E","ARHGAP11A","ARL6IP1","ASF1B","ASPM","ATAD2","AURKA","AURKB","BIRC5","BLM","BRIP1","BUB1","BUB1B","CASP8AP2","CBX5","CCNA2","CCNB1","CCNB2","CCNE2","CDC20","CDC25B","CDC25C","CDC45","CDC6","CDCA2","CDCA3","CDCA5","CDCA7","CDCA8","CDK1","CDKN3","CENPA","CENPE","CENPF","CENPK","CENPM","CENPW","CHAF1B","CKAP2","CKAP2L","CKAP5","CKS1B","CKS2","CLSPN","CTCF","DEK","DHFR","DLGAP5","DNMT1","DSCC1","DSN1","DTL","DTYMK","DUT","E2F8","ECT2","EXO1","EZH2","FABP5","FAM64A","FANCI","FEN1","FOXM1","G2E3","GAS2L3","GINS2","GMNN","GPSM2","GTSE1","H2AFZ","HAT1","HELLS","HIST1H4C","HJURP","HMGB1","HMGB2","HMGB3","HMMR","HN1","KIAA0101","KIF11","KIF20A","KIF20B","KIF22","KIF23","KIF2C","KIF4A","KIFC1","KNSTRN","KPNA2","LBR","LMNB1","MAD2L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MELK","MKI67","MLF1IP","MND1","MSH2","MXD3","MZT1","NASP","NCAPD2","NDC80","NEK2","NUDT1","NUF2","NUSAP1","OIP5","ORC6","PBK","PCNA","PHF19","PKMYT1","PLK1","POLA1","POLA2","POLD3","PRIM1","PSRC1","PTTG1","RACGAP1","RAD51","RAD51AP1","RANBP1","RANGAP1","REEP4","RFC2","RFC3","RFC4","RFC5","RNASEH2A","RPA2","RPL39L","RRM1","RRM2","SAE1","SDF2L1","SHCBP1","SLBP","SMC4","SNRNP25","SPAG5","TACC3","TCF19","TIMELESS","TIPIN","TK1","TMEM106C","TMEM194A","TMPO","TOP2A","TPX2","TROAP","TTK","TUBA1B","TUBA1C","TUBB4B","TUBB6","TUBG1","TYMS","UBE2C","UBE2T","UBR7","UHRF1","UNG","USP1","VRK1","WDR34","WDR76","ZWILCH","ZWINT"))

ggplot(plt, aes(x=logFC, y=-log10(adj.P.Val), col=cell.cycling, label=hugo_symbol)) +
  geom_point(data = plt |>  dplyr::filter(cell.cycling==F),cex=0.5) +
  geom_point(data = plt |>  dplyr::filter(cell.cycling==T)) +
  geom_hline(yintercept = -log10(0.05),col="black", lty=2, size=0.5/2.14) +
  ggrepel::geom_text_repel(data = plt |>  dplyr::filter(cell.cycling==T), col="black",alpha=0.55,
                           nudge_y=-0.25,nudge_x = 0.5,
                           segment.size=0.5 / 2.14,
                           size=7 * (3.88/11)
  ) +
  annotate(geom="text", x=-1.6, y=-log10(0.065), label="Padj = 0.05", color="black",
           size=7 * (3.88/11)
  ) +
  labs(x = "logFC proteomics",
       y = "-log10(adjusted P-value limma)",
       caption = paste0("n=",length(tmp.metadata$Sample_Type),
                        "  (A_IDH: ",sum(tmp.metadata$methylation.sub.diagnosis == "A_IDH"),
                        ", A_IDH_HG: ",sum(tmp.metadata$methylation.sub.diagnosis == "A_IDH_HG")
                        ,")"),
       title = "Differential proteomics: A_IDH  -  A_IDH_HG"
  ) +
  theme_cellpress + # in compliance with most journals
  theme(legend.position = 'bottom', 
        legend.key.height = unit(0, 'cm'))



ggsave("output/figures/vis_DPE_Meth_A_IDH__A_IDH_HG.pdf", width=8.5/2, height=11/3)
rm(plt, tmp.metadata)


## limma - imputed [collagen/c1?] ----


tmp.metadata <- readRDS('cache/meta.proteomics.meth_a_idh__a_idh_hg.Rds')
plt <- readRDS('cache/res.proteomics.meth_a_idh__a_idh_hg.Rds') |> 
  tibble::rownames_to_column('hugo_symbol__proteomics') |> 
  dplyr::left_join(
    expression.proteomics.metadata |> 
      tibble::rownames_to_column('hugo_symbol__proteomics') |> 
      dplyr::select(hugo_symbol__proteomics, hugo_symbol_rna_prot_shared),
    by=c('hugo_symbol__proteomics'='hugo_symbol__proteomics')
  ) |> 
  dplyr::left_join(
    dge.partially.paired.clusters, by=c('hugo_symbol_rna_prot_shared'='gene_name')
  ) |> 
  dplyr::mutate(up.2 = ifelse(is.na(up.2), F, up.2))
stopifnot(nrow(plt) == 3216)


ggplot(plt, aes(x=logFC, y=-log10(adj.P.Val), col=up.2, label=hugo_symbol__proteomics)) +
  geom_point(data = plt |>  dplyr::filter(up.2==F),cex=0.5) +
  geom_point(data = plt |>  dplyr::filter(up.2==T)) +
  geom_hline(yintercept = -log10(0.05),col="black", lty=2, size=0.5/2.14) +
  ggrepel::geom_text_repel(data = plt |>  dplyr::filter(up.2==T), col="black",alpha=0.55,
                           nudge_y=-0.25,nudge_x = 0.5,
                           segment.size=0.5 / 2.14,
                           size=7 * (3.88/11)
  ) +
  annotate(geom="text", x=-1.6, y=-log10(0.065), label="Padj = 0.05", color="black",
           size=7 * (3.88/11)
  ) +
  labs(x = "logFC proteomics",
       y = "-log10(adjusted P-value limma)",
       caption = paste0("n=",length(tmp.metadata$Sample_Type),
                        "  (A_IDH: ",sum(tmp.metadata$methylation.sub.diagnosis == "A_IDH"),
                        ", A_IDH_HG: ",sum(tmp.metadata$methylation.sub.diagnosis == "A_IDH_HG")
                        ,")"),
       title = "Differential proteomics: A_IDH  -  A_IDH_HG"
  ) +
  theme_cellpress + # in compliance with most journals
  theme(legend.position = 'bottom', 
        legend.key.height = unit(0, 'cm'))



ggsave("output/figures/vis_DPE_Meth_A_IDH__A_IDH_HG__c2_ECM.pdf", width=8.5/2, height=11/3)
rm(plt, tmp.metadata)






# DPE: WHO2021: 2 & 3 - WHO2021: 4 ----
## limma - imputed ----


tmp.metadata <- readRDS('cache/meta.proteomics.who21_2_3__4.Rds')
plt <- readRDS("cache/res.proteomics.who21_2_3__4.Rds") |> 
  tibble::rownames_to_column('hugo_symbol') |> 
  dplyr::mutate(cell.cycling = hugo_symbol %in% c("ANLN","ANP32E","ARHGAP11A","ARL6IP1","ASF1B","ASPM","ATAD2","AURKA","AURKB","BIRC5","BLM","BRIP1","BUB1","BUB1B","CASP8AP2","CBX5","CCNA2","CCNB1","CCNB2","CCNE2","CDC20","CDC25B","CDC25C","CDC45","CDC6","CDCA2","CDCA3","CDCA5","CDCA7","CDCA8","CDK1","CDKN3","CENPA","CENPE","CENPF","CENPK","CENPM","CENPW","CHAF1B","CKAP2","CKAP2L","CKAP5","CKS1B","CKS2","CLSPN","CTCF","DEK","DHFR","DLGAP5","DNMT1","DSCC1","DSN1","DTL","DTYMK","DUT","E2F8","ECT2","EXO1","EZH2","FABP5","FAM64A","FANCI","FEN1","FOXM1","G2E3","GAS2L3","GINS2","GMNN","GPSM2","GTSE1","H2AFZ","HAT1","HELLS","HIST1H4C","HJURP","HMGB1","HMGB2","HMGB3","HMMR","HN1","KIAA0101","KIF11","KIF20A","KIF20B","KIF22","KIF23","KIF2C","KIF4A","KIFC1","KNSTRN","KPNA2","LBR","LMNB1","MAD2L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MELK","MKI67","MLF1IP","MND1","MSH2","MXD3","MZT1","NASP","NCAPD2","NDC80","NEK2","NUDT1","NUF2","NUSAP1","OIP5","ORC6","PBK","PCNA","PHF19","PKMYT1","PLK1","POLA1","POLA2","POLD3","PRIM1","PSRC1","PTTG1","RACGAP1","RAD51","RAD51AP1","RANBP1","RANGAP1","REEP4","RFC2","RFC3","RFC4","RFC5","RNASEH2A","RPA2","RPL39L","RRM1","RRM2","SAE1","SDF2L1","SHCBP1","SLBP","SMC4","SNRNP25","SPAG5","TACC3","TCF19","TIMELESS","TIPIN","TK1","TMEM106C","TMEM194A","TMPO","TOP2A","TPX2","TROAP","TTK","TUBA1B","TUBA1C","TUBB4B","TUBB6","TUBG1","TYMS","UBE2C","UBE2T","UBR7","UHRF1","UNG","USP1","VRK1","WDR34","WDR76","ZWILCH","ZWINT"))


update_geom_defaults("text", list(size = 7 * (3.88/11)))
ggplot(plt, aes(x=logFC, y=-log10(adj.P.Val), col=cell.cycling, label=hugo_symbol)) +
  geom_point(data = plt |>  dplyr::filter(cell.cycling==F),cex=0.5) +
  geom_point(data = plt |>  dplyr::filter(cell.cycling==T)) +
  geom_hline(yintercept = -log10(0.05),col="black", lty=2, size=0.5/2.14) +
  ggrepel::geom_text_repel(data = plt |>  dplyr::filter(cell.cycling==T), col="black",alpha=0.55,
                           nudge_y=-0.25,nudge_x = 0.5,
                           segment.size=0.5 / 2.14,
                           size=7 * (3.88/11)
                           ) +
  annotate(geom="text", x=-1.6, y=-log10(0.065), label="Padj = 0.05", color="black",
           size=7 * (3.88/11)
           ) +
  labs(x = "logFC proteomics",
       y = "-log10(adjusted P-value limma)",
       caption = paste0("n=",length(tmp.metadata$Sample_Type),
                        "  (WHO2021 2 & 3: ",sum(tmp.metadata$WHO_Classification2021 == "WHO2021_g23"),
                        ", WHO2021 4: ",sum(tmp.metadata$WHO_Classification2021 == "WHO2021_g4")
                        ,")"),
       title = "Differential proteomics: WHO2021 grade 2 & 3  -  WHO2021 grade 4"
       ) +
  theme_cellpress + # in compliance with most journals
  theme(legend.position = 'bottom', 
        legend.key.height = unit(0, 'cm'))



ggsave("output/figures/vis_DPE_WHO2021.pdf", width=8.5/2, height=11/3)



## limma - imputed [collagen/c1?] ----


tmp.metadata <- readRDS('cache/meta.proteomics.who21_2_3__4.Rds')
plt <- readRDS("cache/res.proteomics.who21_2_3__4.Rds") |> 
  tibble::rownames_to_column('hugo_symbol__proteomics') |> 
  dplyr::left_join(
    expression.proteomics.metadata |> 
      tibble::rownames_to_column('hugo_symbol__proteomics') |> 
      dplyr::select(hugo_symbol__proteomics, hugo_symbol_rna_prot_shared),
    by=c('hugo_symbol__proteomics'='hugo_symbol__proteomics')
  ) |> 
  dplyr::left_join(
    dge.partially.paired.clusters, by=c('hugo_symbol_rna_prot_shared'='gene_name')
  ) |> 
  dplyr::mutate(up.2 = ifelse(is.na(up.2), F, up.2))
stopifnot(nrow(plt) == 3216)


ggplot(plt, aes(x=logFC, y=-log10(adj.P.Val), col=up.2, label=hugo_symbol__proteomics)) +
  geom_point(data = plt |>  dplyr::filter(up.2==F),cex=0.5) +
  geom_point(data = plt |>  dplyr::filter(up.2==T)) +
  geom_hline(yintercept = -log10(0.05),col="black", lty=2, size=0.5/2.14) +
  ggrepel::geom_text_repel(data = plt |>  dplyr::filter(up.2==T), col="black",alpha=0.55,
                           nudge_y=-0.25,nudge_x = 0.5,
                           segment.size=0.5 / 2.14,
                           size=7 * (3.88/11)
  ) +
  annotate(geom="text", x=-1.6, y=-log10(0.065), label="Padj = 0.05", color="black",
           size=7 * (3.88/11)
  ) +
  labs(x = "logFC proteomics",
       y = "-log10(adjusted P-value limma)",
       caption = paste0("n=",length(tmp.metadata$Sample_Type),
                        "  (WHO2021 2 & 3: ",sum(tmp.metadata$WHO_Classification2021 == "WHO2021_g23"),
                        ", WHO2021 4: ",sum(tmp.metadata$WHO_Classification2021 == "WHO2021_g4")
                        ,")"),
       title = "Differential proteomics: WHO2021 grade 2 & 3  -  WHO2021 grade 4"
  ) +
  theme_cellpress + # in compliance with most journals
  theme(legend.position = 'bottom', 
        legend.key.height = unit(0, 'cm'))



ggsave("output/figures/vis_DPE_WHO2021__c2_ECM.pdf", width=8.5/2, height=11/3)






