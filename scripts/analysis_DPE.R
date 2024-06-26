#!usr/bin/env R

# load libs ----

# obtain metadata ---


if(!exists('metadata.glass.per.resection')) {
  source('scripts/load_metadata.R')
}


if(!exists('expression.proteomics.normalised.imputed')) {
  source('scripts/load_proteomics_expression.R')
}


# DPE: primary - recurrence ----
## limma - imputed ----



if(!file.exists("cache/res.proteomics.primary__recurrence.Rds")) {
  tmp.metadata <- metadata.glass.per.resection |>
    dplyr::filter(Sample_Name %in% c(
      
      metadata.glass.per.patient |>
        tidyr::pivot_longer(c('proteomics.sid.I', 'proteomics.sid.R')) |>
        dplyr::filter(!is.na(value)) |>
        dplyr::pull(value)
      
    )) |> 
    dplyr::filter(!is.na(ProtID)) |> 
    dplyr::mutate(Sample_Type = factor(Sample_Type, levels=c('initial','recurrent')))
  stopifnot(nrow(tmp.metadata) == 55)
  
  
  tmp.data <- expression.proteomics.normalised.imputed |> 
    dplyr::select(tmp.metadata$Sample_Name)
  
  
  design <- model.matrix(~Sample_Type, data=tmp.metadata)
  fit.proteomics.primary__recurrence <- limma::lmFit(tmp.data, design)
  res.proteomics.primary__recurrence <- limma::eBayes(fit.proteomics.primary__recurrence, trend = TRUE)
  res.proteomics.primary__recurrence <- limma::topTable(res.proteomics.primary__recurrence, adjust.method="fdr",n=Inf)
  
  saveRDS(res.proteomics.primary__recurrence, "cache/res.proteomics.primary__recurrence.Rds")
  saveRDS(tmp.metadata, 'cache/meta.proteomics.proteomics.primary__recurrence.Rds')
}




## limma - partial ----
# makes barely difference


metadata.glass.per.patient |>
  tidyr::pivot_longer(c('proteomics.sid.I', 'proteomics.sid.R'))

tmp.metadata <- metadata.glass.per.resection |>
  dplyr::filter(!is.na(ProtID)) |> 
  dplyr::mutate(Sample_Type = factor(Sample_Type, levels=c('initial','recurrent')))
stopifnot(nrow(tmp.metadata) == 55)

tmp.data <- expression.proteomics.normalised.partial |> 
  dplyr::select(tmp.metadata$Sample_Name)
stopifnot(ncol(tmp.data) == 55)


design <- model.matrix(~Sample_Type, data=tmp.metadata)
fit <- limma::lmFit(tmp.data, design)
fit2 <- limma::eBayes(fit, trend = TRUE)


# plt <- limma::topTable(fit2, adjust.method="fdr",n=Inf) |> 
#   tibble::rownames_to_column('hugo_symbol') |> 
#   dplyr::mutate(cell.cycling = hugo_symbol %in% c("ANLN","ANP32E","ARHGAP11A","ARL6IP1","ASF1B","ASPM","ATAD2","AURKA","AURKB","BIRC5","BLM","BRIP1","BUB1","BUB1B","CASP8AP2","CBX5","CCNA2","CCNB1","CCNB2","CCNE2","CDC20","CDC25B","CDC25C","CDC45","CDC6","CDCA2","CDCA3","CDCA5","CDCA7","CDCA8","CDK1","CDKN3","CENPA","CENPE","CENPF","CENPK","CENPM","CENPW","CHAF1B","CKAP2","CKAP2L","CKAP5","CKS1B","CKS2","CLSPN","CTCF","DEK","DHFR","DLGAP5","DNMT1","DSCC1","DSN1","DTL","DTYMK","DUT","E2F8","ECT2","EXO1","EZH2","FABP5","FAM64A","FANCI","FEN1","FOXM1","G2E3","GAS2L3","GINS2","GMNN","GPSM2","GTSE1","H2AFZ","HAT1","HELLS","HIST1H4C","HJURP","HMGB1","HMGB2","HMGB3","HMMR","HN1","KIAA0101","KIF11","KIF20A","KIF20B","KIF22","KIF23","KIF2C","KIF4A","KIFC1","KNSTRN","KPNA2","LBR","LMNB1","MAD2L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MELK","MKI67","MLF1IP","MND1","MSH2","MXD3","MZT1","NASP","NCAPD2","NDC80","NEK2","NUDT1","NUF2","NUSAP1","OIP5","ORC6","PBK","PCNA","PHF19","PKMYT1","PLK1","POLA1","POLA2","POLD3","PRIM1","PSRC1","PTTG1","RACGAP1","RAD51","RAD51AP1","RANBP1","RANGAP1","REEP4","RFC2","RFC3","RFC4","RFC5","RNASEH2A","RPA2","RPL39L","RRM1","RRM2","SAE1","SDF2L1","SHCBP1","SLBP","SMC4","SNRNP25","SPAG5","TACC3","TCF19","TIMELESS","TIPIN","TK1","TMEM106C","TMEM194A","TMPO","TOP2A","TPX2","TROAP","TTK","TUBA1B","TUBA1C","TUBB4B","TUBB6","TUBG1","TYMS","UBE2C","UBE2T","UBR7","UHRF1","UNG","USP1","VRK1","WDR34","WDR76","ZWILCH","ZWINT"))
# 
# ggplot(plt, aes(x=logFC, y=-log10(adj.P.Val), col=cell.cycling, label=hugo_symbol)) +
#   geom_point(data = plt |>  dplyr::filter(cell.cycling==F),cex=0.5) +
#   geom_point(data = plt |>  dplyr::filter(cell.cycling==T)) +
#   geom_hline(yintercept = -log10(0.05),col="red", lty=2) +
#   ggrepel::geom_text_repel(data = plt |>  dplyr::filter(cell.cycling==T), col="black",alpha=0.35,nudge_y=-0.25) +
#   theme_bw() +
#   annotate(geom="text", x=-1.0, y=-log10(0.05 + 0.004), label="Padj = 0.05", color="black") +
#   labs(x = "logFC protein [imputed data excluded] (primary - recurrence)",
#        y = "-log10(adjusted P-value limma)",
#        caption = paste0("Differential protein analysis: n=",length(tmp.metadata$Sample_Type),
#                         "  (initial: ",sum(tmp.metadata$Sample_Type == "initial"),
#                         ", recurrent: ",sum(tmp.metadata$Sample_Type == "recurrent")
#                         ,")"))
# 




## wilcox primary recurrence partial ----



tmp.metadata <- metadata.glass.per.resection |>
  dplyr::filter(!is.na(ProtID)) |> 
  dplyr::mutate(Sample_Type = factor(Sample_Type, levels=c('initial','recurrent')))
stopifnot(nrow(tmp.metadata) == 55)


tmp.data <- expression.proteomics.normalised.partial |> 
  dplyr::select(tmp.metadata$Sample_Name)



wtest <- function(expr, mdata) {
  hs <- expr |> 
    dplyr::select(`hugo_symbol`)
  expr <- expr |> 
    dplyr::select(-hugo_symbol)
  
  stopifnot(names(tmp.data[1,]) == names(data))
  
  mdata <- mdata[!is.na(expr)]
  expr <- expr[!is.na(expr)]
  
  vec1 <- expr[mdata == levels(mdata)[1]]
  vec2 <- expr[mdata == levels(mdata)[2]]
  l = length(vec1) + length(vec2)
  
  w = wilcox.test(vec2, vec1)
  
  return(
    data.frame(
               #hugo_symbol = hs,
               p.value = w$p.value,
               statistic = w$statistic,
               stat.norm = w$statistic / ( (l*((l+1)/4))/2)
               )
  )
}


# plt <- tmp.data |> 
#   tibble::rownames_to_column('hugo_symbol') |> 
#   dplyr::rowwise() |> 
#   dplyr::mutate(wilcox = wtest(across(), tmp.metadata |> dplyr::pull(Sample_Type, name=Sample_Name))) |> 
#   dplyr::ungroup() |> 
#   tidyr::unnest(wilcox) |> 
#   dplyr::mutate(p.adj = p.adjust(p.value, method="fdr")) |> 
#   dplyr::mutate(cell.cycling = hugo_symbol %in% c("ANLN","ANP32E","ARHGAP11A","ARL6IP1","ASF1B","ASPM","ATAD2","AURKA","AURKB","BIRC5","BLM","BRIP1","BUB1","BUB1B","CASP8AP2","CBX5","CCNA2","CCNB1","CCNB2","CCNE2","CDC20","CDC25B","CDC25C","CDC45","CDC6","CDCA2","CDCA3","CDCA5","CDCA7","CDCA8","CDK1","CDKN3","CENPA","CENPE","CENPF","CENPK","CENPM","CENPW","CHAF1B","CKAP2","CKAP2L","CKAP5","CKS1B","CKS2","CLSPN","CTCF","DEK","DHFR","DLGAP5","DNMT1","DSCC1","DSN1","DTL","DTYMK","DUT","E2F8","ECT2","EXO1","EZH2","FABP5","FAM64A","FANCI","FEN1","FOXM1","G2E3","GAS2L3","GINS2","GMNN","GPSM2","GTSE1","H2AFZ","HAT1","HELLS","HIST1H4C","HJURP","HMGB1","HMGB2","HMGB3","HMMR","HN1","KIAA0101","KIF11","KIF20A","KIF20B","KIF22","KIF23","KIF2C","KIF4A","KIFC1","KNSTRN","KPNA2","LBR","LMNB1","MAD2L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MELK","MKI67","MLF1IP","MND1","MSH2","MXD3","MZT1","NASP","NCAPD2","NDC80","NEK2","NUDT1","NUF2","NUSAP1","OIP5","ORC6","PBK","PCNA","PHF19","PKMYT1","PLK1","POLA1","POLA2","POLD3","PRIM1","PSRC1","PTTG1","RACGAP1","RAD51","RAD51AP1","RANBP1","RANGAP1","REEP4","RFC2","RFC3","RFC4","RFC5","RNASEH2A","RPA2","RPL39L","RRM1","RRM2","SAE1","SDF2L1","SHCBP1","SLBP","SMC4","SNRNP25","SPAG5","TACC3","TCF19","TIMELESS","TIPIN","TK1","TMEM106C","TMEM194A","TMPO","TOP2A","TPX2","TROAP","TTK","TUBA1B","TUBA1C","TUBB4B","TUBB6","TUBG1","TYMS","UBE2C","UBE2T","UBR7","UHRF1","UNG","USP1","VRK1","WDR34","WDR76","ZWILCH","ZWINT"))
# 
# 
# ggplot(plt, aes(x=stat.norm, y=-log10(p.adj), col=cell.cycling, label=hugo_symbol)) +
#   geom_point(data = plt |>  dplyr::filter(cell.cycling==F),cex=0.5) +
#   geom_point(data = plt |>  dplyr::filter(cell.cycling==T)) +
#   geom_hline(yintercept = -log10(0.05),col="red", lty=2) +
#   ggrepel::geom_text_repel(data = plt |>  dplyr::filter(cell.cycling==T), col="black",alpha=0.35,nudge_y=-0.35) +
#   theme_bw() +
#   annotate(geom="text", x=0.4, y=-log10(0.05 + 0.004), label="Padj = 0.05", color="black") +
#   labs(x = "Wilcox stat [normlised] (primary - recurrence) [excluding imputed values]",
#        y = "-log10(adjusted P-value unpaired wilcox.test)",
#        caption = paste0("Differential protein analysis: n=",length(tmp.metadata$Sample_Type),
#                         "  (initial: ",sum(tmp.metadata$Sample_Type == "initial"),
#                         ", recurrent: ",sum(tmp.metadata$Sample_Type == "recurrent")
#                         ,")"))




## wilcox primary recurrence - imputed ----


tmp.metadata <- metadata.glass.per.resection |>
  dplyr::filter(!is.na(ProtID)) |> 
  dplyr::mutate(resection = factor(ifelse(gsub("^.+(.)$","\\1",Sample_Name) == "1", "primary", "recurrence"), levels=c("primary","recurrence")))

tmp.data <- expression.proteomics.normalised.imputed |> 
  dplyr::select(tmp.metadata$Sample_Name)



wtest <- function(expr, mdata) {
  hs <- expr |> 
    dplyr::select(`hugo_symbol`)
  expr <- expr |> 
    dplyr::select(-hugo_symbol)
  
  stopifnot(names(tmp.data[1,]) == names(data))
  
  mdata <- mdata[!is.na(expr)]
  expr <- expr[!is.na(expr)]
  
  vec1 <- expr[mdata == levels(mdata)[1]]
  vec2 <- expr[mdata == levels(mdata)[2]]
  
  w = wilcox.test(vec2, vec1)
  
  return(
    data.frame(
      #hugo_symbol = hs,
      p.value = w$p.value,
      statistic = w$statistic
    )
  )
}


# plt <- tmp.data |> 
#   tibble::rownames_to_column('hugo_symbol') |> 
#   dplyr::rowwise() |> 
#   dplyr::mutate(wilcox = wtest(across(), tmp.metadata |> dplyr::pull(resection, name=File_Name_Proteomics))) |> 
#   dplyr::ungroup() |> 
#   tidyr::unnest(wilcox) |> 
#   dplyr::mutate(p.adj = p.adjust(p.value, method="fdr")) |> 
#   dplyr::mutate(cell.cycling = hugo_symbol %in% c("ANLN","ANP32E","ARHGAP11A","ARL6IP1","ASF1B","ASPM","ATAD2","AURKA","AURKB","BIRC5","BLM","BRIP1","BUB1","BUB1B","CASP8AP2","CBX5","CCNA2","CCNB1","CCNB2","CCNE2","CDC20","CDC25B","CDC25C","CDC45","CDC6","CDCA2","CDCA3","CDCA5","CDCA7","CDCA8","CDK1","CDKN3","CENPA","CENPE","CENPF","CENPK","CENPM","CENPW","CHAF1B","CKAP2","CKAP2L","CKAP5","CKS1B","CKS2","CLSPN","CTCF","DEK","DHFR","DLGAP5","DNMT1","DSCC1","DSN1","DTL","DTYMK","DUT","E2F8","ECT2","EXO1","EZH2","FABP5","FAM64A","FANCI","FEN1","FOXM1","G2E3","GAS2L3","GINS2","GMNN","GPSM2","GTSE1","H2AFZ","HAT1","HELLS","HIST1H4C","HJURP","HMGB1","HMGB2","HMGB3","HMMR","HN1","KIAA0101","KIF11","KIF20A","KIF20B","KIF22","KIF23","KIF2C","KIF4A","KIFC1","KNSTRN","KPNA2","LBR","LMNB1","MAD2L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MELK","MKI67","MLF1IP","MND1","MSH2","MXD3","MZT1","NASP","NCAPD2","NDC80","NEK2","NUDT1","NUF2","NUSAP1","OIP5","ORC6","PBK","PCNA","PHF19","PKMYT1","PLK1","POLA1","POLA2","POLD3","PRIM1","PSRC1","PTTG1","RACGAP1","RAD51","RAD51AP1","RANBP1","RANGAP1","REEP4","RFC2","RFC3","RFC4","RFC5","RNASEH2A","RPA2","RPL39L","RRM1","RRM2","SAE1","SDF2L1","SHCBP1","SLBP","SMC4","SNRNP25","SPAG5","TACC3","TCF19","TIMELESS","TIPIN","TK1","TMEM106C","TMEM194A","TMPO","TOP2A","TPX2","TROAP","TTK","TUBA1B","TUBA1C","TUBB4B","TUBB6","TUBG1","TYMS","UBE2C","UBE2T","UBR7","UHRF1","UNG","USP1","VRK1","WDR34","WDR76","ZWILCH","ZWINT")) |> 
#   dplyr::mutate(x = statistic + runif(n(), -2.5,2.5))
# 
# 
# ggplot(plt, aes(x=x, y=-log10(p.adj), col=cell.cycling, label=hugo_symbol)) +
#   geom_point(data = plt |>  dplyr::filter(cell.cycling==F),cex=0.5) +
#   geom_point(data = plt |>  dplyr::filter(cell.cycling==T)) +
#   geom_hline(yintercept = -log10(0.05),col="red", lty=2) +
#   ggrepel::geom_text_repel(data = plt |>  dplyr::filter(cell.cycling==T), col="black",alpha=0.35,nudge_y=-0.35) +
#   theme_bw() +
#   annotate(geom="text", x=180, y=-log10(0.05 + 0.004), label="Padj = 0.05", color="black") +
#   labs(x = "Wilcox stat (primary - recurrence) [including imputed values]",
#        y = "-log10(adjusted P-value unpaired wilcox.test)",
#        caption = paste0("Differential protein analysis: n=",length(tmp.metadata$Sample_Type),
#                         "  (initial: ",sum(tmp.metadata$Sample_Type == "initial"),
#                         ", recurrent: ",sum(tmp.metadata$Sample_Type == "recurrent")
#                         ,")"))





# DPE: A_IDH - A_IDH_HG ----
## limma - imputed ----

if(!file.exists("cache/res.proteomics.meth_a_idh__a_idh_hg.Rds")) {
    tmp.metadata <- metadata.glass.per.resection |>
      dplyr::filter(Sample_Name %in% c(
        
        metadata.glass.per.patient |>
          tidyr::pivot_longer(c('proteomics.sid.A_IDH', 'proteomics.sid.A_IDH_HG')) |>
          dplyr::filter(!is.na(value)) |>
          dplyr::pull(value)
        
      )) |> 
      dplyr::filter(!is.na(ProtID)) |> 
      dplyr::mutate(methylation.sub.diagnosis = factor(methylation.sub.diagnosis, levels=c('A_IDH','A_IDH_HG')))
    stopifnot(nrow(tmp.metadata) == 40)


    tmp.data <- expression.proteomics.normalised.imputed |> 
      dplyr::select(tmp.metadata$Sample_Name)
    stopifnot(ncol(tmp.data) == 40)
  
  
  design <- model.matrix(~methylation.sub.diagnosis, data=tmp.metadata)
  fit.proteomics.meth_a_idh__a_idh_hg <- limma::lmFit(tmp.data, design)
  res.proteomics.meth_a_idh__a_idh_hg <- limma::eBayes(fit.proteomics.meth_a_idh__a_idh_hg, trend = TRUE)
  res.proteomics.meth_a_idh__a_idh_hg <- limma::topTable(res.proteomics.meth_a_idh__a_idh_hg, adjust.method="fdr",n=Inf)
  
  saveRDS(res.proteomics.meth_a_idh__a_idh_hg, "cache/res.proteomics.meth_a_idh__a_idh_hg.Rds")
  saveRDS(tmp.metadata, 'cache/meta.proteomics.meth_a_idh__a_idh_hg.Rds')
}




# DPE: WHO2021: 2 & 3 - WHO2021: 4 ----
## limma - imputed ----


if(!file.exists("cache/res.proteomics.who21_2_3__4.Rds")) {
  tmp.metadata <- metadata.glass.per.resection |>
    dplyr::filter(Sample_Name %in% c(
      
      metadata.glass.per.patient |>
        tidyr::pivot_longer(c('proteomics.sid.WHO2021_g23', 'proteomics.sid.WHO2021_g4')) |>
        dplyr::filter(!is.na(value)) |>
        dplyr::pull(value)
      
    )) |> 
    dplyr::filter(!is.na(ProtID)) |> 
    dplyr::mutate(WHO_Classification2021 = case_when(
      WHO_Classification2021 == "Astrocytoma, IDH-mutant, WHO grade 2" ~ "WHO2021_g23",
      WHO_Classification2021 == "Astrocytoma, IDH-mutant, WHO grade 3" ~ "WHO2021_g23",
      WHO_Classification2021 == "Astrocytoma, IDH-mutant, WHO grade 4" ~ "WHO2021_g4"
    )) |> 
  dplyr::mutate(WHO_Classification2021 = factor(WHO_Classification2021, levels=c('WHO2021_g23','WHO2021_g4')))
  stopifnot(nrow(tmp.metadata) == 39+4)
  stopifnot(tmp.metadata$WHO_Classification2021 %in% c("WHO2021_g23","WHO2021_g4"))
  stopifnot(!is.na(tmp.metadata$WHO_Classification2021))
  
  
  tmp.data <- expression.proteomics.normalised.imputed |> 
    dplyr::select(tmp.metadata$Sample_Name)
  stopifnot(ncol(tmp.data) == 39+4)
  
  
  design <- model.matrix(~WHO_Classification2021, data=tmp.metadata)
  fit.proteomics.who21_2_3__4 <- limma::lmFit(tmp.data, design)
  res.proteomics.who21_2_3__4 <- limma::eBayes(fit.proteomics.who21_2_3__4, trend = TRUE)
  res.proteomics.who21_2_3__4 <- limma::topTable(res.proteomics.who21_2_3__4, adjust.method="fdr",n=Inf)
  
  saveRDS(res.proteomics.who21_2_3__4, file="cache/res.proteomics.who21_2_3__4.Rds")
  saveRDS(tmp.metadata, file='cache/meta.proteomics.who21_2_3__4.Rds')
}




