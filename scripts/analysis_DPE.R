#!usr/bin/env R

# load libs ----

# limma

# obtain metadata ---

tmp.metadata <- metadata.glass.per.resection |>
  dplyr::filter(!is.na(proteomics_imputed_id)) |> 
  dplyr::mutate(resection = factor(ifelse(gsub("^.+(.)$","\\1",Sample_Name) == "1", "primary", "recurrence"), levels=c("primary","recurrence")))

stopifnot(nrow(tmp.metadata) == 55)


stopifnot(tmp.metadata |> 
  dplyr::pull(GLASS_ID) |>
  table() |> 
  max() == 2)


# differential protein expression ----

## limma primary recurrence ----

# tmp.data <- expression.proteomics.normalised.imputed |> 
#   dplyr::select(tmp.metadata$File_Name_Proteomics) 
# 
# pca_res <- prcomp(t(tmp.data), scale. = TRUE)
# plot(pca_res)
# autoplot(pca_res, label=T)
# 
# 
# plt <- pca_res$x |> 
#   as.data.frame() |> 
#   tibble::rownames_to_column('tmp') |> 
#   dplyr::mutate(hugo_symbol = tmp) |> 
#   tibble::column_to_rownames('tmp') |> 
#   dplyr::mutate(cell.cycling = hugo_symbol %in% c('PCNA','MCM7','MCM3','TMPO'))
# ggplot(plt, aes(x=PC1,y=PC2,col=cell.cycling,label=hugo_symbol)) +
#   ggrepel::geom_text_repel(data=subset(plt, PC2 < -4.5)) +
#   geom_point()
# 
# 
# 
# 
# plt <- pca_res$x |> 
#   as.data.frame() |> 
#   tibble::rownames_to_column('tmp') |> 
#   dplyr::mutate(sid = tmp) |> 
#   tibble::column_to_rownames('tmp') |> 
#   dplyr::left_join(tmp.metadata, by=c('sid'='File_Name_Proteomics'), suffix=c('',''))
# 
# ggplot(plt, aes(x=PC1,y=PC2,label=Sample_Name, col=CDKN2AB)) +
#   ggrepel::geom_text_repel(col="black",alpha=0.25) +
#   geom_point()
# 
# ggplot(plt, aes(x=PC1,y=PC2,label=Sample_Name, col=methylation.sub.diagnosis)) +
#   ggrepel::geom_text_repel(col="black",alpha=0.25) +
#   geom_point()
# 
# ggplot(plt, aes(x=PC1,y=PC2,label=Sample_Name, col=resection)) +
#   ggrepel::geom_text_repel(col="black",alpha=0.25) +
#   geom_point()
# 
# 
# 
# df <- data.frame()
# sids <- metadata.glass.per.resection |> dplyr::filter(!is.na(proteomics_imputed_id)) |> dplyr::pull(Sample_Name)
# for(sid in sids) {
#   tmp.metadata <- metadata.glass.per.resection |>
#     dplyr::filter(!is.na(proteomics_imputed_id)) |> 
#     dplyr::filter(Sample_Name != sid) |> 
#     dplyr::mutate(resection = factor(ifelse(gsub("^.+(.)$","\\1",Sample_Name) == "1", "primary", "recurrence"), levels=c("primary","recurrence")))
#   
#   tmp.data <- expression.proteomics.normalised.imputed |> 
#     dplyr::select(tmp.metadata$File_Name_Proteomics)
#   
#   design <- model.matrix(~ resection, data=tmp.metadata)
#   fit <- limma::lmFit(tmp.data, design)
#   fit2 <- limma::eBayes(fit, trend = TRUE)
#   
#   o = limma::topTable(fit2, adjust.method="fdr",n=1)$adj.P.Val
#   
#   df <- rbind(df, 
#               data.frame(
#                 sid=sid,
#                 p=o
#               ))
# }
# 
# df |> dplyr::arrange(p)



tmp.metadata <- metadata.glass.per.resection |>
  dplyr::filter(!is.na(proteomics_imputed_id)) |> 
  dplyr::mutate(resection = factor(ifelse(gsub("^.+(.)$","\\1",Sample_Name) == "1", "primary", "recurrence"), levels=c("primary","recurrence")))

tmp.data <- expression.proteomics.normalised.imputed |> 
  dplyr::select(tmp.metadata$File_Name_Proteomics)

# tmp.data <- tmp.data %>%
#   dplyr::mutate(mad= NULL) %>%
#   dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% # old r style
#   #dplyr::filter(mad > 0.35) |> 
#   dplyr::mutate(mad=NULL)



design <- model.matrix(~ resection, data=tmp.metadata)
fit <- limma::lmFit(tmp.data, design)
fit2 <- limma::eBayes(fit, trend = TRUE)


plt <- limma::topTable(fit2, adjust.method="fdr",n=Inf) |> 
  tibble::rownames_to_column('hugo_symbol') |> 
  dplyr::mutate(cell.cycling = hugo_symbol %in% c("ANLN","ANP32E","ARHGAP11A","ARL6IP1","ASF1B","ASPM","ATAD2","AURKA","AURKB","BIRC5","BLM","BRIP1","BUB1","BUB1B","CASP8AP2","CBX5","CCNA2","CCNB1","CCNB2","CCNE2","CDC20","CDC25B","CDC25C","CDC45","CDC6","CDCA2","CDCA3","CDCA5","CDCA7","CDCA8","CDK1","CDKN3","CENPA","CENPE","CENPF","CENPK","CENPM","CENPW","CHAF1B","CKAP2","CKAP2L","CKAP5","CKS1B","CKS2","CLSPN","CTCF","DEK","DHFR","DLGAP5","DNMT1","DSCC1","DSN1","DTL","DTYMK","DUT","E2F8","ECT2","EXO1","EZH2","FABP5","FAM64A","FANCI","FEN1","FOXM1","G2E3","GAS2L3","GINS2","GMNN","GPSM2","GTSE1","H2AFZ","HAT1","HELLS","HIST1H4C","HJURP","HMGB1","HMGB2","HMGB3","HMMR","HN1","KIAA0101","KIF11","KIF20A","KIF20B","KIF22","KIF23","KIF2C","KIF4A","KIFC1","KNSTRN","KPNA2","LBR","LMNB1","MAD2L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MELK","MKI67","MLF1IP","MND1","MSH2","MXD3","MZT1","NASP","NCAPD2","NDC80","NEK2","NUDT1","NUF2","NUSAP1","OIP5","ORC6","PBK","PCNA","PHF19","PKMYT1","PLK1","POLA1","POLA2","POLD3","PRIM1","PSRC1","PTTG1","RACGAP1","RAD51","RAD51AP1","RANBP1","RANGAP1","REEP4","RFC2","RFC3","RFC4","RFC5","RNASEH2A","RPA2","RPL39L","RRM1","RRM2","SAE1","SDF2L1","SHCBP1","SLBP","SMC4","SNRNP25","SPAG5","TACC3","TCF19","TIMELESS","TIPIN","TK1","TMEM106C","TMEM194A","TMPO","TOP2A","TPX2","TROAP","TTK","TUBA1B","TUBA1C","TUBB4B","TUBB6","TUBG1","TYMS","UBE2C","UBE2T","UBR7","UHRF1","UNG","USP1","VRK1","WDR34","WDR76","ZWILCH","ZWINT"))

ggplot(plt, aes(x=logFC, y=-log10(adj.P.Val), col=cell.cycling, label=hugo_symbol)) +
  geom_point(data = plt |>  dplyr::filter(cell.cycling==F),cex=0.5) +
  geom_point(data = plt |>  dplyr::filter(cell.cycling==T)) +
  geom_hline(yintercept = -log10(0.05),col="red", lty=2) +
  ggrepel::geom_text_repel(data = plt |>  dplyr::filter(cell.cycling==T), col="black",alpha=0.35,nudge_y=0.7) +
  theme_bw() +
  annotate(geom="text", x=-1.0, y=-log10(0.05 + 0.004), label="Padj = 0.05",
           color="black")

head(plt)


## wilcox primary recurrence partial ----



tmp.metadata <- metadata.glass.per.resection |>
  dplyr::filter(!is.na(proteomics_imputed_id)) |> 
  dplyr::mutate(resection = factor(ifelse(gsub("^.+(.)$","\\1",Sample_Name) == "1", "primary", "recurrence"), levels=c("primary","recurrence")))

tmp.data <- expression.proteomics.normalised.partial |> 
  dplyr::select(tmp.metadata$File_Name_Proteomics)



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
  
  w = wilcox.test(vec1, vec2)
  
  return(
    data.frame(
               #hugo_symbol = hs,
               p.value = w$p.value,
               statistic = w$statistic
               )
  )
}


plt <- tmp.data |> 
  tibble::rownames_to_column('hugo_symbol') |> 
  dplyr::rowwise() |> 
  dplyr::mutate(wilcox = wtest(across(), tmp.metadata |> dplyr::pull(resection, name=File_Name_Proteomics))) |> 
  dplyr::ungroup() |> 
  tidyr::unnest(wilcox) |> 
  dplyr::mutate(p.adj = p.adjust(p.value, method="fdr")) |> 
  dplyr::mutate(cell.cycling = hugo_symbol %in% c("ANLN","ANP32E","ARHGAP11A","ARL6IP1","ASF1B","ASPM","ATAD2","AURKA","AURKB","BIRC5","BLM","BRIP1","BUB1","BUB1B","CASP8AP2","CBX5","CCNA2","CCNB1","CCNB2","CCNE2","CDC20","CDC25B","CDC25C","CDC45","CDC6","CDCA2","CDCA3","CDCA5","CDCA7","CDCA8","CDK1","CDKN3","CENPA","CENPE","CENPF","CENPK","CENPM","CENPW","CHAF1B","CKAP2","CKAP2L","CKAP5","CKS1B","CKS2","CLSPN","CTCF","DEK","DHFR","DLGAP5","DNMT1","DSCC1","DSN1","DTL","DTYMK","DUT","E2F8","ECT2","EXO1","EZH2","FABP5","FAM64A","FANCI","FEN1","FOXM1","G2E3","GAS2L3","GINS2","GMNN","GPSM2","GTSE1","H2AFZ","HAT1","HELLS","HIST1H4C","HJURP","HMGB1","HMGB2","HMGB3","HMMR","HN1","KIAA0101","KIF11","KIF20A","KIF20B","KIF22","KIF23","KIF2C","KIF4A","KIFC1","KNSTRN","KPNA2","LBR","LMNB1","MAD2L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MELK","MKI67","MLF1IP","MND1","MSH2","MXD3","MZT1","NASP","NCAPD2","NDC80","NEK2","NUDT1","NUF2","NUSAP1","OIP5","ORC6","PBK","PCNA","PHF19","PKMYT1","PLK1","POLA1","POLA2","POLD3","PRIM1","PSRC1","PTTG1","RACGAP1","RAD51","RAD51AP1","RANBP1","RANGAP1","REEP4","RFC2","RFC3","RFC4","RFC5","RNASEH2A","RPA2","RPL39L","RRM1","RRM2","SAE1","SDF2L1","SHCBP1","SLBP","SMC4","SNRNP25","SPAG5","TACC3","TCF19","TIMELESS","TIPIN","TK1","TMEM106C","TMEM194A","TMPO","TOP2A","TPX2","TROAP","TTK","TUBA1B","TUBA1C","TUBB4B","TUBB6","TUBG1","TYMS","UBE2C","UBE2T","UBR7","UHRF1","UNG","USP1","VRK1","WDR34","WDR76","ZWILCH","ZWINT"))


ggplot(plt, aes(x=statistic, y=-log10(p.adj), col=cell.cycling, label=hugo_symbol)) +
  geom_point(data = plt |>  dplyr::filter(cell.cycling==F),cex=0.5) +
  geom_point(data = plt |>  dplyr::filter(cell.cycling==T)) +
  geom_hline(yintercept = -log10(0.05),col="red", lty=2) +
  ggrepel::geom_text_repel(data = plt |>  dplyr::filter(cell.cycling==T), col="black",alpha=0.35,nudge_y=0.7) +
  theme_bw() +
  annotate(geom="text", x=-1.0, y=-log10(0.05 + 0.004), label="Padj = 0.05", color="black")



## wilcox primary recurrence imputed ----


tmp.metadata <- metadata.glass.per.resection |>
  dplyr::filter(!is.na(proteomics_imputed_id)) |> 
  dplyr::mutate(resection = factor(ifelse(gsub("^.+(.)$","\\1",Sample_Name) == "1", "primary", "recurrence"), levels=c("primary","recurrence")))

tmp.data <- expression.proteomics.normalised.imputed |> 
  dplyr::select(tmp.metadata$File_Name_Proteomics)



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


plt <- tmp.data |> 
  tibble::rownames_to_column('hugo_symbol') |> 
  dplyr::rowwise() |> 
  dplyr::mutate(wilcox = wtest(across(), tmp.metadata |> dplyr::pull(resection, name=File_Name_Proteomics))) |> 
  dplyr::ungroup() |> 
  tidyr::unnest(wilcox) |> 
  dplyr::mutate(p.adj = p.adjust(p.value, method="fdr")) |> 
  dplyr::mutate(cell.cycling = hugo_symbol %in% c("ANLN","ANP32E","ARHGAP11A","ARL6IP1","ASF1B","ASPM","ATAD2","AURKA","AURKB","BIRC5","BLM","BRIP1","BUB1","BUB1B","CASP8AP2","CBX5","CCNA2","CCNB1","CCNB2","CCNE2","CDC20","CDC25B","CDC25C","CDC45","CDC6","CDCA2","CDCA3","CDCA5","CDCA7","CDCA8","CDK1","CDKN3","CENPA","CENPE","CENPF","CENPK","CENPM","CENPW","CHAF1B","CKAP2","CKAP2L","CKAP5","CKS1B","CKS2","CLSPN","CTCF","DEK","DHFR","DLGAP5","DNMT1","DSCC1","DSN1","DTL","DTYMK","DUT","E2F8","ECT2","EXO1","EZH2","FABP5","FAM64A","FANCI","FEN1","FOXM1","G2E3","GAS2L3","GINS2","GMNN","GPSM2","GTSE1","H2AFZ","HAT1","HELLS","HIST1H4C","HJURP","HMGB1","HMGB2","HMGB3","HMMR","HN1","KIAA0101","KIF11","KIF20A","KIF20B","KIF22","KIF23","KIF2C","KIF4A","KIFC1","KNSTRN","KPNA2","LBR","LMNB1","MAD2L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MELK","MKI67","MLF1IP","MND1","MSH2","MXD3","MZT1","NASP","NCAPD2","NDC80","NEK2","NUDT1","NUF2","NUSAP1","OIP5","ORC6","PBK","PCNA","PHF19","PKMYT1","PLK1","POLA1","POLA2","POLD3","PRIM1","PSRC1","PTTG1","RACGAP1","RAD51","RAD51AP1","RANBP1","RANGAP1","REEP4","RFC2","RFC3","RFC4","RFC5","RNASEH2A","RPA2","RPL39L","RRM1","RRM2","SAE1","SDF2L1","SHCBP1","SLBP","SMC4","SNRNP25","SPAG5","TACC3","TCF19","TIMELESS","TIPIN","TK1","TMEM106C","TMEM194A","TMPO","TOP2A","TPX2","TROAP","TTK","TUBA1B","TUBA1C","TUBB4B","TUBB6","TUBG1","TYMS","UBE2C","UBE2T","UBR7","UHRF1","UNG","USP1","VRK1","WDR34","WDR76","ZWILCH","ZWINT")) |> 
  dplyr::mutate(x = statistic + runif(n(), -2.5,2.5))


ggplot(plt, aes(x=x, y=-log10(p.adj), col=cell.cycling, label=hugo_symbol)) +
  geom_point(data = plt |>  dplyr::filter(cell.cycling==F),cex=0.5) +
  geom_point(data = plt |>  dplyr::filter(cell.cycling==T)) +
  geom_hline(yintercept = -log10(0.05),col="red", lty=2) +
  ggrepel::geom_text_repel(data = plt |>  dplyr::filter(cell.cycling==T), col="black",alpha=0.35,nudge_y=0.7) +
  theme_bw() +
  annotate(geom="text", x=180, y=-log10(0.05 + 0.004), label="Padj = 0.05", color="black")

## make pairs protein A_IDH A_IDH_HG

a = metadata.glass.per.resection |> 
  dplyr::filter(!is.na(proteomics_imputed_id)) |> 
  dplyr::select(
    genomescan.sid, Sample_Name, GLASS_ID, resection, Sample_Type, Recurrent_Type, Sample_Sex, Exclude.by.Wies.on.complete.pair, Surgery_ID,
    File_Name_Proteomics, proteomics_imputed_id,
    methylation.sub.diagnosis,
    WHO_Classification2021
  )

write.table(a, "/tmp/a.tsv")




