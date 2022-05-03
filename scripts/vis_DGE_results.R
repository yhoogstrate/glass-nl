#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(recursiveCorPlot)


# load data ----

 
source('scripts/R/chrom_sizes.R')
library(recursiveCorPlot)


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


# check DGE info
if("padj.partially.paired.exon" %in% colnames(expression.glass.exon.metadata) == F) {
  warning('DGE analysis results were not loaded')
  
  source('scripts/analysis_DGE.R') 
}



if(!exists('cycling.cell.markers')) {
  warning('cycling cell marker genes were not loaded')
  
  source('scripts/load_cycling_cell_marker_genes.R')
}



# recursiveCorPlot [histone genes + cell types + cycling] ----


sel.hist <- expression.glass.exon.metadata %>% 
  dplyr::filter(padj.partially.paired.exon < 0.01) %>% 
  dplyr::filter(gene_chr == "chr6") %>% 
  dplyr::filter(log2FoldChange.partially.paired.exon > 0.75) %>% 
  dplyr::filter(gene_chr_center_loc > 25000000 & gene_chr_center_loc < 30000000) %>% 
  dplyr::filter(grepl("^H",gene_name)) %>% 
  dplyr::filter(grepl("^HLA-",gene_name) == F) %>% 
  dplyr::pull(gene_name)



sel.tum <- c("SOX4","SOX2","GFAP","OLIG1") # Tum
sel.tam1 <- c("CD74","NEAT1") # TAM.1
sel.tam2 <- c("RBFOX3", "TNNT2", "TMEM130", "GABRG2", "GABRB2") # TAM.2
sel.od <- c("TMEM144", "TMEM125", "MOG", "PLP1") # OD
sel.en <- c("ABCB1","CD34","FLT4","TIE1") # En
sel.pe <- c("RGS5","CD248","PDGFRB")# Pe
sel.tc <- c("CD2","CD3D","TRBC2","TRAC")# TC
sel.cycling <- c("FAM64A","AURKB","TOP2A","TPX2","CDC20") # Cycl

# weird
# "FAM64A" "CD3D"   "TRBC2"  "TRAC" 

sel <- expression.glass.exon.metadata %>%
  dplyr::filter(gene_name %in% c(
    sel.hist,
    sel.tum,
    sel.tam1,
    sel.tam2,
    sel.od,
    sel.en,
    sel.pe,
    sel.tc,
    sel.cycling)
  ) %>%
  dplyr::pull(gene_uid)


plt <- expression.glass.exon.vst %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::filter(gene_uid %in% sel) %>% 
  dplyr::left_join(expression.glass.exon.metadata %>% dplyr::select(gene_uid, gene_name), by=c('gene_uid'='gene_uid')) %>% 
  dplyr::mutate(gene_uid = NULL) %>% 
  tibble::column_to_rownames('gene_name')


metadata <- expression.glass.exon.metadata %>% 
  dplyr::filter(gene_uid %in% sel) %>% 
  dplyr::select(gene_name) %>% 
  dplyr::mutate(pericyte = gene_name %in% sel.pe) %>% 
  dplyr::mutate(tam.1 = gene_name %in% sel.tam1) %>% 
  dplyr::mutate(significant.hist.gene = gene_name %in% sel.hist) %>% 
  dplyr::mutate(cycling.marker = gene_name %in% sel.cycling) %>% 
  dplyr::mutate(tumor = gene_name %in% sel.tum) %>% 
  dplyr::mutate(tam.2 = gene_name %in% sel.tam2) %>% 
  dplyr::mutate(oligodendrocyte = gene_name %in% sel.od) %>% 
  dplyr::mutate(endothelial = gene_name %in% sel.en) %>% 
  tibble::column_to_rownames('gene_name')


recursiveCorPlot(plt, metadata, 7, 1.5)





# recursiveCorPlot [all DGE] ----


sign <- expression.glass.exon.metadata %>% 
  dplyr::filter(padj.partially.paired.exon < 0.01) %>% 
  dplyr::filter(abs(log2FoldChange.partially.paired.exon) >= 0.75)

cp <- expression.glass.exon.vst %>%
  dplyr::select(tmp.metadata$genomescan.sid) %>% 
  dplyr::filter(rownames(.) %in% sign$gene_uid) %>% 
  `rownames<-`(gsub("ENSG00000284906_ARHGAP11B","ENSG00000284906_ARHGAP11B.2",rownames(.),fixed=T)) %>% 
  `rownames<-`(gsub("^ENS.+_","",rownames(.)))


dim(cp)
# find n PCA - 5 features?!
pca <- prcomp(t(cp))
plot(pca)
dev.off()


pcs <- pca$rotation[,1:5] %>%
  as.data.frame() %>% 
  dplyr::mutate(main.PC = {names(.)[max.col(.)]}) %>% 
  dplyr::select(main.PC) %>% 
  tibble::rownames_to_column('gene_name') %>% 
  dplyr::mutate(PC1 = ifelse(main.PC == "PC1", NA , F )) %>% 
  dplyr::mutate(PC2 = ifelse(main.PC == "PC2", NA , F )) %>% 
  dplyr::mutate(PC3 = ifelse(main.PC == "PC3", NA , F )) %>% 
  dplyr::mutate(PC4 = ifelse(main.PC == "PC4", NA , F )) %>% 
  dplyr::mutate(PC5 = ifelse(main.PC == "PC5", NA , F )) %>% 
  dplyr::mutate(main.PC = NULL)



clusters <- data.frame(cluster = cutree(readRDS('cache/h.Rds'), k=5)) %>% 
  dplyr::mutate(cluster = paste0('c',cluster)) %>% 
  dplyr::mutate(val=T) %>% 
  tibble::rownames_to_column('gene_name') %>% 
  dplyr::mutate(marker=T) %>% 
  tidyr::pivot_wider(names_from = cluster,values_from = marker, values_fill=F)


h.up <- dendextend::prune(h,cutree(h,2) %>% purrr::keep(function(x) x == "1") %>%  names )
h.down <- dendextend::prune(h,cutree(h,2) %>% purrr::keep(function(x) x == "2") %>%  names )
#par(mfrow=c(2,2))
#plot(h)
#plot(h.down)
#plot(h.up)
clusters <- 
  rbind(
    data.frame(cluster = cutree(h.down, 2)) %>% 
      dplyr::mutate(cluster = paste0('down.',cluster))
    ,
    data.frame(cluster = cutree(h.up, 3)) %>% 
      dplyr::mutate(cluster = paste0('up.',cluster))
  ) %>%  
  dplyr::mutate(val=T) %>% 
  tibble::rownames_to_column('gene_name') %>% 
  dplyr::mutate(marker=T) %>% 
  tidyr::pivot_wider(names_from = cluster,values_from = marker, values_fill=F) %>% 
  dplyr::mutate(val = NULL)




cpm <- data.frame(gid=rownames(cp)) %>% 
  dplyr::mutate(HOX = grepl("^HOX", gid)) %>% 
  dplyr::mutate(COL = ifelse(grepl("^COL", gid),"red","green")) %>%
  #dplyr::left_join(pcs, by=c('gid'='gene_name')) %>% 
  dplyr::left_join(clusters, by=c('gid'='gene_name')) %>% 
  tibble::column_to_rownames('gid')


#h <- recursiveCorPlot(cp, cpm, 2 ,2, method="ward.D2",T)
recursiveCorPlot(cp, cpm, 2 ,9, method="ward.D2",FALSE)





ggsave("/tmp/glass-supervised.png",height=20 * 1.3,width=30 * 1.3)



"CD248" %in% rownames(cp)


# recursiveCorPlot [histones + cycling + cell type markers] ----

# subsel most prominent cycling markers


tmp.metadata <- metadata.glass.per.patient %>%
  dplyr::select(genomescan.sid.I, genomescan.sid.R) %>%
  tidyr::pivot_longer(cols=c('genomescan.sid.I', 'genomescan.sid.R')) %>%
  tidyr::drop_na(value) %>%
  dplyr::mutate(name=NULL) %>%
  dplyr::rename(genomescan.sid = value) %>%
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(metadata.glass.per.patient %>% dplyr::select(GLASS_ID, patient.correction.id), by=c('GLASS_ID'='GLASS_ID')) %>% 
  dplyr::arrange(Sample_Type, genomescan.sid)


sign <- expression.glass.exon.metadata %>% 
  dplyr::filter(padj.partially.paired.exon < 0.01) %>% 
  dplyr::filter(abs(log2FoldChange.partially.paired.exon) >= 0.75)

sign.hist <- expression.glass.exon.metadata %>% 
  dplyr::filter(padj.partially.paired.exon < 0.01) %>% 
  dplyr::filter(gene_chr == "chr6") %>% 
  dplyr::filter(log2FoldChange.partially.paired.exon > 0.75) %>% 
  dplyr::filter(gene_chr_center_loc > 25000000 & gene_chr_center_loc < 30000000) %>% 
  dplyr::filter(grepl("^H",gene_name)) %>% 
  dplyr::filter(grepl("^HLA-",gene_name) == F) %>% 
  dplyr::pull(gene_uid)

sel.cycling <- cycling.cell.markers %>% 
  tibble::rownames_to_column('gene_name') %>% 
  tibble::column_to_rownames('gene_uid') %>% 
  dplyr::mutate(gene_name = NULL) %>% 
  dplyr::mutate(n = rowSums(.)) %>% 
  dplyr::arrange(-n) %>% 
  dplyr::filter(n >= 2) %>% 
  dplyr::mutate(cycling.melanoma.tirosh = NULL, n = NULL) %>% 
  dplyr::mutate(n = rowSums(.)) %>% 
  dplyr::arrange(-n) %>% 
  dplyr::filter(n == 2) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::filter(gene_uid %in% sign$gene_uid)  # only those considered significant

prog.oligo <- expression.glass.exon.metadata %>% 
  dplyr::filter(gene_name %in% c("OLIG1","NEU4","GPR17","SLC1A1","ATCAY","SIRT2","APOD","MYT1","OLIG2","TMEFF2","OMG","ELMO1","RTKN","HIP1R","TNR","RPSA","MEGF11","EVI2A","OPCML","LHFPL3","RAB33A","GRIA4","SERINC5","NXPH1","BIN1","BMP4","EHD3","GNAI1","CSPG4","DSCAM","GALNT13","ZDHHC9","ABCG1","FKBP1A","LRRN1","ST8SIA3","DNM3","RAPGEF4","CNP","PDGFRA","PTGDS","CHGA","BCAS1","PLXNB3","NFASC","SLC44A1","GNG4","PHLDB1","CD82","PRKCZ")) %>% 
  dplyr::slice_head(n=25) %>% 
  dplyr::pull(gene_uid)

prog.astro <- expression.glass.exon.metadata %>% 
  dplyr::filter(gene_name %in% c("APOE","SPARCL1","VIM","ID4","TIMP3","EDNRB","MLC1","ID3","CLU","TNC","ZFP36L1","ARHGEF26","ATP1B2","AGT","RGMA","JUN","PFKFB3","EZR","SLC1A3","ALDOC","JUNB","ATP1A2","DTNA","ZFP36","SOX9","TRIL","NDRG2","NMB","GFAP","SLC1A2","RFX4","MALAT1","LRIG1","FOS","EGR1","STK17B","FOSB","ATF3","ABCA1","ADCYAP1R1","GLUL","IER2","ZFP36L2","ADHFE1","MSI2","CPE","KLF6","DOCK7","IRF2BP2","SPRY2")) %>% 
  dplyr::slice_head(n=25) %>% 
  dplyr::pull(gene_uid)

prog.stemn <- expression.glass.exon.metadata %>% 
  dplyr::filter(gene_name %in% c("SOX4","DCX","IGFBPL1","SOX11","TCF4","NREP","RND3","CCND2","MIAT","CAMK2N1","STMN4","STMN1","MYT1L","HN1","RNF122","PROX1","KLHDC8A","ELAVL4","NMNAT2","TUBB","ROBO1","NELL2","MLLT11","CELF4","POU3F2","H3F3B","ENC1","GNG2","ACOT7","AKT3","ARL4C","FNBP1L","VOPP1","TOX3","TUBB3","SCG2","TMSB15A","TFDP2","TMSB4X","CDC42","STMN2","KCTD13","RPH3A","KIF5C","NFIX","CALM1","TNPO2","BOC","KLHL13","PGAP1","RBFOX2","TMSB10","DYNLT1","TMSB15B","TCEAL7","PTS","BICD1","UCHL1","COMMD3","MCM7","AMZ2","PDRG1","DDAH2","KLC1","PCSK2","OAZ1","TIMM17A","YWHAG","CBX1","SMS","DGUOK","SNRPG","CDK6","GOLT1B","DUSP10","ATP5J","DYNLRB1","TCP1","GADD45G","SEC31A","CNOT7","DDX39A","SRGAP2","MAST2","PGK1","CELF3","ZFAS1","ENO2","SNRPB","DRG1")) %>% 
  dplyr::slice_head(n=25) %>% 
  dplyr::pull(gene_uid)


cell.type.pericyte <- expression.glass.exon.metadata %>% 
  dplyr::filter(gene_name %in% c('RGS5','CD248','PDGFRB','HEYL','CFH')) %>% 
  dplyr::pull(gene_uid)

cell.type.endothelial <- expression.glass.exon.metadata %>% 
  dplyr::filter(gene_name %in% c('RGS5','NOSTRIN','TIE1','FLT1','CD31')) %>% 
  dplyr::pull(gene_uid)

cell.type.oligodendrocyte <- expression.glass.exon.metadata %>% 
  dplyr::filter(gene_name %in% c('MOG','TMEM144','OPALIN','PLP1')) %>% 
  dplyr::pull(gene_uid)

cell.type.neuron <- expression.glass.exon.metadata %>% 
  dplyr::filter(gene_name %in% c('GABRG2','GABRA1','RBFOX3','GABRB2')) %>% 
  dplyr::pull(gene_uid)

cell.type.astrocyte <- expression.glass.exon.metadata %>% 
  dplyr::filter(gene_name %in% c('GFAP','BMPR1B','CACHD1','GPR37L1')) %>% 
  dplyr::pull(gene_uid)

cell.type.tcell <- expression.glass.exon.metadata %>% 
  dplyr::filter(gene_name %in% c('CD3D','CD2','CD3E','CD8A')) %>% 
  dplyr::pull(gene_uid)

cell.type.tam <- expression.glass.exon.metadata %>% 
  dplyr::filter(gene_name %in% c('C1QA','THEMIS2','CD14','CD163')) %>% 
  dplyr::pull(gene_uid)




sel <- c(
  sign.hist, sel.cycling$gene_uid, prog.astro, prog.oligo, prog.stemn,
  cell.type.pericyte, cell.type.endothelial, cell.type.oligodendrocyte, cell.type.neuron, cell.type.astrocyte, cell.type.tcell, cell.type.tam 
  ) %>% 
  unique()



plt <- expression.glass.exon.vst %>%
  dplyr::select(tmp.metadata$genomescan.sid) %>% 
  dplyr::filter(rownames(.) %in% sel) %>% 
  `rownames<-`(gsub("ENSG00000269226_TMSB15B","ENSG00000269226_TMSB15B.2",rownames(.),fixed=T)) 

#%>% 
  #`rownames<-`(gsub("^ENS.+_","",rownames(.)))


plt.metadata <- data.frame(gene_name = rownames(plt)) %>% 
  dplyr::mutate(`chr6 histon gene locus & DE` = gene_name %in% sign.hist ) %>% 
  dplyr::left_join(cycling.cell.markers , by=c('gene_name'='gene_uid')) %>% 
  dplyr::mutate(prog.oligo = gene_name %in% prog.oligo) %>% 
  dplyr::mutate(prog.astro = gene_name %in% prog.astro) %>% 
  dplyr::mutate(prog.stemn = gene_name %in% prog.stemn) %>% 
  dplyr::mutate(cell.type.pericyte = gene_name %in% cell.type.pericyte) %>% 
  dplyr::mutate(cell.type.endothelial = gene_name %in% cell.type.endothelial) %>% 
  dplyr::mutate(cell.type.oligodendrocyte = gene_name %in% cell.type.oligodendrocyte) %>% 
  dplyr::mutate(cell.type.neuron = gene_name %in% cell.type.neuron) %>% 
  dplyr::mutate(cell.type.astrocyte = gene_name %in% cell.type.astrocyte) %>% 
  dplyr::mutate(cell.type.tcell = gene_name %in% cell.type.tcell) %>% 
  dplyr::mutate(cell.type.tam = gene_name %in% cell.type.tam) %>% 
  tibble::column_to_rownames('gene_name') %>% 
  dplyr::mutate_all(function(arg) { return(ifelse(is.na(arg),FALSE,arg)) }) # verandere overige NA's van niet matchen naar FALSE




recursiveCorPlot(plt, plt.metadata, 5 ,3)




ggsave("/tmp/glass-supervised.png",height=20 * 1.3,width=30 * 1.2)


# recursiveCorPlot [all DGE + hist & cycling labels] ----


tmp.metadata <- metadata.glass.per.patient %>%
  dplyr::select(genomescan.sid.I, genomescan.sid.R) %>%
  tidyr::pivot_longer(cols=c('genomescan.sid.I', 'genomescan.sid.R')) %>%
  tidyr::drop_na(value) %>%
  dplyr::mutate(name=NULL) %>%
  dplyr::rename(genomescan.sid = value) %>%
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid'='genomescan.sid')) %>%
  dplyr::left_join(metadata.glass.per.patient %>% dplyr::select(GLASS_ID, patient.correction.id), by=c('GLASS_ID'='GLASS_ID')) %>% 
  dplyr::arrange(Sample_Type, genomescan.sid)


sign <- expression.glass.exon.metadata %>% 
  dplyr::filter(padj.partially.paired.exon < 0.01) %>% 
  dplyr::filter(abs(log2FoldChange.partially.paired.exon) >= 0.75)

sign.hist <- expression.glass.exon.metadata %>% 
  dplyr::filter(padj.partially.paired.exon < 0.01) %>% 
  dplyr::filter(gene_chr == "chr6") %>% 
  dplyr::filter(log2FoldChange.partially.paired.exon > 0.75) %>% 
  dplyr::filter(gene_chr_center_loc > 25000000 & gene_chr_center_loc < 30000000) %>% 
  dplyr::filter(grepl("^H",gene_name)) %>% 
  dplyr::filter(grepl("^HLA-",gene_name) == F) %>% 
  dplyr::pull(gene_name)




plt <- expression.glass.exon.vst %>%
  dplyr::select(tmp.metadata$genomescan.sid) %>% 
  dplyr::filter(rownames(.) %in% sign$gene_uid) %>% 
  `rownames<-`(gsub("ENSG00000284906_ARHGAP11B","ENSG00000284906_ARHGAP11B.2",rownames(.),fixed=T)) %>% 
  `rownames<-`(gsub("^ENS.+_","",rownames(.)))


plt.metadata <- data.frame(gene_name = rownames(plt)) %>% 
  dplyr::mutate(`chr6 histon gene locus` = gene_name %in% sign.hist ) %>% 
  dplyr::left_join(cycling.cell.markers %>% tibble::rownames_to_column('gene_name'),by=c('gene_name'='gene_name')) %>% 
  tibble::column_to_rownames('gene_name') %>% 
  dplyr::mutate_all(function(arg) { return(ifelse(is.na(arg),FALSE,arg)) }) # verandere overige NA's van niet matchen naar FALSE


stopifnot(sign.hist %in% rownames(plt.metadata)) # should all be in





recursiveCorPlot(plt, plt.metadata, 7 ,13)




ggsave("/tmp/glass-supervised.png",height=20 * 1.3,width=30 * 1.3)





# chromosome plot [all DGE] ----

plt <- expression.glass.exon.metadata %>% 
  dplyr::left_join(chrs_hg38_s, by=c('gene_chr'='chr')) %>% 
  dplyr::mutate(x = gene_chr_center_loc + pos) %>% 
  dplyr::mutate(gene_chr = factor(gene_chr, levels=gtools::mixedsort(unique(as.character(gene_chr))) )) %>% 
  dplyr::mutate(significant = padj.partially.paired.exon < 0.01 & abs(log2FoldChange.partially.paired.exon) > 0.75)



ggplot(plt, aes(x=gene_chr_center_loc / 1000000,y=stat.partially.paired.exon,col=gene_chr)) + 
  facet_grid(cols = vars(gene_chr), scales = "free", space="free") +
  geom_point(pch=19,cex=0.2) +
  geom_point(data = subset(plt, significant==T), pch=21,cex=0.8,col='black',fill=NA) +
  geom_smooth(se=F,col="black", lwd=0.7) +
  youri_gg_theme + labs(x=NULL)



# chromosome plot [all DGE + hist & cycling] ----

plt <- expression.glass.exon.metadata %>% 
  dplyr::left_join(chrs_hg38_s, by=c('gene_chr'='chr')) %>% 
  dplyr::mutate(x = gene_chr_center_loc + pos) %>% 
  dplyr::mutate(gene_chr = factor(gene_chr, levels=gtools::mixedsort(unique(as.character(gene_chr))) )) %>% 
  dplyr::mutate(significant = padj.partially.paired.exon < 0.01 & abs(log2FoldChange.partially.paired.exon) > 0.75)



ggplot(plt, aes(x=gene_chr_center_loc / 1000000,y=stat.partially.paired.exon,col=gene_chr)) + 
  facet_grid(cols = vars(gene_chr), scales = "free", space="free") +
  geom_point(pch=19,cex=0.2) +
  geom_point(data = subset(plt, significant==T), pch=21,cex=0.8,col='black',fill=NA) +
  geom_smooth(se=F,col="black", lwd=0.7) +
  youri_gg_theme + labs(x=NULL)



# rank x basemean? ----

clusters <- 
  rbind(
    data.frame(cluster = cutree(h.down, 2)) %>% 
      dplyr::mutate(cluster = paste0('down.',cluster))
    ,
    data.frame(cluster = cutree(h.up, 3)) %>% 
      dplyr::mutate(cluster = paste0('up.',cluster))
  ) %>% 
  tibble::rownames_to_column('gene_name')



plt <- expression.glass.exon.metadata %>% 
  dplyr::mutate(rank.basemean = rank(baseMean.partially.paired.exon)) %>% 
  dplyr::mutate(rank.lfc = rank(log2FoldChange.partially.paired.exon)) %>% 
  dplyr::mutate(rank.lfcse = rank(lfcSE.partially.paired.exon)) %>% 
  dplyr::left_join(clusters, by=c('gene_name'='gene_name'))




## basemean = ok? ----
ggplot(plt, aes(x=rank.basemean, y=log(1+baseMean.partially.paired.exon),col=cluster)) +
  geom_point(data=subset(plt, is.na(cluster))) +
  geom_point(data=subset(plt, !is.na(cluster) & cluster == "down.1"),cex=0.375) +
  theme_bw()

ggplot(plt, aes(x=rank.basemean, y=log(1+baseMean.partially.paired.exon),col=cluster)) +
  geom_point(data=subset(plt, is.na(cluster))) +
  geom_point(data=subset(plt, !is.na(cluster) & cluster == "down.2"),cex=0.375) +
  theme_bw()


## lfc ----

ggplot(plt, aes(x=rank.lfc, y=log2FoldChange.partially.paired.exon,col=cluster)) +
  geom_point(data=subset(plt, is.na(cluster))) +
  geom_point(data=subset(plt, !is.na(cluster) & cluster == "down.1"),cex=0.375) +
  theme_bw()

ggplot(plt, aes(x=rank.lfc, y=log2FoldChange.partially.paired.exon,col=cluster)) +
  geom_point(data=subset(plt, is.na(cluster))) +
  geom_point(data=subset(plt, !is.na(cluster) & cluster == "down.2"),cex=0.375) +
  theme_bw()


# # lfcSE ----

ggplot(plt, aes(x=rank.lfcse, y=lfcSE.partially.paired.exon,col=cluster)) +
  geom_point(data=subset(plt, is.na(cluster))) +
  geom_point(data=subset(plt, !is.na(cluster) & cluster == "down.1"),cex=0.375) +
  theme_bw()

ggplot(plt, aes(x=rank.lfcse, y=lfcSE.partially.paired.exon,col=cluster)) +
  geom_point(data=subset(plt, is.na(cluster))) +
  geom_point(data=subset(plt, !is.na(cluster) & cluster == "down.2"),cex=0.375) +
  theme_bw()



# OD + AC signature ----







