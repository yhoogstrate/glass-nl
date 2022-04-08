#!/usr/bin/env R


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





## recursiveCorPlot [histone genes + cell types + cycling] ----

#@todo script deze met LFC, chr6 en ^H
sel.hist <- c("H4C1", "H3C2", "H2AC4", "H3C3", "H1-6", "H3C7", "H2BC9", "H2BC11", "H2AC11", "H2BC12", "H2AC12", "H2BC13", "H2AC13", "H3C10", "H2AC14", "H2BC14", "H2AC15", "H2AC16", "H1-5", "H3C11", "H3C12", "H2BC17")

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


sign <- res.paired.a %>% 
  dplyr::filter(padj < 0.01) %>% 
  dplyr::filter(abs(log2FoldChange) >= 0.75)

cp <- expression.glass.vst %>%
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




cpm <- data.frame(gid=rownames(cp)) %>% 
  dplyr::mutate(HOX = grepl("^HOX", gid)) %>% 
  dplyr::mutate(COL = ifelse(grepl("^COL", gid),"red","green")) %>%
  dplyr::left_join(pcs, by=c('gid'='gene_name')) %>% 
  tibble::column_to_rownames('gid')


h <- recursiveCorPlot::recursiveCorPlot(cp, cpm, 2 ,2)

ggsave("/tmp/glass-supervised.png",height=20 * 1.3,width=30 * 1.3)



"CD248" %in% rownames(cp)


