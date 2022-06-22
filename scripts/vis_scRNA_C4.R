#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(Matrix)
library(DropletUtils)
library(Seurat)
#library(infercnv)
library(AnnotationHub)
library(ensembldb)


# clusters data ----


source('scripts/load_hclust.R')
source('scripts/R/lgg.transcriptional.programs.Venteicher.R')


# :: Glimmunology :: ----

# A :: snRNA Sample_Y GBM ----
-



rm(sid, object_1)
gc()


sid <- 'van_Hijfte_Sample_Y'
object_1 <- Read10X(data.dir = "data/Glimmunology_GBM_1/Glioma_Y_and_O/Levi2_Glioma_Y/outs/raw_feature_bc_matrix")
object_1 <- CreateSeuratObject(counts = object_1,
                               min.cells = 3,
                               min.features = 200,
                               project = "glioma_glim")
mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1400,col="red") +
  geom_hline(yintercept = 4500,col="red")

ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 2200,col="red") +
  geom_hline(yintercept = 14000,col="red")
#scale_y_log10()

object_1 <- subset(x = object_1, subset = 
                     nFeature_RNA > 1400 &
                     nFeature_RNA < 4500 & 
                     nCount_RNA > 2200 &
                     nCount_RNA < 14000 &
                     percent.mito < 0.025)

object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 


top10 <- head(VariableFeatures(object_1), 10)


print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))



# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))     
plot1
plot2

## scaling of data ----
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

### PCA plot ---

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

## cluster the cells ----

object_1 <- FindNeighbors(object_1, dims = 1:40)
object_1 <- FindClusters(object_1, resolution = 0.8, algorithm=1)
head(Idents(object_1), 20)


## UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:40)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")

# levels(object_1$seurat_clusters) <- gsub("^21$","\\1.Neuron",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(5|10|13|14)$","Immune + T-cells.\\1",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(4|6|12|19)$","Oligodendrocytes.\\1",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^16$","Endothelial",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^15$","Pericytes",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|9|8|11)$","Tumor.\\1",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^7$","Tumor/Dividing",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^17$","Tumor outlier",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^20$","Astrocyte",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^18$","Tumor/Apoptosis?",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^22$","Immune cell / OD hybrid?",levels(object_1$seurat_clusters))

object_1 <- FindClusters(object_1, resolution = 0.8, algorithm=1)
object_1$class <- as.character(object_1$seurat_clusters)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(21), "NE", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(5,10,13,14), "TAM/MG", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(4,6,12,19), "OD", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(16), "EN", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(15), "PE", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(0,1,2,3,9,8,11), "T", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(7), "T", object_1$class) # dividing
object_1$class <- ifelse(object_1$seurat_clusters %in% c(17), "T ?", object_1$class) #  Outlier?
object_1$class <- ifelse(object_1$seurat_clusters %in% c(20), "AC", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(18), "T ?", object_1$class) # Apoptotic?
object_1$class <- ifelse(object_1$seurat_clusters %in% c(22), "TAM/MG|OD", object_1$class)
object_1$class <- ifelse(object_1@reductions$umap@cell.embeddings[,1] >= 10 &
                           object_1@reductions$umap@cell.embeddings[,1] <= 11 &
                           object_1@reductions$umap@cell.embeddings[,2] >= 1.5 &
                           object_1@reductions$umap@cell.embeddings[,2] <= 3,
                         "TC", object_1$class)
object_1$seurat_clusters <- as.factor(paste0(as.character(object_1$seurat_clusters),". ",object_1$class))

levels(object_1$seurat_clusters)


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T","1. T","2. T","3. T","7. T","8. T","9. T","11. T",
  "17. T ?","18. T ?",
  "20. AC",
  "21. NE",
  "4. OD","6. OD","12. OD","19. OD",
  "16. EN",
  "15. PE",
  "22. TAM/MG|OD" ,
  "5. TAM/MG","10. TAM/MG","13. TAM/MG","14. TAM/MG",
  "14. TC"
))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")  +
  labs(subtitle=sid) +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3)))



#ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_UMAP.pdf"),width=10,height=8)
#ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_UMAP.png"),width=12,height=10)


od.markers <- FindMarkers(object_1, ident.1 = c(4,6,12,19,22))
View(od.markers)



# 
# 
# tmp.17 <- FindMarkers(object_1, ident.1 = 17)
# head(tmp.17,20)
# 
# tmp.22 <- FindMarkers(object_1, ident.1 = 22)
# head(tmp.22,20)

tmp.15 <- FindMarkers(object_1, ident.1 = 15) # PE
View(tmp.15)




## 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR", max.cutoff = 4) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = "PDGFRB") # Tumor/MES


FeaturePlot(object = object_1, features = "EREG") # EREG
FeaturePlot(object = object_1, features = "EGF") # BTC
FeaturePlot(object = object_1, features = "BTC") # BTC

DotPlot(object = object_1, features =list( 'ligands'=c('EGF','EREG','AREG','BTC','EPGN','HBEGF') ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("EGFR and ligands in: ",sid))


FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M


# AC
FeaturePlot(object = object_1, features = c("CST3","S100B","SLC1A3","HEPN1","HOPX","MT3","SPARCL1","MLC1"))
FeaturePlot(object = object_1, features = c("GFAP","FABP7","BCAN","PON2","METTL7B","SPARC","GATM","RAMP1")) # Tumor
FeaturePlot(object = object_1, features = c("PMP2","AQP4","DBI","EDNRB","PTPRZ1","CLU","PMP22","ATP1A2")) # Tumor
FeaturePlot(object = object_1, features = c("S100A16","HEY1","PCDHGC3","TTYH1","NDRG2","PRCP","ATP1B2","AGT","PLTP","GPM6B"))
FeaturePlot(object = object_1, features = c("F3","RAB31","PPAP2B","ANXA5","TSPAN7")) # Tumor

# MES1
FeaturePlot(object = object_1, features = c("CHI3L1","ANXA2","ANXA1","CD44","VIM","MT2A","C1S","NAMPT","EFEMP1","C1R","SOD2")) # Tumor

# MES2
FeaturePlot(object = object_1, features = c("HILPDA","ADM","DDIT3","NDRG1","HERPUD1","DNAJB9","TRIB3","ENO2","AKAP12","SQSTM1","MT1X","ATF3","NAMPT")) # Tumor

# OPC
FeaturePlot(object = object_1, features = c("BCAN","PLP1","GPR17","FIBIN","LHFPL3","OLIG1","PSAT1","SCRG1","OMG","APOD","SIRT2","TNR","THY1","PHYHIPL","SOX2-OT","NKAIN4")) # Tumor

# NPC1
FeaturePlot(object = object_1, features = c("DLL3","DLL1","SOX4","TUBB3","HES6","TAGLN3","NEU4","MARCKSL1","CD24","STMN1","TCF12","BEX1","OLIG1","MAP2","FXYD6","PTPRS","MLLT11","NPPA","BCAN","MEST")) # Tumor

# NPC2
FeaturePlot(object = object_1, features = c("STMN2","CD24","RND3","HMP19","TUBB3","MIAT","DCX","NSG1","ELAVL4","MLLT11","DLX6-AS1","SOX11","NREP","FNBP1L","TAGLN3","STMN4")) # Tumor

# NPC1, NPC2, Neuron
FeaturePlot(object = object_1, features = c("BCAN", "NREP", "RBFOX3")) # Tumor
# BCAN


# GFAP en ANXA1 astro markers?

FeaturePlot(object = object_1, features = "ANXA1") # Tumor / Neuronal? weinig in GFAP-neg tumor cellen?
FeaturePlot(object = object_1, features = "ANXA2") # Tumor

FeaturePlot(object = object_1, features = "SOCS2") # Tumor?
FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # niet alleen in tumor cellen (!)

FeaturePlot(object = object_1, features = "PTPZ1")
FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")

FeaturePlot(object = object_1, features = "OLIG1")

FeaturePlot(object = object_1, features = "NODAL")


## 2. Astrocyte (+) ----

FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "AQP4")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "GJA1")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")



## 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))



## 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


## 3C. Hematopoietic stem cells? ----


FeaturePlot(object = object_1, features = "HBG1") # Tumor
FeaturePlot(object = object_1, features = "HBG2") # Tumor


## 3D. ? Mono/Leukocyte ?? ----


# These are cluster-13 DE genes, of which some at genecards seem related to leukocytes?
FeaturePlot(object = object_1, features = c("LAMP3","IRF4","NCCRP1","CRIP1","SYNPO2","CCR7","EHF","CCL22","VTN","LSP1","CDX2"))



## 4. Neurons (+) ----


tmp <- list('C1'=neuron.genes[neuron.genes %in% NPC2 == F] ,
            'NPC1'=NPC1[NPC1 %in% NPC2 == F] ,
            'NPC1+2' = intersect(NPC1, NPC2),
            'NPC2'=NPC2[NPC2 %in% c(NPC1, neuron.genes) == F],
            'NPC2 + C1' = intersect(neuron.genes, NPC2))

DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C1/NPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C1.pdf"),width=7.5*3, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C1.png"),width=7.5*3, height=4,scale=1.2)




DotPlot(object = object_1, features = c("EGFR", "GFAP","MOG", "PLP1", "TMEM144", 
                                        "RBFOX1", "RBFOX2", "RBFOX3", "CD2",
                                        "CD3D", "P2RY12", "CD163", "ABCB1", "RGS5"
))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))


FeaturePlot(object = object_1, features = c("SOX4", "RBFOX3"))


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")


#FeaturePlot(object = object_1, features = "DCN") # DCN
#FeaturePlot(object = object_1, features = "COL1A2") # DCN
FeaturePlot(object = object_1, features = "ANPEP") # DCN


## 5. Oligodendrocytes (+) ----


# + "PEAR1", "HEYL" , "CFH"
DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C2/OPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C2.pdf"),width=7.5*2, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C2.png"),width=7.5*2, height=4,scale=1.2)



FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


FeaturePlot(object = object_1, features = "ST18") # OD?
FeaturePlot(object = object_1, features = "MBP") # OD?
FeaturePlot(object = object_1, features = "CTNNA3") # OD?
FeaturePlot(object = object_1, features = "SLC24A2") # OD?
FeaturePlot(object = object_1, features = "KIRREL3") # OD?
FeaturePlot(object = object_1, features = "NKAIN2") # OD?
FeaturePlot(object = object_1, features = "MAP7") # OD?
FeaturePlot(object = object_1, features = "RNF220") # OD?
FeaturePlot(object = object_1, features = "PEX5L") # OD?
FeaturePlot(object = object_1, features = "TMEM144") # OD?
FeaturePlot(object = object_1, features = "EDIL3") # OD?
FeaturePlot(object = object_1, features = "DOCK5") # OD?
FeaturePlot(object = object_1, features = "MOBP") # OD?
FeaturePlot(object = object_1, features = "UNC5C") # OD?
FeaturePlot(object = object_1, features = "CLDN11") # OD?
FeaturePlot(object = object_1, features = "SPOCK3") # OD?
FeaturePlot(object = object_1, features = "CNTNAP4") # OD?
FeaturePlot(object = object_1, features = "MAN2A1") # OD?
FeaturePlot(object = object_1, features = "PCSK6") # OD?
FeaturePlot(object = object_1, features = "TTLL7") # OD?

FeaturePlot(object = object_1, features = "OLIG2") # OD?


## 6A. Endothelial (+) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?



## 6B. Pericytes (+) ----

FeaturePlot(object = object_1, features = c("RGS5","PDGFRB","CD248","PEAR1", "HEYL" , "CFH"))

FeaturePlot(object = object_1, features = c("RGS5"))
FeaturePlot(object = object_1, features = c("PDGFRB"))
FeaturePlot(object = object_1, features = c("CD248"))
FeaturePlot(object = object_1, features = c("PEAR1"))
FeaturePlot(object = object_1, features = c("HEYL"))
FeaturePlot(object = object_1, features = c("CFH"))




## 7. Cycling cells (?) ----


FeaturePlot(object = object_1, features = "TOP2A" )
FeaturePlot(object = object_1, features = "AURKA" )
FeaturePlot(object = object_1, features = "AURKB" )
FeaturePlot(object = object_1, features = "BUB1" )
FeaturePlot(object = object_1, features = "BUB1B" )
FeaturePlot(object = object_1, features = "CDC20" )
FeaturePlot(object = object_1, features = "CENPF" )
FeaturePlot(object = object_1, features = "FAM64A" )
FeaturePlot(object = object_1, features = "FOXM1" )
FeaturePlot(object = object_1, features = "TACC3" )
FeaturePlot(object = object_1, features = "TMPO" )
FeaturePlot(object = object_1, features = "TPX2" )
FeaturePlot(object = object_1, features = "TUBA1C" )


# FeaturePlot(object = object_1, features = "RRM2" )
FeaturePlot(object = object_1, features = "PCNA" )
FeaturePlot(object = object_1, features = "KIAA0101" )
# FeaturePlot(object = object_1, features = "HIST1H4C" )
# FeaturePlot(object = object_1, features = "MLF1IP" )
# FeaturePlot(object = object_1, features = "GMNN" )
FeaturePlot(object = object_1, features = "RNASEH2A" )
# FeaturePlot(object = object_1, features = "MELK" )
# FeaturePlot(object = object_1, features = "CENPK" )
# FeaturePlot(object = object_1, features = "TK1" )
FeaturePlot(object = object_1, features = "TMEM106C" )
# FeaturePlot(object = object_1, features = "CDCA5" )
FeaturePlot(object = object_1, features = "CKS1B" )
FeaturePlot(object = object_1, features = "CDC45" )
FeaturePlot(object = object_1, features = "MCM3" )
FeaturePlot(object = object_1, features = "CENPM" )
FeaturePlot(object = object_1, features = "AURKB" )
FeaturePlot(object = object_1, features = "PKMYT1" )
FeaturePlot(object = object_1, features = "MCM4" )
FeaturePlot(object = object_1, features = "ASF1B" )
FeaturePlot(object = object_1, features = "GINS2" )
FeaturePlot(object = object_1, features = "MCM2" )
FeaturePlot(object = object_1, features = "FEN1" )
FeaturePlot(object = object_1, features = "RRM1" )
FeaturePlot(object = object_1, features = "DUT" )
FeaturePlot(object = object_1, features = "RAD51AP1" )
FeaturePlot(object = object_1, features = "MCM7" )
FeaturePlot(object = object_1, features = "CCNE2" )
FeaturePlot(object = object_1, features = "ZWINT" )


FeaturePlot(object = object_1, features = "CCNB1" )
FeaturePlot(object = object_1, features = "CDC20" )
FeaturePlot(object = object_1, features = "CCNB2" )
FeaturePlot(object = object_1, features = "PLK1" )
FeaturePlot(object = object_1, features = "CCNA2" )
FeaturePlot(object = object_1, features = "CKAP2" )
FeaturePlot(object = object_1, features = "KNSTRN" )
FeaturePlot(object = object_1, features = "RACGAP1" )
FeaturePlot(object = object_1, features = "CDCA3" )
FeaturePlot(object = object_1, features = "TROAP" )
FeaturePlot(object = object_1, features = "KIF2C" )
FeaturePlot(object = object_1, features = "KPNA2" )
FeaturePlot(object = object_1, features = "KIF20A" )
FeaturePlot(object = object_1, features = "ECT2" )
FeaturePlot(object = object_1, features = "CDCA8" )
FeaturePlot(object = object_1, features = "TTK" )
FeaturePlot(object = object_1, features = "NCAPD2" )
FeaturePlot(object = object_1, features = "ARL6IP1" )
FeaturePlot(object = object_1, features = "KIF4A" )
FeaturePlot(object = object_1, features = "CKAP2L" )
FeaturePlot(object = object_1, features = "MZT1" )
FeaturePlot(object = object_1, features = "KIFC1" )
FeaturePlot(object = object_1, features = "SPAG5" )
FeaturePlot(object = object_1, features = "ANP32E" )
FeaturePlot(object = object_1, features = "KIF11" )
FeaturePlot(object = object_1, features = "PSRC1" )
FeaturePlot(object = object_1, features = "TUBB4B" )
FeaturePlot(object = object_1, features = "SMC4" )
FeaturePlot(object = object_1, features = "MXD3" )
FeaturePlot(object = object_1, features = "CDC25B" )
FeaturePlot(object = object_1, features = "OIP5" )
FeaturePlot(object = object_1, features = "REEP4" )
FeaturePlot(object = object_1, features = "GPSM2" )
FeaturePlot(object = object_1, features = "HMGB3" )
FeaturePlot(object = object_1, features = "ARHGAP11A" )
FeaturePlot(object = object_1, features = "RANGAP1" )
FeaturePlot(object = object_1, features = "H2AFZ" )

## 8. longitudinal sig. chr6 ----


sig <- c("H4C1", "H3C2", "H2AC4", "H3C3",            "HIST1H4A","HIST1H3B", "HIST1H2AB", "HIST1H3C",
         "H1-6", "H3C7", "H2BC9", "H2BC11",          "HIST1H1T","HIST1H3F","HIST1H2BH","HIST1H2BH",
         "H2AC11", "H2BC12", "H2AC12", "H2BC13",     "HIST1H2AG","HIST1H2BK","HIST1H2AH","HIST1H2BL",
         "H2AC13", "H3C10", "H2AC14", "H2BC14",      "HIST1H2AI","HIST1H3H","HIST1H2AJ","H2BC14",
         "H2AC15", "H2AC16", "H1-5", "H3C11",        "HIST1H2AK","HIST1H2AL","HIST1H1B","HIST1H3I",
         "H3C12", "H2BC17",                          "HIST1H3J","HIST1H2BO"
)
sig <- unique(sig)

sig %in% all.genes
sig <- sig[ sig %in% all.genes]
sig


DotPlot(object_1, features=sig, group.by = "seurat_clusters")



FeaturePlot(object_1, features=sig)

## LGG: up.1 ----

# Cell cycle + histones

set <- dge.partially.paired.clusters %>%
  dplyr::filter(up.1) %>%
  dplyr::pull(gene_name)


set <- c(set, c("H4C1", "H3C2", "H2AC4", "H3C3",            "HIST1H4A","HIST1H3B", "HIST1H2AB", "HIST1H3C",
                       "H1-6", "H3C7", "H2BC9", "H2BC11",          "HIST1H1T","HIST1H3F","HIST1H2BH","HIST1H2BH",
                       "H2AC11", "H2BC12", "H2AC12", "H2BC13",     "HIST1H2AG","HIST1H2BK","HIST1H2AH","HIST1H2BL",
                       "H2AC13", "H3C10", "H2AC14", "H2BC14",      "HIST1H2AI","HIST1H3H","HIST1H2AJ","H2BC14",
                       "H2AC15", "H2AC16", "H1-5", "H3C11",        "HIST1H2AK","HIST1H2AL","HIST1H1B","HIST1H3I",
                       "H3C12", "H2BC17",                          "HIST1H3J","HIST1H2BO"))
set <- unique(set)

DotPlot(object = object_1, features = set, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [LGG.up1/ccycle] in: ",sid))



## LGG: up.2 ----

# COL1A2 + CD248
set <- dge.partially.paired.clusters %>% dplyr::filter(up.2) %>%  dplyr::pull(gene_name)

DotPlot(object = object_1, features = set, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [LGG.up2/col/PE|EN] in: ",sid))



## LGG: up.3 ----

# ??
set <- dge.partially.paired.clusters %>% dplyr::filter(up.3) %>%  dplyr::pull(gene_name)

DotPlot(object = object_1, features = list('UP3'=set, 'NE'=c("RBFOX3")), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [LGG.up3/???] in: ",sid))



## LGG: down.1a ----


set1 <- dge.partially.paired.clusters %>% dplyr::filter(down.1) %>%  dplyr::pull(gene_name)
set2 <- dge.partially.paired.clusters %>% dplyr::filter(down.2) %>%  dplyr::pull(gene_name)
set <- list('down.1'=set1, 'down.2'=set2)

DotPlot(object = object_1, features = set, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [Down] in: ",sid))




# B. LGG_HG sample 101? => Tum ----
# C. LGG_HG sample 103? => Tum ----

# D. LGG combined new dataset ----

## 1/101 ----
# Loading and pre-processing dataset 1
sid <- 'lgg combined'

object_1 <- Read10X(data.dir = "data/Glimmunology_LGG_101P/outs/filtered_feature_bc_matrix/")
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")
mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5)
object_1 <- subset(x = object_1, subset = nFeature_RNA > 700 & nFeature_RNA <7000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 

## 2/102 ----
# Loading and pre-processing dataset 2
object_2 <- Read10X(data.dir = "data/Glimmunology_LGG_102P/outs/filtered_feature_bc_matrix/")
object_2 <- CreateSeuratObject(counts = object_2, min.cells = 3, min.features = 200, project = "glioma")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_2), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_2, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_2, slot = "counts"))
object_2[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5)
object_2 <- subset(x = object_2, subset = nFeature_RNA > 700 & nFeature_RNA <7000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
object_2 <- NormalizeData(object = object_2, normalization.method = "LogNormalize", scale.factor = 1e4)
object_2 <- FindVariableFeatures(object = object_2, selection.method = "vst", nfeatures = 2000)
object_2[["state"]] <- "P2"

## 3 ----
# Loading and pre-processing dataset 3
object_3 <- Read10X(data.dir = "data/Glimmunology_LGG_0/outs/filtered_feature_bc_matrix/")
object_3 <- CreateSeuratObject(counts = object_3, min.cells = 3, min.features = 200, project = "glioma")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_3), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_3, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_3, slot = "counts"))
object_3[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5)
object_3 <- subset(x = object_3, subset = nFeature_RNA > 700 & nFeature_RNA <7000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
object_3 <- NormalizeData(object = object_3, normalization.method = "LogNormalize", scale.factor = 1e4)
object_3 <- FindVariableFeatures(object = object_3, selection.method = "vst", nfeatures = 2000)
object_3[["state"]] <- "P3"

## 4/103 ----
# Loading and pre-processing dataset 4
object_4 <- Read10X(data.dir = "data/Glimmunology_LGG_103P/filtered_feature_bc_matrix/")
object_4 <- CreateSeuratObject(counts = object_4, min.cells = 3, min.features = 200, project = "glioma")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_4), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_4, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_4, slot = "counts"))
object_4[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_4, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5)
object_4 <- subset(x = object_4, subset = nFeature_RNA > 500 & nFeature_RNA <4000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
object_4 <- NormalizeData(object = object_4, normalization.method = "LogNormalize", scale.factor = 1e4)
object_4 <- FindVariableFeatures(object = object_4, selection.method = "vst", nfeatures = 2000)
object_4[["state"]] <- "P4"

## 5/104 ----
# Loading and pre-processing dataset 5
object_5 <- Read10X(data.dir = "data/Glimmunology_LGG_104P/filtered_feature_bc_matrix/")
object_5 <- CreateSeuratObject(counts = object_5, min.cells = 3, min.features = 200, project = "glioma")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_5), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_5, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_5, slot = "counts"))
object_5[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_5, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5)
object_5 <- subset(x = object_5, subset = nFeature_RNA > 300 & nFeature_RNA <1200 & nCount_RNA >200 & nCount_RNA <2000 & percent.mito <0.1)
object_5 <- NormalizeData(object = object_5, normalization.method = "LogNormalize", scale.factor = 1e4)
object_5 <- FindVariableFeatures(object = object_5, selection.method = "vst", nfeatures = 2000)
object_5[["state"]] <- "P5"


## integrate ----

# Identification of integration anchors
reference.list    <- c(object_1, object_2, object_3, object_4, object_5) # Makes a list of the objects you want to merge
glioma.anchors    <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30) # Finds the common sources of variation
glioma.integrated <- IntegrateData(anchorset = glioma.anchors, dims = 1:30) # Integrates the two datasets

rm(reference.list, object_1, object_2, object_3, object_4, object_5)

# Post-processing merged data
DefaultAssay(object = glioma.integrated) <- "integrated"
glioma.integrated <- ScaleData(object = glioma.integrated, verbose = F)
glioma.integrated <- RunPCA(object = glioma.integrated, verbose = F)
ElbowPlot(object = glioma.integrated, ndims = 40)
glioma.integrated <- RunUMAP(object = glioma.integrated, reduction = "pca", dims = 1:30) # Again, you can play around with nrs of dims
glioma.integrated <- FindNeighbors(object = glioma.integrated, dims = 1:30)
glioma.integrated <- FindClusters(object = glioma.integrated, resolution = 1) # Play around with the resolution

DefaultAssay(object = glioma.integrated) = "RNA"


DimPlot(glioma.integrated, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")
#FeaturePlot(object = glioma.integrated, features = "RBFOX3")



levels(glioma.integrated$seurat_clusters) <- gsub("^(16|12|23|18)$",paste0("\\1. NE"),levels(glioma.integrated$seurat_clusters))
levels(glioma.integrated$seurat_clusters) <- gsub("^(4|6)$",paste0("\\1. TAM"),levels(glioma.integrated$seurat_clusters))
levels(glioma.integrated$seurat_clusters) <- gsub("^(1|2|5|8)$",paste0("\\1. T [astr]"),levels(glioma.integrated$seurat_clusters))
levels(glioma.integrated$seurat_clusters) <- gsub("^(13)$",paste0("\\1. T [?]"),levels(glioma.integrated$seurat_clusters))
levels(glioma.integrated$seurat_clusters) <- gsub("^(0|9|20|22)$",paste0("\\1. T [olig]"),levels(glioma.integrated$seurat_clusters))
levels(glioma.integrated$seurat_clusters) <- gsub("^(14)$",paste0("\\1. T+TAM [dupli?]"),levels(glioma.integrated$seurat_clusters))
levels(glioma.integrated$seurat_clusters) <- gsub("^(15)$",paste0("\\1. T+TAM [dupli?]"),levels(glioma.integrated$seurat_clusters))
levels(glioma.integrated$seurat_clusters) <- gsub("^(19)$",paste0("\\1. TAM+? [dupli?]"),levels(glioma.integrated$seurat_clusters))
levels(glioma.integrated$seurat_clusters) <- gsub("^(21)$",paste0("\\1. T+OD [dupli?]"),levels(glioma.integrated$seurat_clusters))
levels(glioma.integrated$seurat_clusters) <- gsub("^(17)$",paste0("\\1. TAM+OD [dupli?]"),levels(glioma.integrated$seurat_clusters))
levels(glioma.integrated$seurat_clusters) <- gsub("^(3|7|10|11)$",paste0("\\1. OD"),levels(glioma.integrated$seurat_clusters))


levels(glioma.integrated$seurat_clusters)



glioma.integrated$seurat_clusters <- factor(glioma.integrated$seurat_clusters, levels=c(
  "1. T [astr]",
  "2. T [astr]",
  "5. T [astr]",
  "8. T [astr]" ,
  "14. T+TAM [dupli?]",
  "15. T+TAM [dupli?]",
  "21. T+OD [dupli?]",
  "13. T [?]",
  "0. T [olig]",
  "9. T [olig]",
  "20. T [olig]",
  "22. T [olig]",
  "3. OD",
  "7. OD",
  "10. OD",
  "11. OD" ,
  "17. TAM+OD [dupli?]" ,
  "19. TAM+? [dupli?]",
  "4. TAM",
  "6. TAM",
  "19. TAM",
  "12. NE",
  "16. NE",
  "18. NE",
  "23. NE",
  "3",
  "7",
  "10",
  "11",
  "24",
  "25",
  "26",
  "27" 
))





## ... ----

FeaturePlot(object = glioma.integrated, features = "MOG")
#FeaturePlot(object = glioma.integrated, features = "OPALIN")
FeaturePlot(object = glioma.integrated, features = "PLP1")
FeaturePlot(object = glioma.integrated, features = "TMEM144")


## LGG: up.1 ----

# Cell cycle + histones

set <- dge.partially.paired.clusters %>%
  dplyr::filter(up.1) %>%
  dplyr::pull(gene_name)


set <- c(set, c("H4C1", "H3C2", "H2AC4", "H3C3",            "HIST1H4A","HIST1H3B", "HIST1H2AB", "HIST1H3C",
                "H1-6", "H3C7", "H2BC9", "H2BC11",          "HIST1H1T","HIST1H3F","HIST1H2BH","HIST1H2BH",
                "H2AC11", "H2BC12", "H2AC12", "H2BC13",     "HIST1H2AG","HIST1H2BK","HIST1H2AH","HIST1H2BL",
                "H2AC13", "H3C10", "H2AC14", "H2BC14",      "HIST1H2AI","HIST1H3H","HIST1H2AJ","H2BC14",
                "H2AC15", "H2AC16", "H1-5", "H3C11",        "HIST1H2AK","HIST1H2AL","HIST1H1B","HIST1H3I",
                "H3C12", "H2BC17",                          "HIST1H3J","HIST1H2BO"))
set <- unique(set)

DotPlot(object = glioma.integrated, features = set, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [LGG.up1/ccycle] in: ",sid))



## LGG: up.2 ----

# COL1A2 + CD248
set <- dge.partially.paired.clusters %>% dplyr::filter(up.2) %>%  dplyr::pull(gene_name)

DotPlot(object = glioma.integrated, features = set, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [LGG.up2/col] in: ",sid))




## LGG: up.3 ----

# ??
set <- dge.partially.paired.clusters %>% dplyr::filter(up.3) %>%  dplyr::pull(gene_name)

DotPlot(object = glioma.integrated, features = list('UP3'=set, 'NE'=c("RBFOX3")), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [LGG.up3/???] in: ",sid))



## LGG: down ----


set1 <- dge.partially.paired.clusters %>% dplyr::filter(down.1) %>%  dplyr::pull(gene_name)
set2 <- dge.partially.paired.clusters %>% dplyr::filter(down.2) %>%  dplyr::pull(gene_name)
set <- list('down.1'=set1, 'down.2'=set2)

DotPlot(object = glioma.integrated, features = set, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [Down] in: ",sid))

## venteicher ----

set <- lgg.transcr.prog.astro
DotPlot(object = glioma.integrated, features = set, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [Vent/ASTR] in: ",sid))


c4d <- dge.partially.paired.clusters %>%  dplyr::filter(down.1 | down.2) %>% dplyr::pull(gene_name)

th <- FindMarkers(glioma.integrated, ident.1 = c(8,2,5,1),
                            ident.2 = c(20,22,0,9,6,4,7,11,3,10,16,23,12,25,27))

th <- th %>% 
  #rownames_to_column('gene_symbol') %>% 
  dplyr::mutate(prog.astro = gene_symbol %in% lgg.transcr.prog.astro) %>% 
  dplyr::mutate(prog.olig = gene_symbol %in% lgg.transcr.prog.oligo) %>% 
  dplyr::mutate(C4.down = gene_symbol %in% 
                  c4d
                )

saveRDS(th,"astro-specific.Rds")



set <- lgg.transcr.prog.oligo
DotPlot(object = glioma.integrated, features = set, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [Vent/OLIG] in: ",sid))

set <- lgg.transcr.prog.stemn
DotPlot(object = glioma.integrated, features = set, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [Vent/STEMN] in: ",sid))






