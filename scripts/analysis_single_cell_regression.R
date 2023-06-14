#!/usr/bin/env R

# load stuff genes ----

source('scripts/load_metadata.R')

source('scripts/load_hclust.R')

if('expression.glass.exon.metadata' %in% colnames(dge.partially.paired.clusters) == F) {
  source("scripts/load_correlation_expression_purity.R")
}



# load Johnson ----

seurat_obj_johnson <- SeuratDisk::LoadH5Seurat("data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_counts.h5")

rename <- read.table('cache/syn25956426_Johnson_cell_identifiers.txt', header=T) |> 
  dplyr::rename(new.name = umi) |> 
  dplyr::mutate(old.name = paste0("Cell",1:dplyr::n())) 
seurat_obj_johnson <- Seurat::RenameCells(seurat_obj_johnson, new.names = rename$new.name)

metadata <- rename |> 
  dplyr::left_join(read.csv("data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_metadata.tsv", sep="\t"), by=c('new.name'='cell_barcode')) |> 
  dplyr::left_join(
    read.csv('data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_syn25880693_clinical_metadata.csv'),
    by=c('case_barcode'='case_barcode')
  )


stopifnot(colnames(seurat_obj_johnson) == metadata$new.name)


for(slot in c("cell_state","case_barcode","case_sex","tumor_location","laterality","driver_mutations","time_point","idh_codel_subtype","who_grade","histological_classification")) {
  seurat_obj_johnson[[slot]] <- metadata[[slot]]
}





seurat_obj_johnson <- Seurat::NormalizeData(object = seurat_obj_johnson, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = F)
seurat_obj_johnson <- Seurat::FindVariableFeatures(object = seurat_obj_johnson, selection.method = "vst", verbose = F)
seurat_obj_johnson <- Seurat::ScaleData(object = seurat_obj_johnson, verbose = T, features = c(Seurat::VariableFeatures(object = seurat_obj_johnson),dge.partially.paired.clusters$gene_name))

seurat_obj_johnson <- Seurat::RunPCA(reduction.key = "PC_", object = seurat_obj_johnson, features = Seurat::VariableFeatures(object = seurat_obj_johnson),verbose = F)
seurat_obj_johnson <- Seurat::FindNeighbors(object = seurat_obj_johnson, dims = 1:35, verbose = F)
seurat_obj_johnson <- Seurat::FindClusters(object = seurat_obj_johnson, resolution = 1.2, algorithm = 1, verbose = F)
seurat_obj_johnson <- Seurat::RunUMAP( object = seurat_obj_johnson, dims = 1:35, verbose = T )




rm(rename, metadata, slot)


# load Bolleboom - Gao ----


# # take processed data from G-SAM cache
# seurat_obj_bolleboom_h243gbm <- readRDS(file="cache/seurat_bolleboom.2.Rds")
# seurat_obj_bolleboom_h243gbm$yaxis <- paste0(seurat_obj_bolleboom_h243gbm$cellt, " ", seurat_obj_bolleboom_h243gbm$seurat_clusters)
# 
# 
# plt.bolleboom <- Seurat::DotPlot(seurat_obj_bolleboom_h243gbm,
#                                   features=dge.partially.paired.clusters |> 
#                                     dplyr::arrange(desc(hclust_rank)) |> 
#                                     dplyr::pull(gene_name),
#                                   col.min = -4.6,
#                                   col.max = 4.6,
#                                   group.by = 'annotated_clusters')
# 
# 
# rm(seurat_obj_bolleboom_h243gbm)
# gc()
# 
# 
# 
# data.bolleboom.ac <- plt.bolleboom$data |> 
#   dplyr::mutate(pct.exp=NULL) |> # not plotted
#   dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
#   dplyr::filter(id %in% c("15. AC")) |> 
#   tibble::remove_rownames() |> 
#   tibble::column_to_rownames('features.plot') |> 
#   dplyr::mutate(id = NULL ) |> 
#   dplyr::rename_with( ~ paste0("Astrocytes [in GBM, Bolleboom] ", .x)) |> 
#   tibble::rownames_to_column('gene_name')
# 
# 
# 
# data.bolleboom.ne <- plt.bolleboom$data |> 
#   dplyr::mutate(pct.exp=NULL) |> # not plotted
#   dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
#   dplyr::filter(id == "4. NE") |> 
#   tibble::remove_rownames() |> 
#   tibble::column_to_rownames('features.plot') |> 
#   dplyr::mutate(id = NULL ) |> 
#   dplyr::rename_with( ~ paste0("Neurons [in GBM, Bolleboom]", .x)) |> 
#   tibble::rownames_to_column('gene_name')
# 
# 




# integrate / merge ----

#metadata_bolleboom_h243gbm <- seurat_obj_bolleboom_h243gbm@meta.data <- data.frame()


#seurat_obj_johnson@meta.data <- data.frame()
#seurat_obj_bolleboom_h243gbm@meta.data <- data.frame()


nrow(seurat_obj_johnson)
nrow(seurat_obj_bolleboom_h243gbm)
length(union(rownames(seurat_obj_johnson),rownames(seurat_obj_bolleboom_h243gbm)))
isct <- intersect(rownames(seurat_obj_johnson),rownames(seurat_obj_bolleboom_h243gbm))
length(isct)


glioma.combined <- merge(x = seurat_obj_johnson[,1:ncol(seurat_obj_johnson)-1], # wow strange error otherwise... https://github.com/satijalab/seurat/issues/6098
                         y = seurat_obj_bolleboom_h243gbm,
                         add.cell.ids = c("j", "b"), project = "glioma")

nrow(glioma.combined)



glioma.combined <- glioma.combined[rownames(glioma.combined) %in% isct,]

stopifnot(nrow(glioma.combined) == length(isct))




table(glioma.combined@meta.data$idh_codel_subtype)


glioma.combined@meta.data <- glioma.combined@meta.data |> 
  dplyr::mutate(idh_codel_subtype = ifelse(is.na(idh_codel_subtype), "IDHwt", idh_codel_subtype))


table(glioma.combined@meta.data$idh_codel_subtype)



glioma.combined <- Seurat::NormalizeData(object = glioma.combined, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = F)
glioma.combined <- Seurat::FindVariableFeatures(object = glioma.combined, selection.method = "vst", verbose = F)
glioma.combined <- Seurat::ScaleData(object = glioma.combined, verbose = T,
                                     features = c(Seurat::VariableFeatures(object = glioma.combined),dge.partially.paired.clusters$gene_name)
                                     )

glioma.combined <- Seurat::RunPCA(reduction.key = "PC_", object = glioma.combined, features = Seurat::VariableFeatures(object = glioma.combined),verbose = F)
glioma.combined <- Seurat::FindNeighbors(object = glioma.combined, dims = 1:35, verbose = F)
glioma.combined <- Seurat::FindClusters(object = glioma.combined, resolution = 1.2, algorithm = 1, verbose = F)
glioma.combined <- Seurat::RunUMAP( object = glioma.combined, dims = 1:35, verbose = T )


glioma.combined@meta.data <- glioma.combined@meta.data |> 
  dplyr::mutate(cell_type_combined = ifelse(is.na(cell_state), cell_type, cell_state)) |> 
  dplyr::mutate(cell_type_combined = dplyr::recode(cell_type_combined, 
                                                   `OD`='oligodendrocyte',
                                                   `Oligodendrocyte`='oligodendrocyte',
                                                   `EN`='endothelial',
                                                   `NE`='neuron',
                                                   `AC`='astrocyte',
                                                   `Pericyte`='pericyte',
                                                   `Fibroblast`='fibroblast',
                                                   `PE`='pericyte',
                                                   `T`='tumor',
                                                   `B cell`='B-cell',
                                                   `Endothelial`='endothelial',
                                                   `T cell`='T-cell',
                                                   `OPC`='OPC',
                                                   `Myeloid` = 'TAM',
                                                   `Dendritic cell` = 'TAM',
                                                   `Granulocyte` = 'TAM',
                                                   `Diff.-like` = 'tumor',
                                                   `Prolif. stem-like` = 'tumor',
                                                   `Stem-like` = 'tumor'
                                                   )) |> 
  dplyr::mutate(tumor_type_str = dplyr::case_when(
    idh_codel_subtype == "IDHwt" ~ "tumor [GBM]",
    idh_codel_subtype == "IDHmut_codel" ~ "tumor [1p/19q]",
    idh_codel_subtype == "IDHmut_noncodel" ~ "tumor [IDH+]"
  )) |> 
  dplyr::mutate(cell_type_combined = ifelse(cell_type_combined == "tumor", tumor_type_str, cell_type_combined))

Seurat::DimPlot(glioma.combined, reduction = "umap", label = TRUE, pt.size = .6, group.by = "cell_type_combined")
Seurat::FeaturePlot(glioma.combined, reduction = "umap", features="TOP2A", label = TRUE)

c39 <- subset(glioma.combined, seurat_clusters == 39)



# Seurat::DimPlot(data.sub, reduction = "umap", label = TRUE, pt.size = .6, group.by = "cell_type_combined")
# Seurat::DimPlot(data.sub, reduction = "umap", label = TRUE, pt.size = .6, group.by = "cell_type_combined")
# 
# 
#  


# regression ----

# quicker subsetting later on
data.sub <- subset(glioma.combined, )[rownames(glioma.combined) %in% dge.partially.paired.clusters$gene_name,]
#data.sub$cell_type_combined <- ifelse(  data.sub$seurat_clusters %in% c(39,26,21,22,38,23),  paste0(data.sub$cell_type_combined ," [cycling]"),  data.sub$cell_type_combined )
data.sub$cell_type_combined[data.sub$cell_type_combined %in% c('fibroblast','pericyte')] <- 'fibroblast & pericyte'
data.sub <- subset(data.sub, cell_type_combined %in% c('tumor [1p/19q]', 'tumor [GBM]',
                                                       
                                               'tumor [1p/19q] [cycling]', 'tumor [GBM] [cycling]',
                                               'endothelial [cycling]',# small cluster. probably wrong
                                               'oligodendrocyte [cycling]',# small cluster. probably wrong
                                               'T-cell [cycling]' # small cluster. probably wrong
) == F)



out.n.estimate <- data.frame()
out.n.stderr <- data.frame()
out.n.tstat <- data.frame()

out.s.estimate <- data.frame()
out.s.stderr <- data.frame()
out.s.tstat <- data.frame()



i = 1
n = length(dge.partially.paired.clusters$gene_name)
#for(gene in c(dge.partially.paired.clusters$gene_name)[1:604 %% 40 == 1]) {
#for(gene in c("TOP2A","COL1A2","SELL","IGSF3","F5")) {
for(gene in dge.partially.paired.clusters$gene_name) {
  print(paste0(gene, " -- " ,i,"/",n , " = ", round(i / n * 100,1), "%"))
  if(gene %in% rownames(data.sub)) {
    
    data <- subset(data.sub, )[rownames(data.sub) == gene,]
    
    if(nrow(data) == 1) {
      df.n <- data.frame(expr.norm = data@assays$RNA@data[1,])
      df.s <- data.frame(expr.norm = data@assays$RNA@scale.data[1,])
      
      
      # for(celltype in unique(data@meta.data$cell_type_combined)) {
      #   df.n[[celltype]] <- as.numeric(data@meta.data$cell_type_combined == celltype)
      #   df.s[[celltype]] <- as.numeric(data@meta.data$cell_type_combined == celltype)
      # }
      for(cluster in levels(data.sub$seurat_clusters)) {
        df.n[[paste0("c",cluster)]] <- as.numeric(data.sub$seurat_clusters == cluster)
        df.s[[paste0("c",cluster)]] <- as.numeric(data.sub$seurat_clusters == cluster)
      }
      
      
      res.n <- lm(expr.norm ~ 0 + ., data=df.n)
      est <- as.data.frame(summary(res.n)$coef) |> 
        dplyr::select('Estimate') |> 
        t() |> 
        as.data.frame() 
      rownames(est) <- gene
      
      
      err <- as.data.frame(summary(res.n)$coef) |> 
        dplyr::select('Std. Error') |> 
        t() |> 
        as.data.frame() 
      rownames(err) <- gene
      
      
      tstat <- as.data.frame(summary(res.n)$coef) |> 
        dplyr::select('t value') |> 
        t() |> 
        as.data.frame() 
      rownames(tstat) <- gene
      
      out.n.estimate <- rbind(out.n.estimate, est)
      out.n.stderr <- rbind(out.n.stderr, err)
      out.n.tstat <- rbind(out.n.tstat, tstat)
      
      rm(est,err,tstat, res.n)
      
      
      
      res.s <- lm(expr.norm ~ 0 + ., data=df.s)
      est <- as.data.frame(summary(res.s)$coef) |> 
        dplyr::select('Estimate') |> 
        t() |> 
        as.data.frame() 
      rownames(est) <- gene
      
      
      err <- as.data.frame(summary(res.s)$coef) |> 
        dplyr::select('Std. Error') |> 
        t() |> 
        as.data.frame() 
      rownames(err) <- gene
      
      
      tstat <- as.data.frame(summary(res.s)$coef) |> 
        dplyr::select('t value') |> 
        t() |> 
        as.data.frame() 
      rownames(tstat) <- gene
      
      out.s.estimate <- rbind(out.s.estimate, est)
      out.s.stderr <- rbind(out.s.stderr, err)
      out.s.tstat <- rbind(out.s.tstat, tstat)
      
      rm(est,err,tstat, res.s)
      
      
  
    } else {
      print(paste0("Skipping: ",gene))
    }
  } else {
    print(paste0("Skipping **: ",gene))
  }
  
  i <- i + 1
  
}



#saveRDS(out.s.estimate, file="cache/out.s.estimate.Rds")
#saveRDS(out.s.stderr, file="cache/out.s.stderr.Rds")
#saveRDS(out.s.tstat , file="cache/out.s.tstat.Rds")



# ## plot ----


out.s.estimate <- readRDS("cache/out.s.estimate.Rds")
out.s.estimate <- readRDS("cache/out.s.stderr.Rds")
out.s.estimate <- readRDS("cache/out.s.tstat.Rds")




df <- data.frame(
  gene_name = dge.partially.paired.clusters |> 
    dplyr::arrange(hclust_rank) |> 
    dplyr::pull(gene_name)
 ) |>
  #dplyr::left_join(out.stderr |> tibble::rownames_to_column('gene_name'),by=c('gene_name'='gene_name'),suffix=c('','')) |> 
  #dplyr::left_join(out.n.estimate |> tibble::rownames_to_column('gene_name'),by=c('gene_name'='gene_name'),suffix=c('','')) |> 
  dplyr::left_join(out.s.estimate |> tibble::rownames_to_column('gene_name'),by=c('gene_name'='gene_name'),suffix=c('','')) |> 
  dplyr::mutate(order.y = 1:dplyr::n()) |> 
  tidyr::pivot_longer(cols=-c(gene_name, order.y))

plt <- rbind(df,
             df |>  dplyr::mutate(value = ifelse(!is.na(value),0,value)))



ggplot(plt, aes(x=value, y=reorder(gene_name, order.y), group = gene_name)) + 
  facet_grid(cols = vars(name)) +
  #geom_hline(yintercept = "COL1A2",col="red",alpha=0.2) +  # COL1A2
  geom_line() +
  labs(y=NULL) + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ) 




data.sub$glial.cluster <- data.sub$seurat_clusters %in% c(18,20,24,25,27,29,32,37, 40, 43, 9) # 37 specifiek
Seurat::DimPlot(data.sub, group.by='glial.cluster')




celltype.per.cluster <- data.frame(
  cluster = paste0('c',data.sub$seurat_clusters),
  celltype = data.sub$cell_type_combined
) |> 
  dplyr::group_by(cluster) |> 
  dplyr::summarize (celltype =names(which.max(table(celltype))))


dataset.per.cluster <- data.frame(
  cluster = paste0('c',data.sub$seurat_clusters),
  dataset = ifelse(is.na(data.sub$cell_state),"Bolleboom","Johnson")
) |> 
  dplyr::group_by(cluster) |> 
  dplyr::summarize (dataset =names(which.max(table(dataset))))


sample.per.cluster <- data.frame(
  cluster = paste0('c',data.sub$seurat_clusters),
  sample = ifelse(is.na(data.sub$case_barcode),"Bolleboom - Gao" ,data.sub$case_barcode)
) |> 
  dplyr::group_by(cluster) |> 
  dplyr::summarize (sample =names(which.max(table(sample))))




df <- data.frame(
  gene_name = dge.partially.paired.clusters |> 
    dplyr::arrange(hclust_rank) |> 
    dplyr::pull(gene_name)
) |>
  dplyr::left_join(out.s.estimate |> tibble::rownames_to_column('gene_name'),by=c('gene_name'='gene_name'),suffix=c('','')) |> 
  dplyr::mutate(order.y = 1:dplyr::n()) |> 
  tidyr::pivot_longer(cols=-c(gene_name, order.y)) |> 
  dplyr::rename(cluster = name) |> 
  dplyr::left_join(celltype.per.cluster, by=c('cluster'='cluster'), suffix=c('','')) |> 
  dplyr::left_join(dataset.per.cluster, by=c('cluster'='cluster'), suffix=c('','')) |> 
  dplyr::left_join(sample.per.cluster, by=c('cluster'='cluster'), suffix=c('','')) |> 
  dplyr::mutate(cluster_rename = paste0(celltype, " [",cluster,"]"))


plt <- rbind(df,
             df |>  dplyr::mutate(value = ifelse(!is.na(value),0,value))) |> 
  dplyr::filter(cluster %in% c('c18','c41','c21','c22','c39','c37')) |>  # ,'c26','c38','c23',
  dplyr::mutate(type = "line")


ggplot(plt, aes(x=value, y=reorder(gene_name, order.y), group = gene_name, col=celltype)) + 
  facet_grid(cols = vars(cluster_rename)) +
  geom_line(data=plt |> dplyr::filter(type == "line")) +
  labs(y=NULL) + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ) 




# data: v Hijfte AC ----


van_hijfte_astro <- readRDS("../scRNA-glioma-workflows/cache/van Hijfte AC__van Hijfte AC - 3PR__[A-IDH]__master.Rds")
Seurat::DefaultAssay(van_hijfte_astro) <- "SCT"

van_hijfte_astro$yaxis <- paste0(van_hijfte_astro$celltype_annotated , 
                                 ifelse(!is.na(van_hijfte_astro$celltype_annotated_venteicher),
                                               paste0(" - ",van_hijfte_astro$celltype_annotated_venteicher)
                                               ,""),
                                 " - C" , van_hijfte_astro$seurat_clusters)




plt.van_hijfte_astros <- Seurat::DotPlot(van_hijfte_astro,
                                         features=dge.partially.paired.clusters |> 
                                           dplyr::arrange(desc(hclust_rank)) |> 
                                           dplyr::pull(gene_name),
                                         col.min = -5.6,
                                         col.max = 5.6,
                                         group.by = 'yaxis')



data.van_hijfte_astros.pericytes.c28 <- plt.van_hijfte_astros |> 
  purrr::pluck('data') |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id %in% c("pericyte - C28")) |> 
  tibble::remove_rownames() |> 
  dplyr::filter(!is.na(features.plot)) |>  # some gene with rna_ prefix is having an NA value
  dplyr::filter(!is.na(avg.exp.scaled)) |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("van Hijfte pericyte #28 ", .x)) |> 
  tibble::rownames_to_column('gene_name')



data.van_hijfte_astros.tumor.c1 <- plt.van_hijfte_astros |> 
  purrr::pluck('data') |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id %in% c("tumor - Remaining - C1")) |> 
  tibble::remove_rownames() |> 
  dplyr::filter(!is.na(features.plot)) |>  # some gene with rna_ prefix is having an NA value
  dplyr::filter(!is.na(avg.exp.scaled)) |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("van Hijfte tumor - Undetermined #1 ", .x)) |> 
  tibble::rownames_to_column('gene_name')



data.van_hijfte_astros.cycling.c23 <- plt.van_hijfte_astros |> 
  purrr::pluck('data') |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id %in% c("tumor - Cycling - C23")) |> 
  tibble::remove_rownames() |> 
  dplyr::filter(!is.na(features.plot)) |>  # some gene with rna_ prefix is having an NA value
  dplyr::filter(!is.na(avg.exp.scaled)) |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("van Hijfte cycling tumor #23 ", .x)) |> 
  tibble::rownames_to_column('gene_name')


data.van_hijfte_astros.aclike.c4 <- plt.van_hijfte_astros |> 
  purrr::pluck('data') |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id %in% c("tumor - Astro. program - C4")) |> 
  tibble::remove_rownames() |> 
  dplyr::filter(!is.na(features.plot)) |>  # some gene with rna_ prefix is having an NA value
  dplyr::filter(!is.na(avg.exp.scaled)) |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("van Hijfte AC-like tumor #4 ", .x)) |> 
  tibble::rownames_to_column('gene_name')



data.van_hijfte_astros.astro.c27 <- plt.van_hijfte_astros |> 
  purrr::pluck('data') |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id %in% c("astrocyte - C27")) |> 
  tibble::remove_rownames() |> 
  dplyr::filter(!is.na(features.plot)) |>  # some gene with rna_ prefix is having an NA value
  dplyr::filter(!is.na(avg.exp.scaled)) |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("van Hijfte astrocytes #27 ", .x)) |> 
  tibble::rownames_to_column('gene_name')






Seurat::DimPlot(van_hijfte_astro)
Seurat::DimPlot(van_hijfte_astro,group.by="yaxis", label=T)
Seurat::FeaturePlot(van_hijfte_astro, features=c('COL1A2'))


# data: Yuan - PJ016 ----
# PJ16 = IDH

yuan_PJ016 <- readRDS("../scRNA-glioma-workflows/cache/Yuan [GSE103224]__PJ016__[A-IDH]__master.Rds")
yuan_PJ016$yaxis <- paste0(yuan_PJ016$celltype.annotated.venteicher, " ", yuan_PJ016$seurat_clusters)



plt.yuan_pj016 <- Seurat::DotPlot(yuan_PJ016,
                features=dge.partially.paired.clusters |> 
                  dplyr::arrange(desc(hclust_rank)) |> 
                  dplyr::pull(gene_name),
                col.min = -5.6,
                col.max = 5.6,
                group.by = 'yaxis')



data.yuan.stem9 <- plt.yuan_pj016$data |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id %in% c("Stemness. program 9")) |> 
  tibble::remove_rownames() |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("PJ016 Stemness 9 ", .x)) |> 
  tibble::rownames_to_column('gene_name')


data.yuan.ac1 <- plt.yuan_pj016$data |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id == "Astro. program 1") |> 
  tibble::remove_rownames() |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("PJ016 AC-like 1 ", .x)) |> 
  tibble::rownames_to_column('gene_name')


# data: Diaz 2019 - SF11136 ----


diaz_2019_SF11136 <- readRDS("../scRNA-glioma-workflows/cache/Diaz 2019 [GSE138794]__SF11136__[A-IDH]__master.Rds")
diaz_2019_SF11136$yaxis <- paste0(diaz_2019_SF11136$celltype_annotated, 
                                
                                ifelse(
                                  !is.na(diaz_2019_SF11136$celltype_annotated_venteicher),
                                  paste0(" ", diaz_2019_SF11136$celltype_annotated_venteicher),
                                  ""
                                ),
                                " ",diaz_2019_SF11136$seurat_clusters
                                )


plt.diaz_SF11136 <- Seurat::DotPlot(diaz_2019_SF11136,
                features=dge.partially.paired.clusters |> 
                  dplyr::arrange(desc(hclust_rank)) |> 
                  dplyr::pull(gene_name),
                col.min = -5.6,
                col.max = 5.6,
                group.by = 'yaxis')



data.diaz.pe <- plt.diaz_SF11136$data |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id == "pericyte 11") |> 
  tibble::remove_rownames() |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("Diaz SF11136 pericyte", .x)) |> 
  tibble::rownames_to_column('gene_name')



data.diaz.ac2 <- plt.diaz_SF11136$data |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id == "tumor Astro. program 2") |> 
  tibble::remove_rownames() |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("Diaz SF11136 AC-like 2", .x)) |> 
  tibble::rownames_to_column('gene_name')




## clusters:



# data: Johnson - only IDH-mut ----


seurat_obj_johnson <- SeuratDisk::LoadH5Seurat("data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_counts.h5")

rename <- read.table('cache/syn25956426_Johnson_cell_identifiers.txt', header=T) |> 
  dplyr::rename(new.name = umi) |> 
  dplyr::mutate(old.name = paste0("Cell",1:dplyr::n())) 
seurat_obj_johnson <- Seurat::RenameCells(seurat_obj_johnson, new.names = rename$new.name)

metadata <- rename |> 
  dplyr::left_join(read.csv("data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_metadata.tsv", sep="\t"), by=c('new.name'='cell_barcode')) |> 
  dplyr::left_join(
    read.csv('data/syn25956426_Johnson/processed_data/analysis_scRNAseq_tumor_syn25880693_clinical_metadata.csv'),
    by=c('case_barcode'='case_barcode')
  )


stopifnot(colnames(seurat_obj_johnson) == metadata$new.name)


for(slot in c("cell_state","case_barcode","case_sex","tumor_location","laterality","driver_mutations","time_point","idh_codel_subtype","who_grade","histological_classification")) {
  seurat_obj_johnson[[slot]] <- metadata[[slot]]
}


rm(rename, metadata, slot)


seurat_obj_johnson <- subset(seurat_obj_johnson, idh_codel_subtype == "IDHmut_noncodel")



seurat_obj_johnson <- Seurat::NormalizeData(object = seurat_obj_johnson, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = F)
seurat_obj_johnson <- Seurat::FindVariableFeatures(object = seurat_obj_johnson, selection.method = "vst", verbose = F)
seurat_obj_johnson <- Seurat::ScaleData(object = seurat_obj_johnson, verbose = T, features = c(Seurat::VariableFeatures(object = seurat_obj_johnson),dge.partially.paired.clusters$gene_name))
seurat_obj_johnson <- Seurat::RunPCA(reduction.key = "PC_", object = seurat_obj_johnson, features = Seurat::VariableFeatures(object = seurat_obj_johnson),verbose = F)
seurat_obj_johnson <- Seurat::FindNeighbors(object = seurat_obj_johnson, dims = 1:35, verbose = F)
seurat_obj_johnson <- Seurat::FindClusters(object = seurat_obj_johnson, resolution = 1.2, algorithm = 1, verbose = F)
seurat_obj_johnson <- Seurat::RunUMAP( object = seurat_obj_johnson, dims = 1:35, verbose = T )


seurat_obj_johnson$annotated_clusters <- paste0(seurat_obj_johnson$seurat_clusters, " ", seurat_obj_johnson$cell_state)


# 
# Seurat::DotPlot(seurat_obj_johnson,
#                 features=dge.partially.paired.clusters |>
#                   dplyr::arrange(desc(hclust_rank)) |>
#                   dplyr::pull(gene_name),
#                 col.min = -4.6,
#                 col.max = 4.6,
#                 group.by = 'annotated_clusters')



plt.johnson <- Seurat::DotPlot(seurat_obj_johnson,
                features=dge.partially.paired.clusters |>
                  dplyr::arrange(desc(hclust_rank)) |>
                  dplyr::pull(gene_name),
                col.min = -5.6,
                col.max = 5.6,
                group.by = 'seurat_clusters')



# 30. Fibroblast / pericyte
# 23. Cycling / tumor
# 19. AC-like tumor
# 11. AC-like tumor



data.johnson.pe30 <- plt.johnson$data |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id == 30) |> 
  tibble::remove_rownames() |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("Fibroblast / pericyte #30 Johnson ", .x)) |> 
  tibble::rownames_to_column('gene_name')


data.johnson.cycling23  <- plt.johnson$data |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id == 23) |> 
  tibble::remove_rownames() |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("cycling cells #23 Johnson ", .x)) |> 
  tibble::rownames_to_column('gene_name')



data.johnson.aclike11  <- plt.johnson$data |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id == 11) |> 
  tibble::remove_rownames() |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("AC-like tumor #11 Johnson ", .x)) |> 
  tibble::rownames_to_column('gene_name')


data.johnson.aclike19  <- plt.johnson$data |> 
  dplyr::mutate(pct.exp=NULL) |> # not plotted
  dplyr::mutate(avg.exp=NULL) |> # not very useful for plotting
  dplyr::filter(id == 19) |> 
  tibble::remove_rownames() |> 
  tibble::column_to_rownames('features.plot') |> 
  dplyr::mutate(id = NULL ) |> 
  dplyr::rename_with( ~ paste0("AC-like tumor #19 Johnson ", .x)) |> 
  tibble::rownames_to_column('gene_name')









# annotate clusters ----
# 
# seurat_obj_johnson$celltype_annotated = seurat_obj_johnson$cell_state
# 
# ## cycling : 23 ----
# 
# 
# Seurat::FeaturePlot(seurat_obj_johnson, features="TOP2A",label=T)
# Seurat::FeaturePlot(seurat_obj_johnson, features="cell_state")
# #Seurat::DimPlot(seurat_obj_johnson, group.by= "cell_state")
# 
# 
# seurat_obj_johnson$celltype_annotated = ifelse(seurat_obj_johnson$seurat_clusters %in% c(23),
#                                              "cycling cells [mostly tumor]",seurat_obj_johnson$celltype_annotated)
# 
# ## tumor ----
# 
# 
# Seurat::FeaturePlot(seurat_obj_johnson,
#                     reduction = "umap",
#                     features = celltype_markers_Venteicher_IDHmut |> dplyr::filter(Astro.program) |>  head(n=8) |> dplyr::pull(gene_symbol),
#                     label=T)
# 
# seurat_obj_johnson$celltype_annotated = ifelse(seurat_obj_johnson$seurat_clusters %in% c(
#   11,8,19,28,
#   7,3,24,
#   14,21,13,
#   16,22
# ),
# "tumor",seurat_obj_johnson$celltype_annotated)
# Seurat::DimPlot(seurat_obj_johnson, group.by= "celltype_annotated")
# 
# 
# 
# 
# Seurat::DimPlot(seurat_obj_johnson, group.by= "seurat_clusters", label=T)
# Seurat::DimPlot(seurat_obj_johnson, group.by= "cell_state", label=T)
# #Seurat::DimPlot(seurat_obj_johnson, group.by= "celltype_annotated", label=T)
# 
# ## OD ----
# 
# seurat_obj_johnson$celltype_annotated = ifelse(seurat_obj_johnson$seurat_clusters %in% c(
#   4,17,15,10,20,1,9
# ),
# "Oligodendrocyte",seurat_obj_johnson$celltype_annotated)
# Seurat::DimPlot(seurat_obj_johnson, group.by= "celltype_annotated")
# 
# 
# Seurat::DimPlot(seurat_obj_johnson, group.by= "seurat_clusters", label=T)
# Seurat::DimPlot(seurat_obj_johnson, group.by= "cell_state", label=T)
# 
# 
# ## left overs 
# 
# 
# seurat_obj_johnson$celltype_annotated = ifelse(seurat_obj_johnson$cell_state %in% c(
#   "Diff.-like", "Stem-like"
# ),
# "tumor",seurat_obj_johnson$celltype_annotated)
# Seurat::DimPlot(seurat_obj_johnson, group.by= "celltype_annotated")
# 
# 
# 
# 
# seurat_obj_johnson$yaxis <- paste0(seurat_obj_johnson$celltype_annotated, " ", seurat_obj_johnson$seurat_clusters)
# 
# 
# 
# 
# 
# 
# Seurat::DotPlot(seurat_obj_johnson,
#                 features=dge.partially.paired.clusters |> 
#                   dplyr::arrange(desc(hclust_rank)) |> 
#                   dplyr::pull(gene_name),
#                 group.by = 'yaxis')
# 
# 
# 
# 
# 
# 
# Seurat::DotPlot(seurat_obj_johnson,
#                 features=
#                   (df |> dplyr::filter(col == "up-3 (fuzzy)") |> dplyr::pull(gene_name))
#                   ,
#                 group.by = 'yaxis')
# 
# 
# 
# 
# 
# table(subset(seurat_obj_johnson, seurat_clusters == 23)$celltype_annotated)
# table(subset(seurat_obj_johnson, seurat_clusters == 31)$celltype_annotated)
# table(subset(seurat_obj_johnson, celltype_annotated == 'tumor')$seurat_clusters)
# 
# table(subset(seurat_obj_johnson, yaxis == 'tumor 23')$cell_state)
# table(subset(seurat_obj_johnson, yaxis == 'tumor 30')$cell_state)
# table(subset(seurat_obj_johnson, seurat_clusters == '30')$cell_state)
# 
# 
# 




# final plot ----


# saveRDS(df, file="/tmp/df.Rds")


df <- dge.partially.paired.clusters |>
       dplyr::arrange(hclust_rank) |>
    dplyr::select(gene_name,  gene_uid, up.1,  up.2 , up.3  ,down, hclust_rank ) |>
    dplyr::mutate(col = case_when(
      up.1 ~ "up-1 (collagen)",
      up.2 ~ "up-2 (cell cycling)",
      up.3 ~ "up-3 (fuzzy)",
      down ~ "down"
    )) |>
  dplyr::mutate(up.1 = NULL)|>
  dplyr::mutate(up.2 = NULL)|>
  dplyr::mutate(up.3 = NULL)|>
  dplyr::mutate(down = NULL) |>
  dplyr::mutate(cluster = 1) |>
  dplyr::left_join(
    expression.glass.exon.metadata |> dplyr::select(c(
      "gene_uid",
      #"cor.dna.wes.VAF_IDH.t.statistic",
      #"cor.dna.wes.VAF_IDH.cor.estimate"            ,
      #"cor.methylation.purity.absolute.t.statistic",
      "cor.methylation.purity.absolute.cor.estimate"
    )
    ), by=c('gene_uid'='gene_uid')
  ) |>
  
  
  dplyr::left_join(data.van_hijfte_astros.aclike.c4 , by=c('gene_name'='gene_name')) |>
  dplyr::left_join(data.van_hijfte_astros.astro.c27 , by=c('gene_name'='gene_name')) |>
  dplyr::left_join(data.van_hijfte_astros.cycling.c23 , by=c('gene_name'='gene_name')) |>
  dplyr::left_join(data.van_hijfte_astros.pericytes.c28 , by=c('gene_name'='gene_name')) |>
  dplyr::left_join(data.van_hijfte_astros.tumor.c1 , by=c('gene_name'='gene_name')) |> 
  dplyr::left_join(data.yuan.ac1, by=c('gene_name'='gene_name')) |>
  dplyr::left_join(data.yuan.stem9, by=c('gene_name'='gene_name')) |> 
  dplyr::left_join(data.diaz.pe, by=c('gene_name'='gene_name'))|>
  dplyr::left_join(data.diaz.ac2, by=c('gene_name'='gene_name')) |> 
  dplyr::left_join(data.johnson.pe30, by=c('gene_name'='gene_name')) |>
  dplyr::left_join(data.johnson.cycling23, by=c('gene_name'='gene_name')) |>
  dplyr::left_join(data.johnson.aclike11, by=c('gene_name'='gene_name')) |>
  dplyr::left_join(data.johnson.aclike19, by=c('gene_name'='gene_name'))


#df <- readRDS("tmp/plt_analysis_single_cell_reg.Rds")




plt <- df |> 
  tidyr::pivot_longer(cols = -c(gene_name, gene_uid, col, hclust_rank), names_to="panel") |> 
  dplyr::mutate(value_c = value) |> 
  dplyr::mutate(panel = dplyr::recode(panel, `cor.methylation.purity.absolute.cor.estimate`='AA corr tumor percentage')) |> 

  dplyr::mutate(panel = dplyr::recode(panel, `Diaz SF11136 AC-like 2avg.exp.scaled`='AC-like tumor #2 SF11137')) |> 
  dplyr::mutate(panel = dplyr::recode(panel, `Diaz SF11136 pericyteavg.exp.scaled`='pericyte #1 SF11137')) |> 
  
  dplyr::mutate(panel = dplyr::recode(panel, `PJ016 AC-like 1 avg.exp.scaled`='AC-like tumor #1 PJ016')) |> 
  dplyr::mutate(panel = dplyr::recode(panel, `PJ016 Stemness 9 avg.exp.scaled`='stemn. #9 PJ016')) |> 
  
  dplyr::mutate(celltype = case_when(
    panel == "AA corr tumor percentage" ~ "correlation bulk",
    
    panel == "AC-like tumor #1 PJ016" ~ "tumor [AC-like]",
    panel == "stemn. #9 PJ016" ~ "tumor [stem. / cycling]",
    
    panel == "AC-like tumor #11 Johnson avg.exp.scaled" ~ "tumor [AC-like]",
    panel == "AC-like tumor #19 Johnson avg.exp.scaled" ~ "tumor [AC-like]",
    panel == "cycling cells #23 Johnson avg.exp.scaled" ~ "tumor [stem. / cycling]",
    panel == "Fibroblast / pericyte #30 Johnson avg.exp.scaled" ~ "pericyte / CAF",
    
    panel == "AC-like tumor #2 SF11137" ~ "tumor [AC-like]",
    panel == "pericyte #1 SF11137" ~ "pericyte / CAF",
    
    panel == "van Hijfte AC-like tumor #4 avg.exp.scaled" ~ "tumor [AC-like]",
    panel == "van Hijfte astrocytes #27 avg.exp.scaled" ~ "astrocyte",
    panel == "van Hijfte cycling tumor #23 avg.exp.scaled" ~ "tumor [stem. / cycling]",
    panel == "van Hijfte pericyte #28 avg.exp.scaled" ~ "pericyte / CAF",
    panel == "van Hijfte tumor - Undetermined #1 avg.exp.scaled" ~ "tumor [undetermined]"
    
  )) |> 
  
  dplyr::mutate(celltype = factor(celltype, levels=
                                    c("correlation bulk",
                                      
                                      "tumor [undetermined]",  "pericyte / CAF",
                                      
                                      "tumor [stem. / cycling]",
                                      
                                      "tumor [AC-like]",         "astrocyte", 
                                      
                                      "neuron"
  ))) |> 

  dplyr::mutate(dataset = case_when(
    grepl("van Hijfte", panel) ~ "van Hijfte AC",
    grepl("Johnson", panel) ~ "Johnson",
    panel == "AA corr tumor percentage" ~ "GLASS-NL [bulk[",
    panel == "AC-like tumor #1 PJ016" ~ "Yuan [PJ016]",
    panel == "AC-like tumor #2 SF11137" ~ "Diaz 2019 [SF11136]",
    panel == "pericyte #1 SF11137" ~ "Diaz 2019 [SF11136]",
    panel == "stemn. #9 PJ016" ~ "Yuan [PJ016]"
  )) |> 
  
  dplyr::mutate(xorder = order(order(celltype, dataset))) |> 
  dplyr::mutate(panel =
                  
                  
                  dplyr::recode(panel, 
                                "AC-like tumor #1 PJ016"="Y-#1",
                                "stemn. #9 PJ016"="Y-#9",
                                
                                "van Hijfte AC-like tumor #4 avg.exp.scaled"="vH-#4",
                                "van Hijfte astrocytes #27 avg.exp.scaled"="vH-#H27",
                                "van Hijfte cycling tumor #23 avg.exp.scaled"="vH-#23",
                                "van Hijfte pericyte #28 avg.exp.scaled"="vH-#28",
                                 "van Hijfte tumor - Undetermined #1 avg.exp.scaled"="vH-#1",
                                
                                "pericyte #1 SF11137" = "D-#1",
                                "AC-like tumor #2 SF11137"="D-#2",
                                
                                "AC-like tumor #11 Johnson avg.exp.scaled"="J-#11",
                                "cycling cells #23 Johnson avg.exp.scaled" = "J-#23",
                                "Fibroblast / pericyte #30 Johnson avg.exp.scaled" = "J-#30",
                                "AC-like tumor #19 Johnson avg.exp.scaled"="J-#19"
                  )
                
                )


unique(plt$panel)



plt.0 <- rbind(plt, plt |>  dplyr::mutate(value=0)) |> 
  dplyr::mutate(col = ifelse(panel == "cluster", col, ".")) |> 
  dplyr::filter(panel == "AA corr tumor percentage")

plt.1 <- rbind(plt, plt |>  dplyr::mutate(value=0)) |> 
  dplyr::mutate(col = ifelse(panel == "cluster", col, ".")) |> 
  dplyr::filter(panel != "cluster") |> 
  dplyr::filter(panel != "AA corr tumor percentage")

plt.2 <- rbind(plt, plt |>  dplyr::mutate(value=0)) |> 
  dplyr::mutate(col = ifelse(panel == "cluster", col, ".")) |> 
  dplyr::filter(panel == "cluster")

plt.3 <- plt |> 
  dplyr::select(panel, celltype, xorder) |> 
  dplyr::filter(panel != "cluster") |> 
  dplyr::filter(panel != "AA corr tumor percentage") |> 
  dplyr::group_by(panel, celltype) |> 
  dplyr::filter(xorder == min(xorder)) |> 
  tidyr::pivot_longer(c(celltype)) |> 
  dplyr::mutate(y= paste0(name, ": ", value)) |> 
  dplyr::mutate(order = case_when(
    value == "correlation bulk" ~ 11,
    value == "tumor" ~ 10,
    value == "astrocyte [GBM]" ~ 9,
    T ~ 0
  ))

plt.4 <- plt |> 
  dplyr::select(panel, dataset, xorder) |> 
  dplyr::filter(panel != "cluster") |> 
  dplyr::filter(panel != "AA corr tumor percentage") |> 
  dplyr::group_by(panel, dataset) |> 
  dplyr::filter(xorder == min(xorder)) |> 
  tidyr::pivot_longer(c(dataset)) |> 
  dplyr::mutate(y= paste0(name, ": ", value)) 





p0 <- ggplot(plt.0, aes(x = value, y = reorder(gene_name, hclust_rank), col=value_c)) +
  facet_grid(cols = vars(panel)) + # , scales = "free_x"
  geom_line() +
  scale_colour_gradient2(
    low = "red",
    mid = "gray",
    high = "blue",
    midpoint = 0
  ) +
  #scale_color_manual(values=c('.'='gray30','up-1 (collagen)' = 'red','up-2 (cell cycling)'='darkgreen','up-3 (fuzzy)'='brown','down'='blue')) +
  labs(y=NULL, x=NULL) + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )  +
  theme(legend.position = 'bottom')


p1 <- ggplot(plt.1, aes(x = value, y = reorder(gene_name, hclust_rank), col=value_c)) +
  facet_grid(cols = vars(reorder(panel,xorder))) + # , scales = "free_x"
  geom_line() +
  scale_colour_gradient2(
    low = "red",
    mid = "gray",
    high = "blue",
    midpoint = 0
  ) +
  #scale_color_manual(values=c('.'='gray30','up-1 (collagen)' = 'red','up-2 (cell cycling)'='darkgreen','up-3 (fuzzy)'='brown','down'='blue')) +
  labs(y=NULL, x=NULL) + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )  +
  theme(legend.position = 'bottom')

p2 <- ggplot(plt.2, aes(x = value, y = reorder(gene_name, hclust_rank), col=col)) +
  geom_line() +
  scale_color_manual(values=c('.'='gray30','up-1 (collagen)' = 'red','up-2 (cell cycling)'='darkgreen','up-3 (fuzzy)'='brown','down'='blue')) +
  labs(y=NULL) + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )  + 
  labs(x=NULL, y=NULL)  +
  theme(legend.position = 'bottom')

p3 <- ggplot(plt.3, aes(x=reorder(panel,xorder), y=reorder(y, order), col=value)) +
  geom_point(size=4) + 
  guides(color = FALSE, size = FALSE) +
  labs(x=NULL, y=NULL) + 
  theme(
    #panel.background = element_rect(fill = 'white', colour = 'white')
        ) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme_bw()


p4 <- ggplot(plt.4, aes(x=reorder(panel, xorder), y=y)) +
  geom_point(size=4, col="gray40") + 
  guides(color = FALSE, size = FALSE) +
  labs(x=NULL, y=NULL) + 
  theme(
    #panel.background = element_rect(fill = 'white', colour = 'white')
    ) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme_bw()


p0 + p1 + p2 + p3 + p4 +  patchwork::plot_layout(widths = c(1,10, 1),heights = c(6,1,1), design= "
ABC
#D#
#E#
")


ggsave("output/figures/single_cell_clusters.pdf", width = 11, height = 8)


