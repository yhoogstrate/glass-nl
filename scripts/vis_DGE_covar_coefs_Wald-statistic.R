#!/usr/bin/env R 

# 1. - load DGE ----

# if(!exists('res.paired.a.exon')) {
#   source('scripts/load_analysis_DGE.R')
# }

dds.paired.a.covar.regression <- readRDS(file = "cache/dds.paired.a.covar.regression.Rds")

# 2. - load hcluster ----

if(!exists('dge.partially.paired.clusters')) {
  source('scripts/load_hclust.R')
}

# merge plot ----

chem <- DESeq2::results(dds.paired.a.covar.regression, contrast=c("status.chemo","chemo","no.chemo"))
rad <- DESeq2::results(dds.paired.a.covar.regression, contrast=c("status.radio","radio","no.radio"))
grad <- DESeq2::results(dds.paired.a.covar.regression, contrast=c("status.grading","Recurrent.High.Grade", "Recurrent.Low.Grade"))
rec <- DESeq2::results(dds.paired.a.covar.regression, contrast=c("Sample_Type","recurrent","initial"))


plt <- rbind(
  chem |> 
    as.data.frame() |> 
    tibble::rownames_to_column('gene_uid') |> 
    dplyr::mutate(covar = "chemo.therapy")
    ,
  rad |>
    as.data.frame() |>
    tibble::rownames_to_column('gene_uid') |> 
    dplyr::mutate(covar = "radio.therapy"),
    #head(n=4),
  grad |>
    as.data.frame() |>
    tibble::rownames_to_column('gene_uid') |> 
    dplyr::mutate(covar = "WHO.grading")
    ,
  rec |>
    as.data.frame() |>
    tibble::rownames_to_column('gene_uid') |> 
    dplyr::mutate(covar = "Recurrence type")
    #head(n=4)
)


plt.wide <- plt |>
  tidyr::pivot_wider(
    id_cols = gene_uid,
    names_from=covar,
    values_from = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  ) |> 
  as.data.frame()


## stat ----

col2 <- grDevices::colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                      "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                      "#4393C3", "#2166AC", "#053061"))

plt2 <-  dge.partially.paired.clusters |> 
  dplyr::left_join(plt.wide, by=c('gene_uid'='gene_uid'), suffix=c('','')) |> 
  dplyr::select(!contains("pvalue") & !contains("baseMean") & !contains("lfcSE") & !contains("log2FoldChange")& !contains("padj")) |> 
  tidyr::pivot_longer(cols=contains("stat"), values_to = "stat",names_to = "covar") |> 
  dplyr::mutate(covar = gsub("stat_","",covar)) |> 
  dplyr::mutate(covar = dplyr::recode(covar,
                                      "radio.therapy" = "Radio therapy [yes]",
                                      "chemo.therapy" = "Chemo therapy [yes]",
                                      "Recurrence type" = "At recurrence",
                                      "WHO.grading" = "WHO grade [4]"
                                      )) |> 
  dplyr::mutate(covar = factor(covar, levels=c("At recurrence", "WHO grade [4]", "Chemo therapy [yes]", "Radio therapy [yes]")))


plt3 <- rbind(
  plt2 |> dplyr::mutate(col=stat),
  plt2 |> dplyr::mutate(col=stat, stat = 0)
) 

ggplot(plt3, aes(x = stat, y = hclust_rank, col=col, group=gene_uid)) +
  facet_grid(cols = vars(covar)) +
  geom_line(lwd=0.75) +
  coord_cartesian(xlim=c(6, -6)) +
  scale_x_continuous(breaks = c(-3, 0, 3)) +
  theme_bw() +
  #geom_vline(xintercept=0,col="black",lwd=0.5,lty=2) +
  ggplot2::scale_color_gradientn(colours = col2(200), 
                                 na.value = "grey50",
                                 limits = c(-3.1, 3.1), 
                                 oob = scales::squish ) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  theme(legend.position = 'bottom')


ggsave("output/figures/vis_DGE_covar_coefs_stat.pdf",width=5,height=11)
ggsave("output/figures/vis_DGE_covar_coefs_stat.svg",width=5,height=11)




## padj ---


plt2 <-  dge.partially.paired.clusters |> 
  dplyr::left_join(plt.wide, by=c('gene_uid'='gene_uid'), suffix=c('','')) |> 
  dplyr::select(!contains("pvalue") & !contains("baseMean") & !contains("lfcSE") & !contains("log2FoldChange")& !contains("stat")) |> 
  tidyr::pivot_longer(cols=contains("padj"), values_to = "padj",names_to = "covar") |> 
  dplyr::mutate(covar = gsub("stat_","",covar))



ggplot(plt2, aes(x = padj, y = hclust_rank, col=down)) +
  facet_grid(cols = vars(covar)) +
  geom_point() +
  xlim(0,1) + 
  geom_vline(xintercept=0,col="black",lwd=0.5,lty=2) +
  geom_vline(xintercept=0.05,col="black",lwd=0.5,lty=2) +
  geom_vline(xintercept=0.01,col="black",lwd=0.5,lty=2)





plt2 <-  dge.partially.paired.clusters |> 
  dplyr::left_join(plt.wide, by=c('gene_uid'='gene_uid'), suffix=c('','')) |> 
  dplyr::select(!contains("pvalue") & !contains("baseMean") & !contains("lfcSE") & !contains("log2FoldChange")& !contains("stat")) |> 
  tidyr::pivot_longer(cols=contains("padj"), values_to = "padj",names_to = "covar") |> 
  dplyr::mutate(covar = gsub("padj_","",covar))



ggplot(plt2, aes(x = padj, y = hclust_rank, col=down)) +
  facet_grid(cols = vars(covar)) +
  geom_point() +
  xlim(0,1) + 
  geom_vline(xintercept=0,col="black",lwd=0.5,lty=2) +
  geom_vline(xintercept=0.05,col="black",lwd=0.5,lty=2) +
  geom_vline(xintercept=0.01,col="black",lwd=0.5,lty=2)






plt2 <-  dge.partially.paired.clusters |> 
  dplyr::left_join(plt.wide, by=c('gene_uid'='gene_uid'), suffix=c('','')) |> 
  dplyr::select(!contains("padj") & !contains("baseMean") & !contains("lfcSE") & !contains("log2FoldChange")& !contains("stat")) |> 
  tidyr::pivot_longer(cols=contains("pvalue"), values_to = "pvalue",names_to = "covar") |> 
  dplyr::mutate(covar = gsub("pvalue_","",covar))



ggplot(plt2, aes(x = pvalue, y = hclust_rank, col=down)) +
  facet_grid(cols = vars(covar)) +
  geom_point() +
  xlim(0,1) + 
  geom_vline(xintercept=0,col="black",lwd=0.5,lty=2) +
  geom_vline(xintercept=0.05,col="black",lwd=0.5,lty=2) +
  geom_vline(xintercept=0.01,col="black",lwd=0.5,lty=2)





# plot chemo ----

stats <- dge.partially.paired.clusters |> 
  dplyr::left_join(
    chem |> 
      as.data.frame() |> 
      tibble::rownames_to_column('gene_uid')
    , by=c('gene_uid'='gene_uid'), suffix=c('','')
  )

ggplot(stats, aes(x=hclust_rank, y=padj, col=up.1)) +
  geom_point()


# plot radio ----

stats <- dge.partially.paired.clusters |> 
  dplyr::left_join(
    rad |> 
      as.data.frame() |> 
      tibble::rownames_to_column('gene_uid')
    , by=c('gene_uid'='gene_uid'), suffix=c('','')
  )


ggplot(stats, aes(x=hclust_rank, y=padj, col=up.1)) +
  geom_point()


## use stats instead ----





