#!/usr/bin/env R

# load covar / coef factors ----


source('scripts/load_metadata.R')
source('scripts/load_hclust.R')



# make plt

plt <- dge.partially.paired.clusters |>
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
  dplyr::mutate(cluster = 1)  |> 
  
  dplyr::left_join(
    expression.glass.exon.metadata |> dplyr::select(
      gene_uid,
      coef.chemotherapy,
      coef.radiotherapy,
      coef.grading,
      coef.resection
    ), by=c('gene_uid'='gene_uid')
  ) |> 
  tidyr::pivot_longer(cols=c(coef.chemotherapy,
                             coef.radiotherapy,
                             coef.grading,
                             coef.resection),
                      values_to="regression_coef",
                      names_to = "covariate"
                      ) |> 
  
  dplyr::mutate(covariate = dplyr::recode(covariate,
                                          "coef.chemotherapy"="Chemotherapy N - Y",
                                          "coef.radiotherapy"="Radiotherapy N - Y",
                                          "coef.grading"="WHO Grade 2 & 3 - 4",
                                          "coef.resection"="Primary - Recurrent"
                                          )) |> 
  dplyr::mutate(covariate = factor(covariate, levels=c(
    "Chemotherapy N - Y", "Radiotherapy N - Y", "WHO Grade 2 & 3 - 4", "Primary - Recurrent"
  )))


plt <- rbind(
  plt,
  plt |> dplyr::mutate(regression_coef = 0)
) 


col2 <- grDevices::colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                      "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                      "#4393C3", "#2166AC", "#053061"))


ggplot(plt, aes(x = -hclust_rank , y = regression_coef , group=gene_name, col=regression_coef)) +
  facet_grid(rows = vars(covariate)) +
  geom_line() +
  theme_bw() +
  labs(y="DESeq2 multi-variate regression coefficient",x="gene") + 
  coord_cartesian(ylim=c(-3.6, 3.6)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() 
  ) +
  ggplot2::scale_color_gradientn(colours = col2(200), na.value = "grey50", limits = c(-1.5, 1.5), 
                                 oob = scales::squish )
  #theme(legend.position = 'bottom', axis.text.x = element_blank())

#saveRDS(plt, "cache/vis_DGE_covar_coefs.Rds")
ggsave("output/figures/vis_DGE_covar_coefs.pdf",width=11,height=5)


