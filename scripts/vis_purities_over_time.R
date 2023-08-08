#!/usr/bin/env R

# load libs ----

library(tidyverse)
source('scripts/R/job_gg_theme.R')
source('scripts/R/youri_gg_theme.R')


# load data ----


source("scripts/load_tumour_purities.R")


# fig test ----


plt <- metadata.glass.per.resection %>% 
  dplyr::select(Sample_Name, resection, methylation.purity.absolute, `dna.shallow.ACE.purity` , `dna.shallow.ACE.purity.below.1`, dna.wes.VAF_IDH ,dna.purity.manual.Erik  ) %>% 
  dplyr::mutate(dna.wes.VAF_IDH = dna.wes.VAF_IDH *2) %>% 
  tidyr::pivot_longer(cols = -c(Sample_Name, resection)) %>% 
  dplyr::rename(Purity.method = name) %>% 
  dplyr::mutate(reliability = ifelse(Purity.method %in% c('methylation.purity.absolute','dna.wes.VAF_IDH'),"most reliable","not most reliable"))


ggplot(plt, aes(x=resection, y=value, fill = resection)) +
  facet_grid(cols = vars(Purity.method),  scales = "free", space="free") +
  geom_violin(alpha=0.1) + 
  ggbeeswarm::geom_beeswarm(pch=21,cex=1.2) +
  theme_cellpress



plt <- metadata.glass.per.resection %>% 
  dplyr::select(Sample_Name, resection, methylation.purity.absolute, `dna.shallow.ACE.purity` , `dna.shallow.ACE.purity.below.1`, dna.wes.VAF_IDH ,dna.purity.manual.Erik  ) %>% 
  dplyr::mutate(dna.wes.VAF_IDH = dna.wes.VAF_IDH *2) %>% 
  tidyr::pivot_longer(cols = -c(Sample_Name, resection)) %>% 
  dplyr::rename(Purity.method = name) %>% 
  dplyr::mutate(reliability = ifelse(Purity.method %in% c('methylation.purity.absolute','dna.wes.VAF_IDH'),"most reliable","not most reliable")) %>% 
  dplyr::filter(Purity.method == "dna.wes.VAF_IDH")


ggplot(plt, aes(x=resection, y=value, fill = resection)) +
  facet_grid(cols = vars(Purity.method),  scales = "free", space="free") +
  geom_violin(alpha=0.1,draw_quantiles = c(0.6), col="black") + 
  ggbeeswarm::geom_beeswarm(pch=21,cex=1.2) +
  youri_gg_theme +
  ylim(0,1.4) +
  ggsignif::geom_signif( comparisons = list(c("S1" , "S2"),c("S2" , "S3"),c("S1" , "S3")), y_position=c(1.125,1.2,1.325),    test="wilcox.test", col="black" )+
  labs(y="fraction tumor cells")





# fig final ----


plt <- metadata.glass.per.resection |> 
  dplyr::select(Sample_Name, 
                GLASS_ID,
                # `dna.shallow.ACE.purity` , `dna.shallow.ACE.purity.below.1`,dna.purity.manual.Erik, # ace is unreliable without manual intervention
                resection, methylation.purity.absolute,
                dna.wes.VAF_IDH ,
                Recurrent_Type) |>  
  tidyr::pivot_longer(cols = -c(Sample_Name, resection, Recurrent_Type, GLASS_ID)) |> 
  dplyr::rename(Purity.method = name) |> 
  dplyr::mutate(facet = dplyr::recode(Purity.method,
                                      `methylation.purity.absolute`='RFpurity Absolute [DNA meth. array]',
                                      `dna.wes.VAF_IDH`='VAF IDH(1/2) [WES DNA]')) |> 
  dplyr::mutate(primary.recurrence = ifelse(is.na(Recurrent_Type), "primary", "recurrence"))


ggplot(plt, aes(x=primary.recurrence, y=value, fill = primary.recurrence)) +
  facet_wrap(~facet,  scales = "free") +
  geom_violin(alpha=0.1) + 
  ggbeeswarm::geom_beeswarm(pch=21,cex=1.2) +
  ggsignif::geom_signif(
    comparisons = list(c("primary" , "recurrence")),  test="t.test", col="black", 
    size=0.5 / 2.14 ,
    textsize = 7 * (3.88/11),
    tip_length = 0
  ) +
  labs(x=NULL, y="Estimated purity or IDH[1/2] VAF") +
  theme_cellpress +
  theme(legend.position = 'bottom', 
        legend.key.height = unit(0, 'cm'))


ggsave("output/figures/vis_purities_over_time.pdf", width=8.5 * 0.9, height=11/ 3) # sizes relative to US-paper

saveRDS(plt, file="cache/plt.vis_purities_over_time.Rds")





