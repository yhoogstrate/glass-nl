#!/usr/bin/env R 


# load libs ----


library(randomForestSRC)# https://cran.r-project.org/web/packages/randomForestSRC/index.html
library(survival)
library(RegParallel)
library(recursiveCorPlot)
library(patchwork)



# load data ----


source('scripts/R/youri_gg_theme.R')


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if(!exists("expression.glass.vst")) {
  source('scripts/load_rna-counts.R')
}


# load per resection metadata ---


metadata <- metadata.glass.per.patient %>% 
  dplyr::left_join(
    metadata.glass.per.resection %>%
      dplyr::select(genomescan.sid, 
                    lts.up1, lts.up2, lts.up3, lts.down , lts.down.a , lts.down.b,
                    methylation.sub.diagnosis, mgmt_status, classification, Epigenetic.subtypes.TCGA, CCND2_stat, CDK4_stat, CDK6_stat, MYC_stat, PTEN_stat, PTCH1_stat, PDGFRA_stat, RB1_stat, TERT_stat, NF2_stat, CDKN2AB
      ) %>% 
      `colnames<-`(paste0(colnames(.),".I")) %>% 
      dplyr::rename(genomescan.sid = genomescan.sid.I)
      ,
    by=c('genomescan.sid.I'='genomescan.sid')) %>% 
  dplyr::left_join(
    metadata.glass.per.resection %>%
      dplyr::select(genomescan.sid, 
                    lts.up1, lts.up2, lts.up3, lts.down , lts.down.a , lts.down.b,
                    methylation.sub.diagnosis, mgmt_status, classification, Epigenetic.subtypes.TCGA, CCND2_stat, CDK4_stat, CDK6_stat, MYC_stat, PTEN_stat, PTCH1_stat, PDGFRA_stat, RB1_stat, TERT_stat, NF2_stat, CDKN2AB
      ) %>% 
      `colnames<-`(paste0(colnames(.),".R")) %>% 
      dplyr::rename(genomescan.sid = genomescan.sid.R),
    by=c('genomescan.sid.R'='genomescan.sid'))


## some test plots ----



ggplot(metadata, aes(x=overall.survival ,y=lts.up2.I )) + 
  geom_point()


p1 <- ggplot(metadata, aes(x=survival.R ,y=lts.up1.R )) + 
  geom_point() + 
  labs(x= "Survival from R2", y="Cell cycling signature @ R2") +
  theme_bw()

p2 <- ggplot(metadata, aes(x=survival.R ,y=lts.up2.R )) + 
  geom_point() + 
  labs(x= "Survival from R2", y="Collagen signature @ R2") +
  theme_bw()

p3 <- ggplot(metadata, aes(x=survival.R ,y=lts.up3.R )) + 
  geom_point() + 
  labs(x= "Survival from R2", y="fuzzy up signature @ R2") +
  theme_bw()

p4 <- ggplot(metadata, aes(x=survival.R ,y=lts.down.R )) + 
  geom_point() + 
  labs(x= "Survival from R2", y="overall down signature @ R2") +
  theme_bw()



(p1 + p2 ) / (p3 + p4)


## paired plot by svvl ----


plt <- metadata %>% 
  dplyr::filter(!is.na( lts.up1.I) & !is.na(survival.R)) %>% 
  dplyr::mutate(order = rank(-survival.R, -survival.I))


plt.expanded <- plt %>% 
  dplyr::select(GLASS_ID, order, survival.I, survival.R, lts.up1.I,  lts.up2.I,  lts.up3.I,  lts.down.I,  lts.up1.R,  lts.up2.R,  lts.up3.R,  lts.down.R) %>% 
  tidyr::pivot_longer(cols = c(lts.up1.I,  lts.up2.I,  lts.up3.I,  lts.down.I, lts.up1.R,  lts.up2.R,  lts.up3.R,  lts.down.R)) %>% 
  dplyr::mutate(facet = gsub("^(.+).[IR]$","\\1",name)) %>% 
  dplyr::mutate(arrow.order = ifelse(gsub("^.+(.)$","\\1",name) == "I",1,2)) %>% 
  dplyr::arrange(GLASS_ID, facet, arrow.order) %>% 
  dplyr::rename(`DGE signature type` = name, `DGE signature contribution` = value)


p1 <- ggplot(plt.expanded, aes(x = reorder(GLASS_ID, order), y = `survival.R` )) +
  geom_segment(aes(x=reorder(GLASS_ID, order), xend=reorder(GLASS_ID, order), y= `survival.R`, yend=0)) +
  geom_point() +
  theme_bw() +
  theme(text = element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25)) +
  labs(y = "Survival time from recurrence")

p2 <- ggplot(plt.expanded, aes(x = order, y = `DGE signature contribution`, group = GLASS_ID )) +
  geom_smooth(method='lm', aes(group=NULL), data = subset(plt.expanded, arrow.order == 2),se=F, col="gray", lwd=1,lty=2) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches"))) +
  facet_grid(rows = vars(facet), scales = "free") +
  theme_bw() +
  theme(text = element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25))




plt.expanded2 <-  plt %>% 
  dplyr::select(GLASS_ID, 
                methylation.sub.diagnosis.I, mgmt_status.I, classification.I, Epigenetic.subtypes.TCGA.I, CCND2_stat.I, CDK4_stat.I, CDK6_stat.I, MYC_stat.I, PTEN_stat.I, PTCH1_stat.I, PDGFRA_stat.I, RB1_stat.I, TERT_stat.I, NF2_stat.I, CDKN2AB.I,
                methylation.sub.diagnosis.R, mgmt_status.R, classification.R, Epigenetic.subtypes.TCGA.R, CCND2_stat.R, CDK4_stat.R, CDK6_stat.R, MYC_stat.R, PTEN_stat.R, PTCH1_stat.R, PDGFRA_stat.R, RB1_stat.R, TERT_stat.R, NF2_stat.R, CDKN2AB.R) %>% 
  tidyr::pivot_longer(cols = c(methylation.sub.diagnosis.I, mgmt_status.I, classification.I, Epigenetic.subtypes.TCGA.I, CCND2_stat.I, CDK4_stat.I, CDK6_stat.I, MYC_stat.I, PTEN_stat.I, PTCH1_stat.I, PDGFRA_stat.I, RB1_stat.I, TERT_stat.I, NF2_stat.I, CDKN2AB.I,
                               methylation.sub.diagnosis.R, mgmt_status.R, classification.R, Epigenetic.subtypes.TCGA.R, CCND2_stat.R, CDK4_stat.R, CDK6_stat.R, MYC_stat.R, PTEN_stat.R, PTCH1_stat.R, PDGFRA_stat.R, RB1_stat.R, TERT_stat.R, NF2_stat.R, CDKN2AB.R)) %>% 
  dplyr::mutate(facet = ifelse(grepl("^.+I$", name), "Initial", "Recurrent")) %>% 
  dplyr::mutate(name = gsub("\\.[IR]$","",name)) %>% 
  dplyr::left_join(plt.expanded %>%  dplyr::select(GLASS_ID, order) %>% dplyr::distinct(), by=c('GLASS_ID'='GLASS_ID'))
  #dplyr::filter(name == "classification")


p3 <- ggplot(plt.expanded2, aes(x=reorder(GLASS_ID, order), y=name, fill=value)) +
  facet_grid(rows = vars(facet), scales = "free") + 
  geom_tile(col="black")



p1 / p2 / p3 + plot_layout(heights = c(1, 3, 2))




## paired plot by sig1 [CCy] ----


rm(plt, plt.expanded, plt.expanded2)

plt <- metadata %>% 
  dplyr::filter(!is.na( lts.up1.I) & !is.na(lts.up1.R)) %>% 
  dplyr::filter(!is.na( survival.I) & !is.na(survival.R)) %>% 
  dplyr::mutate(order = rank(-lts.up1.R, -lts.up1.I))


plt.expanded <- plt %>% 
  dplyr::select(GLASS_ID, order, survival.I, survival.R, lts.up1.I,  lts.up2.I,  lts.up3.I,  lts.down.I,  lts.up1.R,  lts.up2.R,  lts.up3.R,  lts.down.R) %>% 
  tidyr::pivot_longer(cols = c(lts.up1.I,  lts.up2.I,  lts.up3.I,  lts.down.I, lts.up1.R,  lts.up2.R,  lts.up3.R,  lts.down.R)) %>% 
  dplyr::mutate(facet = gsub("^(.+).[IR]$","\\1",name)) %>% 
  dplyr::mutate(arrow.order = ifelse(gsub("^.+(.)$","\\1",name) == "I",1,2)) %>% 
  dplyr::arrange(GLASS_ID, facet, arrow.order) %>% 
  dplyr::rename(`DGE signature type` = name, `DGE signature contribution` = value) %>% 
  dplyr::mutate(facet = factor(facet, levels=c( "lts.up1",  "lts.up2",  "lts.up3", "lts.down")))


p1 <- ggplot(plt.expanded, aes(x = reorder(GLASS_ID, order), y = `survival.R` )) +
  geom_segment(aes(x=reorder(GLASS_ID, order), xend=reorder(GLASS_ID, order), y= `survival.R`, yend=0)) +
  geom_point() +
  theme_bw() +
  theme( axis.text.x  = element_blank()) + # element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25)) +
  labs(y = "Survival time from recurrence")

p2 <- ggplot(plt.expanded, aes(x = reorder(GLASS_ID, order), y = `DGE signature contribution`, group = GLASS_ID )) +
  #geom_smooth(method='lm', aes(group=NULL), data = subset(plt.expanded, arrow.order == 2),se=F, col="gray", lwd=1,lty=2) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches"))) +
  facet_grid(rows = vars(facet)) +
  theme_bw()
theme(text = element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25))




plt.expanded2 <-  plt %>% 
  dplyr::select(GLASS_ID, 
                methylation.sub.diagnosis.I, mgmt_status.I, classification.I, Epigenetic.subtypes.TCGA.I, CCND2_stat.I, CDK4_stat.I, CDK6_stat.I, MYC_stat.I, PTEN_stat.I, PTCH1_stat.I, PDGFRA_stat.I, RB1_stat.I, TERT_stat.I, NF2_stat.I, CDKN2AB.I,
                methylation.sub.diagnosis.R, mgmt_status.R, classification.R, Epigenetic.subtypes.TCGA.R, CCND2_stat.R, CDK4_stat.R, CDK6_stat.R, MYC_stat.R, PTEN_stat.R, PTCH1_stat.R, PDGFRA_stat.R, RB1_stat.R, TERT_stat.R, NF2_stat.R, CDKN2AB.R) %>% 
  tidyr::pivot_longer(cols = c(methylation.sub.diagnosis.I, mgmt_status.I, classification.I, Epigenetic.subtypes.TCGA.I, CCND2_stat.I, CDK4_stat.I, CDK6_stat.I, MYC_stat.I, PTEN_stat.I, PTCH1_stat.I, PDGFRA_stat.I, RB1_stat.I, TERT_stat.I, NF2_stat.I, CDKN2AB.I,
                               methylation.sub.diagnosis.R, mgmt_status.R, classification.R, Epigenetic.subtypes.TCGA.R, CCND2_stat.R, CDK4_stat.R, CDK6_stat.R, MYC_stat.R, PTEN_stat.R, PTCH1_stat.R, PDGFRA_stat.R, RB1_stat.R, TERT_stat.R, NF2_stat.R, CDKN2AB.R)) %>% 
  dplyr::mutate(facet = ifelse(grepl("^.+I$", name), "Initial", "Recurrent")) %>% 
  dplyr::mutate(name = gsub("\\.[IR]$","",name)) %>% 
  dplyr::left_join(plt.expanded %>%  dplyr::select(GLASS_ID, order) %>% dplyr::distinct(), by=c('GLASS_ID'='GLASS_ID')) %>% 
  #dplyr::filter(name %in% c("CDK4_stat","CDK6_stat",  "CCND2_stat", "MYC_stat", "NF2_stat", "PDGFRA_stat", "PTCH1_stat", "PTEN_stat", "RB1_stat", "TERT_stat")) # ok, MYC?
  #dplyr::filter(name == "classification") # ok
  #dplyr::filter(name == "Epigenetic.subtypes.TCGA") # nah?
  #plyr::filter(name == "methylation.sub.diagnosis") # quite well
  #dplyr::filter(name == "mgmt_status") # not at all
  #dplyr::filter(name == "CDKN2AB")
  dplyr::filter(name %in% c("CDKN2AB", "methylation.sub.diagnosis") ) %>%  # no / little 
  dplyr::mutate(value = ifelse(value %in% c("A_IDH", "CONTR_HEMI", "O_IDH", "PLEX_PED_B"), "A_IDH / CONTR_HEMI / O_IDH / PLEX_PED_B", value))


p3 <- ggplot(plt.expanded2, aes(x=reorder(GLASS_ID, order), y=name, fill=value)) +
  facet_grid(rows = vars(facet), scales = "free") + 
  geom_tile(col="black") +
  theme_bw() +
  theme(text = element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25))


p1 / p2 / p3 + plot_layout(heights = c(1, 3, 3)) # quite a bit




## paired plot by sig2 [COL] ----


rm(plt, plt.expanded, plt.expanded2)

plt <- metadata %>% 
  dplyr::filter(!is.na( lts.up1.I) & !is.na(lts.up1.R)) %>% 
  dplyr::filter(!is.na( survival.I) & !is.na(survival.R)) %>% 
  dplyr::mutate(order = rank(-lts.up2.R, -lts.up2.I))


plt.expanded <- plt %>% 
  dplyr::select(GLASS_ID, order, survival.I, survival.R, lts.up1.I,  lts.up2.I,  lts.up3.I,  lts.down.I,  lts.up1.R,  lts.up2.R,  lts.up3.R,  lts.down.R) %>% 
  tidyr::pivot_longer(cols = c(lts.up1.I,  lts.up2.I,  lts.up3.I,  lts.down.I, lts.up1.R,  lts.up2.R,  lts.up3.R,  lts.down.R)) %>% 
  dplyr::mutate(facet = gsub("^(.+).[IR]$","\\1",name)) %>% 
  dplyr::mutate(arrow.order = ifelse(gsub("^.+(.)$","\\1",name) == "I",1,2)) %>% 
  dplyr::arrange(GLASS_ID, facet, arrow.order) %>% 
  dplyr::rename(`DGE signature type` = name, `DGE signature contribution` = value) %>% 
  dplyr::mutate(facet = factor(facet, levels=c( "lts.up1",  "lts.up2",  "lts.up3", "lts.down")))


p1 <- ggplot(plt.expanded, aes(x = reorder(GLASS_ID, order), y = `survival.R` )) +
  geom_segment(aes(x=reorder(GLASS_ID, order), xend=reorder(GLASS_ID, order), y= `survival.R`, yend=0)) +
  geom_point() +
  theme_bw() +
  theme( axis.text.x  = element_blank()) + # element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25)) +
  labs(y = "Survival time from recurrence")

p2 <- ggplot(plt.expanded, aes(x = reorder(GLASS_ID, order), y = `DGE signature contribution`, group = GLASS_ID )) +
  #geom_smooth(method='lm', aes(group=NULL), data = subset(plt.expanded, arrow.order == 2),se=F, col="gray", lwd=1,lty=2) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches"))) +
  facet_grid(rows = vars(facet)) +
  theme_bw()
theme(text = element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25))




plt.expanded2 <-  plt %>% 
  dplyr::select(GLASS_ID, 
                methylation.sub.diagnosis.I, mgmt_status.I, classification.I, Epigenetic.subtypes.TCGA.I, CCND2_stat.I, CDK4_stat.I, CDK6_stat.I, MYC_stat.I, PTEN_stat.I, PTCH1_stat.I, PDGFRA_stat.I, RB1_stat.I, TERT_stat.I, NF2_stat.I, CDKN2AB.I,
                methylation.sub.diagnosis.R, mgmt_status.R, classification.R, Epigenetic.subtypes.TCGA.R, CCND2_stat.R, CDK4_stat.R, CDK6_stat.R, MYC_stat.R, PTEN_stat.R, PTCH1_stat.R, PDGFRA_stat.R, RB1_stat.R, TERT_stat.R, NF2_stat.R, CDKN2AB.R) %>% 
  tidyr::pivot_longer(cols = c(methylation.sub.diagnosis.I, mgmt_status.I, classification.I, Epigenetic.subtypes.TCGA.I, CCND2_stat.I, CDK4_stat.I, CDK6_stat.I, MYC_stat.I, PTEN_stat.I, PTCH1_stat.I, PDGFRA_stat.I, RB1_stat.I, TERT_stat.I, NF2_stat.I, CDKN2AB.I,
                               methylation.sub.diagnosis.R, mgmt_status.R, classification.R, Epigenetic.subtypes.TCGA.R, CCND2_stat.R, CDK4_stat.R, CDK6_stat.R, MYC_stat.R, PTEN_stat.R, PTCH1_stat.R, PDGFRA_stat.R, RB1_stat.R, TERT_stat.R, NF2_stat.R, CDKN2AB.R)) %>% 
  dplyr::mutate(facet = ifelse(grepl("^.+I$", name), "Initial", "Recurrent")) %>% 
  dplyr::mutate(name = gsub("\\.[IR]$","",name)) %>% 
  dplyr::left_join(plt.expanded %>%  dplyr::select(GLASS_ID, order) %>% dplyr::distinct(), by=c('GLASS_ID'='GLASS_ID')) %>% 
  #dplyr::filter(name %in% c("CDK4_stat","CDK6_stat",  "CCND2_stat", "MYC_stat", "NF2_stat", "PDGFRA_stat", "PTCH1_stat", "PTEN_stat", "RB1_stat", "TERT_stat")) # nah
  #dplyr::filter(name == "classification") # nah
  #dplyr::filter(name == "Epigenetic.subtypes.TCGA") # nah
  #dplyr::filter(name == "methylation.sub.diagnosis") # quite well
  #dplyr::filter(name == "mgmt_status") # nah
  #dplyr::filter(name == "CDKN2AB") # nah
  dplyr::filter(name %in% c("CDKN2AB", "methylation.sub.diagnosis") ) %>%  # no / little 
  dplyr::mutate(value = ifelse(value %in% c("A_IDH", "CONTR_HEMI", "O_IDH", "PLEX_PED_B"), "A_IDH / CONTR_HEMI / O_IDH / PLEX_PED_B", value))


p3 <- ggplot(plt.expanded2, aes(x=reorder(GLASS_ID, order), y=name, fill=value)) +
  facet_grid(rows = vars(facet), scales = "free") + 
  geom_tile(col="black") +
  theme_bw() +
  theme(text = element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25))


p1 / p2 / p3 + plot_layout(heights = c(1, 3, 3))




## paired plot by sig3 [COL] ----


rm(plt, plt.expanded, plt.expanded2)

plt <- metadata %>% 
  dplyr::filter(!is.na( lts.up3.I) & !is.na(lts.up3.R)) %>% 
  dplyr::filter(!is.na( survival.I) & !is.na(survival.R)) %>% 
  dplyr::mutate(order = rank(-lts.up3.R, -lts.up3.I))


plt.expanded <- plt %>% 
  dplyr::select(GLASS_ID, order, survival.I, survival.R, lts.up1.I,  lts.up2.I,  lts.up3.I,  lts.down.I,  lts.up1.R,  lts.up2.R,  lts.up3.R,  lts.down.R) %>% 
  tidyr::pivot_longer(cols = c(lts.up1.I,  lts.up2.I,  lts.up3.I,  lts.down.I, lts.up1.R,  lts.up2.R,  lts.up3.R,  lts.down.R)) %>% 
  dplyr::mutate(facet = gsub("^(.+).[IR]$","\\1",name)) %>% 
  dplyr::mutate(arrow.order = ifelse(gsub("^.+(.)$","\\1",name) == "I",1,2)) %>% 
  dplyr::arrange(GLASS_ID, facet, arrow.order) %>% 
  dplyr::rename(`DGE signature type` = name, `DGE signature contribution` = value) %>% 
  dplyr::mutate(facet = factor(facet, levels=c( "lts.up1",  "lts.up2",  "lts.up3", "lts.down")))


p1 <- ggplot(plt.expanded, aes(x = reorder(GLASS_ID, order), y = `survival.R` )) +
  geom_segment(aes(x=reorder(GLASS_ID, order), xend=reorder(GLASS_ID, order), y= `survival.R`, yend=0)) +
  geom_point() +
  theme_bw() +
  theme( axis.text.x  = element_blank()) + # element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25)) +
  labs(y = "Survival time from recurrence")

p2 <- ggplot(plt.expanded, aes(x = reorder(GLASS_ID, order), y = `DGE signature contribution`, group = GLASS_ID )) +
  #geom_smooth(method='lm', aes(group=NULL), data = subset(plt.expanded, arrow.order == 2),se=F, col="gray", lwd=1,lty=2) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches"))) +
  facet_grid(rows = vars(facet)) +
  theme_bw()
theme(text = element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25))




plt.expanded2 <-  plt %>% 
  dplyr::select(GLASS_ID, 
                methylation.sub.diagnosis.I, mgmt_status.I, classification.I, Epigenetic.subtypes.TCGA.I, CCND2_stat.I, CDK4_stat.I, CDK6_stat.I, MYC_stat.I, PTEN_stat.I, PTCH1_stat.I, PDGFRA_stat.I, RB1_stat.I, TERT_stat.I, NF2_stat.I, CDKN2AB.I,
                methylation.sub.diagnosis.R, mgmt_status.R, classification.R, Epigenetic.subtypes.TCGA.R, CCND2_stat.R, CDK4_stat.R, CDK6_stat.R, MYC_stat.R, PTEN_stat.R, PTCH1_stat.R, PDGFRA_stat.R, RB1_stat.R, TERT_stat.R, NF2_stat.R, CDKN2AB.R) %>% 
  tidyr::pivot_longer(cols = c(methylation.sub.diagnosis.I, mgmt_status.I, classification.I, Epigenetic.subtypes.TCGA.I, CCND2_stat.I, CDK4_stat.I, CDK6_stat.I, MYC_stat.I, PTEN_stat.I, PTCH1_stat.I, PDGFRA_stat.I, RB1_stat.I, TERT_stat.I, NF2_stat.I, CDKN2AB.I,
                               methylation.sub.diagnosis.R, mgmt_status.R, classification.R, Epigenetic.subtypes.TCGA.R, CCND2_stat.R, CDK4_stat.R, CDK6_stat.R, MYC_stat.R, PTEN_stat.R, PTCH1_stat.R, PDGFRA_stat.R, RB1_stat.R, TERT_stat.R, NF2_stat.R, CDKN2AB.R)) %>% 
  dplyr::mutate(facet = ifelse(grepl("^.+I$", name), "Initial", "Recurrent")) %>% 
  dplyr::mutate(name = gsub("\\.[IR]$","",name)) %>% 
  dplyr::left_join(plt.expanded %>%  dplyr::select(GLASS_ID, order) %>% dplyr::distinct(), by=c('GLASS_ID'='GLASS_ID')) %>% 
  #dplyr::filter(name %in% c("CDK4_stat","CDK6_stat",  "CCND2_stat", "MYC_stat", "NF2_stat", "PDGFRA_stat", "PTCH1_stat", "PTEN_stat", "RB1_stat", "TERT_stat")) # no
  #dplyr::filter(name == "classification") # barely/nah
  #dplyr::filter(name == "Epigenetic.subtypes.TCGA") # no
  #dplyr::filter(name == "methylation.sub.diagnosis") # a little
  #dplyr::filter(name == "mgmt_status") # no
  #dplyr::filter(name == "CDKN2AB") # a little
  dplyr::filter(name %in% c("CDKN2AB", "methylation.sub.diagnosis") ) %>%  # no / little 
  dplyr::mutate(value = ifelse(value %in% c("A_IDH", "CONTR_HEMI", "O_IDH", "PLEX_PED_B"), "A_IDH / CONTR_HEMI / O_IDH / PLEX_PED_B", value))



p3 <- ggplot(plt.expanded2, aes(x=reorder(GLASS_ID, order), y=name, fill=value)) +
  facet_grid(rows = vars(facet), scales = "free") + 
  geom_tile(col="black") +
  theme_bw() +
  theme(text = element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25))


p1 / p2 / p3 + plot_layout(heights = c(1, 3, 3))




## paired plot by sig4/down [COL] ----


rm(plt, plt.expanded, plt.expanded2)

plt <- metadata %>% 
  dplyr::filter(!is.na( lts.down.I) & !is.na(lts.down.R)) %>% 
  dplyr::filter(!is.na( survival.I) & !is.na(survival.R)) %>% 
  dplyr::mutate(order = rank(-lts.down.R, -lts.down.I))


plt.expanded <- plt %>% 
  dplyr::select(GLASS_ID, order, survival.I, survival.R, lts.up1.I,  lts.up2.I,  lts.up3.I,  lts.down.I,  lts.up1.R,  lts.up2.R,  lts.up3.R,  lts.down.R) %>% 
  tidyr::pivot_longer(cols = c(lts.up1.I,  lts.up2.I,  lts.up3.I,  lts.down.I, lts.up1.R,  lts.up2.R,  lts.up3.R,  lts.down.R)) %>% 
  dplyr::mutate(facet = gsub("^(.+).[IR]$","\\1",name)) %>% 
  dplyr::mutate(arrow.order = ifelse(gsub("^.+(.)$","\\1",name) == "I",1,2)) %>% 
  dplyr::arrange(GLASS_ID, facet, arrow.order) %>% 
  dplyr::rename(`DGE signature type` = name, `DGE signature contribution` = value) %>% 
  dplyr::mutate(facet = factor(facet, levels=c( "lts.up1",  "lts.up2",  "lts.up3", "lts.down")))


p1 <- ggplot(plt.expanded, aes(x = reorder(GLASS_ID, order), y = `survival.R` )) +
  geom_segment(aes(x=reorder(GLASS_ID, order), xend=reorder(GLASS_ID, order), y= `survival.R`, yend=0)) +
  geom_point() +
  theme_bw() +
  theme( axis.text.x  = element_blank()) + # element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25)) +
  labs(y = "Survival time from recurrence")

p2 <- ggplot(plt.expanded, aes(x = reorder(GLASS_ID, order), y = `DGE signature contribution`, group = GLASS_ID )) +
  #geom_smooth(method='lm', aes(group=NULL), data = subset(plt.expanded, arrow.order == 2),se=F, col="gray", lwd=1,lty=2) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches"))) +
  facet_grid(rows = vars(facet)) +
  theme_bw() +
  theme(text = element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25))




plt.expanded2 <-  plt %>% 
  dplyr::select(GLASS_ID, 
                methylation.sub.diagnosis.I, mgmt_status.I, classification.I, Epigenetic.subtypes.TCGA.I, CCND2_stat.I, CDK4_stat.I, CDK6_stat.I, MYC_stat.I, PTEN_stat.I, PTCH1_stat.I, PDGFRA_stat.I, RB1_stat.I, TERT_stat.I, NF2_stat.I, CDKN2AB.I,
                methylation.sub.diagnosis.R, mgmt_status.R, classification.R, Epigenetic.subtypes.TCGA.R, CCND2_stat.R, CDK4_stat.R, CDK6_stat.R, MYC_stat.R, PTEN_stat.R, PTCH1_stat.R, PDGFRA_stat.R, RB1_stat.R, TERT_stat.R, NF2_stat.R, CDKN2AB.R) %>% 
  tidyr::pivot_longer(cols = c(methylation.sub.diagnosis.I, mgmt_status.I, classification.I, Epigenetic.subtypes.TCGA.I, CCND2_stat.I, CDK4_stat.I, CDK6_stat.I, MYC_stat.I, PTEN_stat.I, PTCH1_stat.I, PDGFRA_stat.I, RB1_stat.I, TERT_stat.I, NF2_stat.I, CDKN2AB.I,
                               methylation.sub.diagnosis.R, mgmt_status.R, classification.R, Epigenetic.subtypes.TCGA.R, CCND2_stat.R, CDK4_stat.R, CDK6_stat.R, MYC_stat.R, PTEN_stat.R, PTCH1_stat.R, PDGFRA_stat.R, RB1_stat.R, TERT_stat.R, NF2_stat.R, CDKN2AB.R)) %>% 
  dplyr::mutate(facet = ifelse(grepl("^.+I$", name), "Initial", "Recurrent")) %>% 
  dplyr::mutate(name = gsub("\\.[IR]$","",name)) %>% 
  dplyr::left_join(plt.expanded %>%  dplyr::select(GLASS_ID, order) %>% dplyr::distinct(), by=c('GLASS_ID'='GLASS_ID')) %>% 
  #dplyr::filter(name %in% c("CDK4_stat","CDK6_stat",  "CCND2_stat", "MYC_stat", "NF2_stat", "PDGFRA_stat", "PTCH1_stat", "PTEN_stat", "RB1_stat", "TERT_stat")) # meh, MYC?
  #dplyr::filter(name == "classification") # good
  #dplyr::filter(name == "Epigenetic.subtypes.TCGA") # no?
  #dplyr::filter(name == "methylation.sub.diagnosis") #  good
  # dplyr::filter(name == "mgmt_status") # no
  #dplyr::filter(name == "CDKN2AB") # no / little
  dplyr::filter(name %in% c("CDKN2AB", "methylation.sub.diagnosis") ) %>%  # no / little 
  dplyr::mutate(value = ifelse(value %in% c("A_IDH", "CONTR_HEMI", "O_IDH", "PLEX_PED_B"), "A_IDH / CONTR_HEMI / O_IDH / PLEX_PED_B", value))


p3 <- ggplot(plt.expanded2, aes(x=reorder(GLASS_ID, order), y=name, fill=value)) +
  facet_grid(rows = vars(facet), scales = "free") + 
  geom_tile(col="black") +
  theme_bw() +
  theme(text = element_text(family = 'Helvetica'), axis.text.x = element_text(angle = 90, hjust = 0.25))


p1 / p2 / p3 + plot_layout(heights = c(1, 3, 3))




plot(metadata.glass.per.resection$lts.down, log(metadata.glass.per.resection$A_IDH_HG_cal/metadata.glass.per.resection$A_IDH_cal),xlab="RNA Signature 4 (down) expression",ylab="Heidelberg Classifier log(IDH_LGG_HG score / IDH_LGG score)")

p1 <- ggplot(metadata.glass.per.resection, aes(x = lts.up1, y=log(A_IDH_HG_cal/A_IDH_cal), col=resection)) +
  geom_point() +
  theme_bw() +
  labs(x="RNA Signature 1 / Cell cycling (up) expression",y="Heidelberg Classifier log(IDH_LGG_HG score / IDH_LGG score)")

p2 <- ggplot(metadata.glass.per.resection, aes(x = lts.up2, y=log(A_IDH_HG_cal/A_IDH_cal), col=resection)) +
  geom_point() +
  theme_bw() +
  labs(x="RNA Signature 2 / Collagen (up) expression",y="Heidelberg Classifier log(IDH_LGG_HG score / IDH_LGG score)")

p3 <- ggplot(metadata.glass.per.resection, aes(x = lts.up3, y=log(A_IDH_HG_cal/A_IDH_cal), col=resection)) +
  geom_point() +
  theme_bw() +
  labs(x="RNA Signature 3 / Fuzzy (up) expression",y="Heidelberg Classifier log(IDH_LGG_HG score / IDH_LGG score)")

p4 <- ggplot(metadata.glass.per.resection, aes(x = lts.down, y=log(A_IDH_HG_cal/A_IDH_cal), col=resection)) +
  geom_point() +
  theme_bw() +
  labs(x="RNA Signature 4 / non-malignant astrocytes? (down) expression",y="Heidelberg Classifier log(IDH_LGG_HG score / IDH_LGG score)")






### another plot ----


data = plt.expanded  %>% 
  dplyr::rename(x = order) %>% 
  dplyr::rename(y= `DGE signature contribution`) %>% 
  dplyr::filter(`DGE signature type` == "lts.down.I")

ggplot(data,aes(x=x, y=y, group=GLASS_ID)) +
  geom_smooth(aes(group = NULL),method='lm',se=FALSE) + #, formula= y~x) +
  geom_point()



# example RF model ----

# 
# data(pbc, package = "randomForestSRC")
# 
# 
# pbc.obj <- rfsrc(Surv(days,status) ~ ., pbc, importance = TRUE)
# #find.interaction(pbc.obj, method = "vimp", nvar = 8)
# 
# plot(pbc.obj)
# 

# subset data ----



tmp.metadata <- metadata.glass.per.patient %>%
  tidyr::drop_na(genomescan.sid.R) %>%
  dplyr::select(genomescan.sid.R) %>% 
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid.R'='genomescan.sid')) %>%
  tidyr::drop_na(time.resection.until.last.event) %>% 
  dplyr::rename(genomescan.sid = genomescan.sid.R) %>% 
  dplyr::arrange(Sample_Type, genomescan.sid)





tmp.data.exon <- expression.glass.exon.vst %>%
  dplyr::select(all_of( tmp.metadata$genomescan.sid )) %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad)

#plot(sort(tmp.data$mad, decreasing = T))
#abline(h=1)
k.exon <- sum(tmp.data.exon$mad > 1)

tmp.data.exon <- tmp.data.exon %>% 
  dplyr::mutate(mad = NULL)

stopifnot(colnames(tmp.data.exon) == tmp.metadata$genomescan.sid)





tmp.data.gene <- expression.glass.gene.vst %>%
  dplyr::select(all_of( tmp.metadata$genomescan.sid )) %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad)

# plot(sort(tmp.data.gene$mad, decreasing = T))
# abline(h=1)
k.gene <- sum(tmp.data.gene$mad > 1)

tmp.data.gene <- tmp.data.gene %>% 
  dplyr::mutate(mad = NULL)

stopifnot(colnames(tmp.data.gene) == tmp.metadata$genomescan.sid)





# R2 svvl ----



svvl.R2.exon <- tmp.metadata %>%
  dplyr::rename(svvl = time.resection.until.last.event) %>% 
  dplyr::rename(status = status.resection.until.last.event) %>%
  dplyr::select(genomescan.sid, svvl, status) %>% 
  dplyr::left_join(
    tmp.data.exon %>%
      dplyr::slice_head(n=k.exon) %>% 
      t %>%
      as.data.frame(stringsAsFactors=F) %>% 
      tibble::rownames_to_column('genomescan.sid'),
    by=c('genomescan.sid'='genomescan.sid')
  ) %>% 
  tibble::column_to_rownames('genomescan.sid') %>% 
  dplyr::mutate(svvl = as.numeric(svvl)) %>% 
  `colnames<-`(gsub('-','.',colnames(.),fixed=T))




svvl.R2.gene <- tmp.metadata %>%
  dplyr::rename(svvl = time.resection.until.last.event) %>% 
  dplyr::rename(status = status.resection.until.last.event) %>%
  dplyr::select(genomescan.sid, svvl, status) %>% 
  dplyr::left_join(
    tmp.data.gene %>%
      dplyr::slice_head(n=k.gene) %>% 
      t %>%
      as.data.frame(stringsAsFactors=F) %>% 
      tibble::rownames_to_column('genomescan.sid'),
    by=c('genomescan.sid'='genomescan.sid')
  ) %>% 
  tibble::column_to_rownames('genomescan.sid') %>% 
  dplyr::mutate(svvl = as.numeric(svvl)) %>% 
  `colnames<-`(gsub('-','.',colnames(.),fixed=T))



## RF R2 ----


pbc.R2.exon <- rfsrc(Surv(svvl,status) ~ ., svvl.R2.exon , importance = TRUE,ntree=4000,seed=1234)
plot.rfsrc(pbc.R2.exon)
# plot.variable.rfsrc(pbc.obj)


plt <- data.frame(importance = pbc.R2.exon$importance) %>% 
  dplyr::arrange(abs(importance)) %>% 
  dplyr::top_n(20) %>% 
  dplyr::arrange(-importance) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::mutate(baseline = 0) %>% 
  dplyr::mutate(y = rank(importance)) %>% 
  tidyr::pivot_longer(cols = -c(gene_uid, y)) %>% 
  dplyr::mutate(name=NULL) %>% 
  dplyr::rename(RFSRC.importance.R2.exon = value)


p.exon <- ggplot(plt, aes(x=RFSRC.importance.R2.exon,y=reorder(gene_uid, y), group=gene_uid)) +
  geom_line(lwd=2) +
  youri_gg_theme +
  labs(y=NULL)


rm(plt)




pbc.R2.gene <- rfsrc(Surv(svvl,status) ~ ., svvl.R2.gene , importance = TRUE,ntree=4000,seed=1234)
plot.rfsrc(pbc.R2.gene)
# plot.variable.rfsrc(pbc.obj)


plt <- data.frame(importance = pbc.R2.gene$importance) %>% 
  dplyr::arrange(abs(importance)) %>% 
  dplyr::top_n(20) %>% 
  dplyr::arrange(-importance) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::mutate(baseline = 0) %>% 
  dplyr::mutate(y = rank(importance)) %>% 
  tidyr::pivot_longer(cols = -c(gene_uid, y)) %>% 
  dplyr::mutate(name=NULL) %>% 
  dplyr::rename(RFSRC.importance.R2.gene = value)


p.gene <- ggplot(plt, aes(x=RFSRC.importance.R2.gene,y=reorder(gene_uid, y), group=gene_uid)) +
  geom_line(lwd=2) +
  youri_gg_theme +
  labs(y=NULL)


rm(plt)



p.exon + p.gene


### integrated rank ----



plt <- dplyr::left_join(
  data.frame(importance.exon = pbc.R2.exon$importance) %>% 
    tibble::rownames_to_column('gene_uid'),
  data.frame(importance.gene = pbc.R2.gene$importance) %>% 
    tibble::rownames_to_column('gene_uid'),
  by=c('gene_uid'='gene_uid')
) %>% 
  dplyr::mutate(symbol = gsub("^.+_","",gene_uid)) %>% 
  dplyr::mutate(show.label = importance.exon >= 0.0018)
  

ggplot(plt, aes(x=importance.exon, y=importance.gene, label=symbol)) +
  geom_point() +
  ggrepel::geom_text_repel(data = subset(plt, show.label == T))



### example co-expression top features ----
# ENSG00000232093_DCST1.AS1


plt <- data.frame(pbc.obj.importance = pbc.R2.exon$importance) %>% 
  dplyr::arrange(-abs(pbc.obj.importance)) %>% 
  dplyr::top_n(55) %>% 
  tibble::rownames_to_column("gene_uid") %>% 
  dplyr::mutate(gene_tmp_uid = gsub("_.+$","",gene_uid)) %>% 
  dplyr::mutate(gene_uid = NULL) %>% 
  dplyr::left_join(
    expression.glass.exon.vst %>% 
      tibble::rownames_to_column("gene_uid") %>% 
      dplyr::mutate(gene_tmp_uid = gsub("_.+$","",gene_uid))
    ,
    by=c('gene_tmp_uid'='gene_tmp_uid')
  ) %>% 
  dplyr::mutate(gene_tmp_uid = NULL) %>% 
  tibble::column_to_rownames('gene_uid') 


labels <- data.frame(gid = rownames(plt),'a'=T) %>% 
  tibble::column_to_rownames('gid')


recursiveCorPlot(plt, labels, 7, 7)




plt <- data.frame(pbc.obj.importance = pbc.R2.gene$importance) %>% 
  dplyr::arrange(-abs(pbc.obj.importance)) %>% 
  dplyr::top_n(55) %>% 
  tibble::rownames_to_column("gene_uid") %>% 
  dplyr::mutate(gene_tmp_uid = gsub("_.+$","",gene_uid)) %>% 
  dplyr::mutate(gene_uid = NULL) %>% 
  dplyr::left_join(
    expression.glass.gene.vst %>% 
      tibble::rownames_to_column("gene_uid") %>% 
      dplyr::mutate(gene_tmp_uid = gsub("_.+$","",gene_uid))
    ,
    by=c('gene_tmp_uid'='gene_tmp_uid')
  ) %>% 
  dplyr::mutate(gene_tmp_uid = NULL) %>% 
  tibble::column_to_rownames('gene_uid') 


labels <- data.frame(gid = rownames(plt),'a'=T) %>% 
  tibble::column_to_rownames('gid')


recursiveCorPlot(plt, labels, 7, 7)





## Coxph ----


coxph.exon <- RegParallel(
  data = svvl.R2.exon ,
  formula = 'Surv(svvl, status) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(svvl.R2.exon)[3:ncol(svvl.R2.exon)],
  blocksize = 65,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)

a = coxph.exon %>% dplyr::arrange(LogRank) %>% dplyr::select(Variable, LogRank) %>% 
  dplyr::rename(gene_uid = Variable) %>% 
  dplyr::mutate(pbc.obj.importance = -LogRank)

plt <- a %>% 
  dplyr::arrange(-abs(pbc.obj.importance)) %>% 
  dplyr::top_n(55) %>% 
  dplyr::mutate(gene_tmp_uid = gsub("_.+$","",gene_uid)) %>% 
  dplyr::mutate(gene_uid = NULL) %>% 
  dplyr::left_join(
    expression.glass.exon.vst %>% 
      tibble::rownames_to_column("gene_uid") %>% 
      dplyr::mutate(gene_tmp_uid = gsub("_.+$","",gene_uid))
    ,
    by=c('gene_tmp_uid'='gene_tmp_uid')
  ) %>% 
  dplyr::mutate(gene_tmp_uid = NULL) %>% 
  tibble::column_to_rownames('gene_uid') 


labels <- data.frame(gid = rownames(plt),'a'=T) %>% 
  tibble::column_to_rownames('gid')


recursiveCorPlot(plt, labels, 7, 7)






coxph.gene <- RegParallel(
  data = svvl.R2.gene ,
  formula = 'Surv(svvl, status) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(svvl.R2.gene)[3:ncol(svvl.R2.gene)],
  blocksize = 65,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)

coxph.gene %>% dplyr::arrange(LogRank) %>% head(n=55) %>%  dplyr::pull(Variable)





plt2 <- data.frame(pbc.obj$importance) %>%
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::left_join(
  coxph.res, by=c('gene_uid'='Variable')
)


# Z
plot(log(plt2$pbc.obj.importance), abs(plt2$Z),pch=19,cex=0.25)
plot(log10(plt2$pbc.obj.importance- min(plt2$pbc.obj.importance) + 0.00001), abs(plt2$Z),pch=19,cex=0.25)


# HR
plot(log(plt2$pbc.obj.importance), abs(log(plt2$HR)),pch=19,cex=0.25)

# LR
plot(log(plt2$pbc.obj.importance), -log10(plt2$LogRank),pch=19,cex=0.25)

# W
plot(log(plt2$pbc.obj.importance), -log10(plt2$Wald),pch=19,cex=0.25)


plt2 <- plt2 %>% 
  #dplyr::mutate(log.importance = log(pbc.obj.importance)) %>% 
  dplyr::mutate(log.importance = log10(plt2$pbc.obj.importance- min(plt2$pbc.obj.importance) + 0.00001)) %>% 
  dplyr::mutate(abs.Z = abs(Z)) %>% 
  dplyr::mutate(abs.log.HR = abs(log(HR))) %>% 
  dplyr::mutate(m.10log.LR = -log10(LogRank)) %>% 
  dplyr::mutate(m.10log.W = -log10(Wald)) %>% 
  dplyr::mutate(vis = log.importance > -2.4 | log.importance < -4 | m.10log.LR > 6) %>% 
  dplyr::mutate(label = gsub("^[^_]+_","",gene_uid))



ggplot(plt2, aes(x=log.importance, y=abs.Z, label=label)) +
  geom_point(pch=21) +
  ggrepel::geom_text_repel(data = plt2 %>%  dplyr::filter(vis == T)) +
  youri_gg_theme

ggplot(plt2, aes(x=log.importance, y=m.10log.LR, label=label)) +
  geom_point(pch=21) +
  ggrepel::geom_text_repel(data = plt2 %>%  dplyr::filter(vis == T)) +
  youri_gg_theme

ggplot(plt2, aes(x=log.importance, y=m.10log.LR, label=label)) +
  geom_point(pch=21) +
  ggrepel::geom_text_repel(data = plt2 %>%  dplyr::filter(vis == T)) +
  youri_gg_theme

ggplot(plt2, aes(x=log.importance, y=m.10log.W, label=label)) +
  geom_point(pch=21) +
  ggrepel::geom_text_repel(data = plt2 %>%  dplyr::filter(vis == T)) +
  youri_gg_theme






# R1 svvl ----


tmp.metadata <- metadata.glass.per.patient %>%
  tidyr::drop_na(genomescan.sid.I) %>%
  dplyr::select(genomescan.sid.I) %>% 
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid.I'='genomescan.sid')) %>%
  tidyr::drop_na(time.resection.until.last.event) %>% 
  dplyr::rename(genomescan.sid = genomescan.sid.I) %>% 
  dplyr::arrange(Sample_Type, genomescan.sid)



tmp.data <- expression.glass.exon.vst %>%
  dplyr::select(all_of( tmp.metadata$genomescan.sid )) %>% 
  dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>% 
  dplyr::arrange(-mad) %>% 
  dplyr::mutate(mad=NULL) 



stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)



svvl.obj <- tmp.metadata %>%
  dplyr::rename(svvl = time.resection.until.last.event) %>% 
  dplyr::rename(status = status.resection.until.last.event) %>%
  dplyr::select(genomescan.sid, svvl, status) %>% 
  dplyr::left_join(
    tmp.data %>%
      dplyr::slice_head(n=10000) %>% 
      t %>%
      as.data.frame(stringsAsFactors=F) %>% 
      tibble::rownames_to_column('genomescan.sid'),
    by=c('genomescan.sid'='genomescan.sid')
  ) %>% 
  tibble::column_to_rownames('genomescan.sid') %>% 
  dplyr::mutate(svvl = as.numeric(svvl)) %>% 
  `colnames<-`(gsub('-','.',colnames(.),fixed=T))





## RF (RFSRC) ----


pbc.obj <- rfsrc(Surv(svvl,status) ~ ., svvl.obj , importance = TRUE,ntree=4000,seed=1234)
plot.rfsrc(pbc.obj)
# plot.variable.rfsrc(pbc.obj)


plt <- data.frame(pbc.obj$importance) %>% 
  dplyr::arrange(abs(pbc.obj.importance)) %>% 
  dplyr::top_n(20) %>% 
  dplyr::arrange(-pbc.obj.importance) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::mutate(baseline = 0) %>% 
  dplyr::mutate(y = rank(pbc.obj.importance)) %>% 
  tidyr::pivot_longer(cols = -c(gene_uid, y)) %>% 
  dplyr::mutate(name=NULL) %>% 
  dplyr::rename(RFSRC.importance = value)

ggplot(plt, aes(x=RFSRC.importance,y=reorder(gene_uid, y), group=gene_uid)) +
  geom_line(lwd=2) +
  youri_gg_theme +
  labs(y=NULL)


## Coxph ----


coxph.res <- RegParallel(
  data = svvl.obj ,
  formula = 'Surv(svvl, status) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(svvl.obj)[3:ncol(svvl.obj)],
  blocksize = 65,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)



plt2 <- data.frame(pbc.obj$importance) %>%
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::left_join(
    coxph.res, by=c('gene_uid'='Variable')
  )


# Z
plot(log(plt2$pbc.obj.importance), abs(plt2$Z),pch=19,cex=0.25)
plot(log10(plt2$pbc.obj.importance- min(plt2$pbc.obj.importance) + 0.00001), abs(plt2$Z),pch=19,cex=0.25)


# HR
plot(log(plt2$pbc.obj.importance), abs(log(plt2$HR)),pch=19,cex=0.25)

# LR
plot(log(plt2$pbc.obj.importance), -log10(plt2$LogRank),pch=19,cex=0.25)

# W
plot(log(plt2$pbc.obj.importance), -log10(plt2$Wald),pch=19,cex=0.25)


plt2 <- plt2 %>% 
  #dplyr::mutate(log.importance = log(pbc.obj.importance)) %>% 
  dplyr::mutate(log.importance = log10(plt2$pbc.obj.importance- min(plt2$pbc.obj.importance) + 0.00001)) %>% 
  dplyr::mutate(abs.Z = abs(Z)) %>% 
  dplyr::mutate(abs.log.HR = abs(log(HR))) %>% 
  dplyr::mutate(m.10log.LR = -log10(LogRank)) %>% 
  dplyr::mutate(m.10log.W = -log10(Wald)) %>% 
  dplyr::mutate(vis = log.importance > -2.6 | log.importance < -4 | m.10log.LR > 6) %>% 
  dplyr::mutate(label = gsub("^[^_]+_","",gene_uid))



ggplot(plt2, aes(x=log.importance, y=abs.Z, label=label)) +
  geom_point(pch=21) +
  ggrepel::geom_text_repel(data = plt2 %>%  dplyr::filter(vis == T)) +
  xlim(-3.4,-2.4) +
  youri_gg_theme

ggplot(plt2, aes(x=log.importance, y=m.10log.LR, label=label)) +
  geom_point(pch=21) +
  ggrepel::geom_text_repel(data = plt2 %>%  dplyr::filter(vis == T)) +
  xlim(-3.4,-2.4) +
  youri_gg_theme

ggplot(plt2, aes(x=log.importance, y=m.10log.LR, label=label)) +
  geom_point(pch=21) +
  ggrepel::geom_text_repel(data = plt2 %>%  dplyr::filter(vis == T)) +
  xlim(-3.4,-2.4) +
  youri_gg_theme

ggplot(plt2, aes(x=log.importance, y=m.10log.W, label=label)) +
  geom_point(pch=21) +
  ggrepel::geom_text_repel(data = plt2 %>%  dplyr::filter(vis == T)) +
  xlim(-3.4,-2.4) +
  youri_gg_theme



# cleanup ----


rm(tmp.data, tmp.metadatam, k)



