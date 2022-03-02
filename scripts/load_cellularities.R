#!/usr/bin/env R


library(readxl)
library(tidyverse)
source('scripts/R/job_gg_theme.R')
source('scripts/R/youri_gg_theme.R')



#cnv.cellularities <- read.delim('output/tables/dna-seq/CellularitiesManuallyCurated.xlsx.txt') %>% 
#  dplyr::mutate(names = gsub("-","_",names,fixed=T)) %>% 
#  dplyr::mutate(Sample_ID = gsub("^(.+_.+)_I.+$","\\1",names)) %>% 
#  dplyr::mutate(names = NULL)


glass.cellularities <- readxl::read_excel('data/glass/WES/CombinedDataGLASS/AllDataGLASS.xlsx') %>% 
  dplyr::mutate(samplenames = NULL) %>% 
  dplyr::mutate(names = gsub("-","_",names,fixed=T)) %>% 
  dplyr::mutate(Sample_ID = gsub("^(.+_.+)_I.+$","\\1", names)) %>% 
  dplyr::mutate(VAF_IDH = ifelse(Coverage_IDH == -1 | VAF_IDH == "NA", NA, VAF_IDH)) %>% 
  dplyr::mutate(VAF_IDH = as.numeric(VAF_IDH)) %>% 
  dplyr::mutate(VAF_upperbound = ifelse(Coverage_IDH == -1, NA, VAF_upperbound)) %>% 
  dplyr::mutate(VAF_lowerbound = ifelse(Coverage_IDH == -1, NA, VAF_lowerbound)) %>% 
  dplyr::rename(CNV.purity.shallowseq = Purity) %>% 
  dplyr::mutate(CNV.purity.shallowseq = ifelse(CNV.purity.shallowseq == "NA",NA,CNV.purity.shallowseq)) %>% 
  dplyr::mutate(CNV.purity.shallowseq = as.numeric(CNV.purity.shallowseq )) %>% 
  dplyr::filter(grepl("nDNA", names) == F) %>% 
  dplyr::mutate(CNV_ploidy_IDH = case_when(
    cn_estimate_IDH < 2.5 ~ "n=2",
    cn_estimate_IDH >= 2.5 & cn_estimate_IDH < 3.5 ~ "n=3",
    cn_estimate_IDH >= 3.5 ~ "n=4",
    T ~ "???"
  ))


ggplot(glass.cellularities, aes(x=CNV.purity.shallowseq,y=VAF_IDH * 2, label=names,col=CNV_ploidy_IDH)) +
  geom_abline(intercept = 0, slope = 1, lty=2, col="gray80") +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1.2))  +
  youri_gg_theme + 
  geom_smooth(method="lm", se = FALSE) 
  #ggrepel:: geom_text_repel(size=3, col="gray80")





ggplot(glass.cellularities, aes(x=CNV_ploidy_IDH, y=VAF_IDH)) +
  ggbeeswarm::geom_beeswarm() +
  youri_gg_theme




ggplot(glass.cellularities, aes(x=CNV.purity.shallowseq,y=VAF_IDH * 2, label=names, col=CNV_ploidy_IDH)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1.2))  +
  youri_gg_theme + 
  #geom_smooth(method="lm") +
  ggrepel:: geom_text_repel(size=3, col="gray80")





ggplot(glass.cellularities, aes(x=CNV.purity.shallowseq,y=VAF_IDH * 2, label=names, col=CNV_ploidy_IDH)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1.2))  +
  youri_gg_theme + 
  geom_smooth(method="lm") +
  ggrepel:: geom_text_repel(size=3, col="gray80")


# 1/0.5 * 0.33 * 3

glass.cellularities <- glass.cellularities %>% 
  dplyr::mutate(exp.cn.idh.mut = 1 / glass.cellularities$CNV.purity.shallowseq * glass.cellularities$VAF_IDH * glass.cellularities$Ploidy /2 ) %>% 
  dplyr::mutate(err = abs((2 * VAF_IDH) - CNV.purity.shallowseq))



ggplot(glass.cellularities, aes(x=CNV.purity.shallowseq,y=exp.cn.idh.mut, label=names, col=CNV_ploidy_IDH)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0, 1.2))  +
  youri_gg_theme + 
  #geom_smooth(method="lm") +
  ggrepel:: geom_text_repel(size=3, col="gray60")





ggplot(glass.cellularities, aes(x=CNV.purity.shallowseq,y=VAF_IDH, label=names, col=CNV_ploidy_IDH)) +
  geom_abline(intercept = 0, slope = 0.5, lty=2, col="gray80") +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1.2/2))  +
  youri_gg_theme + 
  #geom_smooth(method="lm") +
  ggrepel:: geom_text_repel(data=filter(glass.cellularities, err > 0.325), size=3, col="gray60")



tmp <- read.table( "output/tables/methylation-array/purities_RFpurity.txt") %>% 
  dplyr::select(Sample_Name, absolute, estimate)


glass.cellularities <- glass.cellularities %>% 
  dplyr::mutate(sid = gsub("^([^_]+_[^_]+).+$","\\1",names)) %>% 
  dplyr::left_join(tmp,
                   by=c('sid'='Sample_Name'))

glass.cellularities$sid


c <- glass.cellularities %>%  dplyr::select(CNV.purity.shallowseq, VAF_IDH, ManualPurity, absolute, estimate )
c <- c %>%  dplyr::filter( !is.na(CNV.purity.shallowseq) &  !is.na(VAF_IDH)  & ! is.na(ManualPurity)  & !is.na(absolute) & ! is.na(estimate))
corrplot(cor(c))



library(patchwork)


p1 <- ggplot(glass.cellularities, aes(x = CNV.purity.shallowseq, y=VAF_IDH)) +
  geom_point() +
  youri_gg_theme

p2 <- ggplot(glass.cellularities, aes(x = CNV.purity.shallowseq, y=absolute)) +
  labs(y="RFpurtity methylation Array absolute-fit") +
  geom_point() +
  youri_gg_theme

p3 <- ggplot(glass.cellularities, aes(x = VAF_IDH, y=absolute)) +
  labs(y="RFpurtity methylation Array absolute-fit") +
  geom_point() +
  youri_gg_theme

p1 + p2 + p3




ggplot(glass.cellularities, aes(y = VAF_IDH * 2, x=absolute, label=names, col=CNV_ploidy_IDH)) +
  labs(x="RFpurtity methylation Array absolute-fit") +
  geom_point() +
  #scale_x_continuous(limits = c(0.2, 0.75)) +
  #scale_y_continuous(limits = c(0, 1)) +
  youri_gg_theme



glass.cellularities <- glass.cellularities %>% 
  dplyr::mutate(resection =  gsub("^.+_([^_]+)_.+$","\\1",names))



p1 <- ggplot(glass.cellularities, aes( y=absolute,  x=resection)) +
  labs(y="RFpurtity methylation Array absolute-fit") +
  ggbeeswarm::geom_beeswarm() +
  youri_gg_theme

p2 <- ggplot(glass.cellularities, aes( y=VAF_IDH ,  x=resection)) +
  ggbeeswarm::geom_beeswarm() +
  youri_gg_theme



# ggplot(glass.cellularities, aes(x = VAF_IDH, y=estimate)) +
#   geom_point()


