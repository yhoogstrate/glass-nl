#!/usr/bin/env R

# load libs ----


library(readxl)
library(tidyverse)

source('scripts/R/job_gg_theme.R')
source('scripts/R/youri_gg_theme.R')


# load data ----


if(!exists("metadata.glass.per.resection")) {
  warning('metadata was not loaded')
  
  source('scripts/load_metadata.R')
}


# DNA-seq purities ----

# results were obtained through the tool ACE
# https://github.com/tgac-vumc/ACE/blob/master/vignettes/ACE_vignette.Rmd
# analysis was performed by the group from Amsterdam

dnaseq.purities <- readxl::read_excel('data/glass/WES/Combined Data WES-METH/AllDataGLASS.xlsx') %>% 
  dplyr::mutate(samplenames = NULL) %>% 
  dplyr::mutate(names = gsub("-","_",names,fixed=T)) %>% 
  dplyr::mutate(VAF_IDH = ifelse(Coverage_IDH == -1 | VAF_IDH == "NA", NA, VAF_IDH)) %>% 
  dplyr::mutate(VAF_IDH = as.numeric(VAF_IDH)) %>% 
  dplyr::mutate(VAF_upperbound = ifelse(Coverage_IDH == -1, NA, VAF_upperbound)) %>% 
  dplyr::mutate(VAF_lowerbound = ifelse(Coverage_IDH == -1, NA, VAF_lowerbound)) %>% 
  dplyr::rename(CNV.purity.shallowseq = Purity) %>% 
  dplyr::mutate(CNV.purity.shallowseq = ifelse(CNV.purity.shallowseq == "NA",NA,CNV.purity.shallowseq)) %>% 
  dplyr::mutate(CNV.purity.shallowseq = as.numeric(CNV.purity.shallowseq )) %>% 
  dplyr::filter(grepl("nDNA", names) == F) %>% 
  dplyr::mutate(dna.wes.ploidy_IDH = case_when(
    cn_estimate_IDH < 2.5 ~ "n=2",
    cn_estimate_IDH >= 2.5 & cn_estimate_IDH < 3.5 ~ "n=3",
    cn_estimate_IDH >= 3.5 ~ "n=4",
    T ~ "???"
  )) %>% 
  dplyr::rename(dna.wes.cn_estimate_IDH = cn_estimate_IDH) %>% 
  dplyr::rename(dna.wes.coverage_IDH = Coverage_IDH) %>% 
  dplyr::rename(dna.wes.VAF_IDH = VAF_IDH) %>% 
  dplyr::rename(dna.wes.VAF_lowerbound = VAF_lowerbound) %>% 
  dplyr::rename(dna.wes.VAF_upperbound = VAF_upperbound) %>% 
  dplyr::mutate(Sample_Name = gsub("^([^_]+_[^_]+).*$","\\1",names)) %>% 
  dplyr::mutate(names = NULL) %>% 
  dplyr::rename(dna.shallow.ACE.heterogeneity = Heterogeneity) %>% 
  dplyr::rename(dna.shallow.ACE.ploidy = Ploidy) %>% 
  dplyr::rename(dna.shallow.ACE.purity = CNV.purity.shallowseq) %>% 
  dplyr::rename(dna.purity.manual.Erik = ManualPurity) %>% 
  dplyr::mutate(dna.wes.cn_estimate_IDH = as.numeric(ifelse(dna.wes.cn_estimate_IDH == "NA", NA , dna.wes.cn_estimate_IDH)) ) %>% 
  dplyr::mutate(dna.wes.coverage_IDH = as.numeric(ifelse(dna.wes.coverage_IDH == "NA", NA , dna.wes.coverage_IDH)) ) %>% 
  dplyr::mutate(dna.wes.VAF_lowerbound = as.numeric(ifelse(dna.wes.VAF_lowerbound == "NA", NA , dna.wes.VAF_lowerbound)) ) %>% 
  dplyr::mutate(dna.wes.VAF_upperbound = as.numeric(ifelse(dna.wes.VAF_upperbound == "NA", NA , dna.wes.VAF_upperbound)) ) %>% 
  dplyr::mutate(dna.shallow.ACE.purity.below.1 = ifelse(dna.shallow.ACE.purity < 1.0,dna.shallow.ACE.purity,NA))




## append to metadata ----



metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(dnaseq.purities, by=c('Sample_Name'='Sample_Name'), keep=F,suffix = c("", "")) # force overwrite





# methylation purities using RFpurity ----


methylation.purities <- read.table("data/glass/Methylation/Analysis/RFpurity/purities_RFpurity.txt") %>% 
  dplyr::mutate(methylation.uid = paste0(Slide, "_", Array)) %>% 
  dplyr::mutate(fn = NULL) %>% 
  dplyr::mutate(Basename = NULL) %>% 
  dplyr::mutate(Slide = NULL) %>% 
  dplyr::mutate(Array = NULL) %>% 
  dplyr::rename(methylation.purity.absolute = absolute) %>% 
  dplyr::rename(methylation.purity.estimate = estimate)

stopifnot(duplicated(methylation.purities$methylation.uid) == F) # no duplicates may exist



tmp.methylation.metdata <- read.csv("data/glass/Methylation/Metadata/Datasheet4.csv") %>% 
  dplyr::mutate(methylation.uid = paste0(Slide, "_", Array)) %>% 
  dplyr::mutate(
    Basename = NULL,
    Array = NULL,
    Slide = NULL,
    Surgery_ID = NULL, 
    GLASS_ID = NULL, 
    Sample_Plate = NULL,       
    Sample_ID = NULL, 
    Sample_Resection = NULL, 
    Sample_Type = NULL, 
    Recurrent_Type = NULL, 
    Sample_Sex = NULL  )

stopifnot(duplicated(tmp.methylation.metdata$methylation.uid) == F) # no duplicates may exist


# I ran all 323 samples I found on the server with RFpurity, only a subset matches the actual glass samples (tmp.2)
stopifnot(tmp.methylation.metdata$methylation.uid %in% methylation.purities$methylation.uid)

methylation.purities <-methylation.purities %>% 
  dplyr::left_join(tmp.methylation.metdata, by=c('methylation.uid'='methylation.uid'))

rm(tmp.methylation.metdata)


## append to metadata ----



metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(methylation.purities, by=c('Sample_Name'='Sample_Name'), keep=F,suffix = c("", "")) # force overwrite



# export for Tobias/Marcel ----

tmp <- metadata.glass.per.resection %>%
  dplyr::select(-contains('fastp')) %>%
  dplyr::select(-contains('featureC')) %>%
  dplyr::select(-contains('idxstats')) %>%
  dplyr::select(-contains('star')) %>%
  dplyr::select(
    -c("resection",
       "Sample_Type",
       "Recurrent_Type",
       "Sample_Sex",
       "Exclude.by.Wies.on.complete.pair", 
       "Surgery_ID",
       "methylation.purity.estimate",
       "institute",
       "rid",
       "assigned.reads.status",
       "excluded.reason",
       "excluded",
       "time.resection.until.last.event", 
       "status.resection.until.last.event"
       )
  ) 

#write.csv(tmp, "/tmp/2022-03-22_purity-export_GLASS_NL_YH.csv")
rm(tmp)



# some basic qc plots ----


plt <- dnaseq.purities %>%
  #dplyr::full_join(methylation.purities, by=c('Sample_Name'='Sample_Name')) %>% 
  dplyr::mutate(status = ifelse(dna.shallow.ACE.purity >= 0.99,"unconfident","regular")) %>% 
  dplyr::mutate(status.manual.fit = ifelse(dna.purity.manual.Erik == dna.shallow.ACE.purity, "ACE", "VAF"))



ggplot(plt, aes(x=dna.shallow.ACE.purity,y=dna.wes.VAF_IDH * 2, label=Sample_Name,col=dna.wes.ploidy_IDH)) +
  geom_abline(intercept = 0, slope = 1, lty=2, col="gray80") +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1.2))  +
  youri_gg_theme +
  geom_smooth(method="lm", se = FALSE)
  #ggrepel:: geom_text_repel(size=3, col="gray80")



rm(plt)

## ACE x WES/VAF ----



plt <- dnaseq.purities %>%
  dplyr::full_join(methylation.purities, by=c('Sample_Name'='Sample_Name')) %>% 
  dplyr::mutate(status = ifelse(dna.shallow.ACE.purity >= 0.99,"unconfident","regular")) %>% 
  dplyr::mutate(status.manual.fit = ifelse(dna.purity.manual.Erik == dna.shallow.ACE.purity, "ACE", "VAF"))



ggplot(plt, aes(x=dna.shallow.ACE.purity,y=dna.wes.VAF_IDH * 2, label=Sample_Name,col=status)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1.2)) +
  scale_color_manual(values=c('unconfident'='red','regular'='black')) +
  youri_gg_theme 


rm(plt)


## ACE x Meth/RF[abs] ----


plt <- dnaseq.purities %>%
  dplyr::full_join(methylation.purities, by=c('Sample_Name'='Sample_Name')) %>% 
  dplyr::mutate(status = ifelse(dna.shallow.ACE.purity >= 0.99,"unconfident","regular")) %>% 
  dplyr::mutate(status.manual.fit = ifelse(dna.purity.manual.Erik == dna.shallow.ACE.purity, "ACE", "VAF"))


ggplot(plt, aes(x=dna.shallow.ACE.purity,y=methylation.purity.absolute, label=Sample_Name,col=status)) +
  geom_point() +
  #scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0.25, 0.75)) +
  scale_color_manual(values=c('unconfident'='red','regular'='black')) +
  youri_gg_theme 

rm(plt)


### log(ACE) x Meth/RF[abs] ----


plt <- dnaseq.purities %>%
  dplyr::full_join(methylation.purities, by=c('Sample_Name'='Sample_Name')) %>% 
  dplyr::mutate(status = ifelse(dna.shallow.ACE.purity >= 0.99,"unconfident","regular")) %>% 
  dplyr::mutate(status.manual.fit = ifelse(dna.purity.manual.Erik == dna.shallow.ACE.purity, "ACE", "VAF"))


ggplot(plt, aes(x=log(dna.shallow.ACE.purity),y=methylation.purity.absolute, label=Sample_Name,col=status)) +
  geom_point() +
  #scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0.25, 0.75)) +
  scale_color_manual(values=c('unconfident'='red','regular'='black')) +
  youri_gg_theme 

rm(plt)

## Manual Fit x Meth/RF[abs] ----


plt <- metadata.glass.per.resection %>%
  dplyr::mutate(status = ifelse(dna.shallow.ACE.purity >= 0.99,"unconfident","regular")) %>% 
  dplyr::mutate(status.manual.fit = ifelse(dna.purity.manual.Erik == dna.shallow.ACE.purity, "ACE", "VAF"))


ggplot(plt, aes(x=dna.purity.manual.Erik,y=methylation.purity.absolute, label=Sample_Name,col=status.manual.fit)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0.25, 0.75))  +
  youri_gg_theme

rm(plt)


## WES/VAF x Meth/RF[abs] ----



plt <- metadata.glass.per.resection %>%
  dplyr::mutate(status = ifelse(dna.shallow.ACE.purity >= 0.99,"unconfident","regular")) %>% 
  dplyr::mutate(status.manual.fit = ifelse(dna.purity.manual.Erik == dna.shallow.ACE.purity, "ACE", "VAF"))


ggplot(plt, aes(x=dna.wes.VAF_IDH * 2,y=methylation.purity.absolute, label=Sample_Name)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0.25, 0.75))  +
  youri_gg_theme 

rm(plt)




## ACE x Meth/RF[est] ----

# ggplot(plt, aes(x=dna.shallow.ACE.purity,y=methylation.purity.estimate, label=Sample_Name)) +
#   geom_point() +
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0, 1.0))  +
#   youri_gg_theme 


## WES/VAF x Meth/RF[est] -- -


# ggplot(plt, aes(x=dna.wes.VAF_IDH * 2,y=methylation.purity.estimate, label=Sample_Name)) +
#   geom_point() +
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0.5, 1.0))  +
#   youri_gg_theme 


## beeswarm ----



plt <- metadata.glass.per.resection %>%
  dplyr::mutate(status = ifelse(dna.shallow.ACE.purity >= 0.99,"unconfident","regular")) %>% 
  dplyr::mutate(status.manual.fit = ifelse(dna.purity.manual.Erik == dna.shallow.ACE.purity, "ACE", "VAF"))


ggplot(plt %>%  dplyr::mutate(resection = gsub("^.+_(.+)$","\\1",Sample_Name)), aes(x=resection, y=dna.wes.VAF_IDH)) +
  ggbeeswarm::geom_beeswarm() +
  youri_gg_theme





# 
# ggplot(plt, aes(x=CNV.purity.shallowseq,y=VAF_IDH * 2, label=names, col=dna.wes.ploidy_IDH)) +
#   geom_point() +
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0, 1.2))  +
#   youri_gg_theme + 
#   #geom_smooth(method="lm") +
#   ggrepel:: geom_text_repel(size=3, col="gray80")
# 
# 
# 
# 
# 
# ggplot(glass.cellularities, aes(x=CNV.purity.shallowseq,y=VAF_IDH * 2, label=names, col=CNV_ploidy_IDH)) +
#   geom_point() +
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0, 1.2))  +
#   youri_gg_theme + 
#   geom_smooth(method="lm") +
#   ggrepel:: geom_text_repel(size=3, col="gray80")
# 

rm(plt)


# 1/0.5 * 0.33 * 3



# 
# plt <- dnaseq.purities %>%
#   dplyr::mutate(status = ifelse(dna.shallow.ACE.purity >= 0.99,"unconfident","regular")) %>% 
#   dplyr::mutate(status.manual.fit = ifelse(dna.purity.manual.Erik == dna.shallow.ACE.purity, "ACE", "VAF"))
# 
# 
# plt <- metadata.glass.per.resection %>% 
#   dplyr::mutate(exp.cn.idh.mut = 1 / glass.cellularities$CNV.purity.shallowseq * glass.cellularities$VAF_IDH * glass.cellularities$Ploidy /2 ) %>% 
#   dplyr::mutate(err = abs((2 * VAF_IDH) - CNV.purity.shallowseq))
# 
# 
# 
# ggplot(glass.cellularities, aes(x=CNV.purity.shallowseq,y=exp.cn.idh.mut, label=names, col=CNV_ploidy_IDH)) +
#   geom_point() +
#   scale_x_continuous(limits = c(0, 1)) +
#   #scale_y_continuous(limits = c(0, 1.2))  +
#   youri_gg_theme + 
#   #geom_smooth(method="lm") +
#   ggrepel:: geom_text_repel(size=3, col="gray60")
# 
# 


# 
# ggplot(glass.cellularities, aes(x=CNV.purity.shallowseq,y=VAF_IDH, label=names, col=CNV_ploidy_IDH)) +
#   geom_abline(intercept = 0, slope = 0.5, lty=2, col="gray80") +
#   geom_point() +
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0, 1.2/2))  +
#   youri_gg_theme + 
#   #geom_smooth(method="lm") +
#   ggrepel:: geom_text_repel(data=filter(glass.cellularities, err > 0.325), size=3, col="gray60")
# 
# 
# 
# 
# glass.cellularities <- glass.cellularities %>% 
#   dplyr::mutate(sid = gsub("^([^_]+_[^_]+).+$","\\1",names)) %>% 
#   dplyr::left_join(tmp,
#                    by=c('sid'='Sample_Name'))
# 
# 
# 
# rm(plt)


# 
# c <- glass.cellularities %>%  dplyr::select(CNV.purity.shallowseq, VAF_IDH, ManualPurity, absolute, estimate )
# c <- c %>%  dplyr::filter( !is.na(CNV.purity.shallowseq) &  !is.na(VAF_IDH)  & ! is.na(ManualPurity)  & !is.na(absolute) & ! is.na(estimate))
# corrplot(cor(c))


# 
# library(patchwork)
# 
# 
# p1 <- ggplot(glass.cellularities, aes(x = CNV.purity.shallowseq, y=VAF_IDH)) +
#   geom_point() +
#   youri_gg_theme
# 
# p2 <- ggplot(glass.cellularities, aes(x = CNV.purity.shallowseq, y=absolute)) +
#   labs(y="RFpurtity methylation Array absolute-fit") +
#   geom_point() +
#   youri_gg_theme
# 
# p3 <- ggplot(glass.cellularities, aes(x = VAF_IDH, y=absolute)) +
#   labs(y="RFpurtity methylation Array absolute-fit") +
#   geom_point() +
#   youri_gg_theme
# 
# p1 + p2 + p3
# 
# 
# 
# 
# ggplot(glass.cellularities, aes(y = VAF_IDH * 2, x=absolute, label=names, col=CNV_ploidy_IDH)) +
#   labs(x="RFpurtity methylation Array absolute-fit") +
#   geom_point() +
#   #scale_x_continuous(limits = c(0.2, 0.75)) +
#   #scale_y_continuous(limits = c(0, 1)) +
#   youri_gg_theme
# 
# 
# 
# glass.cellularities <- glass.cellularities %>% 
#   dplyr::mutate(resection =  gsub("^.+_([^_]+)_.+$","\\1",names))
# 
# 
# 
# p1 <- ggplot(glass.cellularities, aes( y=absolute,  x=resection)) +
#   labs(y="RFpurtity methylation Array absolute-fit") +
#   ggbeeswarm::geom_beeswarm() +
#   youri_gg_theme
# 
# p2 <- ggplot(glass.cellularities, aes( y=VAF_IDH ,  x=resection)) +
#   ggbeeswarm::geom_beeswarm() +
#   youri_gg_theme
# 
# 
# 
# # ggplot(glass.cellularities, aes(x = VAF_IDH, y=estimate)) +
# #   geom_point()
# 

# cleanup ----




