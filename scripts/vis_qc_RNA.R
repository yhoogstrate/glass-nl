#!/usr/bin/env R

# libs & functions ----

library(tidyverse)
library(patchwork)
library(factoextra)


source('scripts/R/youri_gg_theme.R')
source('scripts/R/job_gg_theme.R')


# data ----


source('scripts/metadata.R')


# plots ----

## 1. find those with largest diff in stats ----



# samples separated showed the largest within-patient differences in at least on QC metric:
# 
# low-comp
# 104059-003-009
# 104059-001-065
# 104059-002-090

# %C
# 104059-002-079

# reads passed filtering
# 104059-002-084
# 104059-002-125
# 104059-001-065
# 104059-002-090
# 104059-002-079
# 104059-002-168
# 104059-003-009
# 104059-002-165


# avg.read.trim
# 104059-001-065
# 104059-002-090
# 104059-003-009
# 104059-002-084

# fastp.insert_size.0.ratio
# 104059-002-084
# 104059-002-090




plt  <- data.frame()
for(sid in unique(metadata.glass.per.fastq$genomescan.sid)) {
  sel <- metadata.glass.per.fastq %>% dplyr::filter(genomescan.sid == sid)
  
  plt <- rbind(plt,
                 data.frame(
                   genomescan.sid = sid,
                   
                   fastp.gc.rmse = max(sel$fastp.gc.rmse) - min(sel$fastp.gc.rmse),
                   fastp.percentage.c = max(sel$fastp.percentage.c) - min(sel$fastp.percentage.c),
                   fastp.duplication.rate = max(sel$fastp.duplication.rate) - min(sel$fastp.duplication.rate),
                   fastp.low.complexity.reads = max(sel$fastp.low.complexity.reads) - min(sel$fastp.low.complexity.reads),
                   fastp.ratio.reads.passed.filtering = max(sel$fastp.ratio.reads.passed.filtering) - min(sel$fastp.ratio.reads.passed.filtering),
                   avg.read.trim = max(sel$fastp.avg.read.trim) - min(sel$fastp.avg.read.trim),
                   fastp.insert_size.0.ratio = max(sel$fastp.insert_size.0.ratio) - min(sel$fastp.insert_size.0.ratio)
                 ))
}


plt <- plt %>% 
  dplyr::mutate(order = rank(fastp.insert_size.0.ratio))

p0 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=fastp.gc.rmse)) + 
  geom_point()
p1 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=fastp.percentage.c)) + 
  geom_point()
p2 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=fastp.duplication.rate)) + 
  geom_point()
p3 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=fastp.low.complexity.reads)) + 
  geom_point()
p4 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=fastp.ratio.reads.passed.filtering)) + 
  geom_point()
p5 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=avg.read.trim)) + 
  geom_point()
p6 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=fastp.insert_size.0.ratio)) + 
  geom_point()



p1 / p2 / p3 / p4 / p5 / p6




## 2. fastq : sequencing run ----



plt.top <- metadata.glass.per.fastq %>%
  dplyr::mutate(order = fastp.low.complexity.reads) %>% 
  dplyr::mutate(order = fastp.duplication.rate) %>% 
  dplyr::mutate(facet = case_when(
    genomescan.sid %in% c('104059-001-065') ~ "104059-001-065",
    genomescan.sid %in% c('104059-002-079') ~ "104059-002-079",
    genomescan.sid %in% c('104059-002-084') ~ "104059-002-084",
    genomescan.sid %in% c('104059-002-090') ~ "104059-002-090",
    genomescan.sid %in% c('104059-002-125') ~ "104059-002-125",
    genomescan.sid %in% c('104059-002-165') ~ "104059-002-165",
    genomescan.sid %in% c('104059-002-168') ~ "104059-002-168",
    genomescan.sid %in% c('104059-003-009') ~ "104059-003-009",
    
    T ~ "other"
  ))



plt.pct <- plt.top %>% 
  dplyr::select(fastp.json, genomescan.sid,
                fastp.percentage.a,
                #                fastp.percentage.c,
                #                fastp.percentage.t,
                fastp.percentage.g,
                fastp.gc.rmse,
                order, facet) %>%
  pivot_longer(cols = -c(fastp.json, genomescan.sid, fastp.gc.rmse, order,facet )) %>% 
  dplyr::rename(percentage=value) %>% 
  dplyr::mutate(vis.type = "scatter") 

plt.conf <- plt.pct %>%
  dplyr::select(fastp.json, genomescan.sid, fastp.gc.rmse, order, facet) %>%
  dplyr::mutate(conf.lower = 25 - (fastp.gc.rmse/2)) %>%
  dplyr::mutate(conf.upper = 25 + (fastp.gc.rmse/2)) %>% 
  pivot_longer(cols = -c(fastp.json, genomescan.sid, fastp.gc.rmse, order, facet)) %>% 
  dplyr::rename(percentage=value) %>% 
  dplyr::mutate(vis.type = "lines")


plt <- rbind(plt.pct, plt.conf)


p1 <-ggplot(plt, aes(x = reorder(fastp.json, order ), y=percentage, label=genomescan.sid, color=name)) +
  #geom_line(data = plt %>% dplyr::filter(vis.type == "scatter"), aes(group=fastp.json), col="gray70", alpha=50) + 
  geom_line(data = plt %>% dplyr::filter(vis.type == "lines"), aes(group=fastp.json), col="black", alpha=50) + 
  geom_hline(yintercept=25, lty=2 , lwd=0.5) +
  #geom_point(data = plt %>% dplyr::filter(vis.type == "scatter")) +
  labs(x=NULL) +
  facet_grid(cols = vars(facet), scales = "free", space="free") +
  theme(legend.position = 'bottom', axis.text.x = element_blank())



plt <- rbind(plt.top, plt.top %>% dplyr::mutate(fastp.duplication.rate = 0))
p2 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=fastp.duplication.rate, label=genomescan.sid, group=fastp.json)) +
  geom_line() +
  labs(x=NULL, y='dupl rate') +
  facet_grid(cols = vars(facet), scales = "free", space="free") + 
  theme(legend.position = 'bottom', axis.text.x = element_blank())


plt <- rbind(plt.top, plt.top %>% dplyr::mutate(fastp.low.complexity.reads = 0))
p3 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=fastp.low.complexity.reads, label=genomescan.sid, group=fastp.json)) +
  geom_line() +
  labs(x=NULL, y='low cplxty') +
  facet_grid(cols = vars(facet), scales = "free", space="free") + 
  theme(legend.position = 'bottom', axis.text.x = element_blank())


plt <- rbind(plt.top, plt.top %>% dplyr::mutate(fastp.ratio.reads.not.passed.filtering = 0))
p4 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=fastp.ratio.reads.not.passed.filtering, label=genomescan.sid, group=fastp.json)) +
  geom_line() +
  labs(x=NULL, y='not passed filt') +
  facet_grid(cols = vars(facet), scales = "free", space="free") + 
  theme(legend.position = 'bottom', axis.text.x = element_blank())


plt <- rbind(plt.top, plt.top %>% dplyr::mutate(fastp.avg.read.trim = 0))
p5 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=fastp.avg.read.trim, label=genomescan.sid, group=fastp.json)) +
  geom_line() +
  labs(x=NULL, y='avg read trim') +
  facet_grid(cols = vars(facet), scales = "free", space="free") + 
  theme(legend.position = 'bottom', axis.text.x = element_blank())


plt <- rbind(plt.top, plt.top %>% dplyr::mutate(fastp.insert_size.avg = 0))
p6 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=fastp.insert_size.avg, label=genomescan.sid, group=fastp.json)) +
  geom_line() +
  labs(x=NULL, y='avg ins-size') +
  facet_grid(cols = vars(facet), scales = "free", space="free") +
  theme(legend.position = 'bottom', axis.text.x = element_blank())


plt <- rbind(plt.top, plt.top %>% dplyr::mutate(fastp.insert_size.0.ratio = 0))
p7 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=fastp.insert_size.0.ratio, label=genomescan.sid, group=fastp.json)) +
  geom_line() +
  labs(x=NULL, y='ratio ins-size=0') +
  facet_grid(cols = vars(facet), scales = "free", space="free") + 
  theme(legend.position = 'bottom', axis.text.x = element_blank())


plt <- plt.top %>% 
  dplyr::mutate(y=1)
p8 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=y , col=flow.cell , fill=flow.cell)) +
  geom_tile() +
  facet_grid(cols = vars(facet), scales = "free", space="free")  + 
  labs(x=NULL, y=NULL) +
  theme(legend.position = "none", axis.text.x = element_blank())


p1 / p2 / p3 / p4 / p5 / p7/ p8
# ggsave("output/figures/qc/fastp_qc_metrics.png",width=14,height=7)



## 3. per res pca ----




plt <- metadata.glass.per.resection %>% 
  dplyr::select(
    c("genomescan.sid","fastp.low.complexity.reads",
      "fastp.duplication.rate",
      "fastp.percentage.a",
      "fastp.gc.rmse","fastp.ratio.reads.not.passed.filtering",
      'fastp.insert_size.0.ratio',
      "fastp.avg.read.trim",
      "idxstats.freq.alternate.loci",
      'star.pct.uniquely.mapped.reads')
  ) %>% 
  tibble::column_to_rownames('genomescan.sid')

res.pca <- prcomp(plt, scale = TRUE)

fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = ifelse(metadata.glass.per.resection$featureCounts.Assigned >= 750000, "red","green")
                #col.ind = "#696969"  # Individuals color
                #,label="var" 
)




## per resection idxstats-rainbow ----


plt <- metadata.glass.per.resection %>% 
  dplyr::select('genomescan.sid',contains('idxstats')) %>% 
  dplyr::mutate(idxstats.freq.alternate.loci = NULL, idxstats.total = NULL) %>% 
  tidyr::pivot_longer(cols = -c(genomescan.sid)) %>% 
  dplyr::mutate(name = gsub("idxstats.","",name)) %>% 
  dplyr::mutate(name = gsub("alternate.loci","alternate loci",name)) %>% 
  dplyr::group_by(genomescan.sid) %>% 
  dplyr::mutate(n = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(freq = value / n * 100) 


order <- plt %>% 
  dplyr::select(c('genomescan.sid','name','freq')) %>% 
  dplyr::filter(name %in% c("chr21","alternate loci")) %>% 
  tidyr::pivot_wider(names_from = name, values_from = freq) %>% 
  dplyr::mutate(product = chr21 * `alternate loci`) %>% 
  dplyr::mutate(order = rank(-product)) %>% 
  dplyr::select(genomescan.sid, order)


plt <- plt %>% 
  dplyr::left_join(order, by=c('genomescan.sid'='genomescan.sid'))


ggplot(plt, aes(x = reorder(genomescan.sid,order) ,y = freq, fill=name, label=genomescan.sid)) +
  #coord_flip() + 
  geom_bar(stat = "identity", position = "stack",colour="black") + 
  scale_y_continuous(labels = scales::unit_format(unit = "%")) + 
  theme_bw() + 
  theme( axis.title.y = element_text(size = 11) ,
         axis.text.x = element_text(angle = 90, size = 5 ))

ggsave("output/figures/qc/per_patients__idxstats-rainbox.png",width=15,height=6)



## all / patient sample ----



plt <- metadata.glass.per.resection %>% 
  dplyr::mutate(order = rank(-rank(featureCounts.M.Assigned))) %>% 
  dplyr::mutate(assigned.reads.status = factor(
    ifelse(featureCounts.Assigned > 750000,"PASS","INSUFFICIENT"),
    levels=c("PASS","INSUFFICIENT")))


plt <- rbind(plt, plt %>%
               dplyr::mutate(fastp.insert_size.0.ratio = 0,
                             fastp.avg.read.trim = 0,
                             fastp.gc.rmse = 0,
                             fastp.percentage.c = 23.6,    # qc.stats.collapsed %>% dplyr::filter(order < 120) %>% dplyr::pull(fastp.percentage.c) %>%  median
                             fastp.low.complexity.reads = 0,
                             fastp.duplication.rate = 0,
                             fastp.ratio.reads.not.passed.filtering = 0,
                             
                             featureCounts.M.Assigned = 0,
                             
                             idxstats.freq.alternate.loci = 0)) %>% 
  dplyr::select(genomescan.sid,order,assigned.reads.status,
                fastp.insert_size.0.ratio, featureCounts.M.Assigned,
                fastp.avg.read.trim, fastp.gc.rmse, fastp.percentage.c, fastp.low.complexity.reads,
                fastp.duplication.rate, fastp.ratio.reads.not.passed.filtering, idxstats.freq.alternate.loci) %>% 
  pivot_longer(cols = -c(genomescan.sid,order,assigned.reads.status)) %>% 
  dplyr::mutate(name = gsub( "fastp.avg.read.trim" ,'Avg trim',name)) %>% 
  dplyr::mutate(name = gsub( "fastp.duplication.rate" ,'Dupl rate',name)) %>% 
  dplyr::mutate(name = gsub( "fastp.gc.rmse" ,'ACTG RMSE',name)) %>% 
  dplyr::mutate(name = gsub( "fastp.insert_size.0.ratio" ,'Inst size != 0',name)) %>% 
  dplyr::mutate(name = gsub( "fastp.low.complexity.reads" ,'Low cplxty',name)) %>% 
  dplyr::mutate(name = gsub( "fastp.percentage.c" ,'%C',name)) %>% 
  dplyr::mutate(name = gsub( "fastp.ratio.reads.not.passed.filtering" ,'Not pass filter',name)) %>% 
  dplyr::mutate(name = gsub( "featureCounts.M.Assigned" ,'M Assigned',name)) %>% 
  dplyr::mutate(name = gsub( "idxstats.freq.alternate.loci" ,'rRNA; Alternate loci',name)) %>% 
  dplyr::mutate(name = factor(name, levels=c("M Assigned", "Avg trim", "ACTG RMSE", "%C", "Low cplxty", "Dupl rate","Inst size != 0",  "Not pass filter", "rRNA; Alternate loci")))



p1 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=value, group=genomescan.sid)) +
  geom_line() +
  labs(x=NULL) +
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size=6)) + 
  theme(axis.text.x = element_blank()) + 
  facet_grid(cols = vars(assigned.reads.status), row=vars(name), scale="free", space= "free_x")





plt <- metadata.glass.per.resection %>% 
  dplyr::mutate(order = rank(-rank(featureCounts.M.Assigned))) %>% 
  dplyr::mutate(assigned.reads.status = factor(
    ifelse(featureCounts.Assigned > 750000,"PASS","INSUFFICIENT"),
    levels=c("PASS","INSUFFICIENT")))


p2 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=1, fill=institute)) +
  geom_tile(col="black") +
  facet_grid(cols = vars(assigned.reads.status),space= "free_x",scale="free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size=6)) +
  scale_color_brewer(palette="Accent")



p1 / p2 +   plot_layout(heights = c(10,0.5))



ggsave("output/figures/qc/per_patients__all-stats-overview.png",width=10,height=1.5)





# all / patient combined



plt <- qc.stats.collapsed %>% 
  dplyr::left_join(metadata.glass.per.resection, by = c('genomescan.sid'='GS_ID')) %>%
  dplyr::mutate(col = assigned.reads.status) %>% 
  dplyr::mutate(resection = factor(resection, levels = c('R1','R2','R3','R4','MW1','MW2','MW7') )) %>% 
  dplyr::filter(group == "Glioma")

ggplot(plt, aes(x=pid, y=resection, fill=col)) +
  geom_tile() +
  scale_fill_manual(values=c('TOO LOW'='red','PASSED'='darkgreen')) +
  youri_gg_theme








