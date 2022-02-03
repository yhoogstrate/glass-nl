#!/usr/bin/env R

# visualise all RNA QC stats


library(tidyverse)
library(patchwork)
library(rjson)

source('scripts/R/youri_gg_theme.R')
source('scripts/R/job_gg_theme.R')


# fastp stats ----


parse_fastp_json_files <- function(json_file) {
  out <- fromJSON(file = json_file)
  
  out <- list(
    "a" = sum(out$read1_after_filtering$content_curves$`A`, out$read2_after_filtering$content_curves$`A`),
    "c" = sum(out$read1_after_filtering$content_curves$`C`, out$read2_after_filtering$content_curves$`C`),
    "t" = sum(out$read1_after_filtering$content_curves$`T`, out$read2_after_filtering$content_curves$`T`),
    "g" = sum(out$read1_after_filtering$content_curves$`G`, out$read2_after_filtering$content_curves$`G`),
    "fastp.low.complexity.reads" = out$filtering_result$low_complexity_reads / out$summary$before_filtering$total_reads * 100.0,
    "fastp.duplication.rate" = out$duplication$rate,
    "fastp.avg.duplicates.per.read" = 1.0 / (1.0 - out$duplication$rate),
    "fastp.total_reads" = out$summary$before_filtering$total_reads,
    "fastp.read1_mean_length" = out$summary$after_filtering$read1_mean_length,
    "fastp.read2_mean_length" = out$summary$after_filtering$read2_mean_length,
    'fastp.ratio.reads.passed.filtering' = as.numeric(out$filtering_result$passed_filter_reads) / as.numeric(out$summary$before_filtering$total_reads) * 100,
    'fastp.insert_size.0.ratio' = out$insert_size$histogram[1] / sum(out$insert_size$histogram),
    'fastp.insert_size.unknown.ratio' = as.numeric(out$insert_size$unknown) / as.numeric(sum(c(out$insert_size$unknown, out$insert_size$histogram))) * 100.0,
    'fastp.insert_size.avg' = sum(out$insert_size$histogram * c(0:(length(out$insert_size$histogram)-1)))/sum(out$insert_size$histogram)
    )
  out$n <- out$a + out$c + out$t + out$g
  
  out$fastp.percentage.a = 100.0 * out$a / out$n
  out$fastp.percentage.c = 100.0 * out$c / out$n
  out$fastp.percentage.t = 100.0 * out$t / out$n
  out$fastp.percentage.g = 100.0 * out$g / out$n
  
  out$fastp.gc.rmse = ((
      (26.4 - out$fastp.percentage.a)^2 +
      (23.6 - out$fastp.percentage.c)^2 +
      (26.4 - out$fastp.percentage.t)^2 +
      (23.6 - out$fastp.percentage.g)^2
  ) / 4.0) ^0.5
  
  out$a <- NULL
  out$c <- NULL
  out$t <- NULL
  out$g <- NULL
  out$n <- NULL
  
  return(out)
}


qc.stats <- data.frame(fastp.json = Sys.glob("data/glass/RNAseq/fastq-clean/*.json")) %>% 
  #head(n=25) %>% 
  dplyr::mutate(genomescan.sid = gsub("^.+/[^_]+_([^_]+)_.+$","\\1", fastp.json)) %>% 
  dplyr::mutate(json.stats = lapply(fastp.json, parse_fastp_json_files)) %>% 
  dplyr::mutate(tmp = map(json.stats, ~ data.frame(t(.)))) %>%  unnest(tmp) %>% dplyr::mutate(json.stats = NULL) %>% as.data.frame %>% 
  dplyr::mutate(fastp.avg.duplicates.per.read = as.numeric(fastp.avg.duplicates.per.read)) %>% 
  dplyr::mutate(fastp.percentage.a = as.numeric(fastp.percentage.a)) %>% 
  dplyr::mutate(fastp.percentage.c = as.numeric(fastp.percentage.c)) %>% 
  dplyr::mutate(fastp.percentage.t = as.numeric(fastp.percentage.t)) %>% 
  dplyr::mutate(fastp.percentage.g = as.numeric(fastp.percentage.g)) %>% 
  dplyr::mutate(fastp.duplication.rate = as.numeric(fastp.duplication.rate)) %>% 
  dplyr::mutate(fastp.ratio.reads.passed.filtering = as.numeric(fastp.ratio.reads.passed.filtering)) %>% 
  dplyr::mutate(fastp.ratio.reads.not.passed.filtering = 100 - fastp.ratio.reads.passed.filtering) %>% 
  dplyr::mutate(fastp.insert_size.0.ratio = as.numeric(fastp.insert_size.0.ratio)) %>% 
  dplyr::mutate(fastp.insert_size.unknown.ratio = as.numeric(fastp.insert_size.unknown.ratio)) %>% 
  dplyr::mutate(fastp.insert_size.avg = as.numeric(fastp.insert_size.avg)) %>% 
  dplyr::mutate(fastp.total_reads = as.numeric(fastp.total_reads)) %>% 
  dplyr::mutate(fastp.gc.rmse = as.numeric(fastp.gc.rmse)) %>% 
  dplyr::mutate(fastp.low.complexity.reads = as.numeric(fastp.low.complexity.reads)) %>% 
  dplyr::mutate(fastp.read1_mean_length = as.numeric(fastp.read1_mean_length)) %>% 
  dplyr::mutate(fastp.read2_mean_length = as.numeric(fastp.read2_mean_length)) %>% 
  dplyr::mutate(avg.read.trim = 150 - (fastp.read1_mean_length + fastp.read2_mean_length)/ 2) %>% 
  dplyr::mutate(flow.cell = as.factor(gsub("^.+/([^_]+)_.+$","\\1",fastp.json))) %>% 
  dplyr::mutate(lane = gsub("^.+_(L[0-9]+)_.+$","\\1",fastp.json )) %>% 
  dplyr::mutate(lane.flow.cell = as.factor( paste0(flow.cell , "_", lane )))



qc.stats <- qc.stats %>% 
  dplyr::mutate(order = fastp.low.complexity.reads) %>% 
  dplyr::mutate(order = fastp.duplication.rate) %>% 
  dplyr::mutate(facet = case_when(
    #genomescan.sid %in% c('104059-003-002') ~ "104059-003-002", # rare qc
    #genomescan.sid %in% c('104059-002-099') ~ "104059-002-099",
    #genomescan.sid %in% c('104059-002-087') ~ "104059-002-087",
    
    
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


qc.stats.collapsed <- qc.stats %>% 
  dplyr::mutate(fastp.json = NULL, order=NULL) %>% 
  dplyr::group_by(genomescan.sid) %>% 
  dplyr::summarise(
    avg.read.trim = weighted.mean(avg.read.trim, fastp.total_reads),
    fastp.avg.duplicates.per.read = weighted.mean(fastp.avg.duplicates.per.read, fastp.total_reads),
    fastp.duplication.rate = weighted.mean(fastp.duplication.rate, fastp.total_reads),
    fastp.gc.rmse = weighted.mean(fastp.gc.rmse, fastp.total_reads),
    fastp.insert_size.0.ratio = weighted.mean(fastp.insert_size.0.ratio, fastp.total_reads),
    fastp.insert_size.avg = weighted.mean(fastp.insert_size.avg, fastp.total_reads),
    fastp.insert_size.unknown.ratio = weighted.mean(fastp.insert_size.unknown.ratio, fastp.total_reads),
    fastp.low.complexity.reads = weighted.mean(fastp.low.complexity.reads, fastp.total_reads),
    fastp.percentage.a = weighted.mean(fastp.percentage.a, fastp.total_reads),
    fastp.percentage.c = weighted.mean(fastp.percentage.c, fastp.total_reads),
    fastp.percentage.g = weighted.mean(fastp.percentage.g, fastp.total_reads),
    fastp.percentage.t = weighted.mean(fastp.percentage.t, fastp.total_reads),
    fastp.ratio.reads.not.passed.filtering = weighted.mean(fastp.ratio.reads.not.passed.filtering, fastp.total_reads),
    fastp.ratio.reads.passed.filtering = weighted.mean(fastp.ratio.reads.passed.filtering, fastp.total_reads),
    fastp.read1_mean_length = weighted.mean(fastp.read1_mean_length, fastp.total_reads),
    fastp.read2_mean_length = weighted.mean(fastp.read2_mean_length, fastp.total_reads),
    
    fastp.total_reads = sum(fastp.total_reads)
    )





# idxstats stats  ----


parse_idxstats <- function(idxstats_file) {
  #tmp <- read.delim(idxstats.stats$idxstats[1],header=F) %>% 
  tmp <- read.delim(idxstats_file,header=F) %>% 
    dplyr::mutate(V2 = NULL, V4 = NULL )
  
  out <- tmp$V3
  names(out) <- tmp$V1
  out <- as.list(out)
  
  return(out)
}


idxstats.stats <- data.frame(idxstats = Sys.glob("output/tables/qc/idxstats/*.txt")) %>% 
  dplyr::mutate(genomescan.sid =  gsub("^.+/([^/]+).samtools.idxstats.txt$","\\1", idxstats)   ) %>%
  dplyr::arrange(genomescan.sid) %>% 
  dplyr::mutate(stats = lapply(idxstats, parse_idxstats)) %>% 
  mutate(tmp = map(stats, ~ data.frame(t(.)))) %>%
  tidyr::unnest(tmp) %>%
  dplyr::mutate(stats = NULL) %>%
  as.data.frame %>% 
  dplyr::mutate_if(colnames(.) %in% c('idxstats', 'genomescan.sid') == F , as.numeric)


# STAR alignment metrics ----

parse_star_log_final_out <- function(star_log_final_out) {
  
  # 'data/glass/RNAseq/alignments/alignments-new/104059-001-012/Log.final.out'
  tmp <- read.delim(star_log_final_out,header=F)
  
  out <- list(
    'star.input.reads' = tmp %>% dplyr::filter(grepl('Number of input reads',V1)) %>% dplyr::pull(V2) %>% as.numeric,
    'star.n.uniquely.mapped.reads' = tmp %>% dplyr::filter(grepl('Uniquely mapped reads number',V1)) %>%  dplyr::pull(V2) %>% as.numeric,
    'star.pct.uniquely.mapped.reads' = tmp %>% dplyr::filter(grepl('Uniquely mapped reads %',V1)) %>% dplyr::mutate(V2= gsub('%','',V2)) %>% dplyr::pull(V2) %>%  as.numeric
  )
  return(out)
}




star.log.final.out.stats <- data.frame(star.log.final.out = Sys.glob("data/glass/RNAseq/alignments/alignments-new/*/Log.final.out")) %>% 
  dplyr::mutate(genomescan.sid =  gsub("^.+new/([^/]+)/Log.+$","\\1", star.log.final.out)) %>%
  dplyr::arrange(genomescan.sid) %>% 
  dplyr::mutate(stats = lapply(star.log.final.out, parse_star_log_final_out)) %>% 
  mutate(tmp = map(stats, ~ data.frame(t(.)))) %>%
  tidyr::unnest(tmp) %>%
  dplyr::mutate(stats = NULL) %>%
  as.data.frame %>% 
  dplyr::mutate_if(colnames(.) %in% c('star.log.final.out', 'genomescan.sid') == F , as.numeric) %>% 
  dplyr::mutate(star.log.final.out = NULL)


# featureCounts.stats ----

featureCounts.stats <- read.delim("data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt.summary",skip=0,header=T,check.names=F) %>% 
  tibble::column_to_rownames('Status') %>% 
  t() %>% 
  as.data.frame %>% 
  `colnames<-`(paste0("featureCounts.",colnames(.))) %>% 
  tibble::rownames_to_column('BAM.file') %>% 
  dplyr::mutate(genomescan.sid = gsub("^.+alignments-new/([^/]+)/Aligned.sortedByCoord.+$","\\1",BAM.file))

# merge ----


qc.stats.collapsed <- qc.stats.collapsed %>% 
  dplyr::mutate(freq.alternate.loci = NULL) %>% 
  dplyr::left_join(
    idxstats.stats %>% 
      dplyr::mutate(idxstats = NULL) %>% 
      dplyr::mutate(`Alternate loci` = rowSums(select(., contains("_")))) %>% 
      dplyr::select(!contains("_")) %>% 
      tidyr::pivot_longer(cols = -c(genomescan.sid)) %>% 
      dplyr::filter(name %in% c("chrM","chrEBV","X.") == F) %>% 
      #dplyr::mutate(name = ifelse(grepl("_",name),"Alternate loci",name)) %>% 
      dplyr::group_by(genomescan.sid) %>% 
      dplyr::mutate(n = sum(value)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(freq = value / n * 100) %>%
      dplyr::filter(name == "Alternate loci") %>%
      dplyr::select(genomescan.sid, freq) %>% 
      dplyr::rename(freq.alternate.loci = freq),
    by=c('genomescan.sid'='genomescan.sid')) %>% 
  dplyr::left_join(star.log.final.out.stats , by=c('genomescan.sid'='genomescan.sid')) %>% 
  dplyr::left_join(featureCounts.stats , by=c('genomescan.sid'='genomescan.sid')) %>% 
  dplyr::mutate(featureCounts.M.Assigned = featureCounts.Assigned / 1000000) %>% 
  dplyr::mutate(assigned.reads.status = ifelse(featureCounts.Assigned > 750000,"PASSED","TOO LOW"))



# plots ----

## fastp / sequencing run ----



plt.pct <- qc.stats %>%
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

plt.conf <- qc.stats %>%
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



plt <- rbind(qc.stats, qc.stats %>% dplyr::mutate(fastp.duplication.rate = 0))
p2 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=fastp.duplication.rate, label=genomescan.sid, group=fastp.json)) +
  geom_line() +
  labs(x=NULL, y='dupl rate') +
  facet_grid(cols = vars(facet), scales = "free", space="free") + 
  theme(legend.position = 'bottom', axis.text.x = element_blank())


plt <- rbind(qc.stats, qc.stats %>% dplyr::mutate(fastp.low.complexity.reads = 0))
p3 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=fastp.low.complexity.reads, label=genomescan.sid, group=fastp.json)) +
  geom_line() +
  labs(x=NULL, y='low cplxty') +
  facet_grid(cols = vars(facet), scales = "free", space="free") + 
  theme(legend.position = 'bottom', axis.text.x = element_blank())


plt <- rbind(qc.stats, qc.stats %>% dplyr::mutate(fastp.ratio.reads.not.passed.filtering = 0))
p4 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=fastp.ratio.reads.not.passed.filtering, label=genomescan.sid, group=fastp.json)) +
  geom_line() +
  labs(x=NULL, y='not passed filt') +
  facet_grid(cols = vars(facet), scales = "free", space="free") + 
  theme(legend.position = 'bottom', axis.text.x = element_blank())


plt <- rbind(qc.stats, qc.stats %>% dplyr::mutate(avg.read.trim = 0))
p5 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=avg.read.trim, label=genomescan.sid, group=fastp.json)) +
  geom_line() +
  labs(x=NULL, y='avg read trim') +
  facet_grid(cols = vars(facet), scales = "free", space="free") + 
  theme(legend.position = 'bottom', axis.text.x = element_blank())


# plt <- rbind(qc.stats, qc.stats %>% dplyr::mutate(fastp.insert_size.avg = 0))
# p6 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=fastp.insert_size.avg, label=genomescan.sid, group=fastp.json)) +
#   geom_line() +
#   labs(x=NULL, y='avg ins-size') +
#   facet_grid(cols = vars(facet), scales = "free", space="free") + 
#   theme(legend.position = 'bottom', axis.text.x = element_blank())

plt <- rbind(qc.stats, qc.stats %>% dplyr::mutate(fastp.insert_size.0.ratio = 0))
p7 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=fastp.insert_size.0.ratio, label=genomescan.sid, group=fastp.json)) +
  geom_line() +
  labs(x=NULL, y='ratio ins-size=0') +
  facet_grid(cols = vars(facet), scales = "free", space="free") + 
  theme(legend.position = 'bottom', axis.text.x = element_blank())


plt <- qc.stats %>% 
  dplyr::mutate(y=1)
p8 <- ggplot(plt, aes(x = reorder(fastp.json, order),y=y , col=flow.cell , fill=flow.cell)) +
  geom_tile() +
  facet_grid(cols = vars(facet), scales = "free", space="free")  + 
  labs(x=NULL, y=NULL) +
  theme(legend.position = "none", axis.text.x = element_blank())


p1 / p2 / p3 / p4 / p5 / p7/ p8
# ggsave("output/figures/qc/fastp_qc_metrics.png",width=14,height=7)

## fastp changes / sequencing run ----


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



dists <- data.frame()
for(sid in unique(qc.stats$genomescan.sid)) {
  sel <- qc.stats %>% dplyr::filter(genomescan.sid == sid)
  
  dists <- rbind(dists,
                 data.frame(
                   genomescan.sid = sid,
                   
                   fastp.gc.rmse = max(sel$fastp.gc.rmse) - min(sel$fastp.gc.rmse),
                   fastp.percentage.c = max(sel$fastp.percentage.c) - min(sel$fastp.percentage.c),
                   fastp.duplication.rate = max(sel$fastp.duplication.rate) - min(sel$fastp.duplication.rate),
                   fastp.low.complexity.reads = max(sel$fastp.low.complexity.reads) - min(sel$fastp.low.complexity.reads),
                   fastp.ratio.reads.passed.filtering = max(sel$fastp.ratio.reads.passed.filtering) - min(sel$fastp.ratio.reads.passed.filtering),
                   avg.read.trim = max(sel$avg.read.trim) - min(sel$avg.read.trim),
                   fastp.insert_size.0.ratio = max(sel$fastp.insert_size.0.ratio) - min(sel$fastp.insert_size.0.ratio)
                 ))
}


dists <- dists %>% 
  dplyr::mutate(order = rank(fastp.insert_size.0.ratio))

p0 <- ggplot(dists, aes(x = reorder(genomescan.sid, order),y=fastp.gc.rmse)) + 
  geom_point()
p1 <- ggplot(dists, aes(x = reorder(genomescan.sid, order),y=fastp.percentage.c)) + 
  geom_point()
p2 <- ggplot(dists, aes(x = reorder(genomescan.sid, order),y=fastp.duplication.rate)) + 
  geom_point()
p3 <- ggplot(dists, aes(x = reorder(genomescan.sid, order),y=fastp.low.complexity.reads)) + 
  geom_point()
p4 <- ggplot(dists, aes(x = reorder(genomescan.sid, order),y=fastp.ratio.reads.passed.filtering)) + 
  geom_point()
p5 <- ggplot(dists, aes(x = reorder(genomescan.sid, order),y=avg.read.trim)) + 
  geom_point()
p6 <- ggplot(dists, aes(x = reorder(genomescan.sid, order),y=fastp.insert_size.0.ratio)) + 
  geom_point()


p1 / p2 / p3 / p4 / p5 / p6



library(factoextra)


plt <- qc.stats %>% 
  dplyr::left_join(idxstats.stats %>%
                     dplyr::select(genomescan.sid, freq.alternate.loci),
                   by=c('genomescan.sid'='genomescan.sid') ) %>% 
  dplyr::select(
    c("fastp.json","fastp.low.complexity.reads",
      "fastp.duplication.rate",
      "fastp.percentage.a",
      "fastp.gc.rmse","fastp.ratio.reads.not.passed.filtering",
      'fastp.insert_size.0.ratio',
      "avg.read.trim",
      "freq.alternate.loci")
  ) %>% 
  dplyr::mutate(fastp.json = gsub("data/glass/RNAseq/fastq-clean/","",fastp.json)) %>% 
  dplyr::mutate(fastp.json = gsub("_fastp_log.json","",fastp.json)) %>%
  dplyr::mutate(fastp.json = gsub("[ACTG\\-]{8,}","",fastp.json)) %>% 
  dplyr::mutate(fastp.json = gsub("[_]+","-",fastp.json)) %>% 
  tibble::column_to_rownames('fastp.json')

res.pca <- prcomp(plt, scale = TRUE)

fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
                #,label="var"
)


## idxstats / patient sample - rainbow ----

plt <- idxstats.stats %>% 
  dplyr::mutate(idxstats = NULL) %>% 
  dplyr::mutate(`Alternate loci` = rowSums(select(., contains("_")))) %>% 
  dplyr::select(!contains("_")) %>% 
  tidyr::pivot_longer(cols = -c(genomescan.sid)) %>% 
  dplyr::filter(name %in% c("chrM","chrEBV","X.") == F) %>% 
  #dplyr::mutate(name = ifelse(grepl("_",name),"Alternate loci",name)) %>% 
  dplyr::group_by(genomescan.sid) %>% 
  dplyr::mutate(n = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(freq = value / n * 100) 




plt <- plt %>%
  dplyr::left_join(plt %>% 
                     dplyr::filter(name == "Alternate loci") %>%
                     dplyr::group_by(genomescan.sid) %>% 
                     dplyr::summarise(freq = max(freq)) %>% 
                     dplyr::ungroup() %>% 
                     dplyr::mutate(order = rank(freq, genomescan.sid) ) %>% 
                     dplyr::arrange(order) %>% 
                     dplyr::select(genomescan.sid, order)
                   ,
                   by = c('genomescan.sid'='genomescan.sid')
  )




ggplot(plt, aes(x = reorder(genomescan.sid,order) ,y = freq, fill=name, label=genomescan.sid)) +
  #coord_flip() + 
  geom_bar(stat = "identity", position = "stack",colour="black") + 
  scale_y_continuous(labels = scales::unit_format(unit = "%")) + 
  theme_bw() + 
  theme( axis.title.y = element_text(size = 11) ,
         axis.text.x = element_text(angle = 90, size = 5 ))





## all / patient sample ----



qc.stats.collapsed <- qc.stats.collapsed %>% 
  dplyr::mutate(order = rank(-rank(featureCounts.Assigned))) 





plt <- qc.stats.collapsed %>% 
  dplyr::mutate(featureCounts.Assigned = NULL)
plt <- rbind(plt, plt %>%  dplyr::mutate(featureCounts.M.Assigned = 0))
p1 <- ggplot(plt, aes(x=reorder(genomescan.sid, order),y=featureCounts.M.Assigned,group=genomescan.sid)) +
  geom_line() +
  geom_hline(yintercept = 0.75, col="red",lty=2) +
  geom_hline(yintercept = 20, col="red",lty=2) +
  theme(legend.position = 'bottom', axis.text.x = element_blank()) +
  labs(x=NULL,title=NULL,subtitle=NULL) + 
  facet_grid(cols = vars(assigned.reads.status), scale="free_x", space= "free_x")


plt <- qc.stats.collapsed
plt <- rbind(plt, plt %>%  dplyr::mutate(fastp.gc.rmse = 0))
p2 <- ggplot(plt, aes(x=reorder(genomescan.sid, order),y=fastp.gc.rmse, group=genomescan.sid)) +
  geom_line() +
  theme(legend.position = 'bottom', axis.text.x = element_blank()) +
  labs(x=NULL) + 
  facet_grid(cols = vars(assigned.reads.status), scale="free_x", space= "free_x")


# qc.stats.collapsed %>% dplyr::filter(order < 120) %>% dplyr::pull(fastp.percentage.c) %>%  median
plt <- qc.stats.collapsed
plt <- rbind(plt, plt %>%  dplyr::mutate(fastp.percentage.c = 23.6))
p3 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=fastp.percentage.c, group=genomescan.sid)) +
  geom_hline(yintercept=23.6) +
  geom_line() +
  theme(legend.position = 'bottom', axis.text.x = element_blank()) +
  labs(x=NULL) + 
  facet_grid(cols = vars(assigned.reads.status), scale="free_x", space= "free_x")


plt <- qc.stats.collapsed
plt <- rbind(plt, plt %>%  dplyr::mutate(fastp.duplication.rate = 0))
p4 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=fastp.duplication.rate, group=genomescan.sid)) +
  geom_line() +
  theme(legend.position = 'bottom', axis.text.x = element_blank()) +
  labs(x=NULL) + 
  facet_grid(cols = vars(assigned.reads.status), scale="free_x", space= "free_x")


plt <- qc.stats.collapsed
plt <- rbind(plt, plt %>%  dplyr::mutate(fastp.low.complexity.reads = 0))
p5 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=fastp.low.complexity.reads, group=genomescan.sid)) +
  geom_line() +
  theme(legend.position = 'bottom', axis.text.x = element_blank()) +
  labs(x=NULL) + 
  facet_grid(cols = vars(assigned.reads.status), scale="free_x", space= "free_x")


plt <- qc.stats.collapsed %>% 
  dplyr::mutate(fastp.ratio.reads.failed.filtering = 100 - fastp.ratio.reads.passed.filtering)
plt <- rbind(plt, plt %>%  dplyr::mutate(fastp.ratio.reads.failed.filtering = 0))
p6 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=fastp.ratio.reads.failed.filtering, group=genomescan.sid)) +
  geom_line() +
  theme(legend.position = 'bottom', axis.text.x = element_blank()) +
  labs(x=NULL) + 
  facet_grid(cols = vars(assigned.reads.status), scale="free_x", space= "free_x")


plt <- qc.stats.collapsed
plt <- rbind(plt, plt %>%  dplyr::mutate(avg.read.trim = 0))
p7 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=avg.read.trim, group=genomescan.sid)) +
  geom_line() +
  theme(legend.position = 'bottom', axis.text.x = element_blank()) +
  labs(x=NULL) + 
  facet_grid(cols = vars(assigned.reads.status), scale="free_x", space= "free_x")


plt <- qc.stats.collapsed
plt <- rbind(plt, plt %>%  dplyr::mutate(fastp.insert_size.0.ratio = 0))
p8 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=fastp.insert_size.0.ratio, group=genomescan.sid)) +
  geom_line() +
  theme(legend.position = 'bottom', axis.text.x = element_blank()) +
  labs(x=NULL) + 
  facet_grid(cols = vars(assigned.reads.status), scale="free_x", space= "free_x")


plt <- qc.stats.collapsed
plt <- rbind(plt, plt %>%  dplyr::mutate(freq.alternate.loci = 0))
p9 <- ggplot(plt, aes(x = reorder(genomescan.sid, order),y=freq.alternate.loci, group=genomescan.sid)) +
  geom_line() +
  #theme(legend.position = 'bottom', axis.text.x = element_blank()) +
  labs(x=NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size=6)) + 
  facet_grid(cols = vars(assigned.reads.status), scale="free_x", space= "free_x")



p1 / p2 / p3 / p4 / p5 / p6 / p7 / p8 / p9






plt <- qc.stats.collapsed
plt <- rbind(plt, plt %>%
               dplyr::mutate(fastp.insert_size.0.ratio = 0,
                             featureCounts.M.Assigned = 0,
                             avg.read.trim = 0,
                             fastp.gc.rmse = 0,
                             fastp.percentage.c = 23.6,    # qc.stats.collapsed %>% dplyr::filter(order < 120) %>% dplyr::pull(fastp.percentage.c) %>%  median
                             fastp.low.complexity.reads = 0,
                             fastp.duplication.rate = 0,
                             fastp.ratio.reads.not.passed.filtering = 0,
                             freq.alternate.loci = 0)) %>% 
  dplyr::select(genomescan.sid,order,assigned.reads.status,
                fastp.insert_size.0.ratio, featureCounts.M.Assigned, avg.read.trim, fastp.gc.rmse, fastp.percentage.c, fastp.low.complexity.reads, fastp.duplication.rate, fastp.ratio.reads.not.passed.filtering, freq.alternate.loci) %>% 
  pivot_longer(cols = -c(genomescan.sid,order,assigned.reads.status)) %>% 
  dplyr::mutate(name = gsub( "avg.read.trim" ,'Avg trim',name)) %>% 
  dplyr::mutate(name = gsub( "fastp.duplication.rate" ,'Dupl rate',name)) %>% 
  dplyr::mutate(name = gsub( "fastp.gc.rmse" ,'ACTG RMSE',name)) %>% 
  dplyr::mutate(name = gsub( "fastp.insert_size.0.ratio" ,'Inst size != 0',name)) %>% 
  dplyr::mutate(name = gsub( "fastp.low.complexity.reads" ,'Low cplxty',name)) %>% 
  dplyr::mutate(name = gsub( "fastp.percentage.c" ,'%C',name)) %>% 
  dplyr::mutate(name = gsub( "fastp.ratio.reads.not.passed.filtering" ,'Not pass filter',name)) %>% 
  dplyr::mutate(name = gsub( "featureCounts.M.Assigned" ,'M Assigned',name)) %>% 
  dplyr::mutate(name = gsub( "freq.alternate.loci" ,'rRNA / Alternate loci/',name)) %>% 
  dplyr::mutate(name = factor(name, levels=c("M Assigned", "Avg trim", "ACTG RMSE", "%C", "Low cplxty", "Dupl rate","Inst size != 0",  "Not pass filter", "rRNA / Alternate loci/")))



ggplot(plt, aes(x = reorder(genomescan.sid, order),y=value, group=genomescan.sid)) +
  geom_line() +
  labs(x=NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size=6)) + 
  facet_grid(cols = vars(assigned.reads.status), row=vars(name), scale="free", space= "free_x")





