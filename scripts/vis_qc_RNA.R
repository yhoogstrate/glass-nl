#!/usr/bin/env R

# visualise all RNA QC stats


library(tidyverse)
library(patchwork)
library(rjson)

source('scripts/R/youri_gg_theme.R')
source('scripts/R/job_gg_theme.R')

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



plt <- idxstats.stats %>% 
  dplyr::mutate(idxstats = NULL) %>% 
  tidyr::pivot_longer(cols = -c(genomescan.sid)) %>% 
  dplyr::filter(name %in% c("chrM","chrEBV","X.") == F) %>% 
  dplyr::mutate(name = ifelse(grepl("_",name),"Alternate loci",name)) %>% 
  dplyr::group_by(genomescan.sid) %>% 
  dplyr::mutate(n = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(freq = value / n * 100) 


plt %>%
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
      (25.0 - out$fastp.percentage.a)^2 +
      (25.0 - out$fastp.percentage.c)^2 +
      (25.0 - out$fastp.percentage.t)^2 +
      (25.0 - out$fastp.percentage.g)^2
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
  mutate(tmp = map(json.stats, ~ data.frame(t(.)))) %>%  unnest(tmp) %>% dplyr::mutate(json.stats = NULL) %>% as.data.frame %>% 
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
    genomescan.sid %in% c('104059-003-002') ~ "104059-003-002", # rare qc
    genomescan.sid %in% c('104059-002-099') ~ "104059-002-099",
    genomescan.sid %in% c('104059-002-087') ~ "104059-002-087",
                T ~ "other"
                ))






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



library(factoextra)


plt <- qc.stats %>% 
  dplyr::select(
    c("fastp.json","fastp.low.complexity.reads",
    "fastp.duplication.rate",
    "fastp.percentage.a",
    "fastp.gc.rmse","fastp.ratio.reads.not.passed.filtering",
    'fastp.insert_size.0.ratio',
    "avg.read.trim")
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




