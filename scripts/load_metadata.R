#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(rjson)


# per fastq file ----

## fastp json ----

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



metadata.glass.per.fastq <- data.frame(fastp.json = Sys.glob("data/glass/RNAseq/fastq-clean/*.json")) %>% 
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
  dplyr::mutate(fastp.avg.read.trim = 150 - (fastp.read1_mean_length + fastp.read2_mean_length)/ 2) %>% 
  dplyr::mutate(flow.cell = as.factor(gsub("^.+/([^_]+)_.+$","\\1",fastp.json))) %>% 
  dplyr::mutate(lane = gsub("^.+_(L[0-9]+)_.+$","\\1",fastp.json )) %>% 
  dplyr::mutate(lane.flow.cell = as.factor( paste0(flow.cell , "_", lane ))) 


rm(parse_fastp_json_files)



# per resection ----


# metadata.glass.per.resection <- read.csv('output/tables/metadata/Samplesheet_GLASS.txt',header=T) %>% 
#   dplyr::filter(group == "Glioma") %>% 
#   dplyr::mutate(group=NULL) %>% 
#   dplyr::rename(genomescan.sid = GS_ID) %>% 
#   dplyr::mutate(pid = as.factor(gsub("^([0-9]+)_.+$","g-nl.\\1",GLASS_ID))) %>% 
#   dplyr::mutate(resection = gsub("^.+_(.+)$","\\1",GLASS_ID))

#a = unique(metadata.glass.per.resection$genomescan.sid)
# a = sort(unique(gsub("^.+/","",Sys.glob("data/glass/RNAseq/alignments/alignments-new/10405*"))))
# b = unique(tmp$GS_ID)
# length(a)
# length(b)

# a[a %in% b == F]
# b[b %in% a == F]


metadata.glass.per.resection <- read.csv('data/glass/RNAseq/Metadata/Samplesheet_GLASS_RNAseq__ALL.csv') %>% 
  dplyr::mutate(institute = gsub("^.+_(.+)_.+$","\\1",GLASS_ID)) %>% 
  dplyr::rename(genomescan.sid = GS_ID) %>% 
  dplyr::mutate(rid = paste0(gsub("^(.+_)[^_]+$","\\1",GLASS_ID),Sample_Name)) %>% 
  dplyr::rename(Exclude.by.Wies.on.complete.pair = Exclude)





## aggregated per fastq qc stats ----

tmp <- metadata.glass.per.fastq %>% 
  dplyr::mutate(fastp.json = NULL, order=NULL) %>% 
  dplyr::group_by(genomescan.sid) %>% 
  dplyr::summarise(
    fastp.avg.read.trim = weighted.mean(fastp.avg.read.trim, fastp.total_reads),
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


metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(
    tmp, by=c('genomescan.sid'='genomescan.sid')
  )


stopifnot(is.na(metadata.glass.per.resection$fastp.total_reads) == F) # ENSURE ALL SAMPLE HAVE THIS METADATA

rm(tmp)


## star log stats ----

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




tmp <- data.frame(star.log.final.out = Sys.glob("data/glass/RNAseq/alignments/alignments-new/*/Log.final.out")) %>% 
  dplyr::mutate(genomescan.sid =  gsub("^.+new/([^/]+)/Log.+$","\\1", star.log.final.out)) %>%
  dplyr::arrange(genomescan.sid) %>% 
  dplyr::mutate(stats = lapply(star.log.final.out, parse_star_log_final_out)) %>% 
  mutate(tmp = map(stats, ~ data.frame(t(.)))) %>%
  tidyr::unnest(tmp) %>%
  dplyr::mutate(stats = NULL) %>%
  as.data.frame %>% 
  dplyr::mutate_if(colnames(.) %in% c('star.log.final.out', 'genomescan.sid') == F , as.numeric) %>% 
  dplyr::mutate(star.log.final.out = NULL)



metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(
    tmp, by=c('genomescan.sid'='genomescan.sid')
  )


stopifnot(is.na(metadata.glass.per.resection$star.input.reads) == F)


rm(tmp, parse_star_log_final_out)



## idxstats stats  ----


parse_idxstats <- function(idxstats_file) {
  tmp <- read.delim(idxstats_file,header=F) %>% 
    dplyr::mutate(V2 = NULL, V4 = NULL )
  
  out <- tmp$V3
  names(out) <- tmp$V1
  out <- as.list(out)
  
  return(out)
}


tmp <- data.frame(idxstats = Sys.glob("output/tables/qc/idxstats/*.txt")) %>% 
  dplyr::mutate(genomescan.sid =  gsub("^.+/([^/]+).samtools.idxstats.txt$","\\1", idxstats)   ) %>%
  dplyr::arrange(genomescan.sid) %>% 
  dplyr::mutate(stats = lapply(idxstats, parse_idxstats)) %>% 
  mutate(tmp = map(stats, ~ data.frame(t(.)))) %>%
  tidyr::unnest(tmp) %>%
  dplyr::mutate(stats = NULL) %>%
  as.data.frame %>% 
  dplyr::mutate_if(colnames(.) %in% c('idxstats', 'genomescan.sid') == F , as.numeric) %>% 
  dplyr::mutate(idxstats = NULL, `X.` = NULL)  %>% 
  tibble::column_to_rownames('genomescan.sid') %>% 
  dplyr::mutate(total = rowSums(.))  %>% 
  dplyr::mutate(alternate.loci = rowSums(select(., contains("_")))) %>% 
  dplyr::mutate(freq.alternate.loci = alternate.loci / total ) %>% 
  dplyr::select(!contains("_")) %>% 
  dplyr::select(-c('chrM', 'chrEBV'))  %>% 
  `colnames<-`(paste0("idxstats.",colnames(.))) %>% 
  tibble::rownames_to_column('genomescan.sid')



metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(
    tmp, by=c('genomescan.sid'='genomescan.sid')
  )



stopifnot(is.na(metadata.glass.per.resection$idxstats.alternate.loci) == F)


rm(tmp, parse_idxstats)


##  featurecounts ----


tmp <- read.delim("data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt.summary",skip=0,header=T,check.names=F) %>% 
  tibble::column_to_rownames('Status') %>% 
  t() %>% 
  as.data.frame %>% 
  `colnames<-`(paste0("featureCounts.",colnames(.))) %>% 
  tibble::rownames_to_column('BAM.file') %>% 
  dplyr::mutate(genomescan.sid = gsub("^.+alignments-new/([^/]+)/Aligned.sortedByCoord.+$","\\1",BAM.file)) %>% 
  dplyr::mutate(BAM.file = NULL ,
                featureCounts.Unassigned_Unmapped = NULL,
                featureCounts.Unassigned_Read_Type = NULL,
                featureCounts.Unassigned_Chimera = NULL,
                featureCounts.Unassigned_Singleton = NULL,
                featureCounts.Unassigned_MappingQuality = NULL,
                featureCounts.Unassigned_FragmentLength = NULL,
                featureCounts.Unassigned_Secondary = NULL,
                featureCounts.Unassigned_NonSplit = NULL,
                featureCounts.Unassigned_Overlapping_Length = NULL
                ) %>% 
  dplyr::mutate(featureCounts.M.Assigned = featureCounts.Assigned / 1000000)


metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(
    tmp, by=c('genomescan.sid'='genomescan.sid')
  )



stopifnot(is.na(metadata.glass.per.resection$featureCounts.Assigned) == F)


rm(tmp)


metadata.glass.per.resection <- metadata.glass.per.resection%>% 
  dplyr::mutate(assigned.reads.status = factor(
    ifelse(featureCounts.Assigned > 750000,"PASS","INSUFFICIENT"),
    levels=c("PASS","INSUFFICIENT")))






## classify sample as incl/excl ----



metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::mutate(excluded.reason = case_when(
    assigned.reads.status == "INSUFFICIENT" ~ "featureCounts.Assigned",
    T ~ '-' # keep all those that are not excluded
  )) %>% 
  dplyr::mutate(excluded = ifelse(excluded.reason == "-", FALSE, TRUE))





# per patient ----s

# metadata.glass.per.patient <- 
#   read.csv('data/glass/Clinical data/Cleaned/Survival data_GLASS_12082021.csv') %>% 
#   dplyr::mutate(X=NULL) %>% 
#   dplyr::mutate(data_of_birth = NULL) %>% 
#   dplyr::mutate(Date_of_Death = as.Date(Date_of_Death , format = "%d-%m-%Y")) %>% 
#   dplyr::mutate(Date_of_Diagnosis = as.Date(Date_of_Diagnosis , format = "%d-%m-%Y")) %>% 
#   dplyr::mutate(Date_Last_Followup = as.Date(Date_Last_Followup , format = "%d-%m-%Y")) %>% 
#   dplyr::mutate(overall.survival = difftime(Date_of_Death , Date_of_Diagnosis, units = 'days')) %>% 
#   dplyr::mutate(time.until.last.followup = difftime(Date_Last_Followup, Date_of_Diagnosis, units = 'days'))
# 




