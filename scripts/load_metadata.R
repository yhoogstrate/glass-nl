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



metadata.glass.per.resection <- read.csv('data/glass/Clinical data/Cleaned/metadata_2022/Samplesheet_GLASS_RNAseq__ALL.csv') %>% 
  dplyr::mutate(institute = gsub("^.+_(.+)_.+$","\\1",GLASS_ID)) %>% 
  dplyr::rename(genomescan.sid = GS_ID) %>% 
  dplyr::mutate(rid = paste0(gsub("^(.+_)[^_]+$","\\1",GLASS_ID),Sample_Name)) %>% 
  dplyr::rename(Exclude.by.Wies.on.complete.pair = Exclude) %>% 
  dplyr::mutate(Sample_Type = case_when(Sample_Type == "I" ~ "initial",
                                        Sample_Type == "R" ~ "recurrent",
                                        T ~ "X")) %>% 
  dplyr::mutate(Sample_Type = factor(Sample_Type, levels=c('initial','recurrent','X')))





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
  dplyr::left_join(tmp, by=c('genomescan.sid'='genomescan.sid'))


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
  dplyr::mutate(alternate.loci = rowSums(dplyr::select(., contains("_")))) %>% 
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



## featurecounts ----


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
  dplyr::left_join(tmp, by=c('genomescan.sid'='genomescan.sid'))



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




## attach per-resection survival ----

# find dates of last event

tmp.1 <- read.csv('data/glass/Clinical data/Cleaned/metadata_2022/Surgery data_GLASS RNAseq.csv')
tmp.1 <- rbind(
  tmp.1 %>%
    dplyr::select(`GLASS_ID` | ends_with("_S1")) %>%
    `colnames<-`(gsub('_S[1-4]$','',colnames(.)))
  ,
  tmp.1 %>%
    dplyr::select(`GLASS_ID` | ends_with("_S2")) %>% 
    `colnames<-`(gsub('_S[1-4]$','',colnames(.))),
  tmp.1 %>%
    dplyr::select(`GLASS_ID` | ends_with("_S3")) %>%
    `colnames<-`(gsub('_S[1-4]$','',colnames(.))),
  tmp.1 %>% dplyr::select(`GLASS_ID` | ends_with("_S4")) %>%
    `colnames<-`(gsub('_S[1-4]$','',colnames(.)))
) %>% 
  tidyr::drop_na(Sample_Name) %>% 
  dplyr::arrange(Sample_Name)
  #dplyr::select(Sample_Name, Date_Surgery, GLASS_ID)


tmp.2 <- read.csv('data/glass/Clinical data/Cleaned/metadata_2022/Survival data_GLASS RNAseq__ALL.csv') %>% 
  dplyr::mutate(Date_of_Diagnosis = as.Date(Date_of_Diagnosis , format = "%Y-%m-%d")) %>% 
  dplyr::mutate(Date_of_Diagnosis = as.Date(Date_of_Death , format = "%Y-%m-%d")) %>% 
  dplyr::select(GLASS_ID, Date_of_Death, Date_Last_Followup) %>% 
  dplyr::mutate(Date_Last_Event = ifelse(is.na(Date_of_Death), Date_Last_Followup , Date_of_Death)) %>% 
  dplyr::mutate(Date_Last_Event.status = ifelse(is.na(Date_of_Death), 0 , 1) ) %>% 
  dplyr::mutate(Date_of_Death = NULL,  Date_Last_Followup = NULL)


stopifnot(tmp.1$GLASS_ID %in% tmp.2$GLASS_ID)


tmp <- tmp.1 %>% dplyr::left_join(tmp.2, by=c('GLASS_ID'='GLASS_ID')) %>% 
  dplyr::mutate(time.resection.until.last.event = difftime(Date_Last_Event,Date_Surgery, units = 'days')) %>% 
  dplyr::rename(status.resection.until.last.event = Date_Last_Event.status) %>% 
  dplyr::select(Sample_Name, time.resection.until.last.event, status.resection.until.last.event)


stopifnot(tmp$time.resection.until.last.event > 0)
rm(tmp.1, tmp.2)



metadata.glass.per.resection <- metadata.glass.per.resection %>%
  dplyr::left_join(tmp, by=c('Sample_Name'='Sample_Name'),suffix = c("", ""))



metadata.glass.per.resection %>%  dplyr::filter(is.na(status.resection.until.last.event))


rm(tmp)


## attach growth annotations ----

tmp.1 <- read.csv('data/glass/Imaging/220314_glass_imaging_data_summary.csv')
tmp.2 <- read.csv('data/glass/Imaging/220314_growth_annotations.csv') %>% 
  dplyr::mutate(SUBJECT=NULL,
                `t2_flair_sign`=NULL,
                growth_pattern=NULL)

stopifnot(tmp.2$X %in% tmp.1$Image_ID)
stopifnot(duplicated(tmp.1$Image_ID) == F)
stopifnot(duplicated(tmp.2$X) == F)


tmp <- tmp.1 %>% 
  dplyr::left_join(tmp.2, by=c('Image_ID'='X')) %>% 
  #dplyr::mutate(Date_Surgery = as.Date(Date_Surgery , format = "%d/%m/%Y")) %>% # not needed for current analysis
  #dplyr::mutate(Date = as.Date(Date , format = "%Y-%m-%d")) %>% # not needed for current analysis
  dplyr::mutate(SUBJECT = NULL, Date_Image = NULL, Sample_ID = NULL, Date_Surgery = NULL, Status = NULL, Date = NULL ) %>% 
  dplyr::rename(imaging.volume.cet = Volume_CET) %>% 
  dplyr::rename(imaging.volume.wt = Volume_WT) %>% 
  dplyr::rename(imaging.growth_pattern = growth_pattern) %>% 
  dplyr::mutate(imaging.growth_pattern = ifelse(imaging.growth_pattern %in% c("","Not sure"), NA, imaging.growth_pattern)) %>% 
  dplyr::rename(imaging.t2_flair_sign = t2_flair_sign) %>% 
  dplyr::mutate(imaging.t2_flair_sign = ifelse(imaging.t2_flair_sign %in% c("","Not sure"), NA, imaging.t2_flair_sign)) %>% 
  dplyr::rename(imaging.comments = COMMENTS) %>% 
  dplyr::rename(imaging.rater = RATER) %>% 
  dplyr::rename(imaging.t2_flair_comments = t2_flair_comments) %>% 
  dplyr::rename(imaging.growth_pattern_comment = growth_pattern_comment) %>% 
  dplyr::select(Sample_Name, imaging.volume.cet, imaging.volume.wt, imaging.growth_pattern, imaging.t2_flair_sign)


rm(tmp.1, tmp.2)


#tmp$Sample_Name[tmp$Sample_Name %in% metadata.glass.per.resection$Sample_Name == F] # some are unique for imaging


metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(tmp, by=c('Sample_Name'='Sample_Name'),suffix = c("", ""))


rm(tmp)


## attach longitudinal transcriptional signatures ----


tmp <- readRDS('cache/transcriptional.signatures.Rds')


metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(tmp, by=c('genomescan.sid'='genomescan.sid'),suffix = c("", ""))


rm(tmp)


# per patient ----


metadata.glass.per.patient <- read.csv('data/glass/Clinical data/Cleaned/metadata_2022/Survival data_GLASS RNAseq__ALL.csv') %>% 
  dplyr::mutate(data_of_birth = NULL) %>%
  # dplyr::mutate(Age_at_Diagnosis = NULL) %>%
  dplyr::mutate(Date_of_Diagnosis = as.Date(Date_of_Diagnosis , format = "%Y-%m-%d")) %>%
  dplyr::mutate(Date_of_Death = as.Date(Date_of_Death , format = "%Y-%m-%d")) %>%
  dplyr::mutate(Date_Last_Followup = as.Date(Date_Last_Followup , format = "%Y-%m-%d")) %>%
  dplyr::mutate(overall.survival = difftime(Date_of_Death , Date_of_Diagnosis, units = 'days')) %>%
  dplyr::mutate(time.until.last.followup = difftime(Date_Last_Followup, Date_of_Diagnosis, units = 'days')) %>% 
  dplyr::mutate(deceased = !is.na(Date_of_Death)) %>% 
  dplyr::mutate(Sample_Name.I = NA, Sample_Name.R = NA, genomescan.sid.I = NA, genomescan.sid.R = NA) %>% 
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_EMCR_148", "Tumor" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_EMCR_160", "Tumor" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_EMCR_162", "Tumor" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_EMCR_137", "Tumor" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_EMCR_110", "Unknown" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_EMCR_113", "Unknown" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_EMCR_118", "Unknown" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_EMCR_156", "Unknown" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_EMCR_174", "Unknown" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_AUMC_002", "Unknown" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_AUMC_018", "Unknown" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_AUMC_020", "Unknown" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_AUMC_022", "Unknown" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_AUMC_029", "Unknown" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_AUMC_017", "Other" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  dplyr::mutate(Cause_of_Death = ifelse(GLASS_ID == "GLNL_UMCU_207", "Other" , Cause_of_Death)) %>% # Following e-mail conversation 29-4-2022
  
  dplyr::mutate(Cause_of_Death = ifelse(Cause_of_Death == "", "Unkown" , Cause_of_Death)) %>% # unknown, but deceased

  dplyr::mutate(overall.survival.event = ifelse(is.na(overall.survival) | Cause_of_Death == "Other",0,1),
              overall.survival = ifelse(is.na(overall.survival),time.until.last.followup,overall.survival))
  


# missing metadata for those 3 patients lacking a matching pair
stopifnot(metadata.glass.per.resection$GLASS_ID %in% metadata.glass.per.patient$GLASS_ID)


for(pid in metadata.glass.per.patient$GLASS_ID) {
  slice <- metadata.glass.per.patient %>% 
    dplyr::filter(GLASS_ID == pid)
  
  r.I <- metadata.glass.per.resection %>% 
    dplyr::filter(GLASS_ID == pid & Sample_Type == "initial" & excluded == F) %>% 
    dplyr::arrange(resection) %>% 
    dplyr::slice_head(n=1)
  
  if(nrow(r.I) > 0) {
    metadata.glass.per.patient <- metadata.glass.per.patient %>% 
      dplyr::mutate(Sample_Name.I = ifelse(GLASS_ID == pid, r.I$Sample_Name, Sample_Name.I)) %>% 
      dplyr::mutate(genomescan.sid.I = ifelse(GLASS_ID == pid, r.I$genomescan.sid, genomescan.sid.I))
  }
  
  r.R <- metadata.glass.per.resection %>% 
    dplyr::filter(GLASS_ID == pid & Sample_Type == "recurrent" & excluded == F) %>% 
    dplyr::arrange(resection) %>% 
    dplyr::slice_tail(n=1)

  if(nrow(r.R) > 0) {
    metadata.glass.per.patient <- metadata.glass.per.patient %>% 
      dplyr::mutate(Sample_Name.R = ifelse(GLASS_ID == pid, r.R$Sample_Name, Sample_Name.R) ) %>% 
      dplyr::mutate(genomescan.sid.R = ifelse(GLASS_ID == pid, r.R$genomescan.sid, genomescan.sid.R))
  }
}
rm(r.I, r.R, pid, slice)



metadata.glass.per.patient <- metadata.glass.per.patient %>% 
  dplyr::mutate(pair.status = case_when(
    !is.na(Sample_Name.I) & !is.na(Sample_Name.R) ~ "complete",
    is.na(Sample_Name.I) & !is.na(Sample_Name.R) ~ "missing initial",
    !is.na(Sample_Name.I) & is.na(Sample_Name.R) ~ "missing recurrent",
    T ~ "excluded"
  )) %>% 
  dplyr::mutate(pair.status = factor(pair.status, levels=c("complete","missing recurrent","missing initial" ,"excluded"))) %>% 
  dplyr::mutate(patient.correction.id = ifelse(pair.status == "complete", GLASS_ID, "remaining.individual.samples")) %>% 
  dplyr::mutate(excluded = ifelse(pair.status == "excluded",T,F))



tmp <- read.csv('data/glass/Clinical data/Cleaned/metadata_2022/Surgery data_GLASS RNAseq.csv')
tmp <- rbind(
  tmp %>% dplyr::select(Date_Surgery_S1, Sample_Name_S1, GS_ID_S1) %>%
    dplyr::rename(Date_Surgery = Date_Surgery_S1, Sample_Name = Sample_Name_S1, genomescan.sid = GS_ID_S1),
  tmp %>% dplyr::select(Date_Surgery_S2, Sample_Name_S2, GS_ID_S2) %>% 
    dplyr::rename(Date_Surgery = Date_Surgery_S2, Sample_Name = Sample_Name_S2, genomescan.sid = GS_ID_S2),
  tmp %>% dplyr::select(Date_Surgery_S3, Sample_Name_S3, GS_ID_S3) %>% 
    dplyr::rename(Date_Surgery = Date_Surgery_S3, Sample_Name = Sample_Name_S3, genomescan.sid = GS_ID_S3),
  tmp %>% dplyr::select(Date_Surgery_S4, Sample_Name_S4, GS_ID_S4) %>% 
    dplyr::rename(Date_Surgery = Date_Surgery_S4, Sample_Name = Sample_Name_S4, genomescan.sid = GS_ID_S4)
  ) %>% 
  dplyr::filter(!is.na(Sample_Name) & !is.na(genomescan.sid)) %>% 
  dplyr::mutate(Date_Surgery = as.Date(Date_Surgery , format = "%Y-%m-%d")) %>% 
  dplyr::mutate(Sample_Name = NULL)

metadata.glass.per.patient <- metadata.glass.per.patient %>% 
  dplyr::left_join(tmp %>% dplyr::rename(Date_Surgery.I = Date_Surgery), by=c('genomescan.sid.I' = 'genomescan.sid')) %>% 
  dplyr::left_join(tmp %>% dplyr::rename(Date_Surgery.R = Date_Surgery), by=c('genomescan.sid.R' = 'genomescan.sid')) %>% 
  dplyr::mutate(survival.I = difftime(Date_of_Death , Date_Surgery.I, units = 'days')) %>%
  dplyr::mutate(survival.R = difftime(Date_of_Death , Date_Surgery.R, units = 'days')) 

#  104059-002-121  no date inital
#  104059-002-015  no date recurrent
#  104059-003-014  no date recurrent

rm(tmp)



# remove privacy sensitive information
metadata.glass.per.patient <- metadata.glass.per.patient %>% dplyr::mutate(
  Date_Last_Followup = NULL,
  Date_of_Birth = NULL,
  Date_of_Death = NULL,
  Date_of_Diagnosis = NULL
)


