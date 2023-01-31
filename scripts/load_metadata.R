#!/usr/bin/env R

# load libs ----


library(base)
library(utils)
library(tidyverse)
library(rjson)

library(MASS)
library(fitdistrplus)



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
  
  return(as.data.frame(out)) # as.data.frame for simple unnesting :)
}



metadata.glass.per.fastq <- data.frame(fastp.json = Sys.glob("data/glass/RNAseq/fastq-clean/*.json")) %>% 
  dplyr::mutate(genomescan.sid = gsub("^.+/[^_]+_([^_]+)_.+$","\\1", fastp.json)) %>% 
  dplyr::mutate(json.stats = pbapply::pblapply(fastp.json, parse_fastp_json_files)) %>% 
  tidyr::unnest(json.stats) %>%
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

## load RNA ----

metadata.glass.per.resection <- read.csv('data/glass/Clinical data/Cleaned/metadata_2022/Samplesheet_GLASS_RNAseq__ALL.csv') |> 
  dplyr::mutate(institute = gsub("^.+_(.+)_.+$","\\1",GLASS_ID)) |> 
  dplyr::rename(genomescan.sid = GS_ID) |> 
  dplyr::mutate(rid = paste0(gsub("^(.+_)[^_]+$","\\1",GLASS_ID),Sample_Name)) |> 
  dplyr::rename(Exclude.by.Wies.on.complete.pair = Exclude) |> 
  dplyr::mutate(Sample_Type = case_when(Sample_Type == "I" ~ "initial",
                                        Sample_Type == "R" ~ "recurrent",
                                        T ~ as.character(NA))) |> 
  dplyr::mutate(Sample_Type = factor(Sample_Type, levels=c('initial','recurrent','X'))) |> 
  dplyr::mutate(Customer_ID = NULL) # horribly confusing set of identifiers, some of which match Sample_Name but on different samples


## add proteomics ID's ----


# column does not contain the IDs from the normalised data
tmp <- readxl::read_xlsx('data/glass/Proteomics/Annotation_Reduced_withControls__sample_swap_146_fixed.xlsx') |> 
  dplyr::select(`File_Name_Proteomics`, `Sample_Name`) |> 
  dplyr::filter(grepl("Control",Sample_Name)==F) |> 
  dplyr::left_join(
    read.csv("data/glass/Proteomics/Proteomics_SampleSheet_03112022.csv") |> 
      dplyr::select(Sample_Name, ProtID)
    , by=c('Sample_Name'='Sample_Name'), suffix=c('','')) |> 
  dplyr::filter(!is.na(ProtID)) # not in final results table anyway
# 146_R3 os missing
stopifnot(nrow(tmp) == 55)
stopifnot("146_R1" %in% tmp$Sample_Name == F) # if this file is found in the metadata, the metadata is old, invalid and results in a sample swap
stopifnot("146_R2" %in% tmp$Sample_Name)
stopifnot("146_R3" %in% tmp$Sample_Name)
stopifnot(nrow(tmp) == 55)


metadata.glass.per.resection <- metadata.glass.per.resection |> 
  dplyr::full_join(tmp, by=c('Sample_Name'='Sample_Name'), suffix=c('',''))


metadata.glass.per.resection <- metadata.glass.per.resection |>
  dplyr::mutate(Sample_Type = as.character(Sample_Type)) |>
  dplyr::mutate(Sample_Type = ifelse(is.na(Sample_Type) & grepl("_P$", ProtID), "initial", Sample_Type)) |>
  dplyr::mutate(Sample_Type = ifelse(is.na(`Sample_Type`) & grepl("_R[1-4]$", ProtID), "recurrent", Sample_Type))  |>
  dplyr::mutate(GLASS_ID = ifelse(is.na(GLASS_ID) & !is.na(ProtID), paste0("GLNL_EMCR_",gsub("_.+$","",Sample_Name)), GLASS_ID)) |>
  dplyr::mutate(institute = ifelse(is.na(institute) & !is.na(ProtID), "EMCR", institute))


stopifnot(!is.na(metadata.glass.per.resection$GLASS_ID)) # all must have patient identifier
stopifnot(!is.na(metadata.glass.per.resection |> dplyr::filter(!is.na(ProtID)) |>  dplyr::pull(Sample_Type))) # all protein samples must have initial/recurrent status
stopifnot(metadata.glass.per.resection |> dplyr::filter(!is.na(ProtID)) |>  dplyr::pull(institute) == "EMCR") # all protein samples are from EMC


# ensure no patients with 3 or more resection appear in this data
stopifnot(metadata.glass.per.resection |> 
            dplyr::filter(!is.na(ProtID)) |> 
            dplyr::pull(GLASS_ID) |>
            table() |> 
            max() == 2)




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
  dplyr::left_join(tmp, by=c('genomescan.sid'='genomescan.sid'), suffix = c("", ""))


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

  return(as.data.frame(out))
}




tmp <- data.frame(star.log.final.out = Sys.glob("data/glass/RNAseq/alignments/alignments-new/*/Log.final.out")) |>  
  dplyr::mutate(genomescan.sid =  gsub("^.+new/([^/]+)/Log.+$","\\1", star.log.final.out)) |> 
  dplyr::arrange(genomescan.sid) |> 
  dplyr::mutate(stats = pbapply::pblapply(star.log.final.out, parse_star_log_final_out)) |> 
  tidyr::unnest(stats) |> 
  dplyr::mutate(star.log.final.out = NULL)



metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(tmp, by=c('genomescan.sid'='genomescan.sid'), suffix=c("", ""))


stopifnot(is.na(metadata.glass.per.resection$star.input.reads) == F)


rm(tmp, parse_star_log_final_out)



## idxstats stats  ----


parse_idxstats <- function(idxstats_file) {
  tmp <- read.delim(idxstats_file,header=F) %>% 
    dplyr::mutate(V2 = NULL, V4 = NULL )
  
  out <- tmp$V3
  names(out) <- tmp$V1
  out <- as.list(out)
  
  return(as.data.frame(out))
}


tmp <- data.frame(idxstats = Sys.glob("output/tables/qc/idxstats/*.txt")) %>% 
  dplyr::mutate(genomescan.sid =  gsub("^.+/([^/]+).samtools.idxstats.txt$","\\1", idxstats)   ) %>%
  dplyr::arrange(genomescan.sid) %>% 
  dplyr::mutate(stats = pbapply::pblapply(idxstats, parse_idxstats)) %>% 
  #mutate(tmp = map(stats, ~ data.frame(t(.)))) %>%
  tidyr::unnest(stats) %>%
  #dplyr::mutate(stats = NULL) %>%
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
  dplyr::left_join(    tmp, by=c('genomescan.sid'='genomescan.sid') , suffix=c("", "")  )


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
  dplyr::left_join(tmp, by=c('genomescan.sid'='genomescan.sid'), suffix=c("", ""))



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



#metadata.glass.per.resection %>%  dplyr::filter(is.na(status.resection.until.last.event))


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
  dplyr::left_join(tmp.2, by=c('Image_ID'='X'),suffix=c('','')) %>% 
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

# up.1 = cell cycling
# up.2 = COL/pericyte + CD248


tmp <- readRDS('cache/transcriptional.signatures.Rds')

metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(tmp, by=c('genomescan.sid'='genomescan.sid'),suffix = c("", ""))

rm(tmp)


### transform lts.up1 from gamma to normal ----


fit.data  <- metadata.glass.per.resection |>  
  dplyr::filter(!is.na(lts.up1)) |>  
  dplyr::pull(lts.up1,name='Sample_Name')
fit.data <- fit.data + abs(min(fit.data)) + 0.01

fit.g = fitdistrplus::fitdist(fit.data , "gamma")
fit.ln = fitdistrplus::fitdist(fit.data, "lnorm")
fit.wb = fitdistrplus::fitdist(fit.data, "weibull")

#print(fit.g)

par(mfrow = c(2, 2))
plot.legend <- c("gamma", "lognormal", "weibull")
denscomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)
qqcomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)
cdfcomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)
ppcomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)

dev.off()

# seems gamma fit
fit <- fit.g
rm(fit.g,fit.ln,fit.wb)


metadata.glass.per.resection <- metadata.glass.per.resection |> 
  dplyr::left_join(
    data.frame(lts.up1.norm = qnorm(pgamma(fit.data, shape = fit$estimate['shape'], rate = fit$estimate['rate'] ))) |> 
      tibble::rownames_to_column('Sample_Name'),
    by=c('Sample_Name'='Sample_Name'),suffix = c('',''))


# Transformation, no re-ordering
stopifnot(order(metadata.glass.per.resection$lts.up1) == order(metadata.glass.per.resection$lts.up1.norm))

plot(metadata.glass.per.resection$lts.up1 , metadata.glass.per.resection$lts.up1.norm)


rm(fit, fit.data)



metadata.glass.per.resection <- metadata.glass.per.resection |> 
  dplyr::mutate(lts.up1.norm = ifelse(lts.up1.norm < -3, NA, lts.up1.norm)) #' there is one excessive outlier, troubling our statistics; plot(sort(data$lts.up1.norm)) %>% 


#' there is one excessive outlier, troubling our statistics;
metadata.glass.per.resection <- metadata.glass.per.resection |> 
  dplyr::mutate(lts.up1.norm = ifelse(lts.up1.norm < -3, NA, lts.up1.norm))  

plot(sort(metadata.glass.per.resection$lts.up1.norm))



### transform lts.up2 from gamma to normal ----


fit.data  <- metadata.glass.per.resection |>  
  dplyr::filter(!is.na(lts.up2)) |>  
  dplyr::pull(lts.up2,name='Sample_Name')
fit.data <- fit.data + abs(min(fit.data)) + 0.01

fit.g = fitdistrplus::fitdist(fit.data , "gamma")
fit.ln = fitdistrplus::fitdist(fit.data, "lnorm")
fit.wb = fitdistrplus::fitdist(fit.data, "weibull")

#print(fit.g)

par(mfrow = c(2, 2))
plot.legend <- c("gamma", "lognormal", "weibull")
denscomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)
qqcomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)
cdfcomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)
ppcomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)

dev.off()

# seems gamma fit
fit <- fit.g
rm(fit.g,fit.ln,fit.wb)


metadata.glass.per.resection <- metadata.glass.per.resection |> 
  dplyr::left_join(
    data.frame(lts.up2.norm = qnorm(pgamma(fit.data, shape = fit$estimate['shape'], rate = fit$estimate['rate'] ))) |> 
      tibble::rownames_to_column('Sample_Name'),
    by=c('Sample_Name'='Sample_Name'),suffix = c('',''))


# Transformation, no re-ordering
stopifnot(order(metadata.glass.per.resection$lts.up2) == order(metadata.glass.per.resection$lts.up2.norm))

plot(metadata.glass.per.resection$lts.up1 , metadata.glass.per.resection$lts.up1.norm)


rm(fit, fit.data)


#' there is one excessive outlier, troubling our statistics;
metadata.glass.per.resection <- metadata.glass.per.resection |> 
  dplyr::mutate(lts.up2.norm = ifelse(lts.up2.norm < -3, NA, lts.up2.norm))  

plot(sort(metadata.glass.per.resection$lts.up2.norm))



## attach methylation identifiers ----
# sid mismatch: 204808700074_R07C01 - geen heidelberg folder van ge-upload


tmp.1 <- read.csv('data/glass/Methylation/Metadata/Datasheet4.csv') %>%
  dplyr::mutate(X=NULL) %>% 
  dplyr::mutate(GLASS_ID = NULL,
                Surgery_ID = NULL,
                Sample_Sex = NULL,
                Recurrent_Type = NULL,
                Sample_Type = NULL,
                Sample_Resection = NULL)

  
tmp.2 <- data.frame(Heidelberg.segment.file = Sys.glob("data/glass/Methylation/Heidelberg/Heidelberg_unzip/*/cnvp_v3.0/*.segments.seg")) %>%
  dplyr::mutate(Sample_ID = gsub("^.+_unzip/([^/]+)_Run.+$","\\1",Heidelberg.segment.file))


tmp.3 <- read.csv('data/glass/Clinical data/(Epi)genetic data methylation/(Epi)genetic data_GLASS-NL_01092021.csv') %>% 
  dplyr::mutate(X=NULL)


parse_predictbrain_csv <- function(file, suffix) {
  # file = 'data/glass/Methylation/Heidelberg/Heidelberg_unzip/203989100107_R01C01_Run_78486/predictBrain_v2.1/203989100107_R01C01_scores.csv'
  tmp <- read.csv(file)
  colnames(tmp)[2] <- 'score'
  
  out <- list(
    A_IDH =    tmp %>% dplyr::filter(grepl('^A_IDH$',X))    %>% dplyr::pull(score) %>% as.numeric,
    A_IDH_HG = tmp %>% dplyr::filter(grepl('^A_IDH_HG$',X)) %>% dplyr::pull(score) %>% as.numeric,
    O_IDH =    tmp %>% dplyr::filter(grepl('^O_IDH$',X))    %>% dplyr::pull(score) %>% as.numeric,
    GBM_MES =  tmp %>% dplyr::filter(grepl('^GBM_MES$',X))  %>% dplyr::pull(score) %>% as.numeric
  )
  
  names(out) <- paste0(names(out), suffix)
  
  return(as.data.frame(out))
}

tmp.4 <- data.frame(predictBrain.scores.file = Sys.glob("data/glass/Methylation/Heidelberg/Heidelberg_unzip/*/predictBrain_v2.1/*_scores.csv")) %>%
  dplyr::mutate(Sample_ID = gsub("^.+_v2.1/([^/]+)_scores.csv$","\\1", predictBrain.scores.file)) %>% 
  dplyr::mutate(stats = pbapply::pblapply(predictBrain.scores.file, parse_predictbrain_csv, suffix='')) %>% 
  #mutate(tmp = map(stats, ~ data.frame(t(.)))) %>%
  tidyr::unnest(stats) %>% 
  dplyr::mutate(predictBrain.scores.file = NULL) %>% 
  tibble::column_to_rownames('Sample_ID') %>% 
  dplyr::mutate_all(as.numeric) %>% 
  tibble::rownames_to_column('Sample_ID')


tmp.5 <- data.frame(predictBrain.scores.file = Sys.glob("data/glass/Methylation/Heidelberg/Heidelberg_unzip/*/predictBrain_v2.1/*_scores_cal.csv")) %>%
  dplyr::mutate(Sample_ID = gsub("^.+_v2.1/([^/]+)_scores_cal.csv$","\\1", predictBrain.scores.file)) %>% 
  dplyr::mutate(stats = pbapply::pblapply(predictBrain.scores.file, parse_predictbrain_csv, suffix='_cal')) %>% 
  #mutate(tmp = map(stats, ~ data.frame(t(.)))) %>%
  tidyr::unnest(stats) %>% 
  dplyr::mutate(predictBrain.scores.file = NULL) %>% 
  tibble::column_to_rownames('Sample_ID') %>% 
  dplyr::mutate_all(as.numeric) %>% 
  tibble::rownames_to_column('Sample_ID') |> 
  dplyr::mutate(IDH_HG_IDH_ratio = log(A_IDH_HG_cal/A_IDH_cal))



stopifnot(duplicated(tmp.1$Sample_ID) == FALSE)
stopifnot(duplicated(tmp.2$Sample_ID) == FALSE)
stopifnot(duplicated(tmp.3$Sample_ID) == FALSE)
stopifnot(tmp.1$Sample_ID %in% tmp.3$Sample_ID)
stopifnot(tmp.3$Sample_ID %in% tmp.1$Sample_ID)


tmp <- tmp.1 %>% 
  dplyr::left_join(tmp.2, by=c('Sample_ID'='Sample_ID')) %>% 
  dplyr::left_join(tmp.3, by=c('Sample_ID'='Sample_ID')) %>% 
  dplyr::left_join(tmp.4, by=c('Sample_ID'='Sample_ID')) %>% 
  dplyr::left_join(tmp.5, by=c('Sample_ID'='Sample_ID')) %>% 
  dplyr::rename(methylation.sid = Sample_ID) %>% 
  dplyr::rename(methylation.Sample_Plate = Sample_Plate) %>% 
  dplyr::rename(methylation.Basename = Basename) %>% 
  dplyr::rename(methylation.Array = Array) %>% 
  dplyr::rename(methylation.Slide = Slide) %>% 
  dplyr::rename(methylation.sub.diagnosis = sub.diagnosis)


stopifnot(sum(is.na(tmp$Heidelberg.segment.file)) <= 1) # one file missing so far, that's known


rm(tmp.1,tmp.2,tmp.3,tmp.4,tmp.5, parse_predictbrain_csv)


#dim(tmp)
#dim(metadata.glass.per.resection)

# 
#tmp$Sample_Name %in% metadata.glass.per.resection$Sample_Name
#metadata.glass.per.resection$Sample_Name %in% tmp$Sample_Name
#metadata.glass.per.resection %>%  dplyr::filter(Sample_Name %in% tmp$Sample_Name == F )

# een hoop waarvan wel methylering en metadata is niet in glass metadata tabel
# drie waarvan wel metadata is geen methylaring


metadata.glass.per.resection <- metadata.glass.per.resection %>%
  dplyr::left_join(tmp, by=c('Sample_Name'='Sample_Name'))

rm(tmp)


### transform IDH_HG_IDH_ratio from gamma to normal ----


fit.data  <- metadata.glass.per.resection |>  
  dplyr::filter(!is.na(IDH_HG_IDH_ratio)) |>  
  dplyr::pull(IDH_HG_IDH_ratio,name='Sample_Name')
fit.data <- fit.data + abs(min(fit.data)) + 1

fit.g = fitdistrplus::fitdist(fit.data , "gamma")
fit.ln = fitdistrplus::fitdist(fit.data, "lnorm")
fit.wb = fitdistrplus::fitdist(fit.data, "weibull")

par(mfrow = c(2, 2))
plot.legend <- c("gamma", "lognormal", "weibull")
denscomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)
qqcomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)
cdfcomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)
ppcomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)

dev.off()

# seems gamma fit
fit <- fit.g
rm(fit.g,fit.ln,fit.wb,plot.legend)


metadata.glass.per.resection <- metadata.glass.per.resection |> 
  dplyr::left_join(
    data.frame(IDH_HG_IDH_ratio.norm = qnorm(pgamma(fit.data, shape = fit$estimate['shape'], rate = fit$estimate['rate'] ))) |> 
      tibble::rownames_to_column('Sample_Name'),
    by=c('Sample_Name'='Sample_Name'),suffix = c('',''))


# Transformation, no re-ordering
stopifnot(order(metadata.glass.per.resection$IDH_HG_IDH_ratio) == order(metadata.glass.per.resection$IDH_HG_IDH_ratio.norm))

rm(fit, fit.data)





## attach WHO classification ----

# deze file komt uit de Methylation folder - en is afwezig in Clinical
# de files voor Chemo en Radio in deze folders zijn bijv anders

tmp <- read.csv('data/glass/Methylation/Metadata/WHOclassification_03052022.csv') %>% 
  dplyr::mutate(
    Surgery_ID = NULL,
    Surgery_date = NULL,
    Sample_Resection = NULL,
    Sample_Type = NULL,
    Recurrent_Type = NULL,
    GLASS_ID = NULL
  ) %>% 
  dplyr::filter(!is.na(WHO_Classification2021))


metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(tmp, by=c('Sample_Name'='Sample_Name'),suffix = c("", ""))|>
  dplyr::mutate(WHO_Classification2021 = ifelse(Sample_Name == "153_R2", "Astrocytoma, IDH-mutant, WHO grade 4", WHO_Classification2021))


# ensure R153_R2 is WHO grade IV (not in metadata file, but confirmed by Wies)
stopifnot((metadata.glass.per.resection |> dplyr::filter(Sample_Name == "153_R2") |> dplyr::pull(WHO_Classification2021)) == "Astrocytoma, IDH-mutant, WHO grade 4")


rm(tmp)


## attach Chemo+radio therapy  ----

# volgens Iris zijn behandelings in deze file goed geprocessed:

# alkalating.agent = PCV, CCNU | TMZ

tmp <- read.csv('data/glass/Clinical data/Cleaned/metadata_2022/Chemotherapy data_GLASS RNAseq.csv') %>% 
  dplyr::select(-contains("Date_")) %>% 
  dplyr::select(-contains("KPS_")) %>% 
  dplyr::select(-contains("_Stopped_")) %>% 
  dplyr::select(-contains("WHO_")) %>% 
  dplyr::select(-contains("Multiple_")) %>% 
  dplyr::mutate(Pre_Initial_cat = paste0(Treatment_Window_Pre_Initial , ": ", Chemotherapy_Type_Pre_Initial)) %>% 
  dplyr::mutate(Treatment_Window_Pre_Initial = NULL) %>% 
  dplyr::mutate(Chemotherapy_Type_Pre_Initial = NULL) %>% 
  dplyr::mutate(Post_Initial.1_cat = paste0(Treatment_Window_Post_Initial.1 , ": ", Chemotherapy_Type_Post_Initial.1)) %>% 
  dplyr::mutate(Treatment_Window_Post_Initial.1 = NULL) %>% 
  dplyr::mutate(Chemotherapy_Type_Post_Initial.1 = NULL) %>% 
  dplyr::mutate(Post_Initial.2_cat = paste0(Treatment_Window_Post_Initial.2 , ": ", Chemotherapy_Type_Post_Initial.2)) %>% 
  dplyr::mutate(Treatment_Window_Post_Initial.2 = NULL) %>% 
  dplyr::mutate(Chemotherapy_Type_Post_Initial.2 = NULL) %>% 
  dplyr::mutate(Post_Recurrent1.1_cat = paste0(Treatment_Window_Post_Recurrent1.1 , ": ", Chemotherapy_Type_Post_Recurrent1.1)) %>% 
  dplyr::mutate(Treatment_Window_Post_Recurrent1.1 = NULL) %>% 
  dplyr::mutate(Chemotherapy_Type_Post_Recurrent1.1 = NULL) %>% 
  dplyr::mutate(Post_Recurrent1.2_cat = paste0(Treatment_Window_Post_Recurrent1.2 , ": ", Chemotherapy_Type_Post_Recurrent1.2)) %>% 
  dplyr::mutate(Treatment_Window_Post_Recurrent1.2 = NULL) %>% 
  dplyr::mutate(Chemotherapy_Type_Post_Recurrent1.2 = NULL) %>% 
  dplyr::mutate(Post_Recurrent1.3_cat = paste0(Treatment_Window_Post_Recurrent1.3 , ": ", Chemotherapy_Type_Post_Recurrent1.3)) %>% 
  dplyr::mutate(Treatment_Window_Post_Recurrent1.3 = NULL) %>% 
  dplyr::mutate(Chemotherapy_Type_Post_Recurrent1.3 = NULL) %>% 
  dplyr::mutate(Post_Recurrent1.4_cat = paste0(Treatment_Window_Post_Recurrent1.4 , ": ", Chemotherapy_Type_Post_Recurrent1.4)) %>% 
  dplyr::mutate(Treatment_Window_Post_Recurrent1.4 = NULL) %>% 
  dplyr::mutate(Chemotherapy_Type_Post_Recurrent1.4 = NULL) %>% 
  dplyr::mutate(Post_Recurrent2.1_cat = paste0(Treatment_Window_Post_Recurrent2.1 , ": ", Chemotherapy_Type_Post_Recurrent2.1)) %>% 
  dplyr::mutate(Treatment_Window_Post_Recurrent2.1 = NULL) %>% 
  dplyr::mutate(Chemotherapy_Type_Post_Recurrent2.1 = NULL) %>%
  dplyr::mutate(Post_Recurrent2.2_cat = paste0(Treatment_Window_Post_Recurrent2.2 , ": ", Chemotherapy_Type_Post_Recurrent2.2)) %>% 
  dplyr::mutate(Treatment_Window_Post_Recurrent2.2 = NULL) %>% 
  dplyr::mutate(Chemotherapy_Type_Post_Recurrent2.2 = NULL) %>%
  dplyr::mutate(Post_Recurrent2.3_cat = paste0(Treatment_Window_Post_Recurrent2.3 , ": ", Chemotherapy_Type_Post_Recurrent2.3)) %>% 
  dplyr::mutate(Treatment_Window_Post_Recurrent2.3 = NULL) %>% 
  dplyr::mutate(Chemotherapy_Type_Post_Recurrent2.3 = NULL) %>%
  dplyr::mutate(Post_Recurrent3.1_cat = paste0(Treatment_Window_Post_Recurrent3.1 , ": ", Chemotherapy_Type_Post_Recurrent3.1)) %>% 
  dplyr::mutate(Treatment_Window_Post_Recurrent3.1 = NULL) %>% 
  dplyr::mutate(Chemotherapy_Type_Post_Recurrent3.1 = NULL) %>% 
  tidyr::pivot_longer(!GLASS_ID) %>%
  dplyr::filter(value != 'NA: NA') %>% 
  dplyr::mutate(name = gsub("_cat","",name)) %>% 
  dplyr::mutate(chemotherapy = gsub("^.+: ","", value)) %>% 
  dplyr::mutate(conditional.surgery = gsub("^.+Surgery ([^:]+):.+$","\\1",value)) %>% 
  dplyr::mutate(condition = gsub("^([^ ]+) .+$","\\1",value)) %>% 
  dplyr::mutate(value=NULL) %>% 
  dplyr::left_join(metadata.glass.per.resection %>% dplyr::select(GLASS_ID, Sample_Name), by=c('GLASS_ID'='GLASS_ID')) %>% 
  dplyr::filter(
                (condition == "After" & Sample_Name > conditional.surgery) |
                (condition == "Before" & Sample_Name >= conditional.surgery) 
  ) %>% 
  dplyr::mutate(conditional.surgery = NULL) %>% 
  dplyr::mutate(condition = NULL) %>% 
  dplyr::mutate(name = NULL) %>% 
  dplyr::mutate(GLASS_ID = NULL) %>% 
  dplyr::group_by(Sample_Name) %>% 
  dplyr::summarise(chemotherapy = paste0(chemotherapy,collapse = ",") ) %>% 
  as.data.frame


metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(tmp, by=c('Sample_Name'='Sample_Name'),suffix = c("", "")) 


rm(tmp)


## attach Radio-therapy ----


tmp <- read.csv('data/glass/Clinical data/Cleaned/metadata_2022/Radiotherapy data_GLASS RNAseq.csv') %>% 
  tibble %>% 
  dplyr::select(-contains("Date_")) %>% 
  dplyr::select(-contains("KPS_")) %>% 
  dplyr::select(-contains("_Stopped_")) %>% 
  dplyr::select(-contains("WHO_")) %>% 
  dplyr::select(-contains("Multiple_")) %>% 
  tidyr::pivot_longer(!GLASS_ID) %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::mutate(conditional.surgery = gsub("^.+Surgery ([^:]+)$","\\1",value)) %>% 
  dplyr::mutate(condition = gsub("^([^ ]+) .+$","\\1",value)) %>% 
  dplyr::mutate(value=NULL) %>% 
  dplyr::left_join(metadata.glass.per.resection %>% dplyr::select(GLASS_ID, Sample_Name), by=c('GLASS_ID'='GLASS_ID')) %>% 
  dplyr::filter(
    (condition == "After" & Sample_Name > conditional.surgery) | (condition == "Before" & Sample_Name >= conditional.surgery) 
  ) %>% 
  dplyr::mutate(radiotherapy = "radiotherapy") %>% 
  dplyr::mutate(conditional.surgery = NULL) %>% 
  dplyr::mutate(condition = NULL) %>% 
  dplyr::mutate(name = NULL) %>% 
  dplyr::mutate(GLASS_ID = NULL) %>% 
  dplyr::group_by(Sample_Name) %>% 
  dplyr::summarise(radiotherapy = paste0(radiotherapy, collapse = ",")) %>% 
  as.data.frame


metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::left_join(tmp, by=c('Sample_Name'='Sample_Name'),suffix = c("", "")) 


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
  slice <- metadata.glass.per.patient |> 
    dplyr::filter(GLASS_ID == pid)
  
  
  # define the appropriate primaries and recurrences per patient for RNA
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
  
  
  # define the appropriate primaries and recurrences per patient for protein, for WHO2021
  slice <- slice |> 
    dplyr::filter(!is.na(ProtID)) |> 
    dplyr::mutate(who2021 = case_when(
      WHO_Classification2021 %in% c("Astrocytoma, IDH-mutant, WHO grade 2", "Astrocytoma, IDH-mutant, WHO grade 3") ~ "grade2_3",
      WHO_Classification2021 %in% c("Astrocytoma, IDH-mutant, WHO grade 4") ~ "grade4",
      T ~ "?"
      ))
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




# # zal wel RF score zijn
# a = read.csv('/home/r361003/mnt/neuro-genomic-1-ro/glass/Methylation/Heidelberg/Heidelberg_unzip/203175700013_R01C01_Run_43241/predictBrain_v2.1/203175700013_R01C01_scores.csv')
# sum(a$X203175700013_R01C01)
# 
# # lijkt op soort van probability uit een gefitte density, nooit nul?
# b = read.csv('/home/r361003/mnt/neuro-genomic-1-ro/glass/Methylation/Heidelberg/Heidelberg_unzip/203175700013_R01C01_Run_43241/predictBrain_v2.1/203175700013_R01C01_scores_cal.csv')
# sum(b$X203175700013_R01C01)


# plot(metadata.glass.per.resection$lts.up2, log(metadata.glass.per.resection$A_IDH_HG_cal/metadata.glass.per.resection$A_IDH_cal),xlab="Signature 2 expression",ylab="Heidelberg Classifier log(IDH_LGG_HG score / IDH_LGG score)")
# plot(metadata.glass.per.resection$lts.up2, metadata.glass.per.resection$A_IDH_HG_cal)
# plot(metadata.glass.per.resection$lts.up2, metadata.glass.per.resection$A_IDH_cal)
# 
# 
# plot(metadata.glass.per.resection$lts.up2, metadata.glass.per.resection$A_IDH_HG - metadata.glass.per.resection$A_IDH)
# plot(metadata.glass.per.resection$lts.up2, metadata.glass.per.resection$A_IDH_HG)
# plot(metadata.glass.per.resection$lts.up2, metadata.glass.per.resection$A_IDH)
# 
# 
# plot(metadata.glass.per.resection$lts.down, log(metadata.glass.per.resection$A_IDH_HG_cal/metadata.glass.per.resection$A_IDH_cal),xlab="RNA Signature 4 (down) expression",ylab="Heidelberg Classifier log(IDH_LGG_HG score / IDH_LGG score)")
# plot(metadata.glass.per.resection$lts.up2, metadata.glass.per.resection$A_IDH_HG_cal)
# plot(metadata.glass.per.resection$lts.up2, metadata.glass.per.resection$A_IDH_cal)




