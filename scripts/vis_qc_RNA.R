#!/usr/bin/env R

# visualise all RNA QC stats


library("rjson")


# fastp stats ----


parse_fastp_json_files <- function(json_file) {
  out <- fromJSON(file = json_file)
  
  out <- list(
    "a" = sum(out$read1_after_filtering$content_curves$`A`, out$read2_after_filtering$content_curves$`A`),
    "c" = sum(out$read1_after_filtering$content_curves$`C`, out$read2_after_filtering$content_curves$`C`),
    "t" = sum(out$read1_after_filtering$content_curves$`T`, out$read2_after_filtering$content_curves$`T`),
    "g" = sum(out$read1_after_filtering$content_curves$`G`, out$read2_after_filtering$content_curves$`G`))
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
  dplyr::mutate(genomescan.sid = gsub("^.+/[^_]+_([^_]+)_.+$","\\1", fastp.json)) %>% 
  head(n=5) %>% 
  dplyr::mutate(json.stats = lapply(fastp.json, parse_fastp_json_files)) %>% 
  mutate(tmp = map(json.stats, ~ data.frame(t(.)))) %>%  unnest(tmp) %>% dplyr::mutate(json.stats = NULL) %>% as.data.frame






