#!/usr/bin/env R 


# load data ----


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if(!exists("expression.glass.vst")) {
  source('scripts/load_rna-counts.R')
}

# load libs ----


# https://cran.r-project.org/web/packages/randomForestSRC/index.html
library(randomForestSRC)
library(survival)


# example RF model ----


data(pbc, package = "randomForestSRC")


pbc.obj <- rfsrc(Surv(days,status) ~ ., pbc, importance = TRUE)
#find.interaction(pbc.obj, method = "vimp", nvar = 8)

plot(pbc.obj)


# R2 svvl RF ----

# Van alle resecties tijd tot dood of laatste event nodig



tmp.metadata <- metadata.glass.per.patient %>%
  tidyr::drop_na(genomescan.sid.R) %>%
  dplyr::select(genomescan.sid.R) %>% 
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid.R'='genomescan.sid')) %>%
  tidyr::drop_na(time.resection.until.last.event) %>% 
  dplyr::arrange(Sample_Type, genomescan.sid.R)



tmp.data <- expression.glass %>%
  dplyr::select(all_of( tmp.metadata$genomescan.sid.R ))



stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)



