#!/usr/bin/env R 

#BiocManager::install("minfi")
#BiocManager::install("IlluminaHumanMethylationEPICmanifest")
#BiocManager::install("IlluminaHumanMethylation450kmanifest")
#devtools::install_github('mwsill/RFpurify')


library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation450kmanifest)
library(tidyverse)
library(minfi)
library(RFpurify)
data('IlluminaHumanMethylationEPICmanifest')
data('RFpurify_ABSOLUTE')
data('RFpurify_ESTIMATE')


# obtain idat files ----

# glass/Methylation/20190516_Methylation\ 1\ \(HereisyOUr79Data\)/MET2019-156-014_plate1/MET2019-156-014_plate1/idats/*.idat
# glass/Methylation/20200602_Methylation2_(pTg2@8zRT)/MET2019-186-014/idats/*.idat

# glass/Methylation/MET2019-156-014_plate2/idats/*.idat
# glass/Methylation/MET2019-186-014plate3/MET2019-186-014plate3/idats/*/*.idat
# glass/Methylation/MET2019-186-014_plate4\ redo_kolom\ 1245/idats/*/*.idat
# glass/Methylation/MET2019-186-014_plate4 redo_kolom 3/MET2019-186-014redoplate3kolom3/idats/204086170028/*.idat
# glass/Methylation/MET2019-186-014_plate7/*.idat


c1 <- Sys.glob("data/glass/Methylation/Methylation Array Data/20190516_Methylation 1 (HereisyOUr79Data)/*/*/idats/*/*.idat")
c2 <- Sys.glob("data/glass/Methylation/Methylation Array Data/20200602_Methylation2_(pTg2@8zRT)/MET2019-186-014/idats/*/*.idat")
c3 <- Sys.glob("data/glass/Methylation/Methylation Array Data/MET2019-156-014_plate2/idats/*/*.idat")
c4 <- Sys.glob("data/glass/Methylation/Methylation Array Data/MET2019-186-014plate3/*/idats/*/*.idat")
c5 <- Sys.glob("data/glass/Methylation/Methylation Array Data/MET2019-186-014_plate4\ redo_kolom\ 1245/idats/*/*.idat")
c6 <- Sys.glob("data/glass/Methylation/Methylation Array Data/MET2019-186-014_plate4 redo_kolom 3/MET2019-186-014redoplate3kolom3/idats/204086170028/*.idat")
c7 <- Sys.glob("data/glass/Methylation/Methylation Array Data/MET2019-186-014_plate7/*.idat")

targets <- unique(c(c1,c2,c3,c4,c5,c6,c7))
targets <- targets[grepl("Grn",targets)]

rm(c1,c2,c3,c4,c5,c6,c7)


stopifnot(length(targets) == 329)


targets.df <- data.frame(fn = targets) %>% 
  dplyr::mutate(Basename = gsub("^(.+)_(Grn|Red).idat$","\\1",fn)) %>% 
  dplyr::mutate(Slide = gsub("^.+idats/([^/]+)/.+$","\\1",fn)) %>% 
  dplyr::mutate(Array = gsub("^.+/[^_]+_([^_]+).+$","\\1",fn))




# apply x calc purity from Array ----

out <- data.frame()
for(i in 1:nrow(targets.df)) {
#i <- 323
  print(i)
  
  slice <- targets.df[i,]
  
  RGset <- read.metharray.exp(targets = slice %>% head(n=1))
  MsetEx <- preprocessRaw(RGset)
  
  abs <- predict_purity(MsetEx,method="ABSOLUTE")
  #print(abs)
  
  est <- predict_purity(MsetEx,method="ESTIMATE")
  #print(est)
  
  
  out <- rbind(out, slice %>% dplyr::mutate(absolute=abs,estimate=est))
}




write.table(file, file = "output/tables/methylation-array/purities_RFpurity.txt")



