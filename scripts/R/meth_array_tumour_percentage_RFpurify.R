#!/usr/bin/env R 

library(minfi)
library(RFpurify)

# example 
targets <- read.metharray.sheet('data/glass/Methylation' ,pattern ='t3.csv',recursive=F)
#targets <- read.metharray.sheet('.' ,pattern ='t3.csv',recursive=F, verbose=T)

# 
targets_old <- read.csv('Datasheet3.csv')
#targets_new <- Sys.glob('data/glass/Methylation/*/idats/*/*Grn.idat')
#targets_new <- Sys.glob('data/glass/Methylation/*/idats/*/*/*Grn.idat')


create_targets <- function (fn) {
  #fn <- "data/glass/Methylation/MET2019-156-014_plate2/idats/203430580025/203430580025_R01C01_Grn.idat"
  bbasename <- gsub("^(.+)_(Grn|Red).idat$","\\1",fn)
  sslide <- gsub("^.+idats/([^/]+)/.+$","\\1",fn)
  aarray <- gsub("^.+/[^_]+_([^_]+).+$","\\1",fn)
  
  tmp <- targets %>%
    dplyr::filter(Slide == sslide & Array == aarray) %>% 
    dplyr::mutate(Basename = bbasename)

  return(tmp)
}

# 
# 
# RGset <- read.metharray.exp(targets = targets)
# 
# 
# absolute <- predict_purity(RGset,method="ABSOLUTE")
# estimate <- predict_purity(RGset,method="ESTIMATE")
# 
# 
# 
# t <- create_targets('data/glass/Methylation/MET2019-186-014_plate4 redo_kolom 1245/idats/204088040075/204088040075_R06C01_Grn.idat')
# RGset <- read.metharray.exp(targets = t)
# MsetEx <- preprocessRaw(RGset)
# 
# absolute <- predict_purity(MsetEx,method="ABSOLUTE")
# estimate <- predict_purity(MsetEx,method="ESTIMATE")

# 'data/glass/Methylation/MET2019-186-014_plate4 redo_kolom 1245/idats/204088040075/204088040075_R06C01_Grn.idat'


lin <- read.delim('asd.txt',header=F) %>% 
  dplyr::mutate(V1=gsub("^./","data/glass/Methylation/",V1)) %>% 
  dplyr::pull(V1)


out <- data.frame()
for(a in lin){
  t <- create_targets(a)
  print(nrow(t))
  if(nrow(t) >= 1){
    RGset <- read.metharray.exp(targets = t)
    MsetEx <- preprocessRaw(RGset)
    
    absolute <- predict_purity(MsetEx,method="ABSOLUTE")
    estimate <- predict_purity(MsetEx,method="ESTIMATE")
    
    out <- rbind(out,
                 t %>%  dplyr::mutate(absolute=absolute,estimate=estimate)
                 )
  }
}

file <- out %>% 
  dplyr::filter(!duplicated(Basename)) %>% 
  dplyr::arrange(Sample_Name) %>% 
  dplyr::mutate(X=NULL)


write.table(file, file = "output/tables/methylation-array/purities_RFpurity.txt")



plot(file)


