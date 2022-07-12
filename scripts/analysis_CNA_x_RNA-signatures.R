#!/usr/bin/env R

# see if the CNAs can explain (some part of) the expression of the DGE signatures


# load data ----


if(!exists('cnv2') | !exists('cnv.metadata')) {
  source('scripts/load_genomic_alterations.R')
}


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

# metadata.glass.per.resection$lts.up1
# metadata.glass.per.resection$lts.up2
# metadata.glass.per.resection$lts.up3
# metadata.glass.per.resection$lts.down
# metadata.glass.per.resection$lts.down.a
# metadata.glass.per.resection$lts.down.b



sel <- metadata.glass.per.resection %>% 
  dplyr::filter(excluded == F) %>% 
  dplyr::pull(Sample_Name) %>% 
  intersect(colnames(cnv))

stopifnot(sel %in% metadata.glass.per.resection$Sample_Name)
stopifnot(sel %in% colnames(cnv))


# shape data ----


data <- metadata.glass.per.resection %>% 
  dplyr::filter(Sample_Name %in% sel) %>% 
  dplyr::select(Sample_Name, 
                lts.up1, lts.up2, lts.up3, lts.down, lts.down.a, lts.down.b,
                methylation.sub.diagnosis,mean.DNA.methylation.signature,
                A_IDH_HG_cal,A_IDH_cal,
  ) %>% 
  dplyr::mutate(IDH_HG_IDH_ratio = log(A_IDH_HG_cal/A_IDH_cal)) %>% 
  dplyr::mutate(A_IDH_HG_cal = NULL) %>% 
  dplyr::mutate(A_IDH_cal = NULL) %>%
  dplyr::left_join(
    cnv %>% 
      tibble::rownames_to_column('cnv.segment.id') %>% 
      dplyr::left_join(cnv.metadata %>% 
                         tibble::rownames_to_column('cnv.segment.id')
                       , by=c('cnv.segment.id'='cnv.segment.id')) %>% 
      dplyr::mutate(cnv.segment.id = paste0('cnv.segment.',cnv.segment.id,':',cnv.segments.start,'-',cnv.segments.end)) %>% 
      dplyr::mutate(cnv.segments.chrom = NULL) %>% 
      dplyr::mutate( cnv.segments.start = NULL) %>% 
      dplyr::mutate( cnv.segments.end= NULL) %>% 
      tibble::column_to_rownames('cnv.segment.id') %>% 
      t() %>% 
      as.data.frame %>% 
      tibble::rownames_to_column('Sample_Name')
    , by=c('Sample_Name'='Sample_Name')
  ) %>% 
  tibble::column_to_rownames('Sample_Name')



data.discrete.gain <- metadata.glass.per.resection %>% 
  dplyr::filter(Sample_Name %in% sel) %>% 
  dplyr::select(Sample_Name, 
                lts.up1, lts.up2, lts.up3, lts.down, lts.down.a, lts.down.b,
                methylation.sub.diagnosis,mean.DNA.methylation.signature,
                A_IDH_HG_cal,A_IDH_cal,
  ) %>% 
  dplyr::mutate(IDH_HG_IDH_ratio = log(A_IDH_HG_cal/A_IDH_cal)) %>% 
  dplyr::mutate(A_IDH_HG_cal = NULL) %>% 
  dplyr::mutate(A_IDH_cal = NULL) %>%
  dplyr::left_join(
    cnv %>% 
      tibble::rownames_to_column('cnv.segment.id') %>% 
      dplyr::left_join(cnv.metadata %>% 
                         tibble::rownames_to_column('cnv.segment.id')
                       , by=c('cnv.segment.id'='cnv.segment.id')) %>% 
      dplyr::mutate(cnv.segment.id = paste0('cnv.segment.',cnv.segment.id,':',cnv.segments.start,'-',cnv.segments.end)) %>% 
      dplyr::mutate(cnv.segments.chrom = NULL) %>% 
      dplyr::mutate( cnv.segments.start = NULL) %>% 
      dplyr::mutate( cnv.segments.end= NULL) %>% 
      tibble::column_to_rownames('cnv.segment.id') %>% 
      t() %>% 
      as.data.frame %>% 
      dplyr::mutate_all(function(x){return (ifelse(x > 0,"gain","no-gain"))}) %>% 
      tibble::rownames_to_column('Sample_Name')
    , by=c('Sample_Name'='Sample_Name')
  ) %>% 
  tibble::column_to_rownames('Sample_Name')




data.discrete.loss <- metadata.glass.per.resection %>% 
  dplyr::filter(Sample_Name %in% sel) %>% 
  dplyr::select(Sample_Name, 
                lts.up1, lts.up2, lts.up3, lts.down, lts.down.a, lts.down.b,
                methylation.sub.diagnosis,mean.DNA.methylation.signature,
                A_IDH_HG_cal,A_IDH_cal,
  ) %>% 
  dplyr::mutate(IDH_HG_IDH_ratio = log(A_IDH_HG_cal/A_IDH_cal)) %>% 
  dplyr::mutate(A_IDH_HG_cal = NULL) %>% 
  dplyr::mutate(A_IDH_cal = NULL) %>%
  dplyr::left_join(
    cnv %>% 
      tibble::rownames_to_column('cnv.segment.id') %>% 
      dplyr::left_join(cnv.metadata %>% 
                         tibble::rownames_to_column('cnv.segment.id')
                       , by=c('cnv.segment.id'='cnv.segment.id')) %>% 
      dplyr::mutate(cnv.segment.id = paste0('cnv.segment.',cnv.segment.id,':',cnv.segments.start,'-',cnv.segments.end)) %>% 
      dplyr::mutate(cnv.segments.chrom = NULL) %>% 
      dplyr::mutate( cnv.segments.start = NULL) %>% 
      dplyr::mutate( cnv.segments.end= NULL) %>% 
      tibble::column_to_rownames('cnv.segment.id') %>% 
      t() %>% 
      as.data.frame %>% 
      dplyr::mutate_all(function(x){return (ifelse(x < 0,"loss","no-loss"))}) %>% 
      tibble::rownames_to_column('Sample_Name')
    , by=c('Sample_Name'='Sample_Name')
  ) %>% 
  tibble::column_to_rownames('Sample_Name')

# shape data2 ----


sel <- metadata.glass.per.resection %>% 
  dplyr::filter(excluded == F) %>% 
  dplyr::pull(Sample_Name) %>% 
  intersect(colnames(cnv2))


data2 <- metadata.glass.per.resection %>% 
  dplyr::filter(Sample_Name %in% sel) %>% 
  dplyr::select(Sample_Name, 
                lts.up1, lts.up2, lts.up3, lts.down, lts.down.a, lts.down.b,
                methylation.sub.diagnosis,mean.DNA.methylation.signature,
                A_IDH_HG_cal,A_IDH_cal,
  ) %>% 
  dplyr::mutate(IDH_HG_IDH_ratio = log(A_IDH_HG_cal/A_IDH_cal)) %>% 
  dplyr::mutate(A_IDH_HG_cal = NULL) %>% 
  dplyr::mutate(A_IDH_cal = NULL) %>%
  dplyr::left_join(
    cnv2 %>% 
      t() %>% 
      as.data.frame %>% 
      tibble::rownames_to_column('Sample_Name')
    , by=c('Sample_Name'='Sample_Name')
  ) %>% 
  tibble::column_to_rownames('Sample_Name')


data2.discrete.gain <- metadata.glass.per.resection %>% 
  dplyr::filter(Sample_Name %in% sel) %>% 
  dplyr::select(Sample_Name, 
                lts.up1, lts.up2, lts.up3, lts.down, lts.down.a, lts.down.b,
                methylation.sub.diagnosis,mean.DNA.methylation.signature,
                A_IDH_HG_cal,A_IDH_cal,
  ) %>% 
  dplyr::mutate(IDH_HG_IDH_ratio = log(A_IDH_HG_cal/A_IDH_cal)) %>% 
  dplyr::mutate(A_IDH_HG_cal = NULL) %>% 
  dplyr::mutate(A_IDH_cal = NULL) %>%
  dplyr::left_join(
    cnv2 %>% 
      t() %>% 
      as.data.frame %>% 
      dplyr::mutate_all(function(x){return (ifelse(x > 0,"gain","no-gain"))}) %>% 
      tibble::rownames_to_column('Sample_Name')
    , by=c('Sample_Name'='Sample_Name')
  ) %>% 
  tibble::column_to_rownames('Sample_Name')




data2.discrete.loss <- metadata.glass.per.resection %>% 
  dplyr::filter(Sample_Name %in% sel) %>% 
  dplyr::select(Sample_Name, 
                lts.up1, lts.up2, lts.up3, lts.down, lts.down.a, lts.down.b,
                methylation.sub.diagnosis,mean.DNA.methylation.signature,
                A_IDH_HG_cal,A_IDH_cal,
  ) %>% 
  dplyr::mutate(IDH_HG_IDH_ratio = log(A_IDH_HG_cal/A_IDH_cal)) %>% 
  dplyr::mutate(A_IDH_HG_cal = NULL) %>% 
  dplyr::mutate(A_IDH_cal = NULL) %>%
  dplyr::left_join(
    cnv2 %>% 
      t() %>% 
      as.data.frame %>%
      dplyr::mutate_all(function(x){return (ifelse(x < 0,"loss","no-loss"))}) %>% 
      tibble::rownames_to_column('Sample_Name')
    , by=c('Sample_Name'='Sample_Name')
  ) %>% 
  tibble::column_to_rownames('Sample_Name')






# base stats var ----



# cnv ----

stats <- data.frame(cnv.segment.id = colnames(data)) %>% 
  dplyr::filter(grepl("^cnv.seg", cnv.segment.id))


## lts.up1 [wilcox/gain] ----

stats$lts.up1.gain.wilcox <- NA
for(x in colnames(data.discrete.gain)) {
  if(x %in% stats$cnv.segment.id) {
    #print(x)
    
    tmp <- data.discrete.gain %>%  dplyr::select(lts.up1, c(x))
    stopifnot(names(tmp)[2] == x)
    names(tmp)[1] <- 'expression'
    names(tmp)[2] <- 'segment'
    
    g1 <- tmp %>% dplyr::filter(segment == "gain") %>%  dplyr::pull(expression)
    g2 <- tmp %>% dplyr::filter(segment == "no-gain") %>%  dplyr::pull(expression)
    
    if(length(g1) > 0 & length(g2) > 0)
      w <- wilcox.test(g1, g2)
    stats$lts.up1.gain.wilcox[stats$cnv.segment.id == x] <- w$p.value
  }
}


## lts.up1 [wilcox/loss] ----

stats$lts.up1.loss.wilcox <- NA
for(x in colnames(data.discrete.loss)) {
  if(x %in% stats$cnv.segment.id) {
    #print(x)
    
    tmp <- data.discrete.loss %>%  dplyr::select(lts.up1, c(x))
    stopifnot(names(tmp)[2] == x)
    names(tmp)[1] <- 'expression'
    names(tmp)[2] <- 'segment'
    
    g1 <- tmp %>% dplyr::filter(segment == "loss") %>%  dplyr::pull(expression)
    g2 <- tmp %>% dplyr::filter(segment == "no-loss") %>%  dplyr::pull(expression)
    
    if(length(g1) > 0 & length(g2) > 0)
      w <- wilcox.test(g1, g2)
    stats$lts.up1.loss.wilcox[stats$cnv.segment.id == x] <- w$p.value
  }
}


## lts.up1 [ttest/gain] ----

stats$lts.up1.gain.ttest <- NA
for(x in colnames(data.discrete.gain)) {
  if(x %in% stats$cnv.segment.id) {
    #print(x)
    
    tmp <- data.discrete.gain %>%  dplyr::select(lts.up1, c(x))
    stopifnot(names(tmp)[2] == x)
    names(tmp)[1] <- 'expression'
    names(tmp)[2] <- 'segment'
    
    g1 <- tmp %>% dplyr::filter(segment == "gain") %>%  dplyr::pull(expression)
    g2 <- tmp %>% dplyr::filter(segment == "no-gain") %>%  dplyr::pull(expression)
    
    if(length(g1) > 1 & length(g2) > 1)
      w <- t.test(g1, g2)
    stats$lts.up1.gain.ttest[stats$cnv.segment.id == x] <- w$p.value
  }
}


## lts.up1 [ttest/loss] ----

stats$lts.up1.loss.ttest <- NA
for(x in colnames(data.discrete.loss)) {
  if(x %in% stats$cnv.segment.id) {
    #print(x)
    
    tmp <- data.discrete.loss %>%  dplyr::select(lts.up1, c(x))
    stopifnot(names(tmp)[2] == x)
    names(tmp)[1] <- 'expression'
    names(tmp)[2] <- 'segment'
    
    g1 <- tmp %>% dplyr::filter(segment == "loss") %>%  dplyr::pull(expression)
    g2 <- tmp %>% dplyr::filter(segment == "no-loss") %>%  dplyr::pull(expression)
    
    if(length(g1) > 1 & length(g2) > 1)
    w <- t.test(g1, g2)
    stats$lts.up1.loss.ttest[stats$cnv.segment.id == x] <- w$p.value
  }
}


# cnv2 ----


stats2 <- data.frame(cnv.segment.id = colnames(data2)) %>% 
  dplyr::filter(grepl("^chr", cnv.segment.id))



## lts.up1 [ttest/gain] ----

i <- 0
stats2$lts.up1.gain.ttest <- NA
for(x in colnames(data2.discrete.gain)) {
  if(x %in% stats2$cnv.segment.id) {
    #if(i %% 33 == 1) {
      print(i)
    
      tmp <- data2.discrete.gain %>%  dplyr::select(lts.up1, c(x))
      stopifnot(names(tmp)[2] == x)
      names(tmp)[1] <- 'expression'
      names(tmp)[2] <- 'segment'
      
      g1 <- tmp %>% dplyr::filter(segment == "gain") %>%  dplyr::pull(expression)
      g2 <- tmp %>% dplyr::filter(segment == "no-gain") %>%  dplyr::pull(expression)
      
      if(length(g1) > 1 & length(g2) > 1) {
        w <- t.test(g1, g2)
        stats2$lts.up1.gain.ttest[stats2$cnv.segment.id == x] <- w$p.value
      }
    #}

    i <- i + 1
  }
}


i <- 0
stats2$lts.up1.loss.ttest <- NA
for(x in colnames(data2.discrete.loss)) {
  if(x %in% stats2$cnv.segment.id) {
    #if(i %% 33 == 1) {
    #if(grepl("chr3", x)) {
      print(i)
      
      tmp <- data2.discrete.loss %>%  dplyr::select(lts.up1, c(x))
      stopifnot(names(tmp)[2] == x)
      names(tmp)[1] <- 'expression'
      names(tmp)[2] <- 'segment'
      
      g1 <- tmp %>% dplyr::filter(segment == "loss") %>%  dplyr::pull(expression)
      g2 <- tmp %>% dplyr::filter(segment == "no-loss") %>%  dplyr::pull(expression)
      
      if(length(g1) > 1 & length(g2) > 1) {
        w <- t.test(g1, g2)
        stats2$lts.up1.loss.ttest[stats2$cnv.segment.id == x] <- w$p.value
      }
    #}
    
    i <- i + 1
  }
}


saveRDS(stats2, file="cache/stats2.Rds")


## IDH_IDH_HG p ratio ----

i <- 0
stats2$IDH_HG_IDH_ratio.gain.ttest <- NA
for(x in colnames(data2.discrete.gain)) {
  if(x %in% stats2$cnv.segment.id) {
      if(i %% 33 == 1) {
        print(paste0(i, " -> ", round(i / length(stats2$cnv.segment.id) * 100,2), "%"))
      }
      
      tmp <- data2.discrete.gain %>%  dplyr::select(IDH_HG_IDH_ratio, c(x))
      stopifnot(names(tmp)[2] == x)
      names(tmp)[1] <- 'expression'
      names(tmp)[2] <- 'segment'
      
      g1 <- tmp %>% dplyr::filter(segment == "gain") %>%  dplyr::pull(expression)
      g2 <- tmp %>% dplyr::filter(segment == "no-gain") %>%  dplyr::pull(expression)
      
      if(length(g1) > 1 & length(g2) > 1) {
        w <- t.test(g1, g2)
        stats2$IDH_HG_IDH_ratio.gain.ttest[stats2$cnv.segment.id == x] <- w$p.value
      }

    i <- i + 1
  }
}


i <- 0
stats2$IDH_HG_IDH_ratio.loss.ttest <- NA
for(x in colnames(data2.discrete.loss)) {
  if(x %in% stats2$cnv.segment.id) {
    #if(i %% 133 == 1) {
    if(i %% 33 == 1) {
      print(paste0(i, " -> ", round(i / length(stats2$cnv.segment.id) * 100,2), "%"))
    }
      
    tmp <- data2.discrete.loss %>%  dplyr::select(IDH_HG_IDH_ratio, c(x))
    stopifnot(names(tmp)[2] == x)
    names(tmp)[1] <- 'expression'
    names(tmp)[2] <- 'segment'
    
    g1 <- tmp %>% dplyr::filter(segment == "loss") %>%  dplyr::pull(expression)
    g2 <- tmp %>% dplyr::filter(segment == "no-loss") %>%  dplyr::pull(expression)
    
    if(length(g1) > 1 & length(g2) > 1) {
      w <- t.test(g1, g2)
      stats2$IDH_HG_IDH_ratio.loss.ttest[stats2$cnv.segment.id == x] <- w$p.value
    }
    
    i <- i + 1
  }
}


saveRDS(stats2, file="cache/stats2.Rds")

## mean meth signature ----


i <- 0
stats2$mean.DNA.methylation.signature.gain.ttest <- NA
for(x in colnames(data2.discrete.gain)) {
  if(x %in% stats2$cnv.segment.id) {
    if(i %% 33 == 1) {
      print(paste0(i, " -> ", round(i / length(stats2$cnv.segment.id) * 100,2), "%"))
    }
    
    tmp <- data2.discrete.gain %>%  dplyr::select(mean.DNA.methylation.signature, c(x))
    stopifnot(names(tmp)[2] == x)
    names(tmp)[1] <- 'expression'
    names(tmp)[2] <- 'segment'
    
    g1 <- tmp %>% dplyr::filter(segment == "gain" & !is.na(expression)) %>%  dplyr::pull(expression)
    g2 <- tmp %>% dplyr::filter(segment == "no-gain" & !is.na(expression)) %>%  dplyr::pull(expression)
    
    if(length(g1) > 1 & length(g2) > 1) {
      w <- t.test(g1, g2)
      stats2$mean.DNA.methylation.signature.gain.ttest[stats2$cnv.segment.id == x] <- w$p.value
    }
    
    i <- i + 1
  }
}


i <- 0
stats2$mean.DNA.methylation.signature.loss.ttest <- NA
for(x in colnames(data2.discrete.loss)) {
  if(x %in% stats2$cnv.segment.id) {
    #if(i %% 133 == 1) {
    if(i %% 33 == 1) {
      print(paste0(i, " -> ", round(i / length(stats2$cnv.segment.id) * 100,2), "%"))
    }
    
    tmp <- data2.discrete.loss %>%  dplyr::select(mean.DNA.methylation.signature, c(x))
    stopifnot(names(tmp)[2] == x)
    names(tmp)[1] <- 'expression'
    names(tmp)[2] <- 'segment'
    
    g1 <- tmp %>% dplyr::filter(segment == "loss" & !is.na(expression)) %>%  dplyr::pull(expression)
    g2 <- tmp %>% dplyr::filter(segment == "no-loss" & !is.na(expression)) %>%  dplyr::pull(expression)
    
    #print(g1)
    
    if(length(g1) > 1 & length(g2) > 1) {
      w <- t.test(g1, g2)
      stats2$mean.DNA.methylation.signature.loss.ttest[stats2$cnv.segment.id == x] <- w$p.value
    }
    
    i <- i + 1
  }
}


saveRDS(stats2, file="cache/stats2.Rds")



# plot stats2 ----

source('scripts/R/chrom_sizes.R')

plt <- stats2 %>% 
  dplyr::mutate(chr = gsub(":.+$","",cnv.segment.id)) %>% 
  dplyr::mutate(start = as.numeric(gsub("^.+:([0-9]+)-.+$","\\1",cnv.segment.id))) %>% 
  dplyr::left_join(
    data.frame(chrs_hg19_s) %>% tibble::rownames_to_column('chr'),
    by=c('chr'='chr') 
  ) %>% 
  dplyr::mutate(x = chrs_hg19_s + start) %>% 
  dplyr::mutate(y.gain.lts.up1 = -log10(lts.up1.gain.ttest)) %>% 
  dplyr::mutate(y.loss.lts.up1 = -log10(lts.up1.loss.ttest)) %>% 
  dplyr::mutate(y.gain.IDH_HG_IDH_ratio = -log10(IDH_HG_IDH_ratio.gain.ttest)) %>% 
  dplyr::mutate(y.loss.IDH_HG_IDH_ratio = -log10(IDH_HG_IDH_ratio.loss.ttest)) %>% 
  dplyr::mutate(y.gain.mean.DNA.methylation.signature = -log10(mean.DNA.methylation.signature.gain.ttest)) %>% 
  dplyr::mutate(y.loss.mean.DNA.methylation.signature = -log10(mean.DNA.methylation.signature.loss.ttest))

# https://www.nature.com/articles/s41379-021-00778-x.pdf
# https://www.nature.com/articles/s41379-021-00778-x.pdf

plt.exp <- rbind(plt %>% dplyr::mutate(y =  y.gain.lts.up1, y.loss = NULL, y.gain=NULL, group="gain", facet="lts.up1") ,
                 plt %>% dplyr::mutate(y = -y.loss.lts.up1, y.loss=NULL, y.gain=NULL , group = "loss", facet="lts.up1") ,
                 
                 plt %>% dplyr::mutate(y =  y.gain.IDH_HG_IDH_ratio, y.loss = NULL, y.gain=NULL, group="gain",facet="A_IDH / A_IDH_HG ratio") ,
                 plt %>% dplyr::mutate(y = -y.loss.IDH_HG_IDH_ratio, y.loss=NULL, y.gain=NULL , group = "loss",facet="A_IDH / A_IDH_HG ratio"),
                 
                 plt %>% dplyr::mutate(y =  y.gain.mean.DNA.methylation.signature, y.loss = NULL, y.gain=NULL, group="gain",facet="mean DNA methylation signature") ,
                 plt %>% dplyr::mutate(y = -y.loss.mean.DNA.methylation.signature, y.loss=NULL, y.gain=NULL , group = "loss",facet="mean DNA methylation signature")
                 ) %>% 
  dplyr::mutate(end = as.numeric(gsub("^.+\\-","",cnv.segment.id))) %>% 
  dplyr::filter(!is.na(y)) %>% 
  dplyr::mutate(label = NA) %>% 
  
  
  #dplyr::mutate(label = ifelse(chr == "chr1" &  start < 144427183 & end > 144427222 & group=="gain", "NBPF15", label)) %>% 
  #dplyr::mutate(label = ifelse(chr == "chr1" &  start < 146943172 & end > 146943518 & group=="gain", "NBPF12", label)) %>% 
  dplyr::mutate(label = ifelse(chr == "chr1" &  start < 146586343 & end > 146586805 & group=="gain", "NBPF13P", label)) %>% 
  dplyr::mutate(label = ifelse(chr == "chr1" &  start < 76320465 & end > 76320504 & group=="loss", "MSH4", label)) %>%
  
  dplyr::mutate(label = ifelse(chr == "chr3" &  start < 60494128 & end > 60499364 & group=="loss", "FHIT", label)) %>% 
  dplyr::mutate(label = ifelse(chr == "chr3" &  start < 46355679 & end > 46355718 & group=="loss", "CCR2, CCR3 & CCR5", label)) %>% 
  dplyr::mutate(label = ifelse(chr == "chr4" &  start < 54261320 & end > 54261699 & group=="gain", "PDGFRA", label)) %>% 
  dplyr::mutate(label = ifelse(chr == "chr4" &  start < 183446165 & end > 183446204 & group=="loss", "CDKN2AIP\nchr4:q35.1", label)) %>%
  
  dplyr::mutate(label = ifelse(chr == "chr5" &  start < 178249526 & end > 178249565 & group=="loss", "chr5:q35", label)) %>%
  
  dplyr::mutate(label = ifelse(chr == "chr6" &  start < 21941613 & end > 21949237 & group=="gain", "CASC15", label)) %>%
  dplyr::mutate(label = ifelse(chr == "chr6" &  start < 33319881 & end > 33319920 & group=="gain", "DAXX", label)) %>%
  dplyr::mutate(label = ifelse(chr == "chr6" &  start < 168490294 & end > 168493152 & group=="loss", "SMOC2", label)) %>%
  
  dplyr::mutate(label = ifelse(chr == "chr7" &  start < 14510947 & end > 14512418 & group=="gain", "DGKB", label)) %>%
  dplyr::mutate(label = ifelse(chr == "chr7" &  start < 116796638 & end > 116796677 & group=="gain", "MET", label)) %>%
  #dplyr::mutate(label = ifelse(chr == "chr7" &  start < 92708634 & end > 92709042 & group=="gain", "CDK6", label)) %>%
  
  
  dplyr::mutate(label = ifelse(chr == "chr8" &  start < 128752748 & end > 128752787 & group=="gain", "MYC", label)) %>%
  #dplyr::mutate(label = ifelse(chr == "chr8" &  start < 129955736 & end > 129955775 & group=="gain", "LINC00976", label)) %>%
  
  dplyr::mutate(label = ifelse(chr == "chr9" &  start < 21965752 & end > 21997324 & group=="loss", "CDKN2A/B", label)) %>%
  dplyr::mutate(label = ifelse(chr == "chr9" &  start < 2137095 & end > 2138078 & group=="loss", "SMARCA2", label)) %>%
  
  dplyr::mutate(label = ifelse(chr == "chr10" &  start < 87863980 & end > 87865285 & group=="loss", "PTEN", label))%>%
  dplyr::mutate(label = ifelse(chr == "chr10" &  start < 113033733 & end > 113076506 & group=="loss", "TCF7L2", label))%>%
  dplyr::mutate(label = ifelse(chr == "chr10" &  start < 27700002 & end > 27700001 & group=="gain", "RAB18/MKX", label))%>%
  #dplyr::mutate(label = ifelse(chr == "chr10" &  start < 131271659 & end > 131277709 & group=="loss", "MGMT", label))%>%
  
  dplyr::mutate(label = ifelse(chr == "chr11" &  start < 27709783 & end > 27709822 & group=="loss", "chr11.p", label))%>%
  dplyr::mutate(label = ifelse(chr == "chr11" &  start < 69647835 & end > 69648185 & group=="loss", "CCND1", label))%>%
  
  #dplyr::mutate(label = ifelse(chr == "chr12" &  start < 4286468 & end > 4286714 & group=="gain", "CCND2", label))%>%
  dplyr::mutate(label = ifelse(chr == "chr12" &  start < 57749606 & end > 57749744 & group=="gain", "CDK4", label))%>%
  #dplyr::mutate(label = ifelse(chr == "chr12" &  start < 68827001 & end > 68827079 & group=="gain", "MDM2", label))%>%
  dplyr::mutate(label = ifelse(chr == "chr12" &  start < 25225500 & end > 25225643 & group=="gain", "KRAS", label))%>%
  
  dplyr::mutate(label = ifelse(chr == "chr13" &  start <= 48474808 & end >= 48476242 & group=="loss", "RB1", label)) %>%
  dplyr::mutate(label = ifelse(chr == "chr13" &  start <= 49000001 & end >= 49100000 & group=="loss", "FNDC3A", label)) %>%
  
  dplyr::mutate(label = ifelse(chr == "chr16" &  start < 89689363 & end > 89689402 & group=="gain", "CDK10", label)) %>% 

  dplyr::mutate(label = ifelse(chr == "chr17" &  start < 31366605 & end > 31376825  & group=="loss", "NF1", label)) %>% 
  #dplyr::mutate(label = ifelse(chr == "chr17" &  start < 39702465 & end > 39702516  & group=="gain", "ERBB2", label)) %>%
  
  dplyr::mutate(label = ifelse(chr == "chr19" &  start < 42282635 & end > 42282674  & group=="loss", "CIC", label)) %>% 
  dplyr::mutate(label = ifelse(chr == "chr19" &  start < 29818075 & end > 29818127  & group=="gain", "CCNE1", label)) %>%
  #dplyr::mutate(label = ifelse(chr == "chr19" &  start < 44499806 & end > 44499845  & group=="loss" , "ZNF locus end", label))
  dplyr::mutate(label = ifelse(chr == "chr19" &  start < 39502161 & end > 39502733  & group=="gain", "DLL3", label)) %>% 
  dplyr::mutate(y = ifelse(y > 8, 8 , y))

#unique(plt.exp$label)



#plt.exp <- plt.exp %>%  dplyr::filter(chr %in% c("chr8"))

ggplot(plt.exp, aes(x=x, y=y, col=chr, group=group, label=label)) +
  labs(y = "Association CNA locus with glioma progression",x=NULL) +
  geom_line(col="gray",alpha=0.75) + 
  geom_point(cex=0.84) + 
  theme_bw() + 
  ylim(-8.5,8.5) +
  facet_grid(rows = vars(facet), scales = "free") +
  ggrepel::geom_text_repel(data = subset(plt.exp, y < 0 & !is.na(label)),cex=2.8,col="black",  direction = "y", nudge_y = -5, segment.size=0.25 )+
  ggrepel::geom_text_repel(data = subset(plt.exp, y > 0 & !is.na(label)),cex=2.8,col="black",  direction = "y", nudge_y = 5, segment.size=0.25 )
  # geom_vline(xintercept =  chrs_hg19_s['chr8'] + 131228209)
  #geom_vline(xintercept =  chrs_hg19_s['chr10'] + 32000000)



plt.exp %>%  dplyr::filter(chr == "chr1" & facet == "lts.up1") %>% arrange(-y,start) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr1" & facet == "mean DNA methylation signature") %>% arrange(-y,start) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr1" & facet == "A_IDH / A_IDH_HG ratio") %>% arrange(-y,start) %>% head(n=28)

plt.exp %>%  dplyr::filter(chr == "chr1" & facet == "A_IDH / A_IDH_HG ratio") %>% arrange(y,start) %>% head(n=28)

plt.exp %>%  dplyr::filter(chr == "chr15" & facet == "mean DNA methylation signature") %>% arrange(y,start) %>% head(n=28)


plt.exp %>%  dplyr::filter(chr == "chr2" & facet == "A_IDH / A_IDH_HG ratio") %>% arrange(y,start) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr10" & facet == "A_IDH / A_IDH_HG ratio") %>% arrange(y,start) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr8" & facet == "A_IDH / A_IDH_HG ratio") %>% arrange(-y,start) %>% head(n=28)

plt.exp %>%  dplyr::filter(chr == "chr3") %>% arrange(y,start) %>% head(n=28) # FHIT
plt.exp %>%  dplyr::filter(chr == "chr3" & facet != "lts.up1") %>% arrange(y,start) %>% head(n=28) # FHIT
plt.exp %>%  dplyr::filter(chr == "chr4" & start> 68000000) %>% arrange(-y,start) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr6") %>% arrange(-y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr6") %>% arrange(y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr7") %>% arrange(-y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr9" & start < 17000000) %>% arrange(y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr10") %>% arrange(y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr12" & start < 70000000) %>% arrange(-y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr12" & start > 70000000) %>% arrange(-y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr13") %>% arrange(y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr13" & start < 48226275) %>% arrange(y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr15") %>% arrange(-y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr15") %>% arrange(y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr16") %>% arrange(-y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr16") %>% arrange(y) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr17") %>% arrange(-y,start) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr18") %>% arrange(y,start) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr18") %>% arrange(-y,start) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr19") %>% arrange(y,start) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr19") %>% arrange(-y,start) %>% head(n=28)
plt.exp %>%  dplyr::filter(chr == "chr22") %>% arrange(y,start) %>% head(n=28)



# plt stats1 ----

#stats %>% dplyr::mutate(`min.sig.wilcox` = min(lts.up1.gain.wilcox, lts.up1.loss.wilcox) ) 
plt <- stats %>%
  dplyr::mutate(log.lts.up1.loss.wilcox = -log10(lts.up1.loss.wilcox) )  %>%
  dplyr::mutate(log.lts.up1.gain.wilcox = -log10(lts.up1.gain.wilcox) )  %>% 
  dplyr::mutate(log.lts.up1.loss.ttest = -log10(lts.up1.loss.ttest) )  %>%
  dplyr::mutate(log.lts.up1.gain.ttest = -log10(lts.up1.gain.ttest) )  %>% 
  dplyr::mutate(label = gsub("cnv.segment.chr","",cnv.segment.id))


ggplot(plt , aes(x = log.lts.up1.loss.wilcox, y= log.lts.up1.gain.wilcox, label=label)) + 
  geom_point() +
  geom_text_repel(data = subset(plt, log.lts.up1.gain.wilcox > 4 | log.lts.up1.loss.wilcox > 5),box.padding = 0.5, max.overlaps = Inf)


ggplot(plt , aes(x = log.lts.up1.loss.ttest, y= log.lts.up1.gain.ttest, label=label)) + 
  geom_point()  +
  geom_text_repel(data = subset(plt, log.lts.up1.gain.ttest > 5 | log.lts.up1.loss.ttest > 5),box.padding = 0.5, max.overlaps = Inf)



# regress w/ RNA signatures (linear) ----

## lts.up1 ----

data.lts.up1  <- data %>% 
  dplyr::select(-c("lts.up2", "lts.up3", "lts.down", "lts.down.a", "lts.down.b", "methylation.sub.diagnosis", "IDH_HG_IDH_ratio"))


lm.up1 <- lm(lts.up1 ~ ., data = data.lts.up1)

summary(lm.up1)


summary(lm(lts.up1 ~ `cnv.segment.chr3_26:36600000-37400000`, data = data.lts.up1))
summary(lm(lts.up1 ~ `cnv.segment.chr3_35:47325000-48775000`, data = data.lts.up1))


# regress w/ methylation grading (logistic) ----
# regress w/ methylation grading logratio (linear) ----




