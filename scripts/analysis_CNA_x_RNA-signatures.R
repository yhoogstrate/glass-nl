#!/usr/bin/env R

# see if the CNAs can explain (some part of) the expression of the DGE signatures


# load data ----


if(!exists('cnv') | !exists('cnv.metadata')) {
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
                methylation.sub.diagnosis,
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
                methylation.sub.diagnosis,
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
                methylation.sub.diagnosis,
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





# base stats var ----

stats <- data.frame(cnv.segment.id = colnames(data)) %>% 
  dplyr::filter(grepl("^cnv.seg", cnv.segment.id))



# lts.up1 [wilcox/gain] ----

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


# lts.up1 [wilcox/loss] ----

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


# lts.up1 [ttest/gain] ----

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


# lts.up1 [ttest/loss] ----

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




# plt ----

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




