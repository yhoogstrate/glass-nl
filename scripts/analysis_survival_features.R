#!/usr/bin/env R 


# load libs ----


library(randomForestSRC)# https://cran.r-project.org/web/packages/randomForestSRC/index.html
library(survival)
library(RegParallel)



# load data ----


source('scripts/R/youri_gg_theme.R')


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

if(!exists("expression.glass.vst")) {
  source('scripts/load_rna-counts.R')
}


# example RF model ----


data(pbc, package = "randomForestSRC")


pbc.obj <- rfsrc(Surv(days,status) ~ ., pbc, importance = TRUE)
#find.interaction(pbc.obj, method = "vimp", nvar = 8)

plot(pbc.obj)


# R2 svvl ----


tmp.metadata <- metadata.glass.per.patient %>%
  tidyr::drop_na(genomescan.sid.R) %>%
  dplyr::select(genomescan.sid.R) %>% 
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid.R'='genomescan.sid')) %>%
  tidyr::drop_na(time.resection.until.last.event) %>% 
  dplyr::rename(genomescan.sid = genomescan.sid.R) %>% 
  dplyr::arrange(Sample_Type, genomescan.sid)



tmp.data <- expression.glass.vst %>%
  dplyr::select(all_of( tmp.metadata$genomescan.sid )) %>% 
  dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>% 
  dplyr::arrange(mad) %>% 
  dplyr::mutate(mad=NULL) 



stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)



svvl.obj <- tmp.metadata %>%
  dplyr::rename(svvl = time.resection.until.last.event) %>% 
  dplyr::rename(status = status.resection.until.last.event) %>%
  dplyr::select(genomescan.sid, svvl, status) %>% 
  dplyr::left_join(
    tmp.data %>%
      dplyr::slice_head(n=2000) %>% 
      t %>%
      as.data.frame(stringsAsFactors=F) %>% 
      tibble::rownames_to_column('genomescan.sid'),
    by=c('genomescan.sid'='genomescan.sid')
  ) %>% 
  tibble::column_to_rownames('genomescan.sid') %>% 
  dplyr::mutate(svvl = as.numeric(svvl)) %>% 
  `colnames<-`(gsub('-','.',colnames(.),fixed=T))




## RF (RFSRC) ----


pbc.obj <- rfsrc(Surv(svvl,status) ~ ., svvl.obj , importance = TRUE,ntree=4000,seed=1234)
plot.rfsrc(pbc.obj)
# plot.variable.rfsrc(pbc.obj)


plt <- data.frame(pbc.obj$importance) %>% 
  dplyr::arrange(abs(pbc.obj.importance)) %>% 
  dplyr::top_n(20) %>% 
  dplyr::arrange(-pbc.obj.importance) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::mutate(baseline = 0) %>% 
  dplyr::mutate(y = rank(pbc.obj.importance)) %>% 
  tidyr::pivot_longer(cols = -c(gene_uid, y)) %>% 
  dplyr::mutate(name=NULL) %>% 
  dplyr::rename(RFSRC.importance = value)

ggplot(plt, aes(x=RFSRC.importance,y=reorder(gene_uid, y), group=gene_uid)) +
  geom_line(lwd=2) +
  youri_gg_theme +
  labs(y=NULL)




## Coxph ----


coxph.res <- RegParallel(
  data = svvl.obj ,
  formula = 'Surv(svvl, status) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(svvl.obj)[3:ncol(svvl.obj)],
  blocksize = 65,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)



plt2 <- data.frame(pbc.obj$importance) %>%
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::left_join(
  coxph.res, by=c('gene_uid'='Variable')
)


# Z
plot(log(plt2$pbc.obj.importance), abs(plt2$Z),pch=19,cex=0.25)
plot(log10(plt2$pbc.obj.importance- min(plt2$pbc.obj.importance) + 0.00001), abs(plt2$Z),pch=19,cex=0.25)


# HR
plot(log(plt2$pbc.obj.importance), abs(log(plt2$HR)),pch=19,cex=0.25)

# LR
plot(log(plt2$pbc.obj.importance), -log10(plt2$LogRank),pch=19,cex=0.25)

# W
plot(log(plt2$pbc.obj.importance), -log10(plt2$Wald),pch=19,cex=0.25)


plt2 <- plt2 %>% 
  #dplyr::mutate(log.importance = log(pbc.obj.importance)) %>% 
  dplyr::mutate(log.importance = log10(plt2$pbc.obj.importance- min(plt2$pbc.obj.importance) + 0.00001)) %>% 
  dplyr::mutate(abs.Z = abs(Z)) %>% 
  dplyr::mutate(m.10log.LR = -log10(LogRank)) %>% 
  dplyr::mutate(m.10log.W = -log10(Wald)) %>% 
  dplyr::mutate(vis = log.importance > -2.4 | log.importance < -4 | m.10log.LR > 6) %>% 
  dplyr::mutate(label = gsub("^[^_]+_","",gene_uid))



ggplot(plt2, aes(x=log.importance, y=abs.Z, label=label)) +
  geom_point(pch=21) +
  ggrepel::geom_text_repel(data = plt2 %>%  dplyr::filter(vis == T)) +
  youri_gg_theme

ggplot(plt2, aes(x=log.importance, y=m.10log.LR, label=label)) +
  geom_point(pch=21) +
  ggrepel::geom_text_repel(data = plt2 %>%  dplyr::filter(vis == T)) +
  youri_gg_theme

ggplot(plt2, aes(x=log.importance, y=m.10log.W, label=label)) +
  geom_point(pch=21) +
  ggrepel::geom_text_repel(data = plt2 %>%  dplyr::filter(vis == T)) +
  youri_gg_theme






# R1 svvl ----

# Van alle resecties tijd tot dood of laatste event nodig



tmp.metadata <- metadata.glass.per.patient %>%
  tidyr::drop_na(genomescan.sid.I) %>%
  dplyr::select(genomescan.sid.I) %>% 
  dplyr::left_join(metadata.glass.per.resection, by=c('genomescan.sid.I'='genomescan.sid')) %>%
  tidyr::drop_na(time.resection.until.last.event) %>% 
  dplyr::rename(genomescan.sid = genomescan.sid.I) %>% 
  dplyr::arrange(Sample_Type, genomescan.sid)



tmp.data <- expression.glass.vst %>%
  dplyr::select(all_of( tmp.metadata$genomescan.sid )) %>% 
  dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>% 
  dplyr::arrange(mad) %>% 
  dplyr::mutate(mad=NULL) 



stopifnot(colnames(tmp.data) == tmp.metadata$genomescan.sid)



svvl.obj <- tmp.metadata %>%
  dplyr::rename(svvl = time.resection.until.last.event) %>% 
  dplyr::rename(status = status.resection.until.last.event) %>%
  dplyr::select(genomescan.sid, svvl, status) %>% 
  dplyr::left_join(
    tmp.data %>%
      dplyr::slice_head(n=10000) %>% 
      t %>%
      as.data.frame(stringsAsFactors=F) %>% 
      tibble::rownames_to_column('genomescan.sid'),
    by=c('genomescan.sid'='genomescan.sid')
  ) %>% 
  tibble::column_to_rownames('genomescan.sid') %>% 
  dplyr::mutate(svvl = as.numeric(svvl)) %>% 
  `colnames<-`(gsub('-','.',colnames(.),fixed=T))





## RF (RFSRC) ----


pbc.obj <- rfsrc(Surv(svvl,status) ~ ., svvl.obj , importance = TRUE,ntree=4000,seed=1234)
plot.rfsrc(pbc.obj)
# plot.variable.rfsrc(pbc.obj)


plt <- data.frame(pbc.obj$importance) %>% 
  dplyr::arrange(abs(pbc.obj.importance)) %>% 
  dplyr::top_n(20) %>% 
  dplyr::arrange(-pbc.obj.importance) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::mutate(baseline = 0) %>% 
  dplyr::mutate(y = rank(pbc.obj.importance)) %>% 
  tidyr::pivot_longer(cols = -c(gene_uid, y)) %>% 
  dplyr::mutate(name=NULL) %>% 
  dplyr::rename(RFSRC.importance = value)

ggplot(plt, aes(x=RFSRC.importance,y=reorder(gene_uid, y), group=gene_uid)) +
  geom_line(lwd=2) +
  youri_gg_theme +
  labs(y=NULL)


## Coxph ----


coxph.res <- RegParallel(
  data = svvl.obj ,
  formula = 'Surv(svvl, status) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(svvl.obj)[3:ncol(svvl.obj)],
  blocksize = 65,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)




