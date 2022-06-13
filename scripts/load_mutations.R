#!/usr/bin/env R 


#' De TumorOnly bevat de pipeline geprocessed zonder matched normal (dus ook de samples waar wel een matched normal van is)
#' De MatchedNormal bevat de pipeline geprocessed met matched normal (voor die samples waar dat bekend is).
#' 
#' De files *fun zijn de files geprocessed met de standaard pipeline (gecalled door LoFreq en Mutect)
#' De files *PON zijn de files geprocessed met de standaard pipeline, gevolgd door een filterstap met het panel of normals van HMF.
#' 
#' De map 'RecurrenceOnly' bevat varianten die wel in de recurrence zitten, maar niet in de primary.
#' 
#' Alle files bevatten annotatie van COSMIC - de TumorOnly zal waarschijnlijk daarop gefilterd moeten worden, of op basis van de MatchedNormal data.
#' Helaas bevatten de TumorOnly samples nog te veel varianten.


# 1. file overzicht ----

a = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/intersect/*intersect.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(TumorOnly.intersect = filename)
a1 = rownames(info(header(VariantAnnotation::readVcf(a[1,1], "hg19"))))

b = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/intersect/*intersect_ann.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(TumorOnly.intersect.ann = filename)
b1 = rownames(info(header(VariantAnnotation::readVcf(b[1,1], "hg19"))))

c = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/intersect/*intersect_fun.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(TumorOnly.intersect.fun = filename)
c1 = rownames(info(header(VariantAnnotation::readVcf(c[1,1], "hg19"))))

d = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/intersect/*intersect_PON.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(filename.TumorOnly.intersect.PON = filename)
d1 = rownames(info(header(VariantAnnotation::readVcf(d[1,1], "hg19"))))

e = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/Mutect2/vcf/filtered/*.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(TumorOnly.Mutect = filename)
e1 = rownames(info(header(VariantAnnotation::readVcf(e[1,1], "hg19"))))

f = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/RecurrenceOnly/*Only.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(TumorOnly.RecurrenceOnly = filename)
f1 = rownames(info(header(VariantAnnotation::readVcf(f[1,1], "hg19"))))

g = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/RecurrenceOnly/*OnlyPON.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(TumourOnly.RecurrenceOnly.PON = filename)
g1 = rownames(info(header(VariantAnnotation::readVcf(g[1,1], "hg19"))))  

h = data.frame(filename = Sys.glob("data/glass/WES/Annotated/*_intersect_fun.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(Annotated.intersect.fun = filename)
h1 = rownames(info(header(VariantAnnotation::readVcf(h[1,1], "hg19"))))  

i = data.frame(filename = Sys.glob("data/glass/WES/Annotated/*_Mutect2_passed_fun.vcf")) %>% 
  dplyr::filter(grepl("118_R2_I2_Mutect2_passed_fun.vcf", filename) == F) %>%  # odd file only containing a header
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(Annotated.Mutect.fun = filename)
i1 = rownames(info(header(VariantAnnotation::readVcf(i[1,1], "hg19"))))  

j = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/Annotated/*_ann.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(MatchedNormal.Annotated.ann = filename)
j1 = rownames(info(header(VariantAnnotation::readVcf(j[1,1], "hg19"))))  

k = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/Annotated/*_fun.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(MatchedNormal.Annotated.fun = filename)
k1 = rownames(info(header(VariantAnnotation::readVcf(k[1,1], "hg19"))))  

l = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/Annotated/*_PON.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(MatchedNormal.Annotated.PON = filename)
l1 = rownames(info(header(VariantAnnotation::readVcf(l[1,1], "hg19"))))  

m = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/intersect/*_intersect.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(MatchedNormal.intersect.intersect = filename)
m1 = rownames(info(header(VariantAnnotation::readVcf(m[1,1], "hg19"))))  

n = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/intersect/*_outersect.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(MatchedNormal.intersect.outersect = filename)
n1 = rownames(info(header(VariantAnnotation::readVcf(n[1,1], "hg19"))))  

o = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/RecurrenceOnly/*Onlyfun.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(MatchedNormal.RecurrenceOnly.fun = filename)
o1 = rownames(info(header(VariantAnnotation::readVcf(o[1,1], "hg19"))))  

p = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/RecurrenceOnly/*OnlyPON.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(MatchedNormal.RecurrenceOnly.PON = filename)
p1 = rownames(info(header(VariantAnnotation::readVcf(p[1,1], "hg19"))))  






z <- a %>% 
  dplyr::full_join(b, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(c, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(d, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(e, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(f, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(g, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(h, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(i, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(j, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(k, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(l, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(m, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(n, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(o, by=c('sample_name'='sample_name')) %>% 
  dplyr::full_join(p, by=c('sample_name'='sample_name')) %>% 
  tibble::column_to_rownames('sample_name')




# a <- VariantAnnotation::readVcfAsVRanges(x = c(
#   "data/glass/WES/Tumor_only/intersect/002-R1-I3_intersect_PON.vcf","data/glass/WES/Tumor_only/intersect/005-R1-I3_intersect_PON.vcf","data/glass/WES/Tumor_only/intersect/007_R2_I2_intersect_PON.vcf","data/glass/WES/Tumor_only/intersect/010-R2-I2_intersect_PON.vcf"
# ), genome = 'hg19', use.names = TRUE)



aa <- VariantAnnotation::readVcfAsVRanges(x = c(
  "data/glass/WES/Tumor_only/intersect/002-R1-I3_intersect_PON.vcf"
), genome = 'hg19') # , use.names = TRUE


bb <- VariantAnnotation::readVcf("data/glass/WES/Tumor_only/intersect/002-R1-I3_intersect_PON.vcf", "hg19")




head(rowRanges(b), 15)
head(geno(b), 15)


gr <- rowRanges(b)
gr[seqnames(gr) == "chr2" & 
     #start(gr) > 209100953 &
     end(gr) < 209119806 ] # & strand(gr) == "."


# 2. MatchedNormal ----


idh1f = function(x) grepl("IDH1", x, fixed=TRUE)

filt.idh = FilterRules(list(idh1 = idh1f))


tmp <- data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/Annotated/*_PON.vcf")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(MatchedNormal.Annotated.PON = filename)
#t1 <- rownames(info(header(VariantAnnotation::readVcf(l[1,1], "hg19")))) 
t1 <- VariantAnnotation::readVcf(tmp[1,1], "hg19")
t2 <- VariantAnnotation::readVcf(tmp[2,1], "hg19")
t3 <- VariantAnnotation::readVcf(tmp[3,1], "hg19")
t4 <- VariantAnnotation::readVcf(tmp[4,1], "hg19")

# Geen IDH[1/2],ATRX,TP53
# grep -v '^#' data/glass/WES/MatchedNormal/variant/Annotated/*_PON.vcf | grep -i -E 'IDH|ATRX|TP53'

# wel IDH1/2
# grep -v '^#' data/glass/WES/MatchedNormal/variant/Annotated/*_fun.vcf | grep -i -E 'IDH|ATRX|TP53'



