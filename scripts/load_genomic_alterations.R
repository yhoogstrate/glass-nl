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


# load data ----


if(!exists("metadata.glass.per.resection")) {
  warning('metadata was not loaded')
  
  source('scripts/load_metadata.R')
}



# load libs ----


library(tidyverse)


# 0. file overzicht ----
# 
# a = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/intersect/*intersect.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(TumorOnly.intersect = filename)
# a1 = rownames(info(header(VariantAnnotation::readVcf(a[1,1], "hg19"))))
# 
# b = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/intersect/*intersect_ann.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(TumorOnly.intersect.ann = filename)
# b1 = rownames(info(header(VariantAnnotation::readVcf(b[1,1], "hg19"))))
# 
# c = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/intersect/*intersect_fun.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(TumorOnly.intersect.fun = filename)
# c1 = rownames(info(header(VariantAnnotation::readVcf(c[1,1], "hg19"))))
# 
# d = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/intersect/*intersect_PON.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(filename.TumorOnly.intersect.PON = filename)
# d1 = rownames(info(header(VariantAnnotation::readVcf(d[1,1], "hg19"))))
# 
# e = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/Mutect2/vcf/filtered/*.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(TumorOnly.Mutect = filename)
# e1 = rownames(info(header(VariantAnnotation::readVcf(e[1,1], "hg19"))))
# 
# f = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/RecurrenceOnly/*Only.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(TumorOnly.RecurrenceOnly = filename)
# f1 = rownames(info(header(VariantAnnotation::readVcf(f[1,1], "hg19"))))
# 
# g = data.frame(filename = Sys.glob("data/glass/WES/Tumor_only/RecurrenceOnly/*OnlyPON.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(TumourOnly.RecurrenceOnly.PON = filename)
# g1 = rownames(info(header(VariantAnnotation::readVcf(g[1,1], "hg19"))))  
# 
# h = data.frame(filename = Sys.glob("data/glass/WES/Annotated/*_intersect_fun.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(Annotated.intersect.fun = filename)
# h1 = rownames(info(header(VariantAnnotation::readVcf(h[1,1], "hg19"))))  
# 
# i = data.frame(filename = Sys.glob("data/glass/WES/Annotated/*_Mutect2_passed_fun.vcf")) %>% 
#   dplyr::filter(grepl("118_R2_I2_Mutect2_passed_fun.vcf", filename) == F) %>%  # odd file only containing a header
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(Annotated.Mutect.fun = filename)
# i1 = rownames(info(header(VariantAnnotation::readVcf(i[1,1], "hg19"))))  
# 
# j = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/Annotated/*_ann.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(MatchedNormal.Annotated.ann = filename)
# j1 = rownames(info(header(VariantAnnotation::readVcf(j[1,1], "hg19"))))  
# 
# k = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/Annotated/*_fun.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(MatchedNormal.Annotated.fun = filename)
# k1 = rownames(info(header(VariantAnnotation::readVcf(k[1,1], "hg19"))))  
# 
# l = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/Annotated/*_PON.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(MatchedNormal.Annotated.PON = filename)
# l1 = rownames(info(header(VariantAnnotation::readVcf(l[1,1], "hg19"))))  
# 
# m = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/intersect/*_intersect.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(MatchedNormal.intersect.intersect = filename)
# m1 = rownames(info(header(VariantAnnotation::readVcf(m[1,1], "hg19"))))  
# 
# n = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/intersect/*_outersect.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(MatchedNormal.intersect.outersect = filename)
# n1 = rownames(info(header(VariantAnnotation::readVcf(n[1,1], "hg19"))))  
# 
# o = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/RecurrenceOnly/*Onlyfun.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(MatchedNormal.RecurrenceOnly.fun = filename)
# o1 = rownames(info(header(VariantAnnotation::readVcf(o[1,1], "hg19"))))  
# 
# p = data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/RecurrenceOnly/*OnlyPON.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
#   dplyr::rename(MatchedNormal.RecurrenceOnly.PON = filename)
# p1 = rownames(info(header(VariantAnnotation::readVcf(p[1,1], "hg19"))))  
# 
# 
# 
# 
# 
# 
# z <- a %>% 
#   dplyr::full_join(b, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(c, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(d, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(e, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(f, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(g, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(h, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(i, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(j, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(k, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(l, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(m, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(n, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(o, by=c('Sample_Name'='Sample_Name')) %>% 
#   dplyr::full_join(p, by=c('Sample_Name'='Sample_Name')) %>% 
#   tibble::column_to_rownames('Sample_Name')
# 



# a <- VariantAnnotation::readVcfAsVRanges(x = c(
#   "data/glass/WES/Tumor_only/intersect/002-R1-I3_intersect_PON.vcf","data/glass/WES/Tumor_only/intersect/005-R1-I3_intersect_PON.vcf","data/glass/WES/Tumor_only/intersect/007_R2_I2_intersect_PON.vcf","data/glass/WES/Tumor_only/intersect/010-R2-I2_intersect_PON.vcf"
# ), genome = 'hg19', use.names = TRUE)


# 
# aa <- VariantAnnotation::readVcfAsVRanges(x = c(
#   "data/glass/WES/Tumor_only/intersect/002-R1-I3_intersect_PON.vcf"
# ), genome = 'hg19') # , use.names = TRUE
# 
# 
# bb <- VariantAnnotation::readVcf("data/glass/WES/Tumor_only/intersect/002-R1-I3_intersect_PON.vcf", "hg19")
# 
# 
# 
# 
# head(rowRanges(b), 15)
# head(geno(b), 15)
# 
# 
# gr <- rowRanges(b)
# gr[seqnames(gr) == "chr2" & 
#      #start(gr) > 209100953 &
#      end(gr) < 209119806 ] # & strand(gr) == "."


# 1. outersects ----

# excluded because of some reason

# parse.outersect <- function(fn) {
#   #fn <- 'data/glass/WES/MatchedNormal/variant/intersect/113_R2_I2_outersect.vcf'
#   data <- read.delim(fn, header=F) %>% 
#     dplyr::rename(chr = V1,
#                   start = V2,
#                   s1 = V3,
#                   s2 = V4) %>% 
#     dplyr::mutate(V5 = NULL) %>% 
#     dplyr::mutate(mutations_filename = fn)
#   
#   return(data)
# }
# 
# mutation.outersect.header <- data.frame(outersect = Sys.glob("data/glass/WES/MatchedNormal/variant/intersect/*outersect.vcf")) %>% 
#   dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",outersect))) 
# 
# 
# mutation.outersect <- mutation.outersect.header %>% 
#   dplyr::pull(`outersect`) %>% 
#   pbapply::pblapply(parse.outersect) %>%
#   dplyr::bind_rows() %>% 
#   dplyr::left_join(mutation.outersect.header, by=c('mutations_filename'='outersect')) %>% 
#   dplyr::mutate(mutation.id = paste0(chr, ":",start,"_",s1 ,"/",s2)) %>% 
#   dplyr::group_by(mutation.id) %>% 
#   dplyr::tally(sort=T) %>% 
#   dplyr::rename(outsect.n = n)
# 
# 
# stopifnot(duplicated(mutation.outersect$mutation.id) == F) # ensure uniqueness


# this one was found 6 times while present in filtered other ones
# chr7:100642594_A/C         6




# 2. normals blood ----

# check en integrate met CDC27 - zijn allemaal 'NormalOnly'
# outersect files!?




# 3. MatchedNormal ----


mutation.data.mn.header <- data.frame(filename = Sys.glob("data/glass/WES/vcf_matched-normal/*/Annotated/*_fun.vcf.gz")) %>% 
  dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(MatchedNormal.Annotated.fun = filename)



if(!exists('mutation.data.mn')) {
  if(file.exists('cache/mutation.data.mn.Rds')) {
    
    print("Loading 'mutation.data.mn' from cache.")
    mutation.data.mn <- readRDS('cache/mutation.data.mn.Rds')
    
  } else {
  
    parse.vcf <- function(fn) {
      #fn = "data/glass/WES/vcf_matched-normal/variant/Annotated/113_R2_I2_intersect_fun.vcf.gz"
      tmp <- VariantAnnotation::readVcf(fn, "hg19")
      
      tmp.ref <- data.frame(
        chr = as.character(tmp@rowRanges@seqnames),
        start = tmp@rowRanges@ranges@start,
        ref = as.character(tmp@fixed$REF),
        alt = as.character(unlist(tmp@fixed$ALT))
      )
      
      tmp.info <- tmp %>% 
        VariantAnnotation::info() %>% 
        as.data.frame %>% 
        dplyr::mutate(`FUNCOTATION` = unlist(`FUNCOTATION`)) %>% 
        dplyr::mutate(FUNCOTATION = gsub('^\\[','',FUNCOTATION)) %>% 
        dplyr::mutate(FUNCOTATION = gsub('\\]$','',FUNCOTATION)) %>% 
        dplyr::mutate(FUNCOTATION = iconv(FUNCOTATION, "latin1", "ASCII", sub="")) %>%  # weird non ascii character in linked disease name(s)
        tidyr::separate(FUNCOTATION, sep="\\|", into = c("Gencode_19_hugoSymbol","Gencode_19_ncbiBuild","Gencode_19_chromosome","Gencode_19_start","Gencode_19_end","Gencode_19_variantClassification","Gencode_19_secondaryVariantClassification","Gencode_19_variantType","Gencode_19_refAllele","Gencode_19_tumorSeqAllele1","Gencode_19_tumorSeqAllele2","Gencode_19_genomeChange","Gencode_19_annotationTranscript","Gencode_19_transcriptStrand","Gencode_19_transcriptExon","Gencode_19_transcriptPos","Gencode_19_cDnaChange","Gencode_19_codonChange","Gencode_19_proteinChange","Gencode_19_gcContent","Gencode_19_referenceContext","Gencode_19_otherTranscripts","CGC_Name","CGC_GeneID","CGC_Chr","CGC_Chr_Band","CGC_Cancer_Somatic_Mut","CGC_Cancer_Germline_Mut","CGC_Tumour_Types__(Somatic_Mutations)","CGC_Tumour_Types_(Germline_Mutations)","CGC_Cancer_Syndrome","CGC_Tissue_Type","CGC_Cancer_Molecular_Genetics","CGC_Mutation_Type","CGC_Translocation_Partner","CGC_Other_Germline_Mut","CGC_Other_Syndrome/Disease","ClinVar_VCF_AF_ESP","ClinVar_VCF_AF_EXAC","ClinVar_VCF_AF_TGP","ClinVar_VCF_ALLELEID","ClinVar_VCF_CLNDISDB","ClinVar_VCF_CLNDISDBINCL","ClinVar_VCF_CLNDN","ClinVar_VCF_CLNDNINCL","ClinVar_VCF_CLNHGVS","ClinVar_VCF_CLNREVSTAT","ClinVar_VCF_CLNSIG","ClinVar_VCF_CLNSIGCONF","ClinVar_VCF_CLNSIGINCL","ClinVar_VCF_CLNVC","ClinVar_VCF_CLNVCSO","ClinVar_VCF_CLNVI","ClinVar_VCF_DBVARID","ClinVar_VCF_GENEINFO","ClinVar_VCF_MC","ClinVar_VCF_ORIGIN","ClinVar_VCF_RS","ClinVar_VCF_SSR","ClinVar_VCF_ID","ClinVar_VCF_FILTER","Cosmic_overlapping_mutations","HGNC_HGNC_ID","HGNC_Approved_Name","HGNC_Status","HGNC_Locus_Type","HGNC_Locus_Group","HGNC_Previous_Symbols","HGNC_Previous_Name","HGNC_Synonyms","HGNC_Name_Synonyms","HGNC_Chromosome","HGNC_Date_Modified","HGNC_Date_Symbol_Changed","HGNC_Date_Name_Changed","HGNC_Accession_Numbers","HGNC_Enzyme_IDs","HGNC_Entrez_Gene_ID","HGNC_Ensembl_Gene_ID","HGNC_Pubmed_IDs","HGNC_RefSeq_IDs","HGNC_Gene_Family_ID","HGNC_Gene_Family_Name","HGNC_CCDS_IDs","HGNC_Vega_ID","HGNC_Entrez_Gene_ID(supplied_by_NCBI)","HGNC_OMIM_ID(supplied_by_OMIM)","HGNC_RefSeq(supplied_by_NCBI)","HGNC_UniProt_ID(supplied_by_UniProt)","HGNC_Ensembl_ID(supplied_by_Ensembl)","HGNC_UCSC_ID(supplied_by_UCSC)","dbSNP_ASP","dbSNP_ASS","dbSNP_CAF","dbSNP_CDA","dbSNP_CFL","dbSNP_COMMON","dbSNP_DSS","dbSNP_G5","dbSNP_G5A","dbSNP_GENEINFO","dbSNP_GNO","dbSNP_HD","dbSNP_INT","dbSNP_KGPhase1","dbSNP_KGPhase3","dbSNP_LSD","dbSNP_MTP","dbSNP_MUT","dbSNP_NOC","dbSNP_NOV","dbSNP_NSF","dbSNP_NSM","dbSNP_NSN","dbSNP_OM","dbSNP_OTH","dbSNP_PM","dbSNP_PMC","dbSNP_R3","dbSNP_R5","dbSNP_REF","dbSNP_RS","dbSNP_RSPOS","dbSNP_RV","dbSNP_S3D","dbSNP_SAO","dbSNP_SLO","dbSNP_SSR","dbSNP_SYN","dbSNP_TOPMED","dbSNP_TPA","dbSNP_U3","dbSNP_U5","dbSNP_VC","dbSNP_VLD","dbSNP_VP","dbSNP_WGT","dbSNP_WTD","dbSNP_dbSNPBuildID","dbSNP_ID","dbSNP_FILTER")  ) %>% 
        dplyr::mutate(FUNCOTATION = NULL) %>% 
        dplyr::mutate(mutations_filename = fn) %>% 
        tibble::rownames_to_column('name') %>% 
        dplyr::select(c('name',
                        'AA','AC','AS_SB_TABLE','CDS','CNT','DB','DP','ECNT','GENE','Gencode_19_hugoSymbol','Gencode_19_chromosome',
                        'Gencode_19_start','Gencode_19_end','Gencode_19_variantClassification','Gencode_19_secondaryVariantClassification',
                        'Gencode_19_proteinChange','CGC_Name','CGC_Cancer_Somatic_Mut','CGC_Cancer_Germline_Mut',
                        'CGC_Tumour_Types__(Somatic_Mutations)','CGC_Mutation_Type','ClinVar_VCF_CLNSIG','ClinVar_VCF_MC','dbSNP_ID','MBQ',
                        'MFRL','MMQ','MPOS','NALOD','NLOD','PON','PON_COUNT','POPAF','TLOD','gnomAD_AF','mutations_filename'
        )) %>% 
        dplyr::mutate(PON_LoFreq_2_0_025 = NA, PON_LoFreq_2_0_05 = NA, PON_LoFreq_3_0_025 = NA, PON_LoFreq_3_0_05 = NA)
      
      tmp.geno <- tmp %>% 
        VariantAnnotation::geno()
      
      tmp.geno.ad1 <- tmp.geno$AD[,1] %>%  purrr::map_chr(c(1))
      tmp.geno.ad2 <- tmp.geno$AD[,1] %>%  purrr::map_chr(c(2))
      stopifnot(names(tmp.geno.ad1) == names(tmp.geno.ad2))
      
      tmp.geno.df <- data.frame(FILTER = tmp@fixed$FILTER,
                                DP.ref = as.numeric(tmp.geno.ad1),
                                DP.var = as.numeric(tmp.geno.ad2),
                                AF.geno = as.numeric(tmp.geno$AF[,1])) %>% 
        dplyr::mutate(DP.geno = DP.var + DP.ref) %>% 
        dplyr::mutate(VAF = DP.var / DP.geno)
      
      
      tmp.info <- cbind(tmp.ref, tmp.info, tmp.geno.df) %>% 
        dplyr::mutate(POPAF = as.numeric(POPAF))
      
      return(tmp.info)
    }
    
  
    mutation.data.mn <- mutation.data.mn.header %>% 
      dplyr::pull(`MatchedNormal.Annotated.fun`) %>% 
      pbapply::pblapply(parse.vcf) %>% 
      dplyr::bind_rows() %>% 
      dplyr::left_join(mutation.data.mn.header, by=c('mutations_filename'='MatchedNormal.Annotated.fun'))
    
    
    
    mutation.data.mn %>% 
      dplyr::mutate(pid = gsub("^([0-9]+).+$","\\1",Sample_Name)) %>% 
      dplyr::select(Sample_Name, Gencode_19_hugoSymbol) %>% 
      dplyr::distinct(Sample_Name, Gencode_19_hugoSymbol) %>% 
      count(`Gencode_19_hugoSymbol`) %>%
      dplyr::arrange(-n) %>%
      head(n=25)
    
    
    
    
    mutation.data.mn %>%
      dplyr::filter(Gencode_19_hugoSymbol == "IDH1") %>% 
      dplyr::select(name, FILTER, PON_COUNT, dbSNP_ID, Sample_Name,VAF, Gencode_19_chromosome, Gencode_19_start, Gencode_19_end) 
    
    
    
    mutation.data.mn %>%
      dplyr::filter(Gencode_19_hugoSymbol == "ZNF131") %>%
      dplyr::select(name, FILTER, dbSNP_ID, Sample_Name,VAF, Gencode_19_chromosome, Gencode_19_start, Gencode_19_end)
    
    
    saveRDS(mutation.data.mn, file='cache/mutation.data.mn.Rds')
  }
}


# 4. TumorOnly ----




#mutation.data.to.header <- data.frame(filename = Sys.glob("data/glass/WES/TumorOnly/intersect/*_fun.vcf.gz")) %>% 
mutation.data.to.header <- data.frame(filename = Sys.glob("data/glass/WES/vcf_tumor-only//intersect_with_additional_LoFreq_PON/*_fun_LoFreqPONs.vcf.gz")) %>% 
  dplyr::mutate(Sample_Name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::filter(Sample_Name %in% mutation.data.mn.header$Sample_Name == F) %>% # exclude those present in both
  dplyr::rename(TumorOnly.intersect.fun = filename)

stopifnot(mutation.data.to.header$Sample_Name %in% mutation.data.mn.header$Sample_Name == F)



if(!exists('mutation.data.to')) {
  if(file.exists('cache/mutation.data.to.Rds')) {
    
    print("Loading 'mutation.data.to' from cache.")
    mutation.data.mn <- readRDS('cache/mutation.data.mn.Rds')
    
  } else {
    
    parse.vcf.2  <- function(fn) {
      #fn <- 'data/glass/WES/TumorOnly/intersect/002-R1-I3_intersect_fun.vcf.gz'
      #fn <- 'data/glass/WES/TumorOnly/intersect_with_additional_LoFreq_PON/219-R3-I1_intersect_fun_LoFreqPONs.vcf.gz'
      tmp <- VariantAnnotation::readVcf(fn, "hg19")
      
      tmp.ref <- data.frame(
        chr = as.character(tmp@rowRanges@seqnames),
        start = tmp@rowRanges@ranges@start,
        ref = as.character(tmp@fixed$REF),
        alt = as.character(unlist(tmp@fixed$ALT))
      )
      
      tmp.info <- tmp %>% 
        VariantAnnotation::info() %>% 
        as.data.frame %>% 
        dplyr::mutate(`FUNCOTATION` = unlist(`FUNCOTATION`)) %>% 
        dplyr::mutate(FUNCOTATION = gsub('^\\[','',FUNCOTATION)) %>% 
        dplyr::mutate(FUNCOTATION = gsub('\\]$','',FUNCOTATION)) %>% 
        dplyr::mutate(FUNCOTATION = iconv(FUNCOTATION, "latin1", "ASCII", sub="")) %>%  # weird non ascii character in linked disease name(s)
        tidyr::separate(FUNCOTATION, sep="\\|", into = c("Gencode_19_hugoSymbol","Gencode_19_ncbiBuild","Gencode_19_chromosome","Gencode_19_start","Gencode_19_end","Gencode_19_variantClassification","Gencode_19_secondaryVariantClassification","Gencode_19_variantType","Gencode_19_refAllele","Gencode_19_tumorSeqAllele1","Gencode_19_tumorSeqAllele2","Gencode_19_genomeChange","Gencode_19_annotationTranscript","Gencode_19_transcriptStrand","Gencode_19_transcriptExon","Gencode_19_transcriptPos","Gencode_19_cDnaChange","Gencode_19_codonChange","Gencode_19_proteinChange","Gencode_19_gcContent","Gencode_19_referenceContext","Gencode_19_otherTranscripts","CGC_Name","CGC_GeneID","CGC_Chr","CGC_Chr_Band","CGC_Cancer_Somatic_Mut","CGC_Cancer_Germline_Mut","CGC_Tumour_Types__(Somatic_Mutations)","CGC_Tumour_Types_(Germline_Mutations)","CGC_Cancer_Syndrome","CGC_Tissue_Type","CGC_Cancer_Molecular_Genetics","CGC_Mutation_Type","CGC_Translocation_Partner","CGC_Other_Germline_Mut","CGC_Other_Syndrome/Disease","ClinVar_VCF_AF_ESP","ClinVar_VCF_AF_EXAC","ClinVar_VCF_AF_TGP","ClinVar_VCF_ALLELEID","ClinVar_VCF_CLNDISDB","ClinVar_VCF_CLNDISDBINCL","ClinVar_VCF_CLNDN","ClinVar_VCF_CLNDNINCL","ClinVar_VCF_CLNHGVS","ClinVar_VCF_CLNREVSTAT","ClinVar_VCF_CLNSIG","ClinVar_VCF_CLNSIGCONF","ClinVar_VCF_CLNSIGINCL","ClinVar_VCF_CLNVC","ClinVar_VCF_CLNVCSO","ClinVar_VCF_CLNVI","ClinVar_VCF_DBVARID","ClinVar_VCF_GENEINFO","ClinVar_VCF_MC","ClinVar_VCF_ORIGIN","ClinVar_VCF_RS","ClinVar_VCF_SSR","ClinVar_VCF_ID","ClinVar_VCF_FILTER","Cosmic_overlapping_mutations","HGNC_HGNC_ID","HGNC_Approved_Name","HGNC_Status","HGNC_Locus_Type","HGNC_Locus_Group","HGNC_Previous_Symbols","HGNC_Previous_Name","HGNC_Synonyms","HGNC_Name_Synonyms","HGNC_Chromosome","HGNC_Date_Modified","HGNC_Date_Symbol_Changed","HGNC_Date_Name_Changed","HGNC_Accession_Numbers","HGNC_Enzyme_IDs","HGNC_Entrez_Gene_ID","HGNC_Ensembl_Gene_ID","HGNC_Pubmed_IDs","HGNC_RefSeq_IDs","HGNC_Gene_Family_ID","HGNC_Gene_Family_Name","HGNC_CCDS_IDs","HGNC_Vega_ID","HGNC_Entrez_Gene_ID(supplied_by_NCBI)","HGNC_OMIM_ID(supplied_by_OMIM)","HGNC_RefSeq(supplied_by_NCBI)","HGNC_UniProt_ID(supplied_by_UniProt)","HGNC_Ensembl_ID(supplied_by_Ensembl)","HGNC_UCSC_ID(supplied_by_UCSC)","dbSNP_ASP","dbSNP_ASS","dbSNP_CAF","dbSNP_CDA","dbSNP_CFL","dbSNP_COMMON","dbSNP_DSS","dbSNP_G5","dbSNP_G5A","dbSNP_GENEINFO","dbSNP_GNO","dbSNP_HD","dbSNP_INT","dbSNP_KGPhase1","dbSNP_KGPhase3","dbSNP_LSD","dbSNP_MTP","dbSNP_MUT","dbSNP_NOC","dbSNP_NOV","dbSNP_NSF","dbSNP_NSM","dbSNP_NSN","dbSNP_OM","dbSNP_OTH","dbSNP_PM","dbSNP_PMC","dbSNP_R3","dbSNP_R5","dbSNP_REF","dbSNP_RS","dbSNP_RSPOS","dbSNP_RV","dbSNP_S3D","dbSNP_SAO","dbSNP_SLO","dbSNP_SSR","dbSNP_SYN","dbSNP_TOPMED","dbSNP_TPA","dbSNP_U3","dbSNP_U5","dbSNP_VC","dbSNP_VLD","dbSNP_VP","dbSNP_WGT","dbSNP_WTD","dbSNP_dbSNPBuildID","dbSNP_ID","dbSNP_FILTER")) %>% 
        dplyr::mutate(FUNCOTATION = NULL) %>% 
        dplyr::mutate(mutations_filename = fn) %>% 
        dplyr::mutate(PON_LoFreq_2_0_025 = unlist(ifelse(as.character(PON_LoFreq_2_0_025) == "character(0)", NA, PON_LoFreq_2_0_025))) %>% 
        dplyr::mutate(PON_LoFreq_2_0_05 = unlist(ifelse(as.character(PON_LoFreq_2_0_05) == "character(0)", NA, PON_LoFreq_2_0_05))) %>% 
        dplyr::mutate(PON_LoFreq_3_0_025 = unlist(ifelse(as.character(PON_LoFreq_3_0_025) == "character(0)", NA, PON_LoFreq_2_0_025))) %>% 
        dplyr::mutate(PON_LoFreq_3_0_05 = unlist(ifelse(as.character(PON_LoFreq_3_0_05) == "character(0)", NA, PON_LoFreq_2_0_05))) %>% 
        tibble::rownames_to_column('name') %>% 
        dplyr::select(c('name',
                        'AA','AC','AS_SB_TABLE','CDS','CNT','DB','DP','ECNT','GENE','Gencode_19_hugoSymbol','Gencode_19_chromosome',
                        'Gencode_19_start','Gencode_19_end','Gencode_19_variantClassification','Gencode_19_secondaryVariantClassification',
                        'Gencode_19_proteinChange','CGC_Name','CGC_Cancer_Somatic_Mut','CGC_Cancer_Germline_Mut',
                        'CGC_Tumour_Types__(Somatic_Mutations)','CGC_Mutation_Type','ClinVar_VCF_CLNSIG','ClinVar_VCF_MC','dbSNP_ID','MBQ',
                        'MFRL','MMQ','MPOS','NALOD','NLOD','PON','PON_COUNT','POPAF','TLOD','gnomAD_AF','mutations_filename',
                        'PON_LoFreq_2_0_025', 'PON_LoFreq_2_0_05', 'PON_LoFreq_3_0_025', 'PON_LoFreq_3_0_05'
        ))
      
      tmp.geno <- tmp %>% 
        VariantAnnotation::geno()
      
      tmp.geno.ad1 <- tmp.geno$AD[,1] %>%  purrr::map_chr(c(1))
      tmp.geno.ad2 <- tmp.geno$AD[,1] %>%  purrr::map_chr(c(2))
      stopifnot(names(tmp.geno.ad1) == names(tmp.geno.ad2))
      
      tmp.geno.df <- data.frame(FILTER = tmp@fixed$FILTER,
                                DP.ref = as.numeric(tmp.geno.ad1),
                                DP.var = as.numeric(tmp.geno.ad2),
                                AF.geno = as.numeric(tmp.geno$AF[,1])) %>% 
        dplyr::mutate(DP.geno = DP.var + DP.ref) %>% 
        dplyr::mutate(VAF = DP.var / DP.geno)
      
      
      tmp.info <- cbind(tmp.ref, tmp.info, tmp.geno.df)  %>% 
        dplyr::mutate(POPAF = as.numeric(POPAF))
      
      return(tmp.info)
    }
    
    
    
    mutation.data.to <- mutation.data.to.header %>% 
      dplyr::pull(`TumorOnly.intersect.fun`) %>% 
      pbapply::pblapply(parse.vcf.2) %>%
      dplyr::bind_rows() %>% 
      dplyr::left_join(mutation.data.to.header, by=c('mutations_filename'='TumorOnly.intersect.fun'))
    
    
    saveRDS(mutation.data.to, file='cache/mutation.data.to.Rds')

  }
}




mutation.data.to %>%
  dplyr::filter(Gencode_19_hugoSymbol == "PWRN1") %>% 
  dplyr::select(name, FILTER, dbSNP_ID, Sample_Name, VAF, Gencode_19_chromosome, Gencode_19_start, Gencode_19_end) 



mutatmutation.data.to %>%
  dplyr::filter(Gencode_19_hugoSymbol == "ZNF131") %>%
  dplyr::select(name, FILTER, dbSNP_ID, Sample_Name, VAF, Gencode_19_chromosome, Gencode_19_start, Gencode_19_end)


mutation.data.to %>%
  dplyr::filter(Gencode_19_hugoSymbol == "ATRX") %>%
  dplyr::select(name, FILTER, dbSNP_ID, Sample_Name, VAF, Gencode_19_hugoSymbol, AA, Gencode_19_chromosome, Gencode_19_start, Gencode_19_end) %>% 
  dplyr::filter(FILTER != 'germline') %>% 
  dplyr::filter(grepl('COSM',name) | !is.na(AA))


mutation.data.to %>%
  dplyr::filter(grepl("^H3",Gencode_19_hugoSymbol)) %>%
  dplyr::select(name, FILTER, dbSNP_ID, Sample_Name, VAF, Gencode_19_hugoSymbol, AA,
                Gencode_19_chromosome, Gencode_19_start,
                
                CGC_Cancer_Somatic_Mut,
                CGC_Cancer_Germline_Mut,
                `CGC_Tumour_Types__(Somatic_Mutations)`
                )

# In glioblastomas of young adults and pediatric patients, some studies have identified a
# small percentage of patients who are IDH wild-type and have a loss of ATRX expression (52–54).
# Ebrahimi and colleagues found these patients to have H3F3A G34 or K27 mutations, which is concordant with Ikemura and colleagues’ finding of ATRX-loss glioblastomas in younger patients being most commonly non-hemispheric in location (47, 52)



# 5. merge ----


stopifnot(colnames(mutation.data.to) == colnames(mutation.data.mn))


mutation.data <- rbind(mutation.data.to %>% dplyr::mutate(type = "TumorOnly"),
                       mutation.data.mn %>% dplyr::mutate(type = "MatchingNormal")
                       ) %>% 
  dplyr::filter(FILTER != 'germline') %>% 
  dplyr::mutate(mutation.id = paste0(chr,":",start,"_",ref,"/",alt)) 
  #dplyr::left_join(mutation.outersect, by=c('mutation.id' = 'mutation.id')) # outersect seems to contain practically all IDH1 muts?



## a. find top hit genes ----


mutation.data %>% 
  dplyr::mutate(pid = gsub("^([0-9]+).+$","\\1",Sample_Name)) %>% 
  dplyr::filter(is.na(PON_COUNT) | PON_COUNT < 1) %>% 
  dplyr::filter(Gencode_19_variantClassification != 'SILENT') %>% 
  dplyr::filter(Gencode_19_variantClassification != 'INTRON') %>% 
  dplyr::filter(Gencode_19_variantClassification != 'RNA') %>% 
  dplyr::filter(Gencode_19_variantClassification != "FIVE_PRIME_FLANK") %>% 
  dplyr::filter(Gencode_19_variantClassification != "THREE_PRIME_UTR") %>% 
  dplyr::filter(Gencode_19_variantClassification != "FIVE_PRIME_UTR" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "LINCRNA" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "COULD_NOT_DETERMINE" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "DE_NOVO_START_OUT_FRAME" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "DE_NOVO_START_IN_FRAME" ) %>% 
  dplyr::filter(is.na(PON_LoFreq_2_0_025)) %>%  # PON_LoFreq_3_0_05
  #dplyr::filter(is.na(outsect.n) | outsect.n <= 5) %>% 
  dplyr::select(Sample_Name, Gencode_19_hugoSymbol) %>% 
  dplyr::distinct(Sample_Name, Gencode_19_hugoSymbol) %>% 
  count(`Gencode_19_hugoSymbol`) %>%
  dplyr::arrange(-n) %>%
  head(n=55)





### CDC27 ----

cdc27 <- mutation.data %>% 
  dplyr::mutate(pid = gsub("^([0-9]+).+$","\\1",Sample_Name)) %>% 
  dplyr::filter(is.na(PON_COUNT) | PON_COUNT < 1) %>% 
  dplyr::filter(Gencode_19_variantClassification != 'SILENT') %>% 
  dplyr::filter(Gencode_19_variantClassification != 'INTRON') %>% 
  dplyr::filter(Gencode_19_variantClassification != 'RNA') %>% 
  dplyr::filter(Gencode_19_variantClassification != "FIVE_PRIME_FLANK") %>% 
  dplyr::filter(Gencode_19_variantClassification != "THREE_PRIME_UTR") %>% 
  dplyr::filter(Gencode_19_variantClassification != "FIVE_PRIME_UTR" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "LINCRNA" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "COULD_NOT_DETERMINE" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "DE_NOVO_START_OUT_FRAME" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "DE_NOVO_START_IN_FRAME" ) %>% 
  dplyr::filter(Gencode_19_hugoSymbol == 'CDC27')


cdc27 %>% 
  dplyr::group_by(Gencode_19_start) %>% 
  dplyr::tally(sort=T)


cdc27 %>% 
  dplyr::filter(Gencode_19_start == "45234707") %>% 
  dplyr::select(name, Sample_Name, Gencode_19_chromosome, Gencode_19_start, PON_COUNT, type, Sample_Name, PON_LoFreq_3_0_05)



### FRG1 ----


frg1 <- mutation.data %>% 
  dplyr::mutate(pid = gsub("^([0-9]+).+$","\\1",Sample_Name)) %>% 
  dplyr::filter(is.na(PON_COUNT) | PON_COUNT < 1) %>% 
  dplyr::filter(Gencode_19_variantClassification != 'SILENT') %>% 
  dplyr::filter(Gencode_19_variantClassification != 'INTRON') %>% 
  dplyr::filter(Gencode_19_variantClassification != 'RNA') %>% 
  dplyr::filter(Gencode_19_variantClassification != "FIVE_PRIME_FLANK") %>% 
  dplyr::filter(Gencode_19_variantClassification != "THREE_PRIME_UTR") %>% 
  dplyr::filter(Gencode_19_variantClassification != "FIVE_PRIME_UTR" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "LINCRNA" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "COULD_NOT_DETERMINE" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "DE_NOVO_START_OUT_FRAME" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "DE_NOVO_START_IN_FRAME" ) %>% 
  dplyr::filter(Gencode_19_hugoSymbol == 'FRG1')


frg1 %>% 
  dplyr::group_by(Gencode_19_start) %>% 
  dplyr::tally(sort=T)


frg1 %>% 
  dplyr::filter(Gencode_19_start == "190876263") %>% 
  dplyr::select(name, Sample_Name, Gencode_19_chromosome, Gencode_19_start ,type, PON_COUNT, PON_LoFreq_3_0_05)




### PIK3CA ----



pik3ca <- mutation.data %>% 
  dplyr::mutate(pid = gsub("^([0-9]+).+$","\\1",Sample_Name)) %>% 
  dplyr::filter(is.na(PON_COUNT) | PON_COUNT < 1) %>% 
  dplyr::filter(Gencode_19_variantClassification != 'SILENT') %>% 
  dplyr::filter(Gencode_19_variantClassification != 'INTRON') %>% 
  dplyr::filter(Gencode_19_variantClassification != 'RNA') %>% 
  dplyr::filter(Gencode_19_variantClassification != "FIVE_PRIME_FLANK") %>% 
  dplyr::filter(Gencode_19_variantClassification != "THREE_PRIME_UTR") %>% 
  dplyr::filter(Gencode_19_variantClassification != "FIVE_PRIME_UTR" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "LINCRNA" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "COULD_NOT_DETERMINE" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "DE_NOVO_START_OUT_FRAME" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "DE_NOVO_START_IN_FRAME" ) %>% 
  dplyr::filter(Gencode_19_hugoSymbol == 'PIK3CA')


pik3ca %>% 
  dplyr::group_by(Gencode_19_start) %>% 
  dplyr::tally(sort=T)



pik3ca %>% 
#  dplyr::filter(Gencode_19_start == "22023422") %>% 
  dplyr::select(name, Sample_Name, Gencode_19_chromosome, Gencode_19_start ,type, PON_COUNT, PON_LoFreq_2_0_025, PON_LoFreq_3_0_05)


### MUC12 ----

muc12.old <- muc12


muc12 <- mutation.data %>% 
  dplyr::mutate(pid = gsub("^([0-9]+).+$","\\1",Sample_Name)) %>% 
  dplyr::filter(is.na(PON_COUNT) | PON_COUNT < 1) %>% 
  dplyr::filter(Gencode_19_variantClassification != 'SILENT') %>% 
  dplyr::filter(Gencode_19_variantClassification != 'INTRON') %>% 
  dplyr::filter(Gencode_19_variantClassification != 'RNA') %>% 
  dplyr::filter(Gencode_19_variantClassification != "FIVE_PRIME_FLANK") %>% 
  dplyr::filter(Gencode_19_variantClassification != "THREE_PRIME_UTR") %>% 
  dplyr::filter(Gencode_19_variantClassification != "FIVE_PRIME_UTR" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "LINCRNA" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "COULD_NOT_DETERMINE" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "DE_NOVO_START_OUT_FRAME" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "DE_NOVO_START_IN_FRAME" ) %>% 
  dplyr::filter(Gencode_19_hugoSymbol == 'MUC12') %>% 
  dplyr::filter(is.na(PON_LoFreq_2_0_025))  # PON_LoFreq_3_0_05



a = muc12 %>% dplyr::mutate(fn = basename(mutations_filename )) %>%  dplyr::select(chr, start, ref, alt, fn) %>% dplyr::mutate(id = paste0(fn , "_", start, "_", alt))
b = muc12.old %>% dplyr::mutate(fn = basename(mutations_filename )) %>%  dplyr::select(chr, start, ref, alt, fn) %>% dplyr::mutate(id = paste0(fn , "_", start, "_", alt))

muc12 %>% 
  dplyr::group_by(start) %>% 
  dplyr::tally(sort=T)




### IDH(1/2) ----


idh <- mutation.data %>% 
  dplyr::mutate(pid = gsub("^([0-9]+).+$","\\1",Sample_Name)) %>% 
  dplyr::filter(is.na(PON_COUNT) | PON_COUNT < 1) %>% 
  dplyr::filter(Gencode_19_variantClassification != 'SILENT') %>% 
  dplyr::filter(Gencode_19_variantClassification != 'INTRON') %>% 
  dplyr::filter(Gencode_19_variantClassification != 'RNA') %>% 
  dplyr::filter(Gencode_19_variantClassification != "FIVE_PRIME_FLANK") %>% 
  dplyr::filter(Gencode_19_variantClassification != "THREE_PRIME_UTR") %>% 
  dplyr::filter(Gencode_19_variantClassification != "FIVE_PRIME_UTR" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "LINCRNA" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "COULD_NOT_DETERMINE" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "DE_NOVO_START_OUT_FRAME" ) %>% 
  dplyr::filter(Gencode_19_variantClassification != "DE_NOVO_START_IN_FRAME" ) %>% 
  dplyr::filter(grepl("^IDH",Gencode_19_hugoSymbol))


idh %>% 
  dplyr::select(name, Sample_Name, Gencode_19_chromosome, Gencode_19_start ,type, PON_COUNT, PON_LoFreq_3_0_05)



### :: ----

mutation.data.to %>%
  dplyr::filter(Gencode_19_hugoSymbol == "IDH1")  %>% 
  View


mutation.data.to %>%
  dplyr::filter(Gencode_19_hugoSymbol == "CDC27")  %>% 
  View

mutation.data %>%
  dplyr::filter(Gencode_19_hugoSymbol %in% c("IDH1","IDH2","TP53","ATRX")) %>% 
  View


mutation.data %>% 
  dplyr::filter(is.na(AA) & is.na(Gencode_19_proteinChange))


unique(mutation.data$Gencode_19_variantClassification)
#[1] "MISSENSE"                "SILENT"                  "IN_FRAME_INS"            "INTRON"                 
#[5] "SPLICE_SITE"             "NONSENSE"                "IGR"                     "RNA"                    
#[9] "FIVE_PRIME_FLANK"        "THREE_PRIME_UTR"         "FIVE_PRIME_UTR"          "FRAME_SHIFT_INS"        
#[13] "FRAME_SHIFT_DEL"         "NONSTOP"                 "IN_FRAME_DEL"            "LINCRNA"                
#[17] "COULD_NOT_DETERMINE"     "START_CODON_SNP"         "DE_NOVO_START_OUT_FRAME" "DE_NOVO_START_IN_FRAME" 
#[21] "START_CODON_DEL"         "START_CODON_INS"  




# 6. IDH-mut status (manual curate) ----




tmp <- mutation.data %>% 
  dplyr::mutate(pid = gsub("^([0-9]+).+$","\\1",Sample_Name)) %>% 
  dplyr::mutate(CGC_Name = NULL) %>% 
  dplyr::filter(Gencode_19_hugoSymbol %in% c('IDH1','IDH2')) %>% 
  dplyr::mutate(IDH.mutation.WES = paste0(Gencode_19_hugoSymbol, " ",gsub("^p.","",AA))) %>% 
  dplyr::select(Sample_Name, IDH.mutation.WES) %>% 
  dplyr::mutate(evidence = "VCF")



tmp.manual <- data.frame() %>%
  rbind(data.frame('Sample_Name'='018-R1',IDH.mutation.WES=NA,evidence=NA)) %>% #' @todo x-check bam
  rbind(data.frame('Sample_Name'='020-R3',IDH.mutation.WES=NA,evidence=NA)) %>% #' @todo x-check bam
  rbind(data.frame('Sample_Name'='128-R2',IDH.mutation.WES='IDH1 R132S',evidence='matching R1')) %>%  #' @todo x-check bam
  rbind(data.frame('Sample_Name'='129-R2',IDH.mutation.WES='IDH1 R132H',evidence='matching R1 and R3')) %>% #' @todo x-check bam
  rbind(data.frame('Sample_Name'='135-R1',IDH.mutation.WES='IDH1 R132C',evidence='matching R2')) %>% #' @todo x-check bam
  rbind(data.frame('Sample_Name'='141-R2',IDH.mutation.WES='IDH1 R132H',evidence='matching R1')) %>% #' @todo x-check bam
  rbind(data.frame('Sample_Name'='157-R1',IDH.mutation.WES='IDH1 R132H',evidence='present in TumorOnly?')) %>%
  rbind(data.frame('Sample_Name'='157-R2',IDH.mutation.WES='IDH1 R132H',evidence='present in TumorOnly?')) %>%
  rbind(data.frame('Sample_Name'='158-R1',IDH.mutation.WES='IDH1 R132H',evidence='present in TumorOnly?')) %>%
  rbind(data.frame('Sample_Name'='158-R2',IDH.mutation.WES='IDH1 R132H',evidence='present in TumorOnly?')) %>%
  rbind(data.frame('Sample_Name'='160-R1',IDH.mutation.WES='IDH1 R132H',evidence='present in TumorOnly?')) %>%
  rbind(data.frame('Sample_Name'='160-R2',IDH.mutation.WES='IDH1 R132H',evidence='present in TumorOnly?')) %>%
  rbind(data.frame('Sample_Name'='163-R1',IDH.mutation.WES='IDH1 R132H',evidence='present in TumorOnly?')) %>%
  rbind(data.frame('Sample_Name'='163-R2',IDH.mutation.WES='IDH1 R132H',evidence='present in TumorOnly?')) %>%
  rbind(data.frame('Sample_Name'='163-R4',IDH.mutation.WES='IDH1 R132H',evidence='present in TumorOnly?'))

tmp.manual <- data.frame()


tmp <- rbind(tmp, tmp.manual) %>%
  dplyr::mutate(pid = gsub("\\-.+","",Sample_Name))


# check if there are no patients with inconsistencies throughout time
stopifnot(duplicated(tmp %>% dplyr::mutate(evidence = NULL, Sample_Name = NULL) %>% dplyr::distinct() %>% dplyr::pull(pid)) == F)


# match to pid only
tmp <- tmp %>%
  dplyr::filter(!is.na(IDH.mutation.WES)) %>% 
  dplyr::mutate(evidence = NULL, Sample_Name = NULL) %>%
  dplyr::distinct()



metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::mutate(pid = gsub("\\_.+","",Sample_Name)) %>% 
  dplyr::left_join(tmp, by=c('pid' = 'pid'),suffix = c("", "")) %>% 
  dplyr::mutate(pid = NULL) %>% 
  dplyr::mutate(IDH.mutation = ifelse(is.na(IDH.mutation) & !is.na(IDH.mutation.WES), IDH.mutation.WES, IDH.mutation))



# find discrepancies -- sample 048
metadata.glass.per.resection %>% 
  dplyr::filter(IDH.mutation.WES != IDH.mutation) %>% 
  dplyr::select(Sample_Name, IDH.mutation, IDH.mutation.WES)



metadata.glass.per.resection <- metadata.glass.per.resection %>% 
  dplyr::mutate(IDH.mutation = ifelse(grepl("048",Sample_Name),"IDH1 R132G?/H", IDH.mutation))



# missing? :
# 163-R4-I1_intersect_ann.vcf:chr2        209113112       COSM28746;rs121913500   C       T       .       PASS    AS_FilterStatus=SITE;AS_SB_TABLE=17,17|5,5;DP=46;ECNT=1;GERMQ=40;MBQ=20,20;MFRL=86,97;MMQ=60,60;MPOS=34;POPAF=7.3;ROQ=32;TLOD=18.93;AA=p.R132H;CDS=c.395G>A;CNT=4530;GENE=IDH1;STRAND=-;AC=0;gnomAD_AF=0.00000e+00   GT:AD:AF:DP:F1R2:F2R1:SB        0/1:34,10:0.25:44:20,4:14,6:17,17,5,5


# I. load CNVs [type-1] ----


# segments, seems most meaningful for further stats
#mixedrank = function(x) order(gtools::mixedorder(x)) # https://stackoverflow.com/questions/32378108/using-gtoolsmixedsort-or-alternatives-with-dplyrarrange


cnv <- read.table('data/glass/WES/copynumber_profiles/calls_segments.txt',header=T) %>%
  tibble::remove_rownames() %>% 
  dplyr::rename(cnv.segments.id = breakpoint) %>% 
  dplyr::rename(cnv.segments.chrom = chrom) %>% 
  dplyr::rename(cnv.segments.start = start) %>% 
  dplyr::rename(cnv.segments.end = end) %>% 
  dplyr::mutate(cnv.segments.length = cnv.segments.end - cnv.segments.start) %>% 
  dplyr::arrange(order(gtools::mixedorder(cnv.segments.chrom)), cnv.segments.start, cnv.segments.end)  # rank = order(order())


cnv.metadata <- cnv %>% 
  dplyr::select(cnv.segments.id, cnv.segments.chrom, cnv.segments.start, cnv.segments.end, cnv.segments.length) %>% 
  tibble::column_to_rownames('cnv.segments.id')


# @todo why is 32.7% of all segments NA - which is ?
cnv <- cnv %>% 
  dplyr::mutate(cnv.segments.chrom = NULL ) %>% # first remove all metadata to get a value'ed matrix
  dplyr::mutate(cnv.segments.start = NULL ) %>% 
  dplyr::mutate(cnv.segments.end = NULL ) %>% 
  dplyr::mutate(cnv.segments.length = NULL )  %>% 
  tibble::column_to_rownames('cnv.segments.id') %>% 
  dplyr::mutate(nNA = rowSums(is.na(.))) %>% # count NA's
  dplyr::filter(nNA == 0) %>% # remove NA's
  dplyr::mutate(nNA = NULL) %>%  # remove dummy
  t() %>% # transpose to fix sample names
  as.data.frame() %>% 
  tibble::rownames_to_column('sid') %>% 
  dplyr::mutate(sid = gsub("^X","",sid)) %>% 
  dplyr::mutate(sid = gsub("_I[0-9]+$","",sid)) %>% 
  tibble::column_to_rownames('sid') %>% 
  t() %>%  # transpose back
  as.data.frame


# take same segments subset (no NA's) for the metadata
cnv.metadata <- dplyr::left_join(
    cnv %>%
      tibble::rownames_to_column('cnv.segments.id') %>% 
      dplyr::select(cnv.segments.id), 
    cnv.metadata  %>%
      tibble::rownames_to_column('cnv.segments.id') ,
    by=c('cnv.segments.id'='cnv.segments.id')) %>% 
  tibble::column_to_rownames('cnv.segments.id')
    

# x-check (!)
stopifnot(rownames(cnv) == rownames(cnv.metadata))



# II. load CNVs [type-2] ----

tmp <- read.table('data/glass/WES/copynumber_profiles/100kbp-called_VAFPurity.igv',header=T, sep = "\t") %>% 
  dplyr::mutate(chromosome = paste0('chr', chromosome)) %>% 
  dplyr::mutate(feature = paste0('chr', feature)) %>% 
  `colnames<-`(gsub("^X","",colnames(.))) %>% 
  `colnames<-`(gsub("_I[0-9]$","",colnames(.))) %>% 
  tibble::remove_rownames()

cnv2 <- tmp %>% 
  tibble::column_to_rownames('feature') %>% 
  dplyr::select(-c('chromosome', 'start','end'))

cnv2.metadata <- tmp %>%
  dplyr::select(feature, chromosome, start, end) %>% 
  dplyr::rename(segment.id = feature)


stopifnot(rownames(cnv2) == cnv2.metadata$segment.id)


rm(tmp)





