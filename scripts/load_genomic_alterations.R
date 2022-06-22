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


# 0. file overzicht ----

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


# 1. normals blood ----

# check en integrate met CDC27 - zijn allemaal 'NormalOnly'
# outersect files!?




# 2. MatchedNormal ----


idh1f = function(x) grepl("IDH1", x, fixed=TRUE)

filt.idh = FilterRules(list(idh1 = idh1f))


tmp <- data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/Annotated/*_fun.vcf.gz")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(MatchedNormal.Annotated.fun = filename)
#t1 <- rownames(info(header(VariantAnnotation::readVcf(l[1,1], "hg19")))) 
t1 <- VariantAnnotation::readVcf(tmp[1,1], "hg19")
t2 <- VariantAnnotation::readVcf(tmp[2,1], "hg19")
t3 <- VariantAnnotation::readVcf(tmp[3,1], "hg19")
t4 <- VariantAnnotation::readVcf(tmp[4,1], "hg19")

# Geen IDH[1/2],ATRX,TP53
# grep -v '^#' data/glass/WES/MatchedNormal/variant/Annotated/*_PON.vcf | grep -i -E 'IDH|ATRX|TP53'

# wel IDH1/2
# grep -v '^#' data/glass/WES/MatchedNormal/variant/Annotated/*_fun.vcf | grep -i -E 'IDH|ATRX|TP53'




parse.vcf <- function(fn) {
  tmp <- VariantAnnotation::readVcf(fn, "hg19")
  
  tmp.info <- tmp %>% 
    info() %>% 
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
    ))
  
  tmp.geno <- tmp %>% 
    geno()
  
  tmp.geno.ad1 <- tmp.geno$AD[,1] %>%  purrr::map_chr(c(1))
  tmp.geno.ad2 <- tmp.geno$AD[,1] %>%  purrr::map_chr(c(2))
  stopifnot(names(tmp.geno.ad1) == names(tmp.geno.ad2))
  
  tmp.geno.df <- data.frame(FILTER = tmp@fixed$FILTER,
                            DP.ref = as.numeric(tmp.geno.ad1),
                            DP.var = as.numeric(tmp.geno.ad2),
                            AF.geno = as.numeric(tmp.geno$AF[,1])) %>% 
    dplyr::mutate(DP.geno = DP.var + DP.ref) %>% 
    dplyr::mutate(VAF = DP.var / DP.geno)


  tmp.info <- cbind(tmp.info, tmp.geno.df) %>% 
    dplyr::mutate(POPAF = as.numeric(POPAF))

  return(tmp.info)
}


mutation.data.mn.header <- data.frame(filename = Sys.glob("data/glass/WES/MatchedNormal/*/Annotated/*_fun.vcf.gz")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(MatchedNormal.Annotated.fun = filename)

mutation.data.mn <- mutation.data.header %>% 
  dplyr::pull(`MatchedNormal.Annotated.fun`) %>% 
  pbapply::pblapply(parse.vcf) %>% 
  dplyr::bind_rows() %>% 
  dplyr::left_join(mutation.data.header, by=c('mutations_filename'='MatchedNormal.Annotated.fun'))



mutation.data.mn %>% 
  dplyr::mutate(pid = gsub("^([0-9]+).+$","\\1",sample_name)) %>% 
  dplyr::select(sample_name, Gencode_19_hugoSymbol) %>% 
  dplyr::distinct(sample_name, Gencode_19_hugoSymbol) %>% 
  count(`Gencode_19_hugoSymbol`) %>%
  dplyr::arrange(-n) %>%
  head(n=25)




mutation.data.mn %>%
  dplyr::filter(Gencode_19_hugoSymbol == "IDH1") %>% 
  dplyr::select(name, FILTER, PON_COUNT, dbSNP_ID, sample_name,VAF, Gencode_19_chromosome, Gencode_19_start, Gencode_19_end) 



mutation.data.mn %>%
  dplyr::filter(Gencode_19_hugoSymbol == "ZNF131") %>%
  dplyr::select(name, FILTER, dbSNP_ID, sample_name,VAF, Gencode_19_chromosome, Gencode_19_start, Gencode_19_end)



# 3. TumorOnly ----


parse.vcf.2  <- function(fn) {
  #fn <- 'data/glass/WES/TumorOnly/intersect/002-R1-I3_intersect_fun.vcf.gz'
  tmp <- VariantAnnotation::readVcf(fn, "hg19")
  
  tmp.info <- tmp %>% 
    info() %>% 
    as.data.frame %>% 
    dplyr::mutate(`FUNCOTATION` = unlist(`FUNCOTATION`)) %>% 
    dplyr::mutate(FUNCOTATION = gsub('^\\[','',FUNCOTATION)) %>% 
    dplyr::mutate(FUNCOTATION = gsub('\\]$','',FUNCOTATION)) %>% 
    dplyr::mutate(FUNCOTATION = iconv(FUNCOTATION, "latin1", "ASCII", sub="")) %>%  # weird non ascii character in linked disease name(s)
    tidyr::separate(FUNCOTATION, sep="\\|", into = c("Gencode_19_hugoSymbol","Gencode_19_ncbiBuild","Gencode_19_chromosome","Gencode_19_start","Gencode_19_end","Gencode_19_variantClassification","Gencode_19_secondaryVariantClassification","Gencode_19_variantType","Gencode_19_refAllele","Gencode_19_tumorSeqAllele1","Gencode_19_tumorSeqAllele2","Gencode_19_genomeChange","Gencode_19_annotationTranscript","Gencode_19_transcriptStrand","Gencode_19_transcriptExon","Gencode_19_transcriptPos","Gencode_19_cDnaChange","Gencode_19_codonChange","Gencode_19_proteinChange","Gencode_19_gcContent","Gencode_19_referenceContext","Gencode_19_otherTranscripts","CGC_Name","CGC_GeneID","CGC_Chr","CGC_Chr_Band","CGC_Cancer_Somatic_Mut","CGC_Cancer_Germline_Mut","CGC_Tumour_Types__(Somatic_Mutations)","CGC_Tumour_Types_(Germline_Mutations)","CGC_Cancer_Syndrome","CGC_Tissue_Type","CGC_Cancer_Molecular_Genetics","CGC_Mutation_Type","CGC_Translocation_Partner","CGC_Other_Germline_Mut","CGC_Other_Syndrome/Disease","ClinVar_VCF_AF_ESP","ClinVar_VCF_AF_EXAC","ClinVar_VCF_AF_TGP","ClinVar_VCF_ALLELEID","ClinVar_VCF_CLNDISDB","ClinVar_VCF_CLNDISDBINCL","ClinVar_VCF_CLNDN","ClinVar_VCF_CLNDNINCL","ClinVar_VCF_CLNHGVS","ClinVar_VCF_CLNREVSTAT","ClinVar_VCF_CLNSIG","ClinVar_VCF_CLNSIGCONF","ClinVar_VCF_CLNSIGINCL","ClinVar_VCF_CLNVC","ClinVar_VCF_CLNVCSO","ClinVar_VCF_CLNVI","ClinVar_VCF_DBVARID","ClinVar_VCF_GENEINFO","ClinVar_VCF_MC","ClinVar_VCF_ORIGIN","ClinVar_VCF_RS","ClinVar_VCF_SSR","ClinVar_VCF_ID","ClinVar_VCF_FILTER","Cosmic_overlapping_mutations","HGNC_HGNC_ID","HGNC_Approved_Name","HGNC_Status","HGNC_Locus_Type","HGNC_Locus_Group","HGNC_Previous_Symbols","HGNC_Previous_Name","HGNC_Synonyms","HGNC_Name_Synonyms","HGNC_Chromosome","HGNC_Date_Modified","HGNC_Date_Symbol_Changed","HGNC_Date_Name_Changed","HGNC_Accession_Numbers","HGNC_Enzyme_IDs","HGNC_Entrez_Gene_ID","HGNC_Ensembl_Gene_ID","HGNC_Pubmed_IDs","HGNC_RefSeq_IDs","HGNC_Gene_Family_ID","HGNC_Gene_Family_Name","HGNC_CCDS_IDs","HGNC_Vega_ID","HGNC_Entrez_Gene_ID(supplied_by_NCBI)","HGNC_OMIM_ID(supplied_by_OMIM)","HGNC_RefSeq(supplied_by_NCBI)","HGNC_UniProt_ID(supplied_by_UniProt)","HGNC_Ensembl_ID(supplied_by_Ensembl)","HGNC_UCSC_ID(supplied_by_UCSC)","dbSNP_ASP","dbSNP_ASS","dbSNP_CAF","dbSNP_CDA","dbSNP_CFL","dbSNP_COMMON","dbSNP_DSS","dbSNP_G5","dbSNP_G5A","dbSNP_GENEINFO","dbSNP_GNO","dbSNP_HD","dbSNP_INT","dbSNP_KGPhase1","dbSNP_KGPhase3","dbSNP_LSD","dbSNP_MTP","dbSNP_MUT","dbSNP_NOC","dbSNP_NOV","dbSNP_NSF","dbSNP_NSM","dbSNP_NSN","dbSNP_OM","dbSNP_OTH","dbSNP_PM","dbSNP_PMC","dbSNP_R3","dbSNP_R5","dbSNP_REF","dbSNP_RS","dbSNP_RSPOS","dbSNP_RV","dbSNP_S3D","dbSNP_SAO","dbSNP_SLO","dbSNP_SSR","dbSNP_SYN","dbSNP_TOPMED","dbSNP_TPA","dbSNP_U3","dbSNP_U5","dbSNP_VC","dbSNP_VLD","dbSNP_VP","dbSNP_WGT","dbSNP_WTD","dbSNP_dbSNPBuildID","dbSNP_ID","dbSNP_FILTER")) %>% 
    dplyr::mutate(FUNCOTATION = NULL) %>% 
    dplyr::mutate(mutations_filename = fn) %>% 
    tibble::rownames_to_column('name') %>% 
    dplyr::select(c('name',
                    'AA','AC','AS_SB_TABLE','CDS','CNT','DB','DP','ECNT','GENE','Gencode_19_hugoSymbol','Gencode_19_chromosome',
                    'Gencode_19_start','Gencode_19_end','Gencode_19_variantClassification','Gencode_19_secondaryVariantClassification',
                    'Gencode_19_proteinChange','CGC_Name','CGC_Cancer_Somatic_Mut','CGC_Cancer_Germline_Mut',
                    'CGC_Tumour_Types__(Somatic_Mutations)','CGC_Mutation_Type','ClinVar_VCF_CLNSIG','ClinVar_VCF_MC','dbSNP_ID','MBQ',
                    'MFRL','MMQ','MPOS','NALOD','NLOD','PON','PON_COUNT','POPAF','TLOD','gnomAD_AF','mutations_filename'
    ))
  
  tmp.geno <- tmp %>% 
    geno()
  
  tmp.geno.ad1 <- tmp.geno$AD[,1] %>%  purrr::map_chr(c(1))
  tmp.geno.ad2 <- tmp.geno$AD[,1] %>%  purrr::map_chr(c(2))
  stopifnot(names(tmp.geno.ad1) == names(tmp.geno.ad2))
  
  tmp.geno.df <- data.frame(FILTER = tmp@fixed$FILTER,
                            DP.ref = as.numeric(tmp.geno.ad1),
                            DP.var = as.numeric(tmp.geno.ad2),
                            AF.geno = as.numeric(tmp.geno$AF[,1])) %>% 
    dplyr::mutate(DP.geno = DP.var + DP.ref) %>% 
    dplyr::mutate(VAF = DP.var / DP.geno)
  
  
  tmp.info <- cbind(tmp.info, tmp.geno.df)  %>% 
    dplyr::mutate(POPAF = as.numeric(POPAF))
  
  return(tmp.info)
}


mutation.data.to.header <- data.frame(filename = Sys.glob("data/glass/WES/TumorOnly/intersect/*_fun.vcf.gz")) %>% 
  dplyr::mutate(sample_name = gsub("^.+/([0-9]+[-][^-]+).+$","\\1",gsub("_","-",filename))) %>% 
  dplyr::rename(TumorOnly.intersect.fun = filename)


mutation.data.to <- mutation.data.to.header %>% 
  dplyr::pull(`TumorOnly.intersect.fun`) %>% 
  pbapply::pblapply(parse.vcf.2) %>%
  dplyr::bind_rows() %>% 
  dplyr::left_join(mutation.data.to.header, by=c('mutations_filename'='TumorOnly.intersect.fun'))






mutation.data.to %>%
  dplyr::filter(Gencode_19_hugoSymbol == "PWRN1") %>% 
  dplyr::select(name, FILTER, dbSNP_ID, sample_name, VAF, Gencode_19_chromosome, Gencode_19_start, Gencode_19_end) 



mutation.data.to %>%
  dplyr::filter(Gencode_19_hugoSymbol == "ZNF131") %>%
  dplyr::select(name, FILTER, dbSNP_ID, sample_name, VAF, Gencode_19_chromosome, Gencode_19_start, Gencode_19_end)


mutation.data.to %>%
  dplyr::filter(Gencode_19_hugoSymbol == "ATRX") %>%
  dplyr::select(name, FILTER, dbSNP_ID, sample_name, VAF, Gencode_19_hugoSymbol, AA, Gencode_19_chromosome, Gencode_19_start, Gencode_19_end) %>% 
  dplyr::filter(FILTER != 'germline') %>% 
  dplyr::filter(grepl('COSM',name) | !is.na(AA))


mutation.data.to %>%
  dplyr::filter(grepl("^H3",Gencode_19_hugoSymbol)) %>%
  dplyr::select(name, FILTER, dbSNP_ID, sample_name, VAF, Gencode_19_hugoSymbol, AA,
                Gencode_19_chromosome, Gencode_19_start,
                
                CGC_Cancer_Somatic_Mut,
                CGC_Cancer_Germline_Mut,
                `CGC_Tumour_Types__(Somatic_Mutations)`
                )

# In glioblastomas of young adults and pediatric patients, some studies have identified a
# small percentage of patients who are IDH wild-type and have a loss of ATRX expression (52–54).
# Ebrahimi and colleagues found these patients to have H3F3A G34 or K27 mutations, which is concordant with Ikemura and colleagues’ finding of ATRX-loss glioblastomas in younger patients being most commonly non-hemispheric in location (47, 52)




mutation.data.to %>%
  dplyr::filter(grepl("^IDH",Gencode_19_hugoSymbol)) %>%
  dplyr::select(name, FILTER, dbSNP_ID, sample_name, VAF, Gencode_19_hugoSymbol, AA,
                Gencode_19_chromosome, Gencode_19_start) %>% 
  dplyr::filter(FILTER != "germline")


# 4. merge ----


stopifnot(colnames(mutation.data.to) == colnames(mutation.data.mn))


mutation.data <- rbind(mutation.data.to %>% dplyr::mutate(type = "TumorOnly"),
                       mutation.data.mn %>% dplyr::mutate(type = "MatchingNormal")
                       ) %>% 
  dplyr::filter(FILTER != 'germline')


## a. find top hit genes ----


mutation.data %>% 
  dplyr::mutate(pid = gsub("^([0-9]+).+$","\\1",sample_name)) %>% 
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
  #dplyr::filter(!is.na(AA) | !is.na(Gencode_19_proteinChange)) %>% 
  #dplyr::filter(!is.na(AA)) %>% 
  dplyr::select(sample_name, Gencode_19_hugoSymbol) %>% 
  dplyr::distinct(sample_name, Gencode_19_hugoSymbol) %>% 
  count(`Gencode_19_hugoSymbol`) %>%
  dplyr::arrange(-n) %>%
  head(n=25)




### CDC27 ----

cdc27 <- mutation.data %>% 
  dplyr::mutate(pid = gsub("^([0-9]+).+$","\\1",sample_name)) %>% 
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
  dplyr::select(name, sample_name, Gencode_19_chromosome, Gencode_19_start, PON_COUNT, type, sample_name)



### FRG1 ----


frg1 <- mutation.data %>% 
  dplyr::mutate(pid = gsub("^([0-9]+).+$","\\1",sample_name)) %>% 
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
  dplyr::select(name, sample_name, Gencode_19_chromosome, Gencode_19_start ,type, PON_COUNT)



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





# I. load CNVs ----


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




