#!/usr/bin/env R

# load metadata ----

if(!exists("metadata.glass.per.resection")) {
  warning('metadata was not loaded')
  
  source('scripts/load_metadata.R')
}


# load GTF file (gene annot) ----


expression.glass.gtf <- read.delim('data/glass/RNAseq/gencode.v34.primary_assembly.annotation.gtf',comment.char = "#",sep="\t",header=F) %>% 
  dplyr::filter(V3 == "gene") %>% 
  dplyr::mutate(gene_id = gsub("^.*gene_id[ ]+([^;]+);.+$","\\1", V9)) %>% 
  dplyr::filter(grepl("_PAR_", gene_id) == F) %>%  # these are odd equivalents of chrX positioned at chrY
  dplyr::mutate(ENSID = gsub("\\..+$","",gene_id)) %>% 
  dplyr::mutate(gene_name = gsub("^.*gene_name[ ]+([^;]+);.+$","\\1", V9)) %>% 
  dplyr::mutate(gene_type = gsub("^.*gene_type[ ]+([^;]+);.+$","\\1", V9)) %>% 
  dplyr::mutate(gene_uid = paste0(ENSID , "_", gene_name))


# load read counts ----


expression.glass <- read.delim('data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T) %>% 
  `colnames<-`(gsub("^X.+new.","",colnames(.))) %>% 
  `colnames<-`(gsub(".Aligned.+bam$","",colnames(.))) %>% 
  `colnames<-`(gsub(".","-",colnames(.),fixed=T)) %>% 
  dplyr::filter(grepl("_PAR_", Geneid) == F) %>%  # these are odd equivalents of chrX positioned at chrY
  dplyr::rename(gene_id = Geneid) %>% 
  dplyr::left_join(expression.glass.gtf %>% dplyr::select(gene_id, gene_uid),by=c('gene_id'='gene_id'))


expression.glass %>%
  #tibble::rownames_to_column('gene_name') %>% 
  dplyr::filter(grepl("EGFR",gene_uid)) %>% 
  dplyr::select(`104059-002-042`)



## merge GTF and featureCounts per-gene stats ----


expression.glass.metadata <- expression.glass %>% 
  dplyr::select(gene_id, Chr, Start, End, Strand, Length) %>% 
  dplyr::left_join(expression.glass.gtf, by=c('gene_id' = 'gene_id')) %>% 
  dplyr::rename(gene_chr = V1) %>%  
  dplyr::rename(gene_strand = V7) %>% 
  dplyr::mutate(gene_chr_center_loc = (V4 + V5) /  2) %>% 
  dplyr::mutate(gene_loc = paste0(gene_chr, ":", round(gene_chr_center_loc / 1000000),"M")) 

rm(expression.glass.gtf)



stopifnot(metadata.glass.per.resection$genomescan.sid %in% colnames(expression.glass)) # all metadata included samples must exist expression data


# cleanup counts ----
## remove non count columns ----


expression.glass <- expression.glass %>%
  dplyr::select(-c('gene_id', 'Chr', 'Start', 'End', 'Strand', 'Length')) %>% 
  tibble::column_to_rownames('gene_uid') %>% 
  dplyr::select( # Reorder w/ metadata
    metadata.glass.per.resection %>%
      dplyr::filter(excluded == F) %>% 
      dplyr::pull('genomescan.sid'))


## exclude non relevant gene types ----

# Select protein_coding and lncRNA genes with an average read count >= 3
sel <- expression.glass %>% 
  dplyr::mutate(avg.read.count = rowSums(.) / ncol(.)) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::left_join(expression.glass.metadata %>% dplyr::select('gene_uid','gene_type'), by=c('gene_uid' = 'gene_uid')) %>% 
  dplyr::mutate(keep = (avg.read.count >= 3) & (gene_type %in% c('lncRNA','protein_coding'))) %>% 
  dplyr::pull('keep')



stopifnot(rownames(expression.glass) == expression.glass.metadata$gene_uid)
  
  

expression.glass <- expression.glass %>%
  dplyr::filter(sel)

expression.glass.metadata <- expression.glass.metadata %>%
  dplyr::filter(sel)

rm(sel)


stopifnot(rownames(expression.glass) == expression.glass.metadata$gene_uid)





# VST transform ----


expression.glass.vst <- expression.glass %>% 
  DESeq2::DESeqDataSetFromMatrix( data.frame(cond = as.factor(paste0('c',round(runif(ncol(.)))+1) )), ~cond) %>% 
  DESeq2::vst(blind=T) %>% 
  SummarizedExperiment::assay() %>% 
  as.data.frame(stringsAsFactors=F)




