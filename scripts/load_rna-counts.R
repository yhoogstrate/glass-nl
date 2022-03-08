#!/usr/bin/env R


if(!exists("metadata.glass.per.resection")) {
  warning('metadata was not loaded')
  
  source('scripts/load_metadata.R')
}


expression.glass.gtf <- read.delim('data/gencode.v34.primary_assembly.annotation.gtf',comment.char = "#",sep="\t",header=F) %>% 
  dplyr::filter(V3 == "gene") %>% 
  dplyr::mutate(gene_id = gsub("^.*gene_id[ ]+([^;]+);.+$","\\1", V9)) %>% 
  dplyr::filter(grepl("_PAR_", gene_id) == F) %>%  # these are odd equivalents of chrX positioned at chrY
  dplyr::mutate(ENSID = gsub("\\..+$","",gene_id)) %>% 
  dplyr::mutate(gene_name = gsub("^.*gene_name[ ]+([^;]+);.+$","\\1", V9)) %>% 
  dplyr::mutate(gene_type = gsub("^.*gene_type[ ]+([^;]+);.+$","\\1", V9)) %>% 
  dplyr::mutate(gene_uid = paste0(ENSID , "_", gene_name))




expression.glass <- read.delim('data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T) %>% 
  `colnames<-`(gsub("^X.+new.","",colnames(.))) %>% 
  `colnames<-`(gsub(".Aligned.+bam$","",colnames(.))) %>% 
  `colnames<-`(gsub(".","-",colnames(.),fixed=T)) %>% 
  dplyr::filter(grepl("_PAR_", Geneid) == F) %>%  # these are odd equivalents of chrX positioned at chrY
  dplyr::rename(gene_id = Geneid) %>% 
  dplyr::left_join(expression.glass.gtf %>% dplyr::select(gene_id, gene_uid),by=c('gene_id'='gene_id'))




expression.glass.metadata <- expression.glass %>% 
  dplyr::select(gene_id, Chr, Start, End, Strand, Length) %>% 
  dplyr::left_join(expression.glass.gtf, by=c('gene_id' = 'gene_id')) %>% 
  dplyr::mutate(gene_loc = paste0("chr:", round((V4 + V5) /  2 / 1000000),"M")) %>% 
  dplyr::rename(gene_strand = V7)
rm(expression.glass.gtf)



stopifnot(metadata.glass.per.resection$genomescan.sid %in% colnames(expression.glass)) # all metadata included samples must exist expression data



expression.glass <- expression.glass %>%
  dplyr::select(-c('gene_id', 'Chr', 'Start', 'End', 'Strand', 'Length')) %>% 
  tibble::column_to_rownames('gene_uid') %>% 
  dplyr::select( # Reorder w/ metadata
    metadata.glass.per.resection %>%
      dplyr::filter(excluded == F) %>% 
      dplyr::pull('genomescan.sid'))



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





## VST transform ----


expression.glass.vst <- expression.glass %>% 
  DESeq2::DESeqDataSetFromMatrix( data.frame(cond = as.factor(paste0('c',round(runif(ncol(.)))+1) )), ~cond) %>% 
  DESeq2::vst(blind=F) %>% 
  SummarizedExperiment::assay() %>% 
  as.data.frame(stringsAsFactors=F)







