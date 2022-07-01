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


# expression.glass.gtf %>%
#   dplyr::filter(gene_name == "FCGBP")
# 1x

# expression.glass.gtf %>%
#   dplyr::filter(grepl("ENSG00000275395",ENSID ))
# 1x




# check for XIST
# sum(expression.glass.gtf$gene_name == "XIST")

# load count exonic features ----


expression.glass.exon <- read.delim('data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T) %>% 
  `colnames<-`(gsub("^X.+new.","",colnames(.))) %>% 
  `colnames<-`(gsub(".Aligned.+bam$","",colnames(.))) %>% 
  `colnames<-`(gsub(".","-",colnames(.),fixed=T)) %>% 
  dplyr::filter(grepl("_PAR_", Geneid) == F) %>%  # these are odd equivalents of chrX positioned at chrY
  dplyr::rename(gene_id = Geneid) %>% 
  dplyr::left_join(expression.glass.gtf %>% dplyr::select(gene_id, gene_uid),by=c('gene_id'='gene_id'))


# some sort of qc?
# expression.glass.exon %>%
#   dplyr::filter(grepl("EGFR",gene_uid)) %>% 
#   dplyr::select(`104059-002-042`)


# check for xist
# expression.glass.exon %>%
#   dplyr::filter(grepl("XIST",gene_uid))


stopifnot(metadata.glass.per.resection$genomescan.sid %in% colnames(expression.glass.exon)) # all metadata included samples must exist expression data



## merge GTF and featureCounts per-gene stats ----


expression.glass.exon.metadata <- expression.glass.exon %>% 
  dplyr::select(gene_id, Chr, Start, End, Strand, Length) %>% 
  dplyr::left_join(expression.glass.gtf, by=c('gene_id' = 'gene_id')) %>% 
  dplyr::rename(gene_chr = V1) %>%  
  dplyr::rename(gene_strand = V7) %>% 
  dplyr::mutate(gene_chr_center_loc = (V4 + V5) /  2) %>% 
  dplyr::mutate(gene_loc = paste0(gene_chr, ":", round(gene_chr_center_loc / 1000000),"M")) 


nrow(expression.glass.exon.metadata) == nrow(expression.glass.exon)


# expression.glass.exon.metadata %>% 
#   dplyr::filter(grepl("XIST",gene_name)) %>% 
#   dim




## remove non count columns ----


expression.glass.exon <- expression.glass.exon %>%
  dplyr::select(-c('gene_id', 'Chr', 'Start', 'End', 'Strand', 'Length')) %>% 
  tibble::column_to_rownames('gene_uid') %>% 
  dplyr::select( # Reorder w/ metadata
    metadata.glass.per.resection %>%
      dplyr::filter(excluded == F) %>% 
      dplyr::pull('genomescan.sid'))



## exclude non relevant gene types ----

# Select protein_coding and lncRNA genes with an average read count >= 3
sel <- expression.glass.exon %>% 
  dplyr::mutate(avg.read.count = rowSums(.) / ncol(.)) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::left_join(expression.glass.exon.metadata %>% dplyr::select('gene_uid','gene_type'), by=c('gene_uid' = 'gene_uid')) %>% 
  dplyr::mutate(keep = (avg.read.count >= 3) & (gene_type %in% c('lncRNA','protein_coding'))) %>% 
  dplyr::pull('keep')


nrow(expression.glass.exon.metadata) == nrow(expression.glass.exon)
stopifnot(rownames(expression.glass.exon) == expression.glass.exon.metadata$gene_uid)
  
  

expression.glass.exon <- expression.glass.exon %>%
  dplyr::filter(sel)

expression.glass.exon.metadata <- expression.glass.exon.metadata %>%
  dplyr::filter(sel)

rm(sel)


stopifnot(rownames(expression.glass.exon) == expression.glass.exon.metadata$gene_uid)


# expression.glass.exon %>% 
#   tibble::rownames_to_column('gid') %>% 
#   dplyr::filter(grepl("XIST",gid)) %>% 
#   dim


## add hg19 liftover ----


tmp <- read.table('data/2022-05-27_glass-nl_dge_results_paired_hglft_hg19_genome_120f5_f2660.bed') %>% 
  dplyr::mutate(V5 = NULL) %>% 
  dplyr::rename(chr.hg19 = V1) %>% 
  dplyr::rename(start.hg19 = V2) %>% 
  dplyr::rename(end.hg19 = V3) %>% 
  dplyr::mutate(length = end.hg19 - start.hg19) %>% 
  dplyr::group_by(V4) %>% 
  dplyr::filter(length == max(length)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(length = NULL)


stopifnot(duplicated(tmp$V4) == F)


expression.glass.exon.metadata <- expression.glass.exon.metadata %>% 
  dplyr::left_join(tmp,by=c('gene_id'='V4'),suffix = c("", ""))



## VST transform ----


expression.glass.exon.vst <- expression.glass.exon %>% 
  DESeq2::DESeqDataSetFromMatrix( data.frame(cond = as.factor(paste0('c',round(runif(ncol(.)))+1) )), ~cond) %>% 
  DESeq2::vst(blind=T) %>% 
  SummarizedExperiment::assay() %>% 
  as.data.frame(stringsAsFactors=F)


# load gene body features ----


expression.glass.gene <- read.delim('output/tables/rna-seq/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.per-gene.txt',skip=1,header=T) %>% 
  `colnames<-`(gsub("data.glass.RNAseq.alignments.alignments.new.","",colnames(.))) %>% 
  `colnames<-`(gsub(".Aligned.sortedByCoord.out.markduplicate.bam$","",colnames(.))) %>% 
  `colnames<-`(gsub(".","-",colnames(.),fixed=T)) %>% 
  dplyr::filter(grepl("_PAR_", Geneid) == F) %>%  # these are odd equivalents of chrX positioned at chrY
  dplyr::rename(gene_id = Geneid) %>% 
  dplyr::left_join(expression.glass.gtf %>% dplyr::select(gene_id, gene_uid),by=c('gene_id'='gene_id'))



stopifnot(metadata.glass.per.resection$genomescan.sid %in% colnames(expression.glass.gene)) # all metadata included samples must exist expression data



## merge GTF and featureCounts per-gene stats ----


expression.glass.gene.metadata <- expression.glass.gene %>% 
  dplyr::select(gene_id, Chr, Start, End, Strand, Length) %>% 
  dplyr::left_join(expression.glass.gtf, by=c('gene_id' = 'gene_id')) %>% 
  dplyr::rename(gene_chr = V1) %>%  
  dplyr::rename(gene_strand = V7) %>% 
  dplyr::mutate(gene_chr_center_loc = (V4 + V5) /  2) %>% 
  dplyr::mutate(gene_loc = paste0(gene_chr, ":", round(gene_chr_center_loc / 1000000),"M")) 


## remove non count columns ----


expression.glass.gene <- expression.glass.gene %>%
  dplyr::select(-c('gene_id', 'Chr', 'Start', 'End', 'Strand', 'Length')) %>% 
  tibble::column_to_rownames('gene_uid') %>% 
  dplyr::select( # Reorder w/ metadata
    metadata.glass.per.resection %>%
      dplyr::filter(excluded == F) %>% 
      dplyr::pull('genomescan.sid'))



## exclude non relevant gene types ----
# Select protein_coding and lncRNA genes with an average read count >= 3


sel <- expression.glass.gene %>% 
  dplyr::mutate(avg.read.count = rowSums(.) / ncol(.)) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::left_join(expression.glass.gene.metadata %>% dplyr::select('gene_uid','gene_type'), by=c('gene_uid' = 'gene_uid')) %>% 
  dplyr::mutate(keep = (avg.read.count >= 3) & (gene_type %in% c('lncRNA','protein_coding'))) %>% 
  dplyr::pull('keep')



stopifnot(rownames(expression.glass.gene) == expression.glass.gene.metadata$gene_uid)



expression.glass.gene <- expression.glass.gene %>%
  dplyr::filter(sel)

expression.glass.gene.metadata <- expression.glass.gene.metadata %>%
  dplyr::filter(sel)

rm(sel)


stopifnot(rownames(expression.glass.gene) == expression.glass.gene.metadata$gene_uid)





## VST transform ----


expression.glass.gene.vst <- expression.glass.gene %>% 
  DESeq2::DESeqDataSetFromMatrix( data.frame(cond = as.factor(paste0('c',round(runif(ncol(.)))+1) )), ~cond) %>% 
  DESeq2::vst(blind=T) %>% 
  SummarizedExperiment::assay() %>% 
  as.data.frame(stringsAsFactors=F)



# cleanup ----


rm(expression.glass.gtf)


# power test X/Y plot?




