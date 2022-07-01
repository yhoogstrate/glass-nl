#!/usr/bin/env R


library(tidyverse)


library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)


# load RNA ----


# if metadata was not loaded, load it
if(!exists("metadata.glass.per.resection")) {
  warning('metadata was not loaded')
  
  source('scripts/load_metadata.R')
}



if("padj.partially.paired.exon" %in% colnames(expression.glass.exon.metadata) == F) {
  warning('paired exon-count DGE results were not loaded')
  
  source('scripts/load_analysis_DGE.R')
}



# load meth ----



# plt gviz style ----

lim <- c(20 * 1000000, 40 * 1000000)


ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chr6")
grtrack <- BiomartGeneRegionTrack(genome="hg19", 
                                  #start=20000000,
                                  #end=4000000,
                                  dataset="hsapiens_gene_ensembl",
                                  #chromosome = "chr6", 
                                  name="ENSEMBL genes", stacking="full", showId=T,
                                  symbol="EGFR",
                                  #transcript=c("ENST00000275493","ENST00000459688"),
                                  transcriptAnnotation = "transcript",
                                  rotate.title=T)


plt.dmr <- DMRs_IvR_FDR %>% 
  dplyr::filter(seqnames == "chr6" & end >= lim[1] & start <= lim[]) 

plt.dmr.up <- plt.dmr %>% 
  dplyr::filter(meandiff > 0) 

plt.dmr.down <- plt.dmr %>% 
  dplyr::filter(meandiff < 0)



dmr = DataTrack(GRanges("chr6", IRanges(plt.dmr$start, plt.dmr$end), DMR=plt.dmr$meandiff),
                name = "DMR",type="histogram",groups="TRUE")
dmr.up = DataTrack(GRanges("chr6", IRanges(plt.dmr.up$start, plt.dmr.up$end), DMR=plt.dmr.up$meandiff),
                   name = "DMR [up]",type="histogram")
dmr.down = DataTrack(GRanges("chr6", IRanges(plt.dmr.down$start, plt.dmr.down$end), DMR=plt.dmr.down$meandiff),
                   name = "DMR [down]",type="histogram")


plt.rna <- expression.glass.exon.metadata %>% 
  dplyr::filter(chr.hg19 == "chr6") %>% 
  dplyr::filter(end.hg19 > lim[1]) %>% 
  dplyr::filter(start.hg19 < lim[2]) %>% 
  dplyr::arrange(start.hg19, end.hg19)

rna = DataTrack(GRanges("chr6", IRanges(plt.rna$start.hg19,
                                        plt.rna$end.hg19),
                        RNA=plt.rna$stat.partially.paired.exon  ),
                name = "RNA",type="histogram"
)
plotTracks(list(rna), chromosome="chr6",
           group = c(rep("green",200),rep("red",149)))




plotTracks(list(dmr, rna), chromosome="chr6")#,"boxplot", "a", "g"))


# dunno howto sep colors


aTrack.stacked <- DataTrack(start = c(50, 180, 260, 800, 600, 1240),
                                  width = c(15, 20, 40, 100, 500, 20),
                                  chromosome = "chrX",
                                  strand = "*",
                                  groups = rep(c("Huey", "Dewey", "Louie"), 
                                              c(1, 3, 2)),
                                  genome = "hg19", name = "foo")
plotTracks(aTrack.stacked, groupAnnotation="group")



# plt ggstyle ----

lim <- c(20 * 1000000, 40 * 1000000)

plt.dmr <- DMRs_IvR_FDR %>% 
  dplyr::filter(seqnames == "chr6" & end >= lim[1] & start <= lim[2]) %>% 
  dplyr::rename(intensity = logp) %>% 
  dplyr::rename(chr = seqnames) %>% 
  dplyr::mutate(intensity.type = "<sign> * -log10(Fisher padj)") %>% 
  dplyr::mutate(name="") %>% 
  dplyr::select(name, chr, start, end, intensity, intensity.type)


plt.rna <- expression.glass.exon.metadata %>% 
  dplyr::filter(chr.hg19 == "chr6") %>% 
  dplyr::filter(end.hg19 > lim[1]) %>% 
  dplyr::filter(start.hg19 < lim[2]) %>% 
  dplyr::arrange(start.hg19, end.hg19) %>% 
  dplyr::rename(name  = gene_name) %>% 
  dplyr::rename(start = start.hg19) %>% 
  dplyr::rename(end = end.hg19) %>% 
  dplyr::rename(chr = chr.hg19) %>%
  dplyr::rename(intensity = stat.partially.paired.exon) %>% 
  dplyr::mutate(intensity.type = "DESeq2 Wald statistic") %>% 
  dplyr::select(name, chr, start, end, intensity, intensity.type)


plt <- rbind(plt.dmr, plt.rna) %>% 
  dplyr::mutate(col = grepl("^H[1-5]",name)) %>% 
  dplyr::mutate(center = (end + start) / 2)


ggplot(plt, aes(x = center, y = intensity, fill=col)) +
  facet_grid(rows = vars(intensity.type), scales="free") +
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = intensity)) + 
  geom_point(pch=21) +
  theme_bw()

# col = DAXX (Death Domain Associated Protein) is a Protein Coding gene. Diseases associated with DAXX include Gastric Neuroendocrine Neoplasm and Alpha-Thalassemia. Among its related pathways are Regulation of TP53 Activity and TGF-Beta Pathway. Gene Ontology (GO) annotations related to this gene include protein homodimerization activity and enzyme binding.


