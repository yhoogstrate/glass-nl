#!/usr/bin/env R

source('scripts/R/chrom_sizes.R')



# load cellularities ----

source('scripts/load_cellularities.R')



# load methy identifiers


meth.ids <- read.csv('data/glass/Methylation/Datasheet3.csv') %>% 
  dplyr::select(Sample_Name, Sample_ID)


# load matching meth data

meth.data <- data.frame(Heidelberg.segment.file = Sys.glob("data/glass/Methylation/Heidelberg/Heidelberg_unzip/*/cnvp_v3.0/*.segments.seg")) %>% 
  dplyr::mutate(Sample_ID = gsub("^.+_unzip/([^/]+)_Run.+$","\\1",Heidelberg.segment.file)) %>% 
  dplyr::left_join(meth.ids, by=c('Sample_ID'='Sample_ID')) %>% 
  dplyr::left_join(cnv.cellularities, by=c('Sample_Name'='Sample_ID'))  %>% 
  dplyr::filter(!is.na(cellularity_ManuallyCurated))



#  for each, plot

for(id in unique(meth.data$Sample_ID)) {

  
  #id =  unique(meth.data$Sample_ID)[13]
  
  sel <- meth.data %>% dplyr::filter(`Sample_ID` == id)
  
  plt <- read.delim(sel$Heidelberg.segment.file) %>% 
    dplyr::mutate(id = paste0(chrom, loc.start, "-",loc.end)) %>% 
    dplyr::mutate(chrom.offset = chrs_hg19_s[chrom]) %>% 
    dplyr::mutate(loc.start = loc.start + chrom.offset) %>% 
    dplyr::mutate(loc.end = loc.end + chrom.offset) %>% 
    tidyr::pivot_longer(cols = c(loc.start, loc.end))
  
  p <- sel$cellularity_ManuallyCurated 
  #p <- 0.75
  f1 <- log2( ((1-p) * 2 + p * 1 ) / 2 )  # log2 as checked in conumee: https://github.com/hovestadt/conumee/blob/master/R/process.R#L12 & https://github.com/hovestadt/conumee/blob/master/R/process.R#L89
  f3 <- log2( ((1-p) * 2 + p * 3 ) / 2 ) 
  f4 <- log2( ((1-p) * 2 + p * 4 ) / 2 ) 
  
  ggplot(plt, aes(x = value, y = seg.median, group=id, col=chrom)) +
    geom_vline(xintercept = c(0,chrs_hg19_e),col="gray90") +
    geom_hline(yintercept = f1, col="red") +
    geom_hline(yintercept = f3, col="red") +
    geom_hline(yintercept = f4, col="red",lty=2) +
    geom_line(lwd=3) +
    ylim(-1,1) +
    theme(
      text = element_text(family = 'Helvetica'),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = 'bottom',
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      panel.background = element_rect(fill = 'white', colour = 'white'),
      axis.title = element_text(face = "bold",size = rel(1)),
      axis.title.y = element_text(angle=90,vjust =2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.line = element_line(colour="black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) + labs(x=NULL, title=paste0(sel$Sample_Name," ",id, "    p:",sel$cellularity_ManuallyCurated),
             subtitle="exp seg median: log2(((1-p) * 2 + p * c ) / 2),       where p = purity   &   c = copy number, 1, 3 or 4") 
  
  
  
  
  ggsave(paste0("output/figures/cellularities/",sel$Sample_Name,"__",id,".pdf"),width=10,height=6)
  
}




