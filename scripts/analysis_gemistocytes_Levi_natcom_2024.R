#!/usr/bin/env R


tt <- read.csv("data/glass/Metadata/Samples/Master Datasheet_ALL METHODS_27012023.csv")


metadata <- read.csv("data/Gemistocyte_score.csv") |> 
  dplyr::mutate(gemistocyte_score = factor(gemistocyte_score, levels=c("low","high"))) |> 
  dplyr::mutate(X = NULL) |> 
  assertr::verify(GLASS_ID %in% tt$Sample_Name) |> 
  dplyr::left_join(tt |> dplyr::select(Sample_Name, WES_ID), by=c('GLASS_ID'='Sample_Name'))



# mutations ----




mutation.data <- readRDS("data/glass/WES/VCF Calls/mutation.data.to.Rds") |> 
  assertr::verify(WES_ID %in% tt$WES_ID)


shared <- intersect(metadata |> dplyr::filter(!is.na(WES_ID)) |> dplyr::pull(WES_ID), unique(mutation.data$WES_ID))





plt <- mutation.data |> 
  dplyr::filter(WES_ID %in% shared) |> 
  dplyr::left_join(metadata, by=c("WES_ID"="WES_ID")) |> 
  dplyr::rename(gene = Gencode_19_hugoSymbol) |> 
  
  dplyr::filter(gene != "") |> 
  dplyr::filter(is.na(gnomAD_AF) | gnomAD_AF < 0.05) |> 
  dplyr::filter(VAF >= 0.05) |> 
  
  dplyr::select(name, gene, WES_ID, gemistocyte_score)



stats <- plt |> 
  dplyr::group_by(WES_ID, gene) |> 
  dplyr::summarise(n=dplyr::n(), gemistocyte_score = unique(gemistocyte_score)) |> 
  dplyr::ungroup() |> 
  dplyr::select(WES_ID, gene, n) |> 
  dplyr::mutate(n='y') |> 
  tidyr::pivot_wider(names_from = WES_ID, values_from = n, values_fill='n') |> 
  tibble::column_to_rownames('gene')




nlows <- c()
nhighs <- c()
chisq <- c()
fisher <- c()


samples_low <-  metadata |> dplyr::filter(gemistocyte_score == "low" &  WES_ID %in% shared) |> dplyr::pull(WES_ID)
samples_high <- metadata |> dplyr::filter(gemistocyte_score == "high" & WES_ID %in% shared) |> dplyr::pull(WES_ID)



for(i in 1:nrow(stats)) {
  t <- stats[i,] |> 
    t() |> 
    as.data.frame() |> 
    tibble::rownames_to_column('WES_ID') |> 
    dplyr::left_join(metadata |> dplyr::select(WES_ID, gemistocyte_score),by=c('WES_ID'='WES_ID')) |> 
    tibble::column_to_rownames('WES_ID')
  
  nlow <- t |> 
    tibble::rownames_to_column('WES_ID') |> 
    dplyr::filter(WES_ID %in% samples_low) |> 
    dplyr::rename(count = 2) |> 
    dplyr::mutate(count = count =="y") |> 
    dplyr::pull(count) |> 
    sum()
  
  nhigh <- t |> 
    tibble::rownames_to_column('WES_ID') |> 
    dplyr::filter(WES_ID %in% samples_high) |> 
    dplyr::rename(count = 2) |> 
    dplyr::mutate(count = count =="y") |> 
    dplyr::pull(count) |> 
    sum()
  
  
  nlows[i] <- nlow
  nhighs[i] <- nhigh
  
  
  t <- t |>     table()

  chisq[i] <-  chisq.test(t)$p.value
  fisher[i] <- fisher.test(t)$p.value
}


out <- stats |> 
  dplyr::mutate(`nlowz`=nlows, `nhighz` = nhighs, `chisq` = chisq, `fisher`=fisher)


out |> 
  dplyr::filter(fisher < 0.01 | chisq < 0.01) |>
  dplyr::select(nlowz, nhighz, chisq, fisher, samples_high) |> 
  dplyr::arrange(chisq, fisher) |> 
  View()



write.csv(out, file="output/tables/mutation_gemisto_high_all.csv")
write.csv(out |> 
            dplyr::filter(fisher < 0.01 | chisq < 0.01) |>
            dplyr::select(nlowz, nhighz, chisq, fisher, samples_high) |> 
            dplyr::arrange(chisq, fisher), file="output/tables/mutation_gemisto_high_signi_filtered.csv")



## export ----

# stats table

# associated VCF entries


# CNV ----

cnvs <- read.table('data/glass/WES/Copy Number Calls/copynumber_profiles/100kbp-called_VAFPurity.igv',header=T) |> 
  dplyr::mutate(segment = paste0("chr", feature)) |> 
  dplyr::mutate(chomosome = NULL, start = NULL, end = NULL, feature = NULL, chromosome = NULL) |> 
  tibble::column_to_rownames('segment')


cnvs2 <- cnvs |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('sid') |> 
  dplyr::mutate(sid = gsub("^X", "", sid)) |> 
  dplyr::filter(!grepl("nDNA_", sid)) |> 
  dplyr::left_join(mutation.data |> dplyr::select(Sample_Name, WES_ID) |> dplyr::distinct(),by=c('sid'='WES_ID')) |>
  dplyr::mutate(Sample_Name = ifelse(is.na(Sample_Name), gsub("^([0-9]+)_(R[0-9]+).+$","\\1-\\2", sid), Sample_Name)) |> 
  tibble::column_to_rownames('Sample_Name') |> 
  dplyr::mutate(sid = NULL) |> 
  t() |> 
  as.data.frame()


shared <- intersect(metadata$GLASS_ID, colnames(cnvs2))


cnvs2 <- cnvs2 |> 
  dplyr::select(shared)


samples_low <-  metadata |> dplyr::filter(gemistocyte_score == "low" &  GLASS_ID %in% shared) |> dplyr::pull(GLASS_ID)
samples_high <- metadata |> dplyr::filter(gemistocyte_score == "high" & GLASS_ID %in% shared) |> dplyr::pull(GLASS_ID)





n_losses <- c()
n_losses_gm_low <- c()
n_losses_gm_high <- c()
p_losses_chisq <- c()
p_losses_fisher <- c()
est_losses_fisher <- c()


n_gains <- c()
n_gains_gm_low <- c()
n_gains_gm_high <- c()
p_gains_chisq <- c()
p_gains_fisher <- c()
est_gains_fisher <- c()


for(i in 1:nrow(cnvs2)) {
  t <- cnvs2[i,] |> 
    t() |> 
    as.data.frame() |> 
    tibble::rownames_to_column('Sample_Name') |> 
    dplyr::left_join(metadata |> dplyr::select(GLASS_ID, gemistocyte_score),by=c('Sample_Name'='GLASS_ID')) |> 
    tibble::column_to_rownames('Sample_Name') |> 
    assertr::verify(!is.na(gemistocyte_score)) |> 
    dplyr::rename(cnv_profile = 1) |> 
    dplyr::mutate(losses = ifelse(cnv_profile < 0, 'loss', 'no loss')) |> 
    dplyr::mutate(gains = ifelse(cnv_profile > 0,'gain', 'no gain')) |> 
    dplyr::mutate(losses = factor(losses, levels=c('loss', 'no loss'))) |> 
    dplyr::mutate(gains = factor(gains, levels = c('gain', 'no gain')))
  
  
  n_losses[i] <- sum(t$losses == "loss")
  n_losses_gm_low[i] <- sum(t$losses == "loss" & t$gemistocyte_score == "low")
  n_losses_gm_high[i] <- sum(t$losses == "loss" & t$gemistocyte_score == "high")
  p_losses_chisq[i] <-  t |> dplyr::select(gemistocyte_score, losses) |> table() |> chisq.test()  |> purrr::pluck('p.value')
  p_losses_fisher[i] <- t |> dplyr::select(gemistocyte_score, losses) |> table() |> fisher.test() |> purrr::pluck('p.value')
  est_losses_fisher[i] <- t |> dplyr::select(gemistocyte_score, losses) |> table() |> fisher.test() |> purrr::pluck('estimate')

  n_gains[i] <-        sum(t$gains == "gain")
  n_gains_gm_low[i] <- sum(t$gains == "gain" & t$gemistocyte_score == "low")
  n_gains_gm_high[i] <- sum(t$gains == "gain" & t$gemistocyte_score == "high")
  p_gains_chisq[i] <-  t |> dplyr::select(gemistocyte_score, gains) |> table() |> chisq.test()  |> purrr::pluck('p.value')
  p_gains_fisher[i] <- t |> dplyr::select(gemistocyte_score, gains) |> table() |> fisher.test() |> purrr::pluck('p.value')
  est_gains_fisher[i] <- t |> dplyr::select(gemistocyte_score, gains) |> table() |> fisher.test() |> purrr::pluck('estimate')

}


out_cnv <- data.frame(segment = rownames(cnvs2)[1:length(n_losses)],
                  n_losses = n_losses,
                  n_losses_gm_low = n_losses_gm_low,
                  n_losses_gm_high = n_losses_gm_high,
                  p_losses_chisq = p_losses_chisq,
                  p_losses_fisher = p_losses_fisher,
                  est_losses_fisher = est_losses_fisher,
                  
                  n_gains = n_gains,
                  n_gains_gm_low = n_gains_gm_low,
                  n_gains_gm_high = n_gains_gm_high,
                  p_gains_chisq = p_gains_chisq,
                  p_gains_fisher = p_gains_fisher,
                  est_gains_fisher = est_gains_fisher
                  
                  )



#saveRDS(out_cnv, file="tmp/out_cnv.Rds")



# load 



source('scripts/load_themes.R')

library(ggplot2)
library(patchwork)



## all chrs ----



# losses 

plt <- out_cnv |> 
  dplyr::mutate(chr = gsub("^(.+):(.+)-(.+)$", "\\1", segment)) |> 
  dplyr::mutate(start = as.numeric(gsub("^(.+):(.+)-(.+)$", "\\2", segment))) |> 
  dplyr::mutate(end = as.numeric(gsub("^(.+):(.+)-(.+)$", "\\3", segment))) |> 
  dplyr::mutate(pos = (start + end)/2) |> 
  dplyr::mutate(chr = factor(chr, levels=gtools::mixedsort(unique(as.character(chr))))) |> 
  dplyr::mutate(signi = p_losses_fisher < 0.01) |> 
  dplyr::mutate(x = -log10(p_losses_fisher))


plt <- rbind(
  plt |> dplyr::mutate(x = 0)
  ,
  plt
)




p1 = ggplot(plt , aes(x = x, y=pos / 1000000, group=segment, col=signi)) +
  facet_grid(rows = vars(chr), scale="free", space="free") +
  geom_line(lwd=0.1) +
  
  theme_bw() +
  theme(
    text =          element_text(size = 7, family = "Arial", face = "plain"),
    axis.text =     element_text(size = 7, family = "Arial", face = "plain", color="black"), # , angle=90, vjust =0.5
    axis.title.x =  element_text(size = 7, family = "Arial", face = "plain", color="black"), # , vjust = -0.2
    axis.title.y =  element_text(size = 7, family = "Arial", face = "plain", color="black"),
    axis.line =     element_line(linewidth = theme_cellpress_lwd),
    axis.ticks =    element_line(linewidth = theme_cellpress_lwd),
    
    strip.text =    element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.text.x =  element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.text.y =  element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.background = element_blank(), # clean as possible
    
    legend.title =  element_text(size = 7, family = "Arial", face = "plain", color="black"),
    legend.text =   element_text(size = 7, family = "Arial", face = "plain", color="black"),
    legend.position = 'bottom',
    legend.margin   = margin(t=-2),
    legend.key.size = unit(0.2, 'lines'),
    legend.key = element_blank(), # this should remove the white squares around the legend items, but seems to fail sometimes, probably due to the key.size above?
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    
    plot.title =      element_text(size = 7, family = "Arial", face = "plain", color="black"), # `title` covers both title and subtitle
    plot.subtitle =   element_text(size = 7, family = "Arial", face = "italic", color="darkgray"),
    plot.caption =    element_text(size = 7, family = "Arial", face = "italic", color="black"),
    plot.background = element_blank(),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing =      unit(0.1, "lines"), # facet_grid margin
    #panel.border =       element_blank(), # no sqaure, but two lines instead (axis.line)
    panel.background =   element_blank()
  ) +
  
  
  theme(strip.text.y.left = element_text(angle = 0)) +
  scale_color_manual(values=c('FALSE'=alpha('gray60',0.05),'TRUE'='blue')) +
  scale_x_reverse() +
  labs(x = "-log10(p fisher exact Loss)", col="significant")




# gains

plt <- out_cnv |> 
  dplyr::mutate(chr = gsub("^(.+):(.+)-(.+)$", "\\1", segment)) |> 
  dplyr::mutate(start = as.numeric(gsub("^(.+):(.+)-(.+)$", "\\2", segment))) |> 
  dplyr::mutate(end = as.numeric(gsub("^(.+):(.+)-(.+)$", "\\3", segment))) |> 
  dplyr::mutate(pos = (start + end)/2) |> 
  dplyr::mutate(chr = factor(chr, levels=gtools::mixedsort(unique(as.character(chr))))) |> 
  dplyr::mutate(signi = p_gains_fisher < 0.01) |> 
  dplyr::mutate(x = -log10(p_gains_fisher))


plt <- rbind(
  plt |> dplyr::mutate(x = 0)
  ,
  plt
)


p2 = ggplot(plt , aes(x = x, y=pos / 1000000, group=segment, col=signi)) +
  facet_grid(rows = vars(chr), scale="free", space="free") +
  geom_line(lwd=0.1) +
  theme_bw() +
  theme(
    text =          element_text(size = 7, family = "Arial", face = "plain"),
    axis.text =     element_text(size = 7, family = "Arial", face = "plain", color="black"), # , angle=90, vjust =0.5
    axis.title.x =  element_text(size = 7, family = "Arial", face = "plain", color="black"), # , vjust = -0.2
    axis.title.y =  element_text(size = 7, family = "Arial", face = "plain", color="black"),
    axis.line =     element_line(linewidth = theme_cellpress_lwd),
    axis.ticks =    element_line(linewidth = theme_cellpress_lwd),
    
    strip.text =    element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.text.x =  element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.text.y =  element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.background = element_blank(), # clean as possible
    
    legend.title =  element_text(size = 7, family = "Arial", face = "plain", color="black"),
    legend.text =   element_text(size = 7, family = "Arial", face = "plain", color="black"),
    legend.position = 'bottom',
    legend.margin   = margin(t=-2),
    legend.key.size = unit(0.2, 'lines'),
    legend.key = element_blank(), # this should remove the white squares around the legend items, but seems to fail sometimes, probably due to the key.size above?
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    
    plot.title =      element_text(size = 7, family = "Arial", face = "plain", color="black"), # `title` covers both title and subtitle
    plot.subtitle =   element_text(size = 7, family = "Arial", face = "italic", color="darkgray"),
    plot.caption =    element_text(size = 7, family = "Arial", face = "italic", color="black"),
    plot.background = element_blank(),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing =      unit(0.1, "lines"), # facet_grid margin
    #panel.border =       element_blank(), # no sqaure, but two lines instead (axis.line)
    panel.background =   element_blank()
  ) +
  
  theme(strip.text.y.left = element_text(angle = 0)) +
  scale_color_manual(values=c('FALSE'=alpha('gray60',0.05),'TRUE'='red')) +
  labs(x = "-log10(p fisher exact Gain)", col="significant")




p1 + p2



ggsave("output/figures/gemi_cnv_all.png",width=8.5*0.975/2,height=11*0.975)



## chr6, chr12, chr19 ----



# losses 

plt <- out_cnv |> 
  dplyr::mutate(chr = gsub("^(.+):(.+)-(.+)$", "\\1", segment)) |> 
  dplyr::filter(chr %in% c("chr6", "chr12", "chr19")) |> 
  dplyr::mutate(start = as.numeric(gsub("^(.+):(.+)-(.+)$", "\\2", segment))) |> 
  dplyr::mutate(end = as.numeric(gsub("^(.+):(.+)-(.+)$", "\\3", segment))) |> 
  dplyr::mutate(pos = (start + end)/2) |> 
  dplyr::mutate(chr = factor(chr, levels=gtools::mixedsort(unique(as.character(chr))))) |> 
  dplyr::mutate(signi = p_losses_fisher < 0.01) |> 
  dplyr::mutate(x = -log10(p_losses_fisher))


plt <- rbind(
  plt |> dplyr::mutate(x = 0)
  ,
  plt
)




p1 = ggplot(plt , aes(x = x, y=pos / 1000000, group=segment, col=signi)) +
  facet_grid(rows = vars(chr), scale="free", space="free") +
  geom_line(lwd=0.1) +
  
  theme_bw() +
  theme(
    text =          element_text(size = 7, family = "Arial", face = "plain"),
    axis.text =     element_text(size = 7, family = "Arial", face = "plain", color="black"), # , angle=90, vjust =0.5
    axis.title.x =  element_text(size = 7, family = "Arial", face = "plain", color="black"), # , vjust = -0.2
    axis.title.y =  element_text(size = 7, family = "Arial", face = "plain", color="black"),
    axis.line =     element_line(linewidth = theme_cellpress_lwd),
    axis.ticks =    element_line(linewidth = theme_cellpress_lwd),
    
    strip.text =    element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.text.x =  element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.text.y =  element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.background = element_blank(), # clean as possible
    
    legend.title =  element_text(size = 7, family = "Arial", face = "plain", color="black"),
    legend.text =   element_text(size = 7, family = "Arial", face = "plain", color="black"),
    legend.position = 'bottom',
    legend.margin   = margin(t=-2),
    legend.key.size = unit(0.2, 'lines'),
    legend.key = element_blank(), # this should remove the white squares around the legend items, but seems to fail sometimes, probably due to the key.size above?
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    
    plot.title =      element_text(size = 7, family = "Arial", face = "plain", color="black"), # `title` covers both title and subtitle
    plot.subtitle =   element_text(size = 7, family = "Arial", face = "italic", color="darkgray"),
    plot.caption =    element_text(size = 7, family = "Arial", face = "italic", color="black"),
    plot.background = element_blank(),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing =      unit(0.1, "lines"), # facet_grid margin
    #panel.border =       element_blank(), # no sqaure, but two lines instead (axis.line)
    panel.background =   element_blank()
  ) +
  
  
  theme(strip.text.y.left = element_text(angle = 0)) +
  scale_color_manual(values=c('FALSE'=alpha('gray60',0.05),'TRUE'='blue')) +
  scale_x_reverse() +
  labs(x = "-log10(p fisher exact Loss)", col="significant")




# gains

plt <- out_cnv |> 
  dplyr::mutate(chr = gsub("^(.+):(.+)-(.+)$", "\\1", segment)) |> 
  dplyr::filter(chr %in% c("chr6", "chr12", "chr19")) |> 
  dplyr::mutate(start = as.numeric(gsub("^(.+):(.+)-(.+)$", "\\2", segment))) |> 
  dplyr::mutate(end = as.numeric(gsub("^(.+):(.+)-(.+)$", "\\3", segment))) |> 
  dplyr::mutate(pos = (start + end)/2) |> 
  dplyr::mutate(chr = factor(chr, levels=gtools::mixedsort(unique(as.character(chr))))) |> 
  dplyr::mutate(signi = p_gains_fisher < 0.01) |> 
  dplyr::mutate(x = -log10(p_gains_fisher))


plt <- rbind(
  plt |> dplyr::mutate(x = 0)
  ,
  plt
)


p2 = ggplot(plt , aes(x = x, y=pos / 1000000, group=segment, col=signi)) +
  facet_grid(rows = vars(chr), scale="free", space="free") +
  geom_line(lwd=0.1) +
  theme_bw() +
  theme(
    text =          element_text(size = 7, family = "Arial", face = "plain"),
    axis.text =     element_text(size = 7, family = "Arial", face = "plain", color="black"), # , angle=90, vjust =0.5
    axis.title.x =  element_text(size = 7, family = "Arial", face = "plain", color="black"), # , vjust = -0.2
    axis.title.y =  element_text(size = 7, family = "Arial", face = "plain", color="black"),
    axis.line =     element_line(linewidth = theme_cellpress_lwd),
    axis.ticks =    element_line(linewidth = theme_cellpress_lwd),
    
    strip.text =    element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.text.x =  element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.text.y =  element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.background = element_blank(), # clean as possible
    
    legend.title =  element_text(size = 7, family = "Arial", face = "plain", color="black"),
    legend.text =   element_text(size = 7, family = "Arial", face = "plain", color="black"),
    legend.position = 'bottom',
    legend.margin   = margin(t=-2),
    legend.key.size = unit(0.2, 'lines'),
    legend.key = element_blank(), # this should remove the white squares around the legend items, but seems to fail sometimes, probably due to the key.size above?
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    
    plot.title =      element_text(size = 7, family = "Arial", face = "plain", color="black"), # `title` covers both title and subtitle
    plot.subtitle =   element_text(size = 7, family = "Arial", face = "italic", color="darkgray"),
    plot.caption =    element_text(size = 7, family = "Arial", face = "italic", color="black"),
    plot.background = element_blank(),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing =      unit(0.1, "lines"), # facet_grid margin
    #panel.border =       element_blank(), # no sqaure, but two lines instead (axis.line)
    panel.background =   element_blank()
  ) +
  
  theme(strip.text.y.left = element_text(angle = 0)) +
  scale_color_manual(values=c('FALSE'=alpha('gray60',0.05),'TRUE'='red')) +
  labs(x = "-log10(p fisher exact Gain)", col="significant")




p1 + p2


ggsave("output/figures/gemi_cnv_zoom.png",width=8.5*0.975/2,height=11/3)



# details

out_cnv |> dplyr::filter(grepl("6:",segment) & p_losses_fisher < 0.01)
out_cnv |> dplyr::filter(grepl("12:",segment) & p_gains_fisher < 0.01)



