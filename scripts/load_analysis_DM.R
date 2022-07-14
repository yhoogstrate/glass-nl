#!/usr/bin/env R

# These are the differential methylated regions and probes from the main analysis from Wies

# all these data are hg19(!)


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}




#DMPs_IvR_FDR <- read.csv('data/DM_regions_probes_Wies/DMPs_IvR_FDR 0.05_10052022.csv')

# Differentially methylated regions (DMR) primary recurrency ----

DMRs_IvR_FDR <- read.csv('data/DM_regions_probes_Wies/DMRs_IvR_FDR 0.05_20052022.csv') |> 
  dplyr::mutate(sign = ifelse(meandiff < 9, -1, 1)) |> 
  dplyr::mutate(logp = sign * -log10(Fisher)) |> 
  dplyr::arrange(start, end)



# (De)-methylation signature ----


tmp <- read.csv('data/DM_regions_probes_Wies/meanDNAm_DMPs IvR (FDR 1e-9 Delta 1)_11062022.csv') |> 
  dplyr::select(Sample_Name, mean.DNAm) |> 
  dplyr::rename(mean.DNA.methylation.signature = mean.DNAm)


metadata.glass.per.resection <- metadata.glass.per.resection |> 
  dplyr::left_join( tmp , by=c('Sample_Name'='Sample_Name') ,suffix = c("", "") )

rm(tmp)


# ggplot(metadata.glass.per.resection, 
#        aes(x=mean.DNA.methylation.signature,y=lts.down,
#            shape=methylation.sub.diagnosis,
#            col=CDKN2AB)) + 
#   # geom_smooth() + 
#   geom_point(size=1.8) +
#   theme_bw()
# 
# 
# plt <- metadata.glass.per.resection |> 
#   dplyr::mutate(IDH_HG_IDH_ratio = log(A_IDH_HG_cal/A_IDH_cal)) |> 
#   dplyr::select(lts.up1, lts.up2, lts.up3, lts.down, mean.DNA.methylation.signature, IDH_HG_IDH_ratio ) |> 
#   dplyr::filter(!is.na(lts.up1)) |> 
#   dplyr::filter(!is.na(lts.up2)) |> 
#   dplyr::filter(!is.na( lts.up3)) |> 
#   dplyr::filter(!is.na( lts.down)) |> 
#   dplyr::filter(!is.na(mean.DNA.methylation.signature)) |> 
#   dplyr::filter(!is.na(IDH_HG_IDH_ratio))
# 
# corrplot::corrplot.mixed(cor(plt), upper = 'pie',lower='ellipse')
# 
# corrplot::corrplot(cor(plt))




### transform mean.DNA.methylation.signature from gamma to normal ----


fit.data  <- metadata.glass.per.resection |>  
  dplyr::filter(!is.na(mean.DNA.methylation.signature)) |>  
  dplyr::pull(mean.DNA.methylation.signature,name='Sample_Name')
fit.data <- - fit.data + max(fit.data) + 0.01

fit.g = fitdistrplus::fitdist(fit.data , "gamma")
fit.ln = fitdistrplus::fitdist(fit.data, "lnorm")
fit.wb = fitdistrplus::fitdist(fit.data, "weibull")

par(mfrow = c(2, 2))
plot.legend <- c("gamma", "lognormal", "weibull")
denscomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)
qqcomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)
cdfcomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)
ppcomp(list(fit.g, fit.ln, fit.wb), legendtext = plot.legend)

dev.off()

# seems gamma fit
fit <- fit.g
rm(fit.g,fit.ln,fit.wb)


metadata.glass.per.resection <- metadata.glass.per.resection |> 
  dplyr::left_join(
    data.frame(mean.DNA.methylation.signature.norm = - qnorm(pgamma(fit.data, shape = fit$estimate['shape'], rate = fit$estimate['rate'] ))) |> 
      tibble::rownames_to_column('Sample_Name'),
    by=c('Sample_Name'='Sample_Name'),suffix = c('',''))


plot(metadata.glass.per.resection$mean.DNA.methylation.signature, metadata.glass.per.resection$mean.DNA.methylation.signature.norm)


# Transformation, no re-ordering
stopifnot(order(metadata.glass.per.resection$mean.DNA.methylation.signature) == order(metadata.glass.per.resection$mean.DNA.methylation.signature.norm))


rm(fit, fit.data)



