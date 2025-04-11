#______________________________________________________________________________#
# DNA methylation QC & Preprocessing----
#______________________________________________________________________________#
## Load libraries ----
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#______________________________________________________________________________#
## Raw data files ----
targets <- read.metharray.sheet(base=basedir, pattern="Datasheet.csv$")
RGSet <- read.metharray.exp(targets=targets, force=T)

#______________________________________________________________________________#
## Preprocessing & Filtering----
proc <- preprocessNoob(RGSet, offset=0, dyeCorr=T, verbose=T, dyeMethod="reference")

detP <- detectionP(RGSet, type="m+u") 
keep<- rowSums(detP <0.01) == ncol(proc)
proc <- proc[keep,]
failed <- colMeans(detP) > 0.05 # Samples with >5% failed positions

maskdata <- read.delim("EPIC.hg38.manifest.tsv.gz") #(Zhou et al. Nucleic Acids Research, 2017)
maskdata <- subset(maskdata, maskdata$MASK_general == F)
keep<- rownames(proc) %in% maskdata$probeID
proc <- proc[keep,]

annoEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
sexProbes <- rownames(annoEPIC)[which(annoEPIC[,"chr"] %in% c("chrX","chrY"))]
keep<- !rownames(proc) %in% sexProbes
proc <- proc[keep,]

#______________________________________________________________________________#
## Get Beta and M-values----
mVals<- getM(proc)
mVals<- na.omit(mVals)
# betas<- getBeta(proc)
# betas<- na.omit(betas)

#______________________________________________________________________________#
# Differentially Methylated Positions----
#______________________________________________________________________________#
## Load libraries ----
library(limma)

#______________________________________________________________________________#
## Prep data----
initial<- subset(targets, targets$Sample_Type == "I")
recurrent<- subset(targets, targets$Recurrent_Select == "TRUE")
identical(initial$GLASS_ID, recurrent$GLASS_ID)
SampUsed<- rbind(initial, recurrent)

MvalUsed <- mVals[,colnames(mVals) %in% SampUsed$Sample_ID] 
identical(colnames(MvalUsed), SampUsed$Sample_ID)
col_order<- SampUsed$Sample_ID
MvalUsed<- MvalUsed[, col_order]
identical(colnames(MvalUsed), SampUsed$Sample_ID)

#______________________________________________________________________________#
## Linear model----
resection<- factor(SampUsed$Sample_Type) # "I", "R"
patient<- factor(SampUsed$GLASS_ID) 

design <- model.matrix(~0 +resection + patient)
colnames(design) <- c(levels(resection),levels(patient)[-1])
fit <- lmFit(MvalUsed, design)
contMatrix <- makeContrasts(R-I, levels=design)

fit2<- contrasts.fit(fit, contMatrix)
fit2<- eBayes(fit2)

testfit<- as.data.frame(decideTests(fit2, method = "separate", adjust.method = "fdr", p.value = 0.05))
setDT(testfit, keep.rownames = "Probe_ID")
DMPs<- topTable(fit2, number = Inf, genelist = annoEPIC, coef = "R - I", sort.by = "B")
setDT(DMPs, keep.rownames = "Probe_ID")

DMPs<- merge(testfit, DMPs, by = "Probe_ID")

DMPs<- DMPs[order(-DMPs$B), ]
DMPs$log10P<- log10(DMPs$adj.P.Val)

#______________________________________________________________________________#

#______________________________________________________________________________#
## Signature Selection Probes ----
DMPsUsed<- DMPs[abs(DMPs$logFC) > 1,] 
DMPsUsed<- DMPsUsed[DMPsUsed$adj.P.Val < 1e-9,] #1389

Signature<- mVals[rownames(mVals) %in% DMPsUsed$Probe_ID,]
Signature<- as.data.frame(t(Signature))

# Calculate mean DNA methylation per sample
Signature$GMSscore = apply(Signature[,1:1389], 1, median, na.rm = TRUE)
Signature<- setDT(Signature, keep.rownames = "Sample_ID")

Signature$GMSclass<- ifelse(Signature$GMSscore < 0, "hypo", 
                             ifelse(Signature$GMSscore > 1, "hyper", "intermediate"))
#______________________________________________________________________________#

#______________________________________________________________________________#