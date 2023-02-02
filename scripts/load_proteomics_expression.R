#!/usr/bin/env R


if(!exists("metadata.glass.per.resection")) {
  source('scripts/load_metadata.R')
}

# needed to match id's to RNA Data
if(!exists('expression.glass.exon.metadata')) {
  source('scripts/load_rna-counts.R')
}


# raw data, useful for gene identifiers and 'Filtered' values (i.e. imputed in normalised data)

expression.proteomics.raw <- read.table("data/glass/Proteomics/2022-03-31_data_update/20210729_084405_GLASSNL_LGG_dia_ProteinReport.tsv", sep = "\t", quote = "", header = T) |> 
  dplyr::rename_with( ~ gsub('^X','',.x)) |> 
  dplyr::filter(grepl("Y-FGCZCont",PG.ProteinAccessions) == F) |> 
  dplyr::filter(PG.ProteinAccessions != 'P01892') |>  # one of the four duplicate genes
  dplyr::filter(PG.ProteinAccessions != 'Q04826') |> # one of the four duplicate genes
  dplyr::filter(PG.ProteinAccessions != 'Q29974') |> # one of the four duplicate genes
  dplyr::filter(PG.ProteinAccessions != 'P42167') |> # one of the four duplicate genes
  dplyr::mutate(PG.Gene.first = gsub(";.*?$","",PG.Genes)) |> 
  tibble::column_to_rownames('PG.Gene.first')

expression.proteomics.metadata <- expression.proteomics.raw |>
  dplyr::select('PG.ProteinAccessions', "PG.Genes", "PG.CellularComponent", "PG.BiologicalProcess", "PG.MolecularFunction")

stopifnot(nrow(expression.proteomics.raw) == 4828)
stopifnot(nrow(expression.proteomics.metadata) == 4828)



# match identifiers with RNA -- all except 3 are missing in the gtf
expression.proteomics.metadata <- expression.proteomics.metadata |> 
  dplyr::mutate(hugo_symbol_rna_prot_shared = gsub(";.*?$","",PG.Genes)) |> 
  dplyr::mutate(hugo_symbol_rna_prot_shared = 
                  dplyr::recode(hugo_symbol_rna_prot_shared,
                                'AARS' = 'AARS1',
                                'ADPRHL2' = 'ADPRS',
                                'ADSS' = 'ADSS2',
                                'ASNA1' = 'GET3',
                                'CARS' = 'CARS1',
                                'CRAD' = 'CRACD',
                                'DARS' = 'DARS1',
                                'EPRS' = 'EPRS1',
                                'FAM129B' = 'NIBAN2',
                                'FAM49A' = 'CYRIA',
                                'FAM49B' = 'CYRIB',
                                'GARS' = 'GARS1',
                                'H1F0' = 'H1-0',
                                'H1FX' = 'H1-10',
                                'H2AFY2' = 'MACROH2A2',
                                'H2AFY' = 'MACROH2A1',
                                'H2AFZ' = 'H2AZ1',
                                'H3F3A' = 'H3-3A',
                                'HARS' = 'HARS1',
                                'HIST1H1B' = 'H1-5',
                                'HIST1H1C' = 'H1-2',
                                'HIST1H1D' = 'H1-3',
                                'HIST1H1E' = 'H1-4',
                                'HIST1H2AB' = 'H2AC4',
                                'HIST1H2AG' = 'H2AC11',
                                'HIST1H2BA' = 'H2BC1',
                                'HIST1H2BK' = 'H2BC12',
                                'HIST1H2BO' = 'H2BC17',
                                'HIST1H4A' = 'H4C1',
                                'HIST2H2AB' = 'H2AC21',
                                'HIST2H2AC' = 'H2AC20',
                                'IARS' = 'IARS1',
                                'KARS' = 'KARS1',
                                'KIF1BP' = 'KIFBP',
                                'LARS' = 'LARS1',
                                'MARS' = 'MARS1',
                                'NARS' = 'NARS1',
                                'QARS' = 'QARS1',
                                'RARS' = 'RARS1',
                                'RYDEN' = 'SHFL',
                                'SARS' = 'SARS1',
                                'TARSL2' = 'TARS3',
                                'TARS' = 'TARS1',
                                'VARS' = 'VARS1',
                                'WARS' = 'WARS1',
                                'YARS' = 'YARS1',
                                'ICK' = 'CILK1',
                                'KIAA1107' = 'BTBD8',
                                'AKAP2' = 'PALM2AKAP2',
                                'MARC2' = 'MTARC2',
                                'BAP18' = 'C17orf49',
                                'C6orf203' = 'MTRES1',
                                'MARC1' = 'MTARC1',
                                'CNK3/IPCEF1' = 'IPCEF1',
                                'FAM129A' = 'NIBAN1',
                                'FAM192A'='PSME3IP1',
                                'FAM45BP' = 'DENND10P1',
                                'GSTT1' = '',
                                'LCHN' = 'DENND11',
                                'NUPL2' = 'NUP42',
                                'PALM2' = ''
                         )
                  )


stopifnot((expression.proteomics.metadata |> 
  dplyr::filter(duplicated(hugo_symbol_rna_prot_shared)) |> 
  dplyr::pull(hugo_symbol_rna_prot_shared) |> 
  unique()) == "")



# only 77 / 99 of the samples have metadata - 3 controls make sense, the other 19 unclear, most likely replicates
# 20210720_010_GLASSNL_LGG7_S307431.d.PG.Quantity
# 20210720_012_GLASSNL_LGG9_S307433.d.PG.Quantity
# 20210720_017_2_GLASSNL_LGG14_S307438.d.PG.Quantity
# 20210720_019_GLASSNL_LGG16_S307440.d.PG.Quantity
# 20210720_027_2_GLASSNL_LGG24_S307448.d.PG.Quantity
# 20210720_033_GLASSNL_LGG30_S307454.d.PG.Quantity
# 20210720_037_GLASSNL_LGG32_S307457.d.PG.Quantity
# 20210720_045_GLASSNL_LGG40_S307465.d.PG.Quantity
# 20210720_056_GLASSNL_LGG51_S307476.d.PG.Quantity
# 20210720_059_GLASSNL_LGG54_S307479.d.PG.Quantity
# 20210720_062_GLASSNL_LGG57_S307482.d.PG.Quantity
# 20210720_064_GLASSNL_LGG59_S307484.d.PG.Quantity
# 20210720_076_GLASSNL_LGG69_S307495.d.PG.Quantity
# 20210720_080_2_GLASSNL_LGG73_S307499.d.PG.Quantity
expression.proteomics.raw <- expression.proteomics.raw |>
  dplyr::rename(!!! ( # rename to names as used in the raw data
    metadata.glass.per.resection |> 
      dplyr::filter(!is.na(File_Name_Proteomics)) |> 
      dplyr::pull(File_Name_Proteomics, name=Sample_Name)
  )) |> 
  dplyr::select(metadata.glass.per.resection |> dplyr::filter(!is.na(ProtID)) |>  dplyr::pull(Sample_Name) |> sort())

stopifnot(sum(colnames(expression.proteomics.raw) %in% metadata.glass.per.resection$Sample_Name) == 55) # from 77 to 55
stopifnot(sum(colnames(expression.proteomics.raw) %in% metadata.glass.per.resection$Sample_Name == F) == 0)
stopifnot((metadata.glass.per.resection |> dplyr::filter(!is.na(ProtID)) |> dplyr::pull(Sample_Name)) %in% colnames(expression.proteomics.raw))




# normalised data contains 55 samples, not 99 or 77
expression.proteomics.normalised.imputed <- read.csv('data/glass/Proteomics/ProteinMatrix_30percentNA_cutoff_75percent_proteincutoff_MADnorm_MixedImputed_correct annotations_fixed-quotes_fixedspaces.csv',header=T) |> 
  dplyr::filter(X %in% c("HLA-B.1") == F) |>  # duplicated
  tibble::column_to_rownames('X') |> 
  dplyr::rename(!!! ( # rename to names as used in the raw data
    metadata.glass.per.resection |> 
      dplyr::filter(!is.na(ProtID)) |> 
      dplyr::pull(ProtID, name=Sample_Name)
  )) |> 
  dplyr::select(metadata.glass.per.resection |> dplyr::filter(!is.na(ProtID)) |>  dplyr::pull(Sample_Name) |> sort())


# ensure identifiers between raw and imputed table exist
stopifnot(nrow(expression.proteomics.normalised.imputed) == 3247)
stopifnot(ncol(expression.proteomics.normalised.imputed) == 55)
stopifnot(rownames(expression.proteomics.normalised.imputed) %in% rownames(expression.proteomics.metadata))
stopifnot(colnames(expression.proteomics.normalised.imputed) %in% (metadata.glass.per.resection |> dplyr::filter(!is.na(ProtID)) |> dplyr::pull(Sample_Name)))
stopifnot((metadata.glass.per.resection |> dplyr::filter(!is.na(ProtID)) |> dplyr::pull(Sample_Name)) %in% colnames(expression.proteomics.normalised.imputed))




# create the same matrix but with the raw values, to obtain 'Filtered' values
# this matrix consists of 1's and NA's
tmp <- expression.proteomics.raw |> 
  dplyr::mutate_all(function(arg) { return (ifelse(arg == "Filtered", as.numeric(NA), arg)) }) |> 
  dplyr::select(colnames(expression.proteomics.normalised.imputed)) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(rownames(expression.proteomics.normalised.imputed)) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::mutate_all(function(arg) { return (ifelse(is.na(arg), as.numeric(NA), 1)) })


stopifnot(colnames(tmp) == colnames(expression.proteomics.normalised.imputed))
stopifnot(rownames(tmp) == rownames(expression.proteomics.normalised.imputed))



# multipying the imputed matrix with the 1's and NA's will set the imputed values to NA and leave the rest as is
expression.proteomics.normalised.partial <- expression.proteomics.normalised.imputed * tmp
rm(tmp)




