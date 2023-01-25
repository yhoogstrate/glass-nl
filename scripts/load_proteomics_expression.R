#!/usr/bin/env R


# raw data, useful for gene identifiers and 'Filtered' values (i.e. imputed in normalised data)

expression.proteomics.raw <- read.table("data/glass/Proteomics/2022-03-31_data_update/20210729_084405_GLASSNL_LGG_dia_ProteinReport.tsv", sep = "\t", quote = "", header = T) |> 
  dplyr::rename_with( ~ gsub('^X','',.x)) |> 
  dplyr::filter(grepl("^Y-FGCZCont",PG.ProteinAccessions) == F) |> 
  dplyr::filter(PG.ProteinAccessions != 'P01892') |>  # one of the four duplicate genes
  dplyr::filter(PG.ProteinAccessions != 'Q04826') |> # one of the four duplicate genes
  dplyr::filter(PG.ProteinAccessions != 'Q29974') |> # one of the four duplicate genes
  dplyr::filter(PG.ProteinAccessions != 'P42167') |> # one of the four duplicate genes
  dplyr::mutate(PG.Gene.first = gsub(";.*?$","",PG.Genes)) |> 
  tibble::column_to_rownames('PG.Gene.first')

expression.proteomics.metadata <- expression.proteomics.raw |>
  dplyr::select('PG.ProteinAccessions', "PG.Genes", "PG.CellularComponent", "PG.BiologicalProcess", "PG.MolecularFunction")


stopifnot(nrow(expression.proteomics.raw) == 4881)
stopifnot(nrow(expression.proteomics.metadata) == 4881)


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
  dplyr::select(metadata.glass.per.resection |> dplyr::filter(!is.na(File_Name_Proteomics)) |>  dplyr::pull(File_Name_Proteomics))

stopifnot(sum(colnames(expression.proteomics.raw) %in% metadata.glass.per.resection$File_Name_Proteomics) == 77)
stopifnot(sum(colnames(expression.proteomics.raw) %in% metadata.glass.per.resection$File_Name_Proteomics == F) == 0)




# normalised data contains 55 samples, not 99 or 77
expression.proteomics.normalised.imputed <- read.csv('data/glass/Proteomics/ProteinMatrix_30percentNA_cutoff_75percent_proteincutoff_MADnorm_MixedImputed_correct annotations_fixed-quotes_fixedspaces.csv',header=T) |> 
  dplyr::filter(X %in% c("HLA-B.1") == F) |>  # duplicated
  dplyr::rename_with( ~ gsub('^[A-Z]+_[A-Z]+_','',.x)) |> # do not use sample identifiers for metadata
  tibble::column_to_rownames('X') |> 
  dplyr::rename(!!! ( # rename to names as used in the raw data
    metadata.glass.per.resection |> 
      dplyr::filter(!is.na(proteomics_imputed_id)) |> 
      dplyr::select(File_Name_Proteomics, proteomics_imputed_id) |> 
      dplyr::pull(proteomics_imputed_id, name=File_Name_Proteomics)
  )) 


# ensure identifiers between raw and imputed table exist
stopifnot(nrow(expression.proteomics.normalised.imputed) == 3247)
stopifnot(ncol(expression.proteomics.normalised.imputed) == 55)
stopifnot(rownames(expression.proteomics.normalised.imputed) %in% rownames(expression.proteomics.metadata))



# create the same matrix but with the raw values, to obtain 'Filtered' values
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



expression.proteomics.normalised.partial <- expression.proteomics.normalised.imputed * tmp

rm(tmp)




