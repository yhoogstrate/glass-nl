#!/usr/bin/env R


# raw data, useful for gene identifiers and 'Filtered' values (i.e. imputed in normalised data)

expression.proteomics.raw <- read.table("data/glass/Proteomics/2022-03-31_data_update/20210729_084405_GLASSNL_LGG_dia_ProteinReport.tsv", sep = "\t", quote = "", header = T) |> 
  dplyr::rename_with( ~ gsub('^X','',.x)) |> 
  dplyr::mutate(protein.id = gsub(';.+$','',PG.ProteinAccessions)) |> 
  tibble::column_to_rownames('protein.id')

expression.proteomics.metadata <- expression.proteomics.raw |>
  dplyr::select('PG.ProteinAccessions', "PG.Genes", "PG.CellularComponent", "PG.BiologicalProcess", "PG.MolecularFunction")


# only 77 / 99 of the samples have metadata - 3 controls make sense, the other 19 unclear
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
expression.proteomics.imputed <- read.csv('data/glass/Proteomics/ProteinMatrix_30percentNA_cutoff_75percent_proteincutoff_MADnorm_MixedImputed_correct annotations_fixed-quotes_fixedspaces.csv',header=T) |> 
  dplyr::filter(X %in% c("HLA-B.1") == F) |>  # duplicated
  dplyr::rename_with( ~ gsub('^[A-Z]+_[A-Z]+_','',.x)) |> # do not use sample identifiers for metadata
  tibble::column_to_rownames('X')


#sum(colnames(expression.proteomics.imputed) %in% == 55)

