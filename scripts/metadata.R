#!/usr/bin/env R

# ---- load libs ----


library(tidyverse)


# ---- per resection ----


#metadata.glass.per.resection <- 
  


# ---- per patient ----s

metadata.glass.per.patient <- 
  read.csv('data/glass/Clinical data/Cleaned/Survival data_GLASS_12082021.csv') %>% 
  dplyr::mutate(X=NULL) %>% 
  dplyr::mutate(data_of_birth = NULL) %>% 
  dplyr::mutate(Date_of_Death = as.Date(Date_of_Death , format = "%d-%m-%Y")) %>% 
  dplyr::mutate(Date_of_Diagnosis = as.Date(Date_of_Diagnosis , format = "%d-%m-%Y")) %>% 
  dplyr::mutate(Date_Last_Followup = as.Date(Date_Last_Followup , format = "%d-%m-%Y")) %>% 
  dplyr::mutate(overall.survival = difftime(Date_of_Death , Date_of_Diagnosis, units = 'days')) %>% 
  dplyr::mutate(time.until.last.followup = difftime(Date_Last_Followup, Date_of_Diagnosis, units = 'days'))



  
# reshape?!
