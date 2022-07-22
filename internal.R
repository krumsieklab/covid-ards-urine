#### internal data processing script
# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath")))
# set working directory to location of source code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### init ----
library(maplet)
library(tidyverse)
library(magrittr)

## define all files
# clinical
file_clin <- file.path(data.makepath("COVID19/Choi_COVID19_data/ICU/clinical_lab_annotations/ICU_urine_annotations.xlsx"))
file_clin_sheet <- "annotations"
#file_clin_checksum <- "56e2eb1ef52a30c2b2e67d45bd7b3fb7"

# urine metabolomics
file_urine_metabo <- file.path(data.makepath("COVID19/Choi_COVID19_data/ICU/metabolomics/urine_metabolite_intensity.xlsx"))
#file_urine_metabo_checksum <- "269691828d3e9159970a83b27c58d920"
# urine proteomics
file_urine_prot <- file.path(data.makepath("COVID19/Choi_COVID19_data/ICU/olink/urine_proteomics.xlsx"))
#file_urine_prot_checksum <- "7e8b2c62599fc69b972af6c49852240b"
file_urine_prot_anno <- file.path(data.makepath("COVID19/Choi_COVID19_data/ICU/olink/proteomics_sampleinfo.xlsx"))
#file_urine_prot_anno_checksum <- "f1c8dbd752a95b3ca3140bb8b01c3d41"
# For loading any of the ICU datasets
# @returns D: a summarized experiment object
data_loader <- function(dataset) {
  
  if(!(dataset %in% c("urine_metabo","urine_proteo")))
    stop("dataset can only be urine_metabo or urine_proteo")
  if (dataset=="urine_metabo"){
    # load data
    D <- mt_load_metabolon_v1(file=file_urine_metabo, sheet = "OrigScale") 
    # load sample annotations
    D <- D %>%
      # verify data and annotations are unchanged
      #mt_files_checksum(file = file_urine_metabo, checksum = file_urine_metabo_checksum) %>%
      #mt_files_checksum(file = file_uclin, checksum = file_clin_checksum) %>%
      # load clinical data
      mt_anno_xls(file=file_clin, sheet = file_clin_sheet, anno_type = "samples", anno_id_col = "Sample_ID", data_id_col= "CLIENT SAMPLE ID")%>%
      # adjust the age variable with ">" sign
      mt_anno_mutate(anno_type = "samples", col_name='age', term=case_when(age%in%">90" ~ "90", TRUE~age)) %>%
      # make sure age is numeric
      mt_anno_mutate(anno_type = "samples", col_name = "age", term = as.numeric(as.matrix(age))) %>%
      # same names as in plasma dataset
      mt_anno_mutate(anno_type = "samples", col_name = "Group", term = Group %>% recode(`Bacterial sepsis ARDS`="Bact-Seps", `Covid-19 ARDS`="Co19-ARDS",`Non-sepsis control`="Cont") ) %>%
      # create COVID- vs non-COVID disease
      mt_anno_mutate(anno_type = "samples", col_name = "Disease", term = Group %>% recode(`Bact-Seps`="non-COVID", `Co19-ARDS`="COVID", `Cont`="non-COVID") ) %>%
      
      {.}
    
  } else if(dataset=="urine_proteo"){
    # inpute data
    anno_sheet <- "ID matching urine"
    # load data
    D <- mt_load_olink(file = file_urine_prot)
    
    # load sample annotations
    D <- D %>%
      # verify data and annotations are unchanged
      # mt_files_checksum(file = file_urine_prot, checksum = file_urine_prot_checksum) %>%
      #mt_files_checksum(file = file_uclin, checksum = file_clin_checksum ) %>%
      #mt_files_checksum(file = file_urine_prot_anno, checksum = file_urine_prot_anno_checksum) %>%
      # load clinical data
      mt_anno_xls(file = file_urine_prot_anno, sheet = anno_sheet, anno_type = "samples", anno_id_col = "sample id WCM-Q", data_id_col = "sample_id")%>%
      mt_anno_xls(file=file_clin, sheet = file_clin_sheet, anno_type = "samples", anno_id_col = "Prot_ID", data_id_col= "subject.id.WCM.NY")%>%
      # adjust the age variable with ">" sign
      mt_anno_mutate(anno_type = "samples", col_name='age', term=case_when(age%in%">90" ~ "90", TRUE~age)) %>%
      # make sure age is numeric
      mt_anno_mutate(anno_type = "samples", col_name = "age", term = as.numeric(as.matrix(age))) %>%
      # filter out samples without group information
      mt_modify_filter_samples(filter = !is.na(Group)) %>%
      mt_modify_filter_samples(filter = endsWith(sample.id.WCM.Q, "_2")) %>%
      # same names as in plasma dataset
      mt_anno_mutate(anno_type = "samples", col_name = "Group", term = Group %>% recode(`sepsis_ARDS`="Bact-Seps", `COVID19_ARDS`="Co19-ARDS", `CNTRL`="Cont")) %>%
      # create COVID- vs non-COVID disease
      mt_anno_mutate(anno_type = "samples", col_name = "Disease", term = Group %>% recode(`Bact-Seps`="non-COVID", `Co19-ARDS`="COVID", `Cont`="non-COVID") ) %>%
      
      {.}
  }
  
  # return 
  D
}

# run preprocessing for datasets
data_preprocessing <- function(D, dataset) {
  
  D <- D %>%
    mt_pre_filter_missingness(samp_max=0.5) %>%
    mt_pre_filter_missingness(feat_max=0.25) %>%
    {.}
  # normalization
  if(grepl("_proteo", dataset)){
    # undo log transformation of data for quot norm
    D %<>% mt_pre_trans_exp()
  } 
  # normalization, log, imputation
  D  %<>%
    mt_pre_norm_quot(feat_max = 0) %>%
    # transformation & imputation
    mt_pre_trans_log() %>% 
    mt_pre_impute_knn() %>% 
    {.}
  
  # average duplicate molecules --> specifically needed for proteins
  if (grepl("proteo", dataset)){
    D %<>% mt_modify_avg_features(group_col = 'name')  
  }
  # final dataset info
  D %<>%  # dataset info
    mt_reporting_data()
  
}
# make sure result directory exists
dir.create("results/", showWarnings = F, recursive = T)
datasets <- c( "urine_metabo", "urine_proteo") #

##### start analysis ----

# loop over datasets ----
for (dataset in datasets) {
  # load and preprocess
  D <-
    # load data
    data_loader(dataset=dataset) %>%
    # preprocess
    data_preprocessing(dataset=dataset)
  # select columns needed
  Y <- D %>% colData() %>% data.frame() %>% select(Subject_ID, sex, age, aki, platelet, pf, death, Group)
  D1 <- SummarizedExperiment(assay= assay(D), rowData=rowData(D), colData=Y)
  # save processed data 
  D1 %<>% mt_write_se_xls(file = paste0("input/", dataset, '_processed.xlsx'))
  
} # closing dataset loop

## finished
print("Done! per omics data processing!") 
print("Generated excel files with processed data in input folder!") 
