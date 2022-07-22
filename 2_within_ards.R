# Script to perform per omics per clinical manifestation comparison within ARDS etiologies
#    Generates supplementary files 4, 5

#### SET UP --------
### clear workspace and set directories ------
rm(list = ls())
# set working directory to location of source code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# make sure result directory exists
dir.create("results/", showWarnings = F, recursive = T)

### input variables ------
# datasets
datasets <- c("urine_metabo", "urine_proteo")
# clinical manifestations
outcomes <- data.frame(outcome = c('aki', 'platelet', 'pf', 'death'),
                       outcome_type = c(rep('numeric', 3), 'binary'),
                       outcome_mode = c(rep('numeric', 3), 'character'))
# groups to compare
complst <- list(comps = list(c("Group","Co19-ARDS","Bact-Seps")), 
                check_groups = c("Bact-Seps", "Co19-ARDS"))
# adjusted p-value cutoff
pcut <- 0.05
sfile_num <- 4

### libraries ------
library(tidyverse)
source('custom_functions.R')

#### MAIN --------
### load, prepare, and perform within-ards analysis on each dataset ------
for (dataset in datasets) {
  intermediate_file <- sprintf("results/tmp_%s_within_stats.xlsx", dataset)
  formatted_file <- sprintf("results/supplementary_table_%d_%s_within_ards_clinical_manifestation_stats.xlsx", sfile_num, dataset)
  
  # load processed data 
  D <- mt_load_se_xls(file=paste0('input/', dataset, '_processed.xlsx')) %>%
    # flag that data is logged
    mt_load_flag_logged() %>%
    # format clinical manifestation columns as per data type (numeric/factor)
    outcome_type_conversion(outcomes)
  
  # get pathway annotations if there is an UniProt column (used for proteomics)
  if ("UniprotID" %in% (D %>% rowData() %>% colnames())) {
    D %<>%  
      mt_anno_pathways_uniprot(
        in_col = "UniprotID", 
        out_col = "kegg_db") %>% 
      mt_anno_pathways_remove_redundant(feat_col = "UniprotID", pw_col = "kegg_db")
  }
  
  ## perform analysis per outcome per group ----
  for (i in 1:nrow(outcomes)) {
    for (group in complst$check_groups) {
      comp_name <- sprintf("%s_%s", outcomes$outcome[[i]], group)
      D %<>%
        within_ards_analysis(
          outcome_info = outcomes[i,],
          comp_name = comp_name,
          keep = group, pcut =pcut)
    }
  }
  # save intermediate output
  write_stats(D, out_file = intermediate_file)
  
  ## format the dataset output excel files ----
  # read all sheets
  tmp_stats <- intermediate_file %>%
    excel_sheets() %>%
    purrr::set_names() %>%
    map(read_excel, path = intermediate_file)
  
  # for each dataset, select desired columns from the statistical result for mortality and
  # other clinical outcomes; rename columns from each where necessary
  if(dataset=='urine_proteo'){ 
    for(i in 1:length(tmp_stats)){
      if(grepl('death', names(tmp_stats)[i])){
        tmp_stats[[i]] %<>% select(name, outcome, estimate, std_error, statistic, fold_change, p_value, adj_p,
                                   effect_high_in, OlinkID, UniprotID, Panel, Panel_Version, kegg_db) %>%
          dplyr::rename(KEGG_Pathway_IDs=kegg_db)
      } else {
        tmp_stats[[i]] %<>% select(name, outcome, estimate, std_error, statistic, p_value, adj_p,
                                   OlinkID, UniprotID, Panel, Panel_Version, kegg_db) %>%
          dplyr::rename(KEGG_Pathway_IDs=kegg_db)
        
      }
    }
  } else if(dataset=='urine_metabo'){ 
    for(i in 1:length(tmp_stats)){
      if(grepl('death', names(tmp_stats)[i])){
        tmp_stats[[i]] %<>% select(name, outcome,	estimate, std_error, 
                                   statistic, fold_change, p_value, adj_p,	effect_high_in,
                                   SUPER_PATHWAY, SUB_PATHWAY, COMP_ID, PUBCHEM, CAS, KEGG, HMDb)%>%
          dplyr::rename(HMDB=HMDb)
      } else {
        tmp_stats[[i]] %<>% select(name, outcome,	estimate, std_error, 
                                   statistic, p_value, adj_p,	
                                   SUPER_PATHWAY, SUB_PATHWAY, COMP_ID, PUBCHEM, CAS, KEGG, HMDb)%>%
          dplyr::rename(HMDB=HMDb)
      }
    }
  } else{
    break('unidentified dataset!')
  }
  ## write formatted output ----
  write_dataframes (df_list=tmp_stats, out_file=formatted_file)
  sfile_num <- sfile_num + 1
} # END dataset loop


### finished ----
print("Done! per omics per clinical manifestation comparison within ARDS etiologies!") 
print("Generated excel files with supplementary tables 4, 5 in results folder!") 