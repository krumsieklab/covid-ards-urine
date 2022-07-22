# Script to compare two ARDS etiologies per omic
#    Generates supplementary files 2 and 3

### SET UP --------
## clear workspace and set directories ----
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir.create("results/", showWarnings = F, recursive = T)

# input variables ----
datasets <- c("urine_metabo","urine_proteo")
pwgroups <- list(urine_metabo='SUB_PATHWAY',
                 urine_proteo='kegg_db')
# groups to compare
complst <- list(comps = list(c("Group","Co19-ARDS","Bact-Seps")),
                check_groups = c("Bact-Seps", "Co19-ARDS"))
# adjusted p-value cutoff
pcut <- 0.05

## libraries ----
library(tidyverse)
source('custom_functions.R')

#### MAIN --------
stats_list <- path_list <- list()

### load, prepare, and perform between-ards analysis for each dataset ----
for (dataset in datasets) {
  tmp_stats_file <- sprintf("results/tmp_%s_between_stats.xlsx", dataset)
  tmp_path_file <- sprintf("results/tmp_%s_between_pathstats.xlsx", dataset)
  
  # load processed data
  D <- mt_load_se_xls(file=paste0('input/', dataset, '_processed.xlsx')) %>%
    # flag that data is logged
    mt_load_flag_logged()
  
  # get pathway annotations if there is an UniProt column (used for proteomics)
  if ("UniprotID" %in% (D %>% rowData() %>% colnames())) {
    D %<>%
      mt_anno_pathways_uniprot(
        in_col = "UniprotID",
        out_col = "kegg_db") %>%
      mt_anno_pathways_remove_redundant(feat_col = "UniprotID", pw_col = "kegg_db") %>%
      mt_write_pathways(pw_col = "kegg_db", file = 'results/supplementary_table_2_sheet3_proteins_kegg_pathway_annotations.xlsx')
  }
  
  ## perform analysis per comparison ----
  for (comp in complst$comps) {
    compnamebase <- paste0(comp[-1], collapse = "_")
    D %<>% between_ards_comparison(comp_name = sprintf("%s", compnamebase),
                                   comp_info = comp,
                                   pwgroup = pwgroups[[dataset]],
                                   p_adj_cut=pcut,
                                   path_outfile=tmp_path_file)
  }
  write_stats(D, out_file = tmp_stats_file)
  
  ## format the dataset output excel files ----
  # read intermediate outputs
  tmp_stats <- read.xlsx(tmp_stats_file)
  tmp_path <- read.xlsx(tmp_path_file, sheet='IndividualResults')
  
  # for each dataset, select desired columns from statistical result and pathway annotation outputs
  # rename columns from each where necessary
  if(dataset=='urine_proteo'){
    tmp_stats %<>% select(name, outcome, estimate, std_error, statistic, fold_change, p_value, adj_p, effect_high_in,
                          OlinkID, UniprotID, Panel, Panel_Version, kegg_db) %>%
      dplyr::rename(KEGG_Pathway_IDs=kegg_db)
    tmp_path  %<>% select(name, pathway, pathway_id, estimate, std.error, statistic,
                          fc, p.value, p.adj) %>%
      dplyr::rename(KEGG_ID=pathway_id, std_error=std.error, fold_change=fc,
                    p_value=p.value, adj_p=p.adj)
  } else if (dataset=='urine_metabo') {
    tmp_stats %<>% select(name, outcome,	estimate, std_error,
                          statistic, fold_change, p_value, adj_p, effect_high_in,
                          SUPER_PATHWAY, SUB_PATHWAY, COMP_ID, PUBCHEM, CAS, KEGG, HMDb)%>%
      dplyr::rename(HMDB=HMDb)
    tmp_path  %<>% select(name, pathway, color, estimate,
                          std.error, statistic, fc, p.value, p.adj) %>%
      dplyr::rename(SUB_PATHWAY=pathway,
                    SUPER_PATHWAY=color,
                    std_error=std.error,
                    fold_change=fc,
                    p_value=p.value,
                    adj_p=p.adj)
  }else{
    break('unidentified dataset!')
  }
  path_list[[dataset]] <- tmp_path
  stats_list[[dataset]] <- tmp_stats
} # END dataset loop
names(path_list) <- names(stats_list) <- c('Metabolomics', 'Proteomics')

### write supplementary files ----
write_dataframes (df_list=stats_list, out_file='results/supplementary_table_2_between_ards_stats.xlsx')
write_dataframes (df_list=path_list, out_file='results/supplementary_table_3_between_ards_pathway_annotations.xlsx')

### finished ----
print("Done! Finished per omics comparison of ARDS etiologies!")
print("Generated excel files with supplementary tables 2 and 3 in results folder!")
