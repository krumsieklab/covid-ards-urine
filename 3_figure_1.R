# Script to generate panels of Figure 1
# Once the PDFs are generated the figure was assembled in Adobe Illustrator

#### SET UP --------
### clear workspace and set directories ------
rm(list =ls())
# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### libraries ------
library(gridExtra)
library(ggh4x)
library(cowplot)
library(RColorBrewer)
source('custom_functions.R')
### input variables ------
# output directory
results_dir <- 'results/'
# number of molecules in volcano plots
num_mol <- 5
# significance threshold 
pcut <- 0.05
# annotation colors
annocols <- brewer.pal(n = 9, name = "Pastel1")
annocols[9] <- '#999999' 
plot_cols <- annocols[c(2, 1, 9)]
# datasets
datasets <- c("urine_metabo", "urine_proteo")
# ards groups
grp1 <- 'Bact-Seps'
grp2 <- 'Co19-ARDS'



#### MAIN --------
### compile 'between' stat result files ------
compiled_stats <- get_compiled_stats(
  sub_set=paste0('results/tmp_', datasets, '_between_stats.xlsx'), stat_type='between')

### create volcano and pathway plots for each dataset ------
volcano_list <- pathway_list <- list() 
for(dataset in datasets){
  ## volcano plots ----
  plot_mat <- read.xlsx(sprintf("results/tmp_%s_between_stats.xlsx", dataset))
  # change statistic direction to fit 
  plot_mat %<>% mutate(fold_change = -1*fold_change)
  volcano_list[[dataset]] <- get_volcano_plots (plot_mat=plot_mat, 
                                                pcut=pcut,
                                                grp1=grp1,
                                                grp2=grp2,
                                                num_mol=num_mol,
                                                plot_cols=plot_cols)
  
  ## pathway plots ----
  plot_mat <- read.xlsx(sprintf("results/tmp_%s_between_pathstats.xlsx", dataset), sheet='IndividualResults')
  outfile=NULL
  # select pathway groups per omic
  if(grepl('metabo', dataset)){
    plot_mat <- plot_mat [grep('Metabolism', plot_mat$pathway), ]  
  } else if (grepl('proteo', dataset)){
    # remove disease pathways
    plot_mat %<>% filter(!pathway %in%c( # kegg cancer pathways
      'Pathways in cancer', 'Breast cancer',
      'Proteoglycans in cancer',
      'Transcriptional misregulation in cancer',
      # parasitic diseases
      'Amoebiasis',  'Malaria',
      'Toxoplasmosis',  'Chagas disease',
      # infectious diseases
      'Pertussis', 'Tuberculosis',
      'Influenza A', 'Measles',
      'Legionellosis', 'Hepatitis B',
      'Epithelial cell signaling in Helicobacter pylori infection',
      'Herpes simplex virus 1 infection',
      'Human cytomegalovirus infection', 
      'Salmonella infection',
      'Pathogenic Escherichia coli infection',
      'Human papillomavirus infection',
      'Kaposi sarcoma-associated herpesvirus infection',
      # immune diseases 
      'Inflammatory bowel disease','Cushing syndrome',
      'Allograft rejection', 'Rheumatoid arthritis',
      # metabolic and endocrine diseases
      'AGE-RAGE signaling pathway in diabetic complications',
      'Non-alcoholic fatty liver disease'
    )) 
  } else {
    plot_mat %<>%  filter(pathway!='Total')
    outfile <- sprintf('%s/Figure1b_%s_case_control_path_break.pdf', results_dir, dataset)
  }
  pathway_list[[dataset]] <- get_pathway_plot(plot_mat=plot_mat,
                                              annocols=annocols,
                                              pcut=pcut,
                                              outfile=outfile, analyte_thresh=3)
}

### print volcano plots ------
pdf(sprintf('%s/Figure1a_urine_case_control_volcano.pdf', results_dir), 
    height=3, width=6, useDingbats = F)
grid.arrange(grobs=volcano_list,ncol=2)
dev.off()

### print pathway plots ------
pdf(sprintf('%s/Figure1b_urine_case_control_pathway.pdf', results_dir), 
    height=8, width=5)
grid.arrange(grobs=pathway_list,ncol=1)
dev.off()

### finished ------
print("Done!") 
print("Generated three PDF files with panels of Figure 1 in results folder!") 