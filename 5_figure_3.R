# Script to generate panels of Figure 3
# Once the PDFs are generated the figure was assembled in Adobe Illustrator

#### SET UP --------
### clear workspace and set directories ------
rm(list =ls())

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
### libraries ------
library(tidyverse) # tidy
library(magrittr) # %<>%
library(RColorBrewer) # colors 
library(readxl) # read excel sheets at once
library(gridExtra) # plots in columns
source('custom_functions.R') # customized functions


### input variables ------
# make sure result directory exists
results_dir <- 'results/'
# my color palette
annocols <- c('#E6EE9C', '#C5E1A5', '#80CBC4', '#999999')
# adjusted p-value cutoff
pcut <- 0.05
# datasets
datasets <- c("urine_metabo", "urine_proteo")
# outcomes == clinical manifestations
outcomes <- data.frame(outcome = c('aki', 'platelet', 'pf', 'death'),
                       outcome_type = c(rep('numeric', 3), 'binary'),
                       outcome_mode = c(rep('numeric', 3), 'character'))
# groups to compare
complst <- list(comps = list(c("Group","Co19-ARDS","Bact-Seps")), 
                check_groups = c("Bact-Seps", "Co19-ARDS"))

# compile within stat files
compiled_stats <- get_compiled_stats(
  sub_set=paste0('results/tmp_', datasets, '_within_stats.xlsx'))

# comparing covid-19 death variability in urine and plasma ----
# fluids
fluids <- c('urine', 'plasma'); 
# plot colors
plot_cols <- annocols[c(2, 4, 3)]
# group lables
grp1 <- 'Survivors'; grp2 <- 'Non-survivors'; 
# number of molecules named in volcano plots
num_mol <- 5
# initialize lists for plots
volcano_list <- pathway_list <- list() 
# Figure 3a loop over fluids -----
for(fluid in fluids){
  # input file
  x <- sprintf('results/tmp_%s_proteo_within_stats.xlsx', fluid)
  # all sheets of input file read
  this_mat <- x %>%
    excel_sheets() %>%
    purrr::set_names() %>%
    map(read_excel, path = x)
  # choose mortality sheet for covid ----
  plot_mat <- this_mat[['death_Co19-ARDS']]
  # rename fc column name to foldchange
  if(length(grep('^fc$', names(plot_mat)))>0){
    plot_mat %<>% dplyr::rename(fold_change=fc)
  }
  # get volcano plot
  volcano_list[[fluid]] <- get_volcano_plots (plot_mat=plot_mat, pcut=pcut,
                                                grp1=grp2, grp2=grp1,
                                                num_mol=num_mol, 
                                              plot_cols=plot_cols,
                                              ymax=2.5
                                              )
  # choose mortality sheet for bacterial ards ----
  plot_mat <- this_mat[['death_Bact-Seps']]
  # rename fc column name to foldchange
  if(length(grep('^fc$', names(plot_mat)))>0){
    plot_mat %<>% dplyr::rename(fold_change=fc)
  }
  # get volcano plot
  volcano_list[[sprintf('%s_bs', fluid)]] <- get_volcano_plots (plot_mat=plot_mat, pcut=pcut,
                                              grp1=grp2, grp2=grp1,
                                              num_mol=num_mol, 
                                              plot_cols=plot_cols,ymax=2.5
                                              )
} # fluid loop ends
# print volcano plots
pdf(sprintf('%s/Figure3a_death_urine_plasma_volcano.pdf', results_dir), 
    height=6, width=6, useDingbats = F)
grid.arrange(grobs=volcano_list,ncol=2, nrow=2)
dev.off()

# Figure 3b box plots urine top death proteins ----

# read proteomics processed data
Dp <- mt_load_xls(file='input/urine_proteo_processed.xlsx', 
                 samples_in_rows = F,
                 is_mt_write_se_xls_output = T)

# proteins with fold change >2 and adj_p <-0.05 selected from results
pois <- c('FABP4', '4E-BP1', 'GH', 'NT-proBNP', 'IGFBP-1', 'CXCL16', 'CXCL9', 'IGFBP-2',
          'HO-1', 'CCL3', 'MB', 'FGF-21', 'HSP 27', 'VEGFA')
# initialize empty list for plots
box_list <- list()
# loop over proteins of interest
for(poi in pois){
  # extract protein expression from the data for covid ards
  u_df <- extract_mol_profiles(D=Dp, poi=poi) %>% 
    filter(Group=='Co19-ARDS') %>% select(-Group)
  # format data for plotting
  plot_mat <- u_df %>% reshape2::melt(id='death') %>%
    mutate(death=case_when(death=='Yes' ~ 'C19-NS', death=='No' ~ 'C19-S', TRUE~death))
  plot_mat$death <- factor(plot_mat$death, levels=rev(c('C19-NS', 'C19-S')))
 
   # plot
  p <- ggplot(plot_mat, aes(x=death, y=value, fill=death))+
    geom_boxplot () + geom_jitter(alpha=0.4, size=0.8) + 
    theme_bw() + ggtitle(poi) +
    xlab('') + ylab('Log2 abundance (a.u.)')+ 
    theme(text = element_text(size=12)) +
    theme(legend.position="none")+
    scale_fill_manual(values=plot_cols[c(1, 3)])
  
  # save plot in list
  box_list[[poi]]<- p
} # proteins of interest loop ends

# print in pdf
pdf('results/Figure3b_urine_death_boxplots_2fc.pdf', height = 6 , width=10)
grid.arrange(grobs=box_list[1:length(box_list)], ncol=5)
dev.off()
### finished ------
print("Done!") 
print("Generated PDF files with panels of Figure 3 in results folder!")
