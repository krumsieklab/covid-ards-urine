# Script to generate panels of Figure 2
# Once the networks are generated in Cytoscape the figure was assembled in Adobe Illustrator
# !! NOTE: Cytoscape needs to be running and allow incoming connections to
#          generate GGM networks (see also RCy3 package). !!
#
# Once the networks are printed in cytoscape save it using File --> save as --> 'supplementary_file_1.cys'
# Then remove self-loops: Edit --> remove self-loops


#### SET UP --------
### clear workspace and set directories ------
rm(list =ls())
# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
### libraries ------
library(stringi) # string match
library(RColorBrewer) # colors
source('custom_functions.R') # customized functions

### input variables ------
# my color palette
annocols <- brewer.pal(n = 9, name = "Pastel1")
annocols[9] <- '#999999'; annocols[7] <- '#80DEEA'
annocols[3:4] <- c('#80CBC4', '#C5E1A5'); annocols[5] <- '#E6EE9C' 
# datasets
datasets <- c("urine_metabo", "urine_proteo")

# groups to compare
complst <- list(comps = list(c("Group","Co19-ARDS","Bact-Seps")), 
                check_groups = c("Bact-Seps", "Co19-ARDS"))
group_names <- sub('-', '_', complst$check_groups)
#### MAIN --------
### compile 'between' stat result files ------
compiled_stats <- get_compiled_stats(
  sub_set=paste0('results/tmp_', datasets, '_between_stats.xlsx'), stat_type='between')

### get ggm edge list ------
ggm_edges <- get_multiomics_network()

### build a node attribute data frame ------
node_attributes <- data.frame(node_name=unique(c(as.matrix(ggm_edges$source), as.matrix(ggm_edges$target))), 
                              node_type="Metabolite") %>% 
  dplyr::mutate(node_type = case_when(node_name%in%(compiled_stats %>% 
                                                      filter(mol_type=='Proteins') %>% 
                                                      pull(mol_uname)) ~ "Protein", 
                                      TRUE~"Metabolite")) %>% 
  dplyr::left_join(compiled_stats, by=c("node_name"= "mol_uname"))

### create Figure 2a ------

elist_to_cytoscape(node_attributes=node_attributes,
                   edge_list=ggm_edges,
                   outcome_name='overall', 
                   group_names=group_names,
                   collection_name='ggm_bh_05',
                   whole_net=T)

### create Figure 2b ------
# subgraph with 2 neighbour
get_subnetwork_of_interest(ggm_edges, seed_node=c('GP6', 'tiglyl.carnitine..C5.'), neigh_order=2, 
                     node_attributes=node_attributes, 
                     group_names=group_names,outcome_name='between-ards',
                     collection_name='met_path_subnet', column_name='Co19_ARDS_Bact_Seps_pval_sign') 

### finished ------
print("Done!") 
print("Generated panels of Figure 2 in Cytoscape interface!") 