# Customized functions used in one or more number-prefixed scripts.
#
# MAIN FUNCTIONS ---
#  Used by one or more number-prefixed scripts. Number(s) in parentheses is the number prefix(es) 
#   of the script(s) it is used in (e.g. (s: 1) indicates the function is used in 
#   script 1_bewtween_ards)
# (1) between_ards_comparison (s: 1) - compare between two ards groups
# (2) within_ards_analysis (s: 2) - compare within an ards group
# (3) write_stats (s: 1, 2) - format the output of comparisons
# (4) write_dataframes (s: 1, 2) - write a list of dataframes in a customized format
# (5) outcome_type_conversion (s: 2) - convert clinical manifestations to according to 
#       datatype (numeric/factor)
# (6) get_volcano_plots (s: 3) - generate volcano barplots for between group analysis
# (7) get_pathway_plot (s: 3) - generate stacked barplots for between group analysis
# (8) get_compiled_stats (s: 3, 4, 5) - compile all statistics within/between group comparison in one dataframe
# (9) within_ards_barplots (s: 4) - stack barplots per clinical manifestation all omics, all ards group
# (10) get_multiomics_network (s: 5) - generate a multiomic ggm
# (11) elist_to_cytoscape - helper to print networks to cytoscape using an edgelist
# (12) get_manifestation_subnetwork (s: 5) - print subnetwork to cytoscape for given a clinical 
#       manifestation
#
# HELPER FUNCTIONS ---
#   Used internally by main functions
# (13) add_self_edges - simple helper function to add self edges for nodes without any connections
# (14) edge_paste - simple helper paste to be used for interactions list

# libraries ----
library(tidyverse)
library(magrittr)
library(GeneNet)
library(readxl)
library(igraph)
library(RCy3)
library(maplet)
library(glue)
library(openxlsx)

### MAIN FUNCTIONS ------
# main function for between ards analysis
#   create glm comparing two sample groups and generate stats pathway bar plot
between_ards_comparison <- function(D, # SE object
                                    comp_name, # name for statistic result
                                    comp_info, # vector containing: (1) formula, (2) ards group 1, (3) ards group 2
                                    pwgroup, # column name containing pathway groups
                                    p_adj_cut, # adjusted p-value cutoff
                                    path_outfile # file to write bar plot data to
) {
  # formula provided
  symgroup <- sym(comp_info[1])
  
  ##create linear model ----
  D %<>%
    mt_stats_univ_lm(
      samp_filter = (!!symgroup %in% comp_info[2:3]), 
      formula      = as.formula(glue("~{comp_info[1]}")),
      stat_name = comp_name) 
  
  ## p-value correction ----
  D <- D %>% 
    mt_post_fold_change(stat_name = comp_name) %>%
    mt_post_multtest(stat_name = comp_name, method = "BH") %>%
    {.}
  
  ## generate stats pathway bar plots ----
  if ("SUB_PATHWAY" %in% pwgroup) {
    # if metabolon data, create super and sub pathways bar plots
    D %<>%
      mt_plots_stats_pathway_bar(stat_list = comp_name,  
                                 feat_filter = p.adj < !!enquo(p_adj_cut),
                                 group_col = "SUB_PATHWAY", color_col = "SUPER_PATHWAY",
                                 y_scale = "count", assoc_sign_col = "statistic", 
                                 outfile=path_outfile)
  } else {
    # else, create pathway group bar plots
    D %<>%
      mt_plots_stats_pathway_bar(stat_list = comp_name,  
                                 feat_filter = p.adj < !!enquo(p_adj_cut),
                                 group_col = pwgroup,
                                 y_scale = "count", 
                                 assoc_sign_col = "statistic", outfile=path_outfile)
  }
  D
}

# main function for within ards analysis
#   create glm to analyze an outcome within a group of interest
within_ards_analysis <- function(D, # SE object
                                 outcome_info, 
                                 comp_name,
                                 keep,
                                 pcut
) {
  
  ## create linear model for binary and numeric outcomes ----
  if (outcome_info$outcome_type %in% c("binary", "numeric")) {
    if (!is.na(keep)) {
      # filter down to group of interest
      D %<>%
        mt_stats_univ_lm(
          formula     = as.formula(glue("~{outcome_info$outcome}")),  
          samp_filter = (Group %in% keep),
          stat_name   = comp_name)
      
    } else {
      # use all samples
      D %<>%
        mt_stats_univ_lm(
          formula     = as.formula(glue("~{outcome_info$outcome}")),  
          stat_name   = comp_name
        )
    }
  } else {
    stop("only binary and numeric outcomes implemented so far")
  }
  
  ## fold changes and multiple testing correction ----
  D <- D %>%
    {if (outcome_info$outcome_type=="binary"){mt_post_fold_change(., stat_name = comp_name)}else{.}} %>%
    mt_post_multtest(stat_name = comp_name, method = "BH")
  D %<>% 
    mt_reporting_stats(stat_name = comp_name, stat_filter =  p.adj < !!enquo(pcut))
  
  D
}

# write out stats in excel files
#   extract all statistical results table from an SE object, apply desired formatting, and write
#   each table to a sheet in an excel file
write_stats <- function(D, # SE object
                        out_file # file to write stats to
){
  ### extract all stats entries ------
  S <- D %>% maplet:::mtm_res_get_entries("stats")
  allcomps <- S %>% purrr::map("output") %>% purrr::map("name") %>% unlist()
  
  ### format and output all stats entries ------
  wb <- openxlsx::createWorkbook()
  for (i in 1:length(S)) {
    df <- S[[i]]$output$table
    
    ## add direction ----
    if ("groups" %in% names(S[[i]]$output) && length(S[[i]]$output$groups)==2) {
      # generate vector of directions for indexing
      inds <- as.numeric(df$estimate>0)+1
      # translate into names
      df$effect_high_in <- S[[i]]$output$groups[inds]
    }
    
    ## rename and arrange columns ----
    df%<>% dplyr::rename(var=var, outcome=term, covariates=formula,	
                         std_error=std.error, p_value=p.value, adj_p=p.adj)
    # got fc?
    fc_flag <- grep('^fc$', names(df))
    if(length(fc_flag)>0){df%<>% dplyr::rename(fold_change=fc)}
    # add rowdata
    df %<>% bind_cols(rowData(D) %>% data.frame())
    #got kegg?
    kegg_flag <- grep('kegg_db', names(df))
    if(length(kegg_flag)>0){
      # convert it to string
      df$kegg_db <- apply(df, 1, FUN=function(x)toString(x$kegg_db))
    }
    df %<>% select(name, outcome, covariates, p_value, adj_p, everything())
    df %<>% arrange(adj_p)
    
    ## create and write to worksheet ----
    name <- S[[i]]$output$name
    ws=openxlsx::addWorksheet(wb,sheetName=name)
    openxlsx::writeData(wb=wb, sheet=name, x=df)
    
    ## create and add a style to the column headers and the body ----
    headerStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center', textDecoration = 'bold')
    bodyStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center')
    addStyle(wb, sheet = name, bodyStyle, rows = 1:(nrow(df)+1), cols = 1:ncol(df), gridExpand = TRUE)
    addStyle(wb, sheet = name, headerStyle, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  }
  ### write excel file ------
  openxlsx::saveWorkbook(wb, file=out_file, overwrite=T)
}

# write a list of dataframes in a customized format
#   for a user-provided list of dataframes, apply desired formatting and write each dataframe
#   to a sheet in an excel file
write_dataframes <- function(df_list, # list of dataframes
                             out_file # file to write dataframes to
){
  ### format and output all dataframes ------
  wb <- openxlsx::createWorkbook()
  for (i in 1:length(df_list)) {
    df <- df_list[[i]]
    
    ## create and write to worksheet ----
    name <- names(df_list)[i]
    ws=openxlsx::addWorksheet(wb,sheetName=name)
    openxlsx::writeData(wb=wb, sheet=name, x=df)
    
    ## create and add a style to the column headers and body ----
    headerStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center', textDecoration = 'bold')
    bodyStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center')
    addStyle(wb, sheet = name, bodyStyle, rows = 1:(nrow(df)+1), cols = 1:ncol(df), gridExpand = TRUE)
    addStyle(wb, sheet = name, headerStyle, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  }
  ### write excel file ------
  openxlsx::saveWorkbook(wb, file=out_file, overwrite=T)
}

# convert outcomes to right data types
#   for each outcome, use a dataframe that includes the outcome column name and data type to
#   convert outcome to assigned data type
outcome_type_conversion <- function(D, # SE object
                                    outcomes # a dataframe with columns: (1) outcome, (2) outcome_type, and (3) outcome_mode
) {
  # ensure numeric outcomes are numeric
  for (outcome in (outcomes %>% filter(outcome_mode=="numeric") %>% .$outcome)) {
    D %<>% mt_anno_mutate(anno_type = "samples", col_name = outcome, term = as.numeric(!!sym(outcome)))
  }
  # ensure binary outcomes are factor
  for (outcome in (outcomes %>% filter(outcome_type=="binary") %>% .$outcome)) {
    D %<>% mt_anno_mutate(anno_type = "samples", col_name = outcome, term = as.factor(!!sym(outcome)))
  }
  D
}

# generate volcano plots for between ards analysis
#   create volcano plot highlighting the significant molecules within each ards group
get_volcano_plots <- function(plot_mat, # plot data matrix
                              pcut=0.05, # significance threshold
                              grp1, # ards group 1
                              grp2, # ards group 2
                              num_mol, # number of molecules to annotate
                              plot_cols, # point colors
                              ymax=NULL
){
  
  ### create y-axis value and group ------
  plot_mat %<>% dplyr::mutate(logp = -1 * log10(adj_p), 
                              point_col = case_when(plot_mat$adj_p>pcut  ~ 'none', TRUE ~ effect_high_in))
  if((plot_mat %>% pull(point_col) %>% unique())=='none'){
    plot_cols <- "#999999"
  }
  
  ### create volcano plot ------
  p <- plot_mat %>% 
    ggplot(aes(x = fold_change, y = logp)) +
    # horizontal at p-value threshold
    geom_hline(yintercept = min(plot_mat$logp[which(plot_mat$adj_p<=pcut)]), 
               linetype='dashed', color='grey') +
    annotate(geom="text", x=min(plot_mat$statistic)+2, 
             y=min(plot_mat$logp[which(plot_mat$adj_p <=pcut)])+0.2, 
             label="FDR X%", color="grey") +
    geom_point(aes(color=point_col)) + xlab('Log2 fold change') + 
    ylab("- Log10 (adj. p)") + theme_bw() +
    theme(legend.position='none')+
    # add group names 
    annotate(geom="text", x=min(plot_mat$statistic)+3, y=0, 
             label=sprintf('High in %s', grp2), color="grey") +
    annotate(geom="text", x=max(plot_mat$statistic)-3, y=0, 
             label=sprintf('High in %s', grp1), color="grey") +
    scale_color_manual(values=plot_cols)
  
  ### add names of the top molecules to plot ------
  data_annotate <- plot_mat %>% .[order(plot_mat$adj_p, decreasing = F), ] %>% .[1:num_mol, ]
  p <- p + ggrepel::geom_text_repel(data = data_annotate,
                                    aes(label = name),
                                    size=3)+
    theme(text = element_text(size=12))
  if(!missing(ymax)){
    
    p <- p+ coord_cartesian(ylim = c(0, ymax))
  }
  
  p
}

# generate pathway bar plots for between ards analysis
#   created stacked bar plot with number of significant analytes within each pathway for each
#   ards group
get_pathway_plot <- function(plot_mat, # plot data matrix
                             annocols, # bar group colors
                             pcut, # significance threshold
                             outfile, # pdf file to save to
                             analyte_thresh
){
  
  y_lab <- sprintf('# analytes with padj < %.2f', pcut)
  
  ### prepare data for ggplot ------
  plot_mat %<>% dplyr::rename(label=pathway) %>% unique()
  # count total number of significant molecules
  plot_mat %<>% dplyr::left_join(., plot_mat %>% dplyr::group_by(label) %>%
                                   dplyr::count(name) %>% dplyr::group_by(label) %>% 
                                   dplyr::count(label), by='label') %>% 
    dplyr::rename(path_wt=n)  %>% filter(path_wt>=analyte_thresh)
  # determine subset of significant molecules belonging to each group
  plot_mat %<>% dplyr::left_join(plot_mat %>% dplyr::group_by(label) %>% 
                                   dplyr::count(statistic<0) %>% 
                                   dplyr::filter(`statistic < 0` ==TRUE) %>% 
                                   dplyr::select(- `statistic < 0`), by='label') %>% 
    dplyr::rename(high_grp1=n) %>%
    dplyr::left_join(plot_mat %>% dplyr::group_by(label) %>% 
                       dplyr::count(statistic>0) %>% 
                       dplyr::filter(`statistic > 0` ==TRUE) %>% 
                       dplyr::select(- `statistic > 0`), by='label') %>% 
    dplyr::rename(high_grp2=n) %>% 
    dplyr::select(label, path_wt, high_grp1, high_grp2) %>% unique()
  # sort pathways by total number molecules - max at the top
  path_order <- plot_mat$label[order(plot_mat$path_wt, decreasing = F)]
  max_mol <- max(plot_mat$path_wt)
  # reshape data for plotting
  plot_mat <- reshape2::melt(plot_mat, id=c('label')) %>% 
    dplyr::filter(variable!='path_wt')%>% 
    dplyr::mutate(value=as.numeric(as.matrix(value)),
                  value=case_when(is.na(value)~0, TRUE~value))
  # order pathways based on path_wt
  plot_mat$label <- factor(plot_mat$label, levels=path_order)
  
  ### create path bar plot ------
  p <- ggplot(plot_mat, aes(fill=variable, x=label, y=value)) +
    geom_bar(stat='identity') + theme_bw()+
    xlab('') + ylab(y_lab) + ggtitle('')+
    theme(text=element_text(size=8))+
    theme (axis.ticks = element_blank())+
    coord_flip(clip = 'off') + theme(legend.position='none')+
    scale_fill_manual(values=annocols[c(2, 1)]) +
    scale_y_continuous(breaks = seq(0, max_mol, 2))
  
  return(p)
} 

# compile stat results from multiple sheets and omics
#   combine stat results from a list of either 'within_ards' or 'between_ards' excel files
get_compiled_stats <- function(sub_set, # list of files
                               stat_type='within' # analysis type / file substring
) {
  ### read in each input file ------
  stat_list<- lapply(sub_set, FUN=function(x){
    res <- x %>%
      excel_sheets() %>%
      purrr::set_names() %>%
      map(read_excel, path = x)
    res
  })
  # add dataset names to list
  names(stat_list) <- gsub('..*/', '', sub_set) %>% gsub('tmp_', '', .) %>%
    strsplit(., sprintf('_%s_stats.xlsx', stat_type)) %>% unlist()
  
  ### join all omics using selected columns ------
  stat_data <- lapply(1:length(stat_list[[1]]), FUN=function(i){
    res <- do.call(plyr::rbind.fill, list(stat_list[['urine_metabo']][[i]] %>% 
                                            dplyr::select(name, p_value, adj_p, statistic, SUB_PATHWAY) %>%
                                            dplyr::mutate(mol_type='Metabolites'), 
                                          
                                          stat_list[['urine_proteo']][[i]] %>% 
                                            dplyr::select(name, p_value, adj_p, statistic, UniprotID) %>%
                                            dplyr::mutate(mol_type='Proteins')))
    
    ## create new columns ----
    res %<>% dplyr::mutate (
      pval_score = -1*log10(adj_p), # log10 adjusted p-value
      asso_type = case_when(statistic > 0 ~ 1, statistic < 0 ~ -1), # effect direction
      sig_fdr5 = case_when(pval_score >= 1.3 ~ "sig", pval_score < 1.3 ~ "nsig"), # significance at 0.05
      sig_fdr10 = case_when(pval_score >= 1 ~ "sig", pval_score < 1 ~ "nsig"), # significance at 0.10
      pval_sign = pval_score * asso_type # directed pscore
    )
    
    res
  })
  # add group names
  names(stat_data) <- names(stat_list[[1]])
  
  ### combine all data into one large data frame ------
  compiled_stats <- lapply(1:length(stat_data), FUN=function(x){
    # create unique names and then order by the names
    stat_data[[x]] %<>% dplyr::mutate(mol_uname=make.names(name, unique = T))
    stat_data[[x]]  <- stat_data[[x]][order(stat_data[[x]]$mol_uname), ]
    # remove name columns (because they will duplicate otherwise)
    res <- stat_data[[x]] %>% dplyr::select(-mol_uname, -name, -mol_type)
    # add the clinical manifestation and ards group name in the column names
    names(res) <- paste0(gsub("-", "_", names(stat_data)[x]), "_", names(res))
    # add the names back
    res <- dplyr::bind_cols(res, stat_data[[x]] %>% select(mol_uname, name, mol_type))
    
    res
  }) %>% plyr::join_all(by=c('mol_uname', 'name', 'mol_type'), type='left')
  # rename 'name' column to 'mol_name'
  compiled_stats %<>% dplyr::rename(mol_name=name)
  
  compiled_stats
  
}

# plot within ards analysis results in barplots
#   for each omic, generate a stacked barplot for each outcome with the number of significant
#   molecules; colored based on whether the molecules are significant in ards group1, 
#   ards group2, or both
within_ards_barplots <- function(compiled_stats, # df of all omics results combined
                                 outcomes, # a dataframe with columns: (1) outcome, (2) outcome_type, and (3) outcome_mode
                                 group_names, # vector of groups names
                                 sig_col_suf='fdr5', # suffix of significance column
                                 annocols # bar group colors
){
  
  tmp_stats <- compiled_stats
  tmp_list <- list()
  
  ### for each omic, collect information to be plotted ------
  for(this_omics in unique(tmp_stats$mol_type)){
    # filter for this omics
    compiled_stats <- tmp_stats %>% filter(mol_type%in%this_omics)
    # empty data frame
    tmp <- data.frame(regulation=c('both', 'c', 'none', 'b'))
    
    ## for each outcome, record stats ----
    for(outcome in outcomes$outcome){
      grp1_col <- sprintf('%s_%s_sig_%s', outcome, group_names[1], sig_col_suf)
      grp2_col <- sprintf('%s_%s_sig_%s', outcome, group_names[2], sig_col_suf)
      if(grp1_col %in% names(compiled_stats) && grp2_col %in%names(compiled_stats )){
        # create column to record the stats
        tmp %<>% mutate(!!outcome := c(compiled_stats  %>% #sig in both
                                         filter(!!sym(grp1_col)=='sig' & !!sym(grp2_col)=='sig') %>% nrow(), 
                                       compiled_stats  %>% # sig only in CO
                                         filter(!!sym(grp1_col)=='sig' & !!sym(grp2_col)=='nsig') %>% nrow(), 
                                       compiled_stats  %>% # sig in none
                                         filter(!!sym(grp1_col)=='nsig' & !!sym(grp2_col)=='nsig') %>% nrow(), 
                                       compiled_stats  %>% # sig only in BS
                                         filter(!!sym(grp1_col)=='nsig' & !!sym(grp2_col)=='sig') %>% nrow()))
        
      }
    }
    tmp_list[[this_omics]] <- tmp
  }
  
  ### plot results for each outcome ------
  plot_list <- list()
  for(attribute in outcomes$outcome){
    ## reshape and order data ----
    plot_mat <- reshape2::melt(tmp_list) %>% filter(variable == attribute) %>% 
      filter(!regulation =='none')
    plot_mat$L1 <- factor(plot_mat$L1, levels=c("Proteins", "Metabolites"))
    
    ## generate plots ----
    # handle 'no results' case
    if(sum(plot_mat$value)==0){
      plot_list[[attribute]] <- ggplot() +
        geom_text(aes(x=0,y=0, label="No significant results"), size=10) +
        xlab("omics") + ylab("# analytes p < threshold") + theme_classic() +
        theme(text = element_text(size=15))+
        ggtitle (attribute) + theme(legend.position='none') +
        coord_flip(clip = 'off')
    } else {
      # create plot bar
      plot_mat$regulation <- factor(plot_mat$regulation, levels=c("c", "both", "b", "none"))
      plot_list[[attribute]] <- 
        ggplot(plot_mat, aes(x=L1, y=value, fill=regulation, label=value)) + 
        geom_bar(stat="identity", position = 'stack') +
        geom_text(size = 3, position = position_stack(vjust = 0.5))+
        xlab("omics") + ylab("# analytes p < threshold") + theme_classic() +
        theme(text = element_text(size=15))+
        ggtitle (attribute) + theme(legend.position='none') +
        scale_fill_manual(values=annocols, name='significant in', labels=c("Both", "CO-19", 'BS'))+
        coord_flip(clip = 'off')
    }
  }
  
  plot_list
}

# function to extract edge list from multiomics networks
#   generates a multiomics network and extracts the edge list; returns an edge list consisting 
#   of all significant edges (including 'island edges' - 'edges' from single nodes)
get_multiomics_network <- function(datasets = c("urine_metabo", "urine_proteo"),   # datasets
                                   formula_conf = "sex + age + Group",   # confounder formula
                                   ggm_thresh=0.05 # pvalue threshold;
                                   ){
 
  outcomes <- list(data.frame(outcome="sex", outcome_type="binary", outcome_mode="character"),
                   data.frame(outcome="age", outcome_type="numeric", outcome_mode="numeric"),
                   data.frame(outcome="Group", outcome_type="ordinal", outcome_mode="character")
  ) %>% do.call(rbind,.)
  
  
  ### input and preprocess data ------
  Ds <- list()
  for (dataset in datasets){
    # load data in a SE
    D <- mt_load_se_xls(file=paste0('input/', dataset, '_processed.xlsx')) %>%
      # flag that data is logged
      mt_load_flag_logged() %>%
      # adjust the age variable with ">" sign
      mt_anno_mutate(anno_type = "samples", col_name='age', term=case_when(age%in%">90" ~ "90", TRUE~age)) %>%
      # convert confounders to correct data types
      outcome_type_conversion(outcomes) %>%
      # correct for confouders and groups
      mt_pre_confounding_correction(as.formula(glue("~ {formula_conf}"))) %>%
      # for coding convenience (empty statement to terminate %>% pipe)
      {.}
    Ds[[dataset]] <- D
  }
  names(Ds) <- c('Metabolites','Proteins')
  
  ### assemble necessary data ------
  # group molecules by type
  mol_types <- lapply(Ds, FUN=function(d) {
    res <- d %>% rowData() %>% as_tibble() %>% select(name)
    res %<>% mutate(unique_name=make.names(name, unique=T))
    res
  })
  # get samples common to all 3 omics
  com_samples <- lapply(Ds, FUN=function(D) D %>% colData() %>% 
                          as_tibble() %>% select(Subject_ID) %>% 
                          unlist()) %>% Reduce(intersect, .)
  # create list of dataframes, one for each molecule type, consisting of common samples
  data_mat <- lapply(Ds, FUN=function(D) {
    res <- D %>% assay() %>% t() %>% as_tibble %>% filter(colData(D)$Subject_ID%in%com_samples) %>% data.frame()
    rownames(res) <- D %>% colData() %>% as_tibble %>% filter(colData(D)$Subject_ID%in%com_samples) %>% 
      select(Subject_ID) %>% unlist()
    colnames(res) <- D %>% rowData() %>% as_tibble() %>% select(name) %>% unlist() %>% make.names(unique = T)
    res <- res[order(rownames(res)), ]
    res
  })
  
  ### create Triomics GGMs ------
  ## format network data ----
  this_mat <- bind_cols(data_mat$Proteins, data_mat$Metabolites)
  pcor_mat <- ggm.estimate.pcor(as.matrix(this_mat), method = "dynamic", verbose = F)
  pval_mat <- network.test.edges(pcor_mat, plot = F, verbose = F)
  pval_mat$p.adj.bh <- p.adjust(pval_mat$pval, method="BH")
  pval_mat$p.adj.bon <- p.adjust(pval_mat$pval, method="bonferroni")
  
  ## generate network from pcor matrix ----
  tmp <- pcor_mat %>% 
    graph_from_adjacency_matrix(mode='undirected', weighted = T) %>% 
    igraph::simplify()
  
  ### extract edge list from network ------
  ggm_edges <- cbind.data.frame(get.edgelist(tmp), edge_attr(tmp)$weight)
  names(ggm_edges) <- c("source", "target", "pcor_val")
  
  ## filter edges based on pvalues ----
  ggm_edges %<>% dplyr::filter(abs(pcor_val)>=min(abs(pval_mat$pcor[pval_mat$p.adj.bh<=ggm_thresh])))
  
  ## add 'island edges' ----
  # all nodes possible in the network
  all_nodes <- c(as.matrix(mol_types$Proteins$unique_name),
                 as.matrix(mol_types$Metabolites$unique_name))
  # nodes with at least one edge
  con_nodes <- c(as.matrix(ggm_edges%>% pull(source)), 
                 as.matrix(ggm_edges%>% pull(target)))
  # nodes with no edge
  island_nodes <- setdiff(all_nodes, con_nodes)
  # island edges with the node itself
  island_edges <- data.frame(source=island_nodes, target=island_nodes, pcor_val=1) 
  # add nodes with no sig edge to edge list
  ggm_edges <- dplyr::bind_rows(ggm_edges, island_edges)
  
  ## add edge_type based on sign of pcor ----
  ggm_edges %<>% dplyr::mutate(edge_type = case_when(pcor_val > 0 ~ "pos", 
                                                     pcor_val < 0 ~ "neg"))%>%
    data.frame(.,stringsAsFactors=FALSE) %>% 
    mutate(edge_id=purrr::map2_chr(source, target, edge_paste))
  
  ggm_edges
}

# function to print networks to cytoscape
#    node_attributes must have following columns:
#       node_name, mol_name, node_type, 'outcome_name'_'group_name'_pval_sign
#    edge_list must have columns named source , target, edge_type
elist_to_cytoscape <- function(node_attributes, # df containing nodes and all associated stat res
                               edge_list, # df of all significant edges
                               outcome_name='test', # outcome name to be used from node_attributes 
                               group_names=NULL, # group name to be used from node_attributes 
                               collection_name='test', # name of the cytoscape collection s
                               whole_net=F, # generate whole network?
                               column_name =NULL # to be used to color nodes
){
  
  ### prepare network data ------
  # generate numeric node ids
  node_attributes %<>% mutate(id=as.character(1:nrow(node_attributes)))
  # convert source and target to numeric values based on ids of nodes
  edge_list <- node_attributes %>% dplyr::select(id, node_name) %>% 
    dplyr::left_join(edge_list,.,by=c("source"= "node_name")) %>% 
    dplyr::rename(from=source, source=id)
  edge_list <- node_attributes %>% dplyr::select(id, node_name) %>% 
    dplyr::left_join(edge_list,.,by=c("target"= "node_name")) %>% 
    dplyr::rename(to=target, target=id)
  # data frame specification needed for cytoscape
  edge_list <- data.frame(edge_list, stringsAsFactors=FALSE)
  node_attributes <- data.frame(node_attributes, row.names = node_attributes$id, stringsAsFactors = F)
  # name styles of cytoscape networks based on the outcome_name
  this_style <- this_net <- outcome_name 
  
  ### load the network to cytoscape ------
  RCy3::createNetworkFromDataFrames(edges=edge_list[,c('source', 'target')], title=this_net, collection=collection_name)
  # workaround for a bug in cytoscape
  edge_list %<>% dplyr::mutate(key=paste(source, "(interacts with)", target))
  # get SUID for edges and match with edge attributes
  cy_edges <- getTableColumns(table = 'edge')
  cy_edges <- cy_edges[order(cy_edges$name), ]
  edge_list <- edge_list[order(edge_list$key), ]
  edge_list$cpSUID <- cy_edges$SUID
  
  ### add edges and nodes ------
  # add edge attributes
  RCy3::loadTableData(subset(edge_list, select = -c(source, target)), 
                      data.key.column = 'cpSUID', table.key.column = 'SUID',
                      table = 'edge')
  
  # add node attributes
  RCy3::loadTableData(node_attributes, data.key.column = 'id', table = 'node')
  
  ### apply styling to network ------
  # prepare style variables
  style_name <- this_style
  nodeLabels <- RCy3::mapVisualProperty('node label','mol_name','p')
  nodeShapes <- RCy3::mapVisualProperty('node shape','node_type','d',
                                        c('Protein',  'Metabolite'),
                                        c ('DIAMOND', 'ELLIPSE'))
  edgeStyles <- RCy3::mapVisualProperty('edge line type', 'edge_type', 'd', 
                                        c("neg", "pos"), c("LONG_DASH", "SOLID"))
  # create style
  RCy3::createVisualStyle(style.name=style_name, base.url = 'http://localhost:1234/v1',
                          mappings = list(nodeLabels,nodeShapes,edgeStyles)) 
  # apply style
  RCy3::setVisualStyle(style.name=style_name, base.url = 'http://localhost:1234/v1',
                       network = this_net)
  if(length(group_names)>1){
    if(whole_net){
      RCy3::setNodeColorMapping(
        table.column = "node_type",
        table.column.values = c('Protein',  'Metabolite'),
        mapping.type = "d",
        colors= c('#75B771' ,'#DB6769', '#2C86B5'),
        style.name=style_name, base.url = 'http://localhost:1234/v1',
        network = this_net)
      # setting the layout
      RCy3::layoutNetwork(layout.name = 'cose', base.url = 'http://localhost:1234/v1')
    } else if(!missing(column_name)){
      max_val <- max(abs(node_attributes[column_name]), na.rm=T)
      RCy3::setNodeColorMapping(
        table.column = column_name,
        table.column.values = c(-1*max_val, 0, max_val),
        mapping.type = "c",
        colors= c('#E74C3C' ,'#F0F3F4', '#3498DB'),
        style.name=style_name, base.url = 'http://localhost:1234/v1',
        network = this_net)
      # setting the layout
      RCy3::layoutNetwork(layout.name = 'cose', base.url = 'http://localhost:1234/v1')
      
     }
    else{
      RCy3::setNodeCustomBarChart(columns=c(sprintf('%s_%s_pval_score', outcome_name, group_names[1]),
                                            sprintf('%s_%s_pval_score', outcome_name, group_names[2])),
                                  type = "GROUPED", colors = list("#80CBC4", "#E6EE9C"), range = NULL,
                                  orientation = "VERTICAL", colAxis = FALSE, rangeAxis = FALSE,
                                  zeroLine = FALSE, axisWidth = 0.25, axisColor = "#000000",
                                  axisFontSize = 1, separation = 0, slot = 1,
                                  style.name = style_name, base.url = 'http://localhost:1234/v1')
      setNodeColorDefault(new.color='#F8F9F9',
                          style.name = style_name, base.url = 'http://localhost:1234/v1')
    }
  }
}

# print network to cytoscape given an attribute/clinical manifestation
#    node_attributes must have following columns:
#       node_name, mol_name, node_type, 'outcome_name'_'group_name'_pval_sign, and
#       'outcome_name'_'group_name'_pval_score
#    edge_list must have columns named source , target, edge_type
get_subnetwork_of_interest <-  function(ggm_edges, # network
                                        seed_node=NULL, # if only this node should be used for plotting
                                        neigh_order=0, # neighbor order to be used
                                        neigh_nodes='all', # if neighbors of all sig nodes be used
                                        node_attributes, # node attributes of g
                                        outcome_name, # outcome name to be used from node_attributes 
                                        group_names, # group name to be used from node_attributes 
                                        collection_name, # name of the cytoscape collection s
                                        sig_thresh=1.3, # suffix of significant column
                                        column_name =NULL # to be used to color nodes
) 
{
  ### create graph from edgelist and extract basic object ------
  tmp_g <- graph_from_edgelist(ggm_edges %>% select(source, target) %>% as.matrix, 
                               directed = FALSE) %>% 
    set_edge_attr(., 'weight', value = ggm_edges$pcor_val) %>%
    set_edge_attr(., 'edge_type', value = ggm_edges$edge_type) %>%
    set_vertex_attr(., 'node_name', value = V(.)$name)
  
  ## extract dataframe of the graph ----
  df <- igraph::as_data_frame(tmp_g, 'both')
  # add selected node attributes to it
  col_ids <- grep(column_name, names(node_attributes))
  # add columns to vertices
  df$vertices <- df$vertices %>% 
    dplyr::left_join(node_attributes %>% 
                       select(node_name, node_type, names(node_attributes)[col_ids]), 
                     by='node_name')
  
  ### create new graph object with node attributes from dataframe ------
  g <- graph_from_data_frame(df$edges,directed = F, vertices = df$vertices) %>% 
    igraph::simplify() # remove self loops
  # significant nodes in outcome of interest
  
  this_sig_nodes <- lapply(group_names, FUN=function(group_name) {
    which(abs(vertex_attr(g)[[column_name]])>=sig_thresh)}) %>% 
    unlist()
  
  # should plot only selected significant nodes? 
  if(!missing(seed_node)){
    seed_node <- which(V(g)$name%in%seed_node)
    tmp <- this_sig_nodes[which(this_sig_nodes %in%seed_node)]
    if(length(tmp)>0){
      this_sig_nodes <- tmp
    } else{
      warning('Could not find the seed node, returning results with all sig_nodes')
    }
  }
  # should have neighbours?
  if (neigh_order >0) {
    # should we neighbours of all nodes?
    if(neigh_nodes=='all'){
      this_sig_nodes <- neighborhood(g, nodes=this_sig_nodes, order=neigh_order) %>%
        lapply(., FUN=function(x) which(V(g)$name%in%names(x))) %>% 
        unlist()
    } else{ # if not, select neighbours of only provided nodes like IL6...
      this_sig_nodes <- c(this_sig_nodes, 
                          neighborhood(g, nodes=which(V(g)$name%in%make.names(neigh_nodes, unique = T)), order=neigh_order) %>%
                            lapply(., FUN=function(x) which(V(g)$name%in%names(x))) %>% unlist())
    }
  }
  
  ### induce subgraph of those sig nodes ------
  this_net <- induced_subgraph(g, this_sig_nodes)
  # the edge list of that subgraph
  this_elist <- get.edgelist(this_net) %>% data.frame()
  names(this_elist) <- c('source', 'target')
  # add edge ids to that edge list
  this_elist %<>% dplyr::mutate(edge_id=purrr::map2_chr(source, target, edge_paste)) %>%
    dplyr::left_join(ggm_edges %>% select(edge_id, edge_type), by='edge_id')
  # nodes with at least one edge
  conn_nodes <- unique(c(as.matrix(this_elist%>% pull(source)), 
                         as.matrix(this_elist%>% pull(target))))
  # are there island nodes?
  if(length(conn_nodes) < length(this_sig_nodes)){
    # add the nodes which are islands
    this_elist <- add_self_edges(all_nodes = V(g)$name[this_sig_nodes], edge_list = this_elist)
  }
  
  ### print to cytoscape ------
  elist_to_cytoscape (node_attributes=node_attributes, 
                      edge_list=this_elist, 
                      outcome_name=outcome_name, 
                      group_names=group_names,
                      collection_name=collection_name, column_name = column_name)
}
### HELPER FUNCTIONS ------

# function to add self edges for nodes without any connections
#    edge_list must have columns named source , target
add_self_edges <- function(all_nodes, edge_list){
  # nodes with at least one edge
  con_nodes <- c(as.matrix(edge_list%>% pull(source)), 
                 as.matrix(edge_list%>% pull(target)))
  # nodes with no edge
  island_nodes <- setdiff(all_nodes, con_nodes)
  if(length(island_nodes)>0){
    # island edges with the node itself
    island_edges <- data.frame(source=island_nodes, target=island_nodes, pcor_val=1) 
    # add nodes with no sig edge to edge list?
    edge_list <- dplyr::bind_rows(edge_list, island_edges)
  }
  return(edge_list)
}

# paste to be used for interactions
edge_paste <- function(x, y){
  paste(sort(c(as.matrix(x), as.matrix(y))), collapse='(interacts with)')}

# define segments in y-Axis for 'ggplot2'
gg_gap <- function(plot,ylim,segments,tick_width,rel_heights,vjust=0,margin=c(top=1,right=2,bottom=1,left=1),...){
  #check whether segments is list
  if (!is.list(segments)){
    segments=list(segments)
  }
  #get and check y limits
  if (all(missing(ylim),is.null(plot$coordinates$limits$y))){
    stop('ylim is undefined')
  }else if(ylim[1]==ylim[2]){
    stop('ylim should not be the same number')
  }else if (missing(ylim)){
    ylim=plot$coordinates$limits$y
  }
  #check segments in order from small to large or from large to small
  for (j in 1:length(segments)) {
    seg1=segments[[j]][1]
    seg2=segments[[j]][2]
    if (seg1 > seg2){
      if (ylim[1]<ylim[2]){ #y-axis is from small to large
        msg=paste0('No.',j,' segment: c(',seg1,',',seg2,') is wrong. It should be ','c(',seg2,',',seg1,')')
        stop(msg)
      }
    }else if(seg1 < seg2){
      if (ylim[1]>ylim[2]){ #y-axis is from large to small
        msg=paste0('No.',j,' segment: c(',seg1,',',seg2,') is wrong. It should be ','c(',seg2,',',seg1,')')
        stop(msg)
      }
    }else if(seg1==seg2){
      msg=paste0('No.',j,' segment: c(',seg1,',',seg2,') is wrong. tick_width should not be equal')
      stop(msg)
    }
  }
  # check segments vectors sequence from large to small or from small to large
  if (length(segments)>=2){
    if (ylim[1] < ylim[2]){
      for (k in 2:length(segments)) {
        pre.2=segments[[k-1]][2]
        suf.1=segments[[k]][1]
        if (pre.2 > suf.1){
          pre=paste0('c(',segments[[k-1]][1],',',segments[[k-1]][2],')')
          suf=paste0('c(',segments[[k]][1],',',segments[[k]][2],')')
          msg=paste0('Segments ',k-1,' and ',k,': ',pre,',',suf,' are wrong. They should be ',suf,',',pre)
          stop(msg)
        }
      }
    }else if (ylim[1] > ylim[2]){
      for (k in 2:length(segments)) {
        pre.2=segments[[k-1]][2]
        suf.1=segments[[k]][1]
        if (pre.2 < suf.1){
          pre=paste0('c(',segments[[k-1]][1],',',segments[[k-1]][2],')')
          suf=paste0('c(',segments[[k]][1],',',segments[[k]][2],')')
          msg=paste0('Segments ',k-1,' and ',k,': ',pre,',',suf,' are wrong. They should be ',suf,',',pre)
          stop(msg)
        }
      }
    }
  }
  if (ylim[1] < ylim[2]){
    #check the minimum of segments must be more than min of ylim
    if (min(unlist(segments)) <= ylim[1]) stop('the minimum of segments must be more than the minium of ylim')
    #check the maximum of segments must be lower than maximum of ylim
    if (max(unlist(segments)) > ylim[2]) stop('the maximum of segments must be lower than maximum of ylim')
  }else if (ylim[1] > ylim[2]){
    #check the minimum of segments must be more than min of ylim
    if (min(unlist(segments)) <= ylim[2]) stop('the minimum of segments must be more than the minium of ylim')
    #check the maximum of segments must be lower than maximum of ylim
    if (max(unlist(segments)) > ylim[1]) stop('the maximum of segments must be lower than maximum of ylim')
  }
  #auto add tick_width if missing
  if (missing(tick_width)){
    tick_width=rep(abs(ylim[2]-ylim[1])/10,(length(segments)+1))
  }
  #check and add tick_width
  if ((length(tick_width)-length(segments)) < 1){
    int_len=length(tick_width)
    for (m in (int_len+1):(length(segments)+1)) {
      tick_width[m]=tick_width[int_len]
    }
  }
  seg_heights=0
  y_heights=1
  #check and add seg_heights
  if (length(seg_heights)<length(segments)){
    seg_heights_len=length(seg_heights)
    for (m in (seg_heights_len+1):length(segments)) {
      seg_heights[m]=seg_heights[seg_heights_len]
    }
  }
  #check and add y_heights
  if (length(y_heights)<(length(segments)+1)){
    y_heights_len=length(y_heights)
    for (m in (y_heights_len+1):(length(segments)+1)) {
      y_heights[m]=y_heights[y_heights_len]
    }
  }
  ### -------- plot -------- ###
  #get elements
  ##trans
  if (length(plot$scales$scales)==0){
    trans="identity"
  }else if ('trans' %in% names(plot$scales$scales[[1]])){
    trans=plot$scales$scales[[1]]$trans
  }else{
    trans="identity"
  }
  if ('reverse' %in% trans){
    if (ylim[1] < ylim[2]){
      msg=paste0('ylim: ','c(',ylim[1],',',ylim[2],')',' is wrong. It should be ','c(',ylim[2],',',ylim[1],')')
      stop(msg)
    }
  }
  if ('identity' %in% trans){
    if (ylim[1] > ylim[2]){
      msg=paste0('ylim: ','c(',ylim[1],',',ylim[2],')',' is wrong. It should be ','c(',ylim[2],',',ylim[1],')')
      stop(msg)
    }
  }
  #loop to plot 3 parts
  #the lowest, median and the toppest part by segments
  for (i in 1:length(segments)) {
    gap=unlist(segments[i])
    if (i==1){
      #plot the lowest part
      if (ylim[1] < ylim[2]){
        breaks=seq(ylim[1],gap[1],by=tick_width[i])
      }else if (ylim[1] > ylim[2]){
        breaks=seq(gap[1],ylim[1],by=tick_width[i])
      }
      p_segment.i<-plot+coord_cartesian(ylim=c(ylim[1],gap[1]))+
        theme(panel.border = element_blank())+
        theme(axis.line.y=element_line(),
              axis.line.x.bottom = element_line(),
              plot.title = element_blank(),
              legend.position = "none",
              strip.text.x = element_blank())+
        scale_y_continuous(expand = c(0,0),
                           trans = trans,
                           breaks = breaks)+
        ylab(label=NULL)
      p_segment=list(p_segment.i)
      names(p_segment)[length(p_segment)]=i
      rel_heigh=c(y_heights[i],seg_heights[i])
    }else{
      #plot the median part
      if (ylim[1] < ylim[2]){
        breaks=seq(ylim[1],gap[1],by=tick_width[i])
      }else if (ylim[1] > ylim[2]){
        breaks=seq(gap[1],ylim[1],by=tick_width[i])
      }
      p_segment.i<-plot+
        coord_cartesian(ylim=c(unlist(segments[i-1])[2],
                               gap[1]))+
        theme(panel.border = element_blank())+
        theme(axis.line.y=element_line(),
              #axis.line.x.bottom = element_line(),
              legend.position = "none",
              axis.text.x=element_blank(),
              axis.ticks.x =element_blank(),
              title = element_blank(),
              axis.title.x=element_blank(),
              strip.text.x = element_blank())+
        scale_y_continuous(expand = c(0,0),
                           breaks = breaks,
                           trans = trans)+
        ylab(label=NULL)
      #add y label in the middle median part
      p_segment=c(p_segment,list(NULL),list(p_segment.i))
      names(p_segment)[length(p_segment)]=i
      rel_heigh=c(rel_heigh,y_heights[i],seg_heights[i])
    }
    #plot the toppest part in the end
    if (i==length(segments)){
      if (ylim[1]<ylim[2]){
        breaks=seq(gap[2],ylim[2],by=tick_width[i+1])
      }else if (ylim[1]>ylim[2]){
        breaks=seq(ylim[2],gap[2],by=tick_width[i+1])
      }
      
      p_segment.i<-plot+
        coord_cartesian(ylim=c(gap[2],ylim[2]))+
        theme(panel.border = element_blank())+
        theme(axis.line.y=element_line(),
              axis.line.x.top = element_line(),
              #axis.line.x.bottom = element_line(),
              legend.position = "none",
              axis.text.x=element_blank(),
              axis.ticks.x =element_blank(),
              axis.title.x=element_blank())+
        scale_y_continuous(expand = c(0,0),
                           breaks = breaks,
                           trans = trans)+
        ylab(label=NULL)
      p_segment=c(p_segment,list(NULL),list(p_segment.i))
      names(p_segment)[length(p_segment)]=i+1
      rel_heigh=c(rel_heigh,y_heights[i])
    }
  }
  #reverse order
  p_segment=rev(p_segment)
  return(p_segment)
}
extract_mol_profiles <- function(D, poi){
  # rowids of proteins of interest
  poidx <- which(D %>% rowData()%>% data.frame()%>% pull(name) %in% poi)
  # df with proteins as columns and samples as rows
  X <- D %>% assay() %>% t() %>% data.frame() %>% .[, poidx]
  # column /protein names
  names(X) <- D %>% rowData()%>% data.frame()%>% pull(name) %>% .[poidx]
  # add annotation columns
  X <- bind_cols(X, D%>% colData() %>% data.frame() %>% 
                   select(Group, death)) %>% dplyr::rename(!!as.name(poi):='...1')
  return(X)
}
