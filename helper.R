# Functions for the analysis of C. elegans Wormtracker data set

installPackages <- function(packageList=c(), biocPackageList=c(), loadAll=F){
  if ( length(setdiff(packageList, rownames(installed.packages()))) > 0 ) {
    install.packages(setdiff(packageList, rownames(installed.packages())))
  }
  if( length(setdiff(biocPackageList, rownames(installed.packages())))>0 ){source("http://bioconductor.org/biocLite.R")}
  for(package in setdiff(biocPackageList, rownames(installed.packages())) ){
    biocLite(package)
  }
  if(loadAll){
    for(package in packageList){library(package, character.only = TRUE)}
    for(package in biocPackageList){library(package, character.only = TRUE)}
  }
  
}
installPackages(packageList = CRAN_packageList,biocPackageList = biocPackageList,loadAll = TRUE)


# Import experiment
import_file_directly <- function(path = data_path,experiment_name = "20181212-glp-atf5", file_name = "report.csv"){
  file <- paste0(path,"/",experiment_name,"/",file_name)
  print(paste("file: ", file))
  if(grepl(pattern = ".csv",x = file)) table <- read.csv(file = file,header = FALSE,stringsAsFactors = FALSE)
  if(grepl(pattern = ".xlsx",x = file)) table <- xlsx::read.xlsx(file = file,sheetIndex = 1,header = FALSE, stringsAsFactors = FALSE)
  return(table)
}


import_subtables <- function(table,tables_to_extract = "Well Activity"){
  start_idx <- grep(pattern = tables_to_extract,x = table[,1],useBytes = TRUE)
  add_rows_to_start <- grep(pattern = "<<<<<<",x = table[(start_idx + 1) :nrow(table),1],useBytes = TRUE)[1]
  end_idx <- start_idx + add_rows_to_start
  if(is.na(end_idx)) end_idx <- nrow(table)
  subtable <- table[start_idx : end_idx,]
  return(subtable)
}

import_quality_check <- function(table,max_na_ratio = 0.5){
  # Remove NA rows
  na_ratio <- apply(table, 1, function(x) sum(is.na(x) / length(x)))
  keep_row <- ifelse(na_ratio <=  max_na_ratio, TRUE, FALSE)
  table_qc <- table[keep_row,]
  return(table_qc)
}

import_tidy_table <- function(table){
  table <- as_tibble(table,rownames = NULL)
  if(colnames(table)[1] == "V1"){
    names <- unname(table[1,])
    colnames(table) <- as.character(names)
    table <- table[-1,]
  }
  return(table)
}


import_metainfo <- function(table, tables_to_extract = "Color"){
  start <- grep(pattern = tables_to_extract,x = table[,1])
  sub <- table[start : (start + 8),2:ncol(table)]
  colnames(sub) <- sub[1,]
  sub <- as_tibble(sub[-1,]) %>% 
    gather(key = "Col",value = "Content",2:13) %>%
    mutate(ID = paste0(Row,Col))
  sub <- sub[,c("ID","Row","Col","Content")]
  colnames(sub)[4] <- tables_to_extract
  return(sub)
}

import_combine <- function(import_raw,import_meta){
  tab_wellspergroup <- import_raw %>% 
    import_subtables(tables_to_extract = "Number of wells/group") %>%
    import_quality_check() %>%
    import_tidy_table()
  
  tab_samplesperinterval <- import_raw %>% 
    import_subtables(tables_to_extract = "Number of Samples per Acquisition Interval") %>%
    import_quality_check() %>%
    import_tidy_table()
  
  tab_wellactivity <- import_raw %>% 
    import_subtables(tables_to_extract = "Well Activity") %>%
    import_quality_check() %>%
    import_tidy_table()
  
  tab_channelactivity <- import_raw %>% 
    import_subtables(tables_to_extract = "Channel Activity") %>%
    import_quality_check() %>%
    import_tidy_table() %>%
    separate(Well, c("Well", "Well_sub")) %>%
    dplyr::select(-Well_sub)
  
  ############ Meta data ############ 
  Color = import_meta %>% import_metainfo(tables_to_extract = "Color") %>% dplyr::select(1,4) 
  Linetype = import_meta %>% import_metainfo(tables_to_extract = "Linetype") %>% dplyr::select(1,4) %>% mutate(Linetype = ifelse(is.na(Linetype),"No Linetype",Linetype))
  Stressor = import_meta %>% import_metainfo(tables_to_extract = "Stressor") %>% dplyr::select(1,4)
  meta_raw <- cbind(Color, Linetype, Stressor)
  ID_cols <- c(3,5,7)
  meta <- meta_raw[,-ID_cols] %>% as_tibble()
  
  
  ############ create master table ############ 
  dat <- full_join(x = meta, y = tab_wellactivity, by = c("ID" = "Well"))
  dat <- gather(dat, key =  "time", value = "activity", 7:ncol(dat))
  dat$activity <- as.numeric(as.character(dat$activity))
  dat$time <- as.numeric(as.character(dat$time))
  base_level <- unique(dat$Stressor)[grep(pattern = "M9",x = unique(dat$Stressor))]
  stress_levels <- unique(dat$Stressor)[grep(pattern = "M9",x = unique(dat$Stressor),invert = TRUE)]
  dat$Stressor <- factor(dat$Stressor,levels = c(base_level,stress_levels))
  dat <- wrangle_qc(dat)
  return(dat)
}

# Wrangle imported datasets
wrangle_qc <- function(data,required_cols = c("Color","Linetype","Stressor", "activity")){
  na_idx <- apply(data[,required_cols],1,function(x) any(is.na(x)))
  data <- data[!na_idx,]
  return(data)
}

wrangle_remove_wells <- function(data,remove_wells = NULL){
  remove_idx <- data$ID %in% remove_wells
  data <- data[!remove_idx,]
  return(data)
}

wrangle_summarize_well <- function(dat, activity_padding, time_padding){
  # add Statistics information for every WELL
  summary <- dat %>%
    dplyr::group_by(ID) %>%
    summarize(n = n(), 
              Activ_Mean_perWell = mean(activity,na.rm = TRUE),
              Activ_SD_perWell = sd(activity,na.rm = TRUE),
              Activ_STE_perWell = stat_ste(activity),
              HAT_perWell = stat_activitydecay(activity = activity,time = time, activity_padding = activity_padding, time_padding = time_padding)$time,
              HAA_perWell = stat_activitydecay(activity = activity,time = time, activity_padding = activity_padding, time_padding = time_padding)$activity,
              Min_perWell = min(activity,na.rm = TRUE),
              Q1_perWell = quantile(activity,na.rm = TRUE)["25%"],
              Q2_perWell = median(activity,na.rm = TRUE),
              Q3_perWell = quantile(activity,na.rm = TRUE)["75%"],
              Max_perWell = max(activity,na.rm = TRUE)
    ) %>%
    ungroup()

  return(summary)
}

wrangle_summarize_groups <- function(df){
  for_plotting <- df %>% 
    group_by(Color,Linetype,Stressor) %>%
    summarize(`Half-activity time AVG` = mean(HAT_perWell,na.rm = TRUE),
              `Half-activity activity AVG` = mean(HAA_perWell,na.rm = TRUE)
              )
  return(for_plotting)
}

wrangle_summarize_stressor <- function(df){
  summary_stressor <- df %>% 
    group_by(Color,Linetype) %>%
    summarize(stressor_pval = stat_wilcox(label = Stressor,value = activity)$pval,
              stressor_Low95Conf = stat_wilcox(label = Stressor,value = activity)$conf_low,
              stressor_High95Conf = stat_wilcox(label = Stressor,value = activity)$conf_high) %>%
    mutate(stressor_pval = format(stressor_pval, scientific = 5, digits = 3)) %>%
    mutate(Stressor = factor(levels(df$Stressor)[2],levels = levels(df$Stressor)))
  return(summary_stressor)
}

# 
# wrangle_summarize_control <- function(df,ctr_color, ctr_cond){
#   color_ctr_idx <- df$Color %in% ctr_color
#   color_cond_idx <- df$Linetype %in% ctr_cond
#   M9 <- levels(df$Stressor)[1]
#   Arsenite <- levels(df$Stressor)[2]
#   
#   # compute control values without error metric for the control Linetype
#   df_ctr <- df %>%
#     filter(color_ctr_idx & color_cond_idx)  %>%
#     group_by(Stressor) %>%
#     dplyr::select(Stressor,Activ_Mean_perWell,HAT_perWell,HAA_perWell) %>%
#     summarise(Activ_Mean_ctr = mean(Activ_Mean_perWell,na.rm = TRUE),
#            HAT_perWell_ctr = mean(HAT_perWell,na.rm = TRUE),
#            HAA_perWell_ctr = mean(HAA_perWell,na.rm = TRUE)
#            )
#   
#   # compute the per well differences to the control, once in M9 and once in Arsenite
#   summary_ctr_M9_well <- df %>% 
#     filter(Stressor == M9) %>%
#     dplyr::select(ID,Color,Linetype,Stressor,Activ_Mean_perWell,HAT_perWell,HAA_perWell,AUC_ID) %>%
#     ungroup() %>%
#     distinct() %>%
#     group_by(ID) %>%
#     mutate(Activ_reltoctr_well = (Activ_Mean_perWell - df_ctr[1,]$Activ_Mean_ctr)/df_ctr[1,]$Activ_Mean_ctr * 100,
#            HAT_reltoctr_well = (HAT_perWell - df_ctr[1,]$HAT_perWell_ctr)/df_ctr[1,]$HAT_perWell_ctr * 100,
#            HAA_reltoctr_well = (HAA_perWell - df_ctr[1,]$HAA_perWell_ctr)/df_ctr[1,]$HAA_perWell_ctr * 100) %>%
#     ungroup()
#   
#   summary_ctr_Ars_well <- df %>% 
#     filter(Stressor == Arsenite) %>%
#     dplyr::select(ID,Color,Linetype,Stressor,Activ_Mean_perWell,HAT_perWell,HAA_perWell,AUC_ID) %>%
#     ungroup() %>%
#     distinct() %>%
#     group_by(ID) %>%
#     mutate(Activ_reltoctr_well = (Activ_Mean_perWell - df_ctr[2,]$Activ_Mean_ctr)/df_ctr[2,]$Activ_Mean_ctr * 100,
#            HAT_reltoctr_well = (HAT_perWell - df_ctr[2,]$HAT_perWell_ctr)/df_ctr[2,]$HAT_perWell_ctr * 100,
#            HAA_reltoctr_well = (HAA_perWell - df_ctr[2,]$HAA_perWell_ctr)/df_ctr[2,]$HAA_perWell_ctr * 100
#     ) %>%
#     ungroup()
#   
#   summary_ctr_well <- bind_rows(summary_ctr_M9_well,summary_ctr_Ars_well)
#   
#   
#   #######  NEW SCRIPT 191007
#   summary_ctr_strain_cond <- summary_ctr_well # All control vs. experiment information is already ENCODED IN THE FACTOR LEVELS
#   #######  NEW SCRIPT 191007
#   
# 
#   # summaruze the well results to a single value with error metric
#   summary_ctr_well <- summary_ctr_well %>% mutate(CONTROL = Color == ctr_color & Linetype == ctr_cond)
#   
#   summary_ctrM9_group <- summary_ctr_well %>% filter(Stressor == M9)
#   summary_ctrM9_group <- summary_ctrM9_group %>%
#     group_by(Color,Linetype) %>%
#     summarize(Stressor = M9,
#               Activ_reltoctr_group_mean = mean(Activ_reltoctr_well,na.rm = TRUE),
#               Activ_reltoctr_group_sd = sd(Activ_reltoctr_well,na.rm = TRUE),
#               Activ_reltoctr_group_ste = stat_ste(Activ_reltoctr_well),
#               Activ_reltoctr_group_pval = stat_student_t_pval(df = summary_ctrM9_group,experiment_vals = Activ_reltoctr_well, variable_name = "Activ_reltoctr_well"),
#               HAT_reltoctr_group_mean = mean(HAT_reltoctr_well,na.rm = TRUE),
#               HAT_reltoctr_group_sd = sd(HAT_reltoctr_well,na.rm = TRUE),
#               HAT_reltoctr_group_ste = stat_ste(HAT_reltoctr_well),
#               HAT_reltoctr_group_pval = stat_student_t_pval(df = summary_ctrM9_group,experiment_vals = HAT_reltoctr_well, variable_name = "HAT_reltoctr_well"),
#               HAA_reltoctr_group_mean = mean(HAA_reltoctr_well,na.rm = TRUE),
#               HAA_reltoctr_group_sd = sd(HAA_reltoctr_well,na.rm = TRUE),
#               HAA_reltoctr_group_ste = stat_ste(HAA_reltoctr_well),
#               HAA_reltoctr_group_pval = stat_student_t_pval(df = summary_ctrM9_group,experiment_vals = HAA_reltoctr_well, variable_name = "HAA_reltoctr_well")
#     ) %>%
#     ungroup()
#   
#   summary_ctrArs_group <- summary_ctr_well %>% filter(Stressor == Arsenite)
#   summary_ctrArs_group <- summary_ctrArs_group %>%
#     group_by(Color,Linetype) %>%
#     summarize(Stressor = Arsenite,
#               Activ_reltoctr_group_mean = mean(Activ_reltoctr_well,na.rm = TRUE),
#               Activ_reltoctr_group_sd = sd(Activ_reltoctr_well,na.rm = TRUE),
#               Activ_reltoctr_group_ste = stat_ste(Activ_reltoctr_well),
#               Activ_reltoctr_group_pval = stat_student_t_pval(df = summary_ctrArs_group,experiment_vals = Activ_reltoctr_well, variable_name = "Activ_reltoctr_well"),
#               HAT_reltoctr_group_mean = mean(HAT_reltoctr_well,na.rm = TRUE),
#               HAT_reltoctr_group_sd = sd(HAT_reltoctr_well,na.rm = TRUE),
#               HAT_reltoctr_group_ste = stat_ste(HAT_reltoctr_well),
#               HAT_reltoctr_group_pval = stat_student_t_pval(df = summary_ctrArs_group,experiment_vals = HAT_reltoctr_well, variable_name = "HAT_reltoctr_well"),
#               HAA_reltoctr_group_mean = mean(HAA_reltoctr_well,na.rm = TRUE),
#               HAA_reltoctr_group_sd = sd(HAA_reltoctr_well,na.rm = TRUE),
#               HAA_reltoctr_group_ste = stat_ste(HAA_reltoctr_well),
#               HAA_reltoctr_group_pval = stat_student_t_pval(df = summary_ctrArs_group,experiment_vals = HAA_reltoctr_well, variable_name = "HAA_reltoctr_well")
#     ) %>%
#     ungroup()
#   
#   summary_ctr_group <- bind_rows(summary_ctrM9_group,summary_ctrArs_group)
#   
#   # reorder table
#   out <- summary_ctr_group %>% arrange(desc(Stressor),Linetype,Color)
#   row_bool <- out %>% dplyr::select(Color) %>% unlist() %in% ctr_color
#   top_idx <- which(row_bool)
#   rest_idx <- which(!row_bool)
#   out <- out[c(top_idx,rest_idx),]
# 
#   return(out)
# }

wrangle_grouped_df_to_base_level <- function(df, col_label, col_value){
  df <- df %>% filter(!duplicated(ID))
  col_value <- enquo(col_value)
  ctr_level <- df %>% pull(!!col_label) %>% levels() %>% .[[1]]
  idx_ctr <- df %>% pull(!!col_label) == ctr_level
  ctr_population <- df[idx_ctr,]
  out <- tibble(idx = levels(df %>% pull(!!col_label))[1],Pvalue_WC = NA, conf_low_WC = NA, conf_high_WC = NA,Mean = mean(ctr_population %>% pull(!!col_value), na.rm = T), 
                Mean_ste = stat_ste(ctr_population %>% pull(!!col_value)), Mean_diff_abs = c(NA), Mean_diff_rel = c(NA), Individual_AUCs = ctr_population %>% pull(!!col_value) %>% paste(.,collapse = ", "))
  for(col_labl_level in levels(df %>% pull(!!col_label))[-1]){
    add_row <- tibble(idx = col_labl_level,Pvalue_WC = NA, conf_low_WC = NA, conf_high_WC = NA, Mean_diff_abs = c(NA), Mean_diff_rel = c(NA), Individual_AUCs = "")
    idx_exp <- df %>% pull(!!col_label) == col_labl_level
    exp_pop <- df[idx_exp,] %>% pull(!!col_value)
    sub_comp <- df[idx_ctr | idx_exp,]
    add_row[["Pvalue_WC"]] <- stat_wilcox(label = sub_comp %>% pull(!!col_label),value = sub_comp %>% pull(!!col_value))$pval
    add_row[["conf_low_WC"]] <- stat_wilcox(label = sub_comp %>% pull(!!col_label),value = sub_comp %>% pull(!!col_value))$conf_low
    add_row[["conf_high_WC"]] <- stat_wilcox(label = sub_comp %>% pull(!!col_label),value = sub_comp %>% pull(!!col_value))$conf_high
    add_row[["Mean"]] <-  mean(df[idx_exp,] %>% pull(!!col_value) %>% unlist() %>% unname())
    add_row[["Mean_ste"]] <- stat_ste(df[idx_exp,] %>% pull(!!col_value) %>% unlist() %>% unname())
    add_row[["Mean_diff_abs"]] <- stat_change_two_levels(label = sub_comp %>% pull(!!col_label),value = sub_comp %>% pull(!!col_value))$chg_abs
    add_row[["Mean_diff_rel"]] <- stat_change_two_levels(label = sub_comp %>% pull(!!col_label),value = sub_comp %>% pull(!!col_value))$chg_rel
    add_row[["Individual_AUCs"]] <- df[idx_exp,] %>% pull(!!col_value) %>% paste(.,collapse = ", ")
    out <- bind_rows(out,add_row)
  }
  # original levels
  out <- out %>% mutate(idx = factor(idx,levels = levels(df %>% pull(!!col_label))))
  #rename
  column_name <- quo_name(col_label)
  out <- out %>% rename(!!column_name := idx)
  return(out)
}

wrangle_normalize_to_first_row <- function(df,colnames, colname_suffix = "_norm"){
  for (col in colnames) {
    first <- df[[col]][[1]]
    normalized <- df[[col]] / first
    df[[paste0(col,colname_suffix)]] <- normalized
  }
  return(df)
}

wrangle_normalize_to_first_value_in_single_row <- function(df,normalize_to_first_value_in_this_row = "Mean", columns_to_normalize = c("Mean", "Mean_ste"), colname_suffix = "_nrm"){
  normalize_by_constant <- df[1,normalize_to_first_value_in_this_row] %>% unlist() %>% unname()
  for (col in columns_to_normalize) {
    normalized <- df[[col]] / normalize_by_constant
    df[[paste0(col,colname_suffix,"_by_global",normalize_to_first_value_in_this_row)]] <- normalized
  }
  return(df)
}

wrangle_summarize_control <- function(df,col_label, ...){
  col_label <- enquo(col_label)
  groupCol <- enquos(...)
 tab_summary <- df %>% 
    group_by(!!!groupCol) %>%
    group_modify(~wrangle_grouped_df_to_base_level(df = .x,col_label = col_label, col_value = AUC_ID)) %>%
    wrangle_normalize_to_first_value_in_single_row(df = .,normalize_to_first_value_in_this_row = "Mean",columns_to_normalize = c("Mean","Mean_ste"),colname_suffix = "_norm") %>%
   ungroup() 
 tab_summary <- tab_summary %>% mutate(`Activity [a.u.]` = paste0(round(Mean_norm_by_globalMean,3)," ± ",round(Mean_ste_norm_by_globalMean,3)), 
                        `Increase [%]` = paste0(ifelse(Mean_diff_rel >= 0,"+ ","- "),round((Mean_diff_rel - 1) * 100,3)) %>% ifelse(. == "NANA",NA,.),
                        `P-value [Wilcox]` = round(Pvalue_WC,6),
                        `Significance` = case_when(is.na(Pvalue_WC) ~ NA_character_,Pvalue_WC > 0.05 ~ "n.s.", Pvalue_WC > 0.01 ~ "*", Pvalue_WC > 0.001 ~ "**", Pvalue_WC > 0.0001 ~ "***", TRUE ~ "***"),
                        `Individual AUC values [a.u. raw]` = Individual_AUCs
                        )
 fcts <- tab_summary %>% map_lgl(.,is.factor) %>% which()
 tab_summary <- tab_summary %>% select(fcts,`Activity [a.u.]`,`Increase [%]`,`P-value [Wilcox]`,`Significance`,`Individual AUC values [a.u. raw]`)
  return(list(tab_summary = tab_summary))
}




wrangle_summarize_pvalmatrix <- function(df){
  tab_AUC <- df %>% 
    filter(!duplicated(ID)) %>%
    select(ID,Color,Linetype,Stressor,AUC_ID) %>%
    arrange(Color, Linetype, Stressor) %>%
    unite(Color,Linetype,sep = ", ",col = "Label",remove = FALSE) %>%
    unite(Color,Linetype, Stressor, sep = ", ",col = "Label_Stressor") %>%
    mutate(Label = factor(Label,levels = unique(Label)),
           Label_Stressor = factor(Label_Stressor,levels = unique(Label_Stressor)))
  # No stressor 
  uniq_all <- levels(tab_AUC$Label) 
  out_outer <- list()
  for (uniq_outer in uniq_all) {
    sub_outer <- tab_AUC %>% filter(Label == uniq_outer)
    out_inner <- list()
    for (uniq_inner in uniq_all) {
      sub_inner <- tab_AUC %>% filter(Label == uniq_inner)
      merge <- bind_rows(sub_outer,sub_inner) %>% droplevels()
      out_inner[[uniq_inner]] <- stat_wilcox(label = merge$Label,value = merge$AUC_ID)$pval %>% unlist()
    }
    out_outer[[uniq_outer]] <- out_inner
  }
  out <- data.table::rbindlist(out_outer) %>% as.data.frame()
  rowlab <- uniq_all %>% enframe(value = "Condition") %>% select(Condition)
  pval_mat_vals <- bind_cols(rowlab,out)
  pval_mat_symb <- apply(out, 2, function(x){
    case_when(is.na(x) ~ "",
              x >= 0.05 ~ "n.s.",
              x >= 0.01 ~ "*",
              x >= 0.001 ~ "**",
              x >= 0.0001 ~ "***",
              TRUE ~ "****")
  }) %>% as.data.frame()
  pval_mat_symb <- bind_cols(rowlab,pval_mat_symb)
  
  # With stressor 
  uniq_all <- levels(tab_AUC$Label_Stressor) 
  out_outer <- list()
  for (uniq_outer in uniq_all) {
    sub_outer <- tab_AUC %>% filter(Label_Stressor == uniq_outer)
    out_inner <- list()
    for (uniq_inner in uniq_all) {
      sub_inner <- tab_AUC %>% filter(Label_Stressor == uniq_inner)
      merge <- bind_rows(sub_outer,sub_inner) %>% droplevels()
      out_inner[[uniq_inner]] <- stat_wilcox(label = merge$Label_Stressor,value = merge$AUC_ID)$pval %>% unlist()
    }
    out_outer[[uniq_outer]] <- out_inner
  }
  out <- data.table::rbindlist(out_outer) %>% as.data.frame()
  rowlab <- uniq_all %>% enframe(value = "Condition") %>% select(Condition)
  pval_mat_vals_stressor <- bind_cols(rowlab,out)
  pval_mat_symb_stressor <- apply(out, 2, function(x){
    case_when(is.na(x) ~ "",
              x >= 0.05 ~ "n.s.",
              x >= 0.01 ~ "*",
              x >= 0.001 ~ "**",
              x >= 0.0001 ~ "***",
              TRUE ~ "****")
  }) %>% as.data.frame()
  pval_mat_symb_stressor <- bind_cols(rowlab,pval_mat_symb_stressor)
  return(list(pval_mat_vals = pval_mat_vals, pval_mat_symb = pval_mat_symb, pval_mat_vals_stressor = pval_mat_vals_stressor, pval_mat_symb_stressor = pval_mat_symb_stressor))
}



wrangle_threshold_set_0_activity <- function(rollsum, activity,rollmean_activity_threshold){
  first_zero <- which(rollsum <= rollmean_activity_threshold)[1]
  if(is.na(first_zero)) return(activity)
  to_replace <- first_zero:length(rollsum)
  activity[to_replace] <- 0
  activity
}

wrangle_threshold_filter <- function(rollsum,rollmean_activity_threshold){
  bool <- rep_len(TRUE,length(rollsum))
  first_zero <- which(rollsum <= rollmean_activity_threshold)[1]
  if(is.na(first_zero)) return(bool)
  to_replace <- first_zero:length(rollsum)
  bool[to_replace] <- FALSE
  bool
}

wrangle_rescale <- function(dat,rescale_to_value,use_timepoint_to_scale){
  use_timepoint_to_scale <- as.numeric(as.character(use_timepoint_to_scale))
  rescale_to_value <- as.numeric(as.character(rescale_to_value))
  dat <- dat %>%
    group_by(ID) %>%
    arrange(time) %>%
    mutate(activity = activity/median(activity[time == use_timepoint_to_scale],na.rm = TRUE)) %>%
    mutate(activity = activity * rescale_to_value) %>%
    ungroup()
  return(dat)
}

wrangle_convert_time <- function(time, prev, change){
  if(prev == change) out <- time
  if(prev == "minutes" & change == "hours")    out <- time / 60
  if(prev == "minutes" & change == "days")     out <- time / 1440
  if(prev == "hours"   & change == "minutes")  out <- time * 60
  if(prev == "hours"   & change == "days")     out <- time / 24
  if(prev == "days"    & change == "minutes")  out <- time * 1440
  if(prev == "days"    & change == "hours")    out <- time * 24
  if(length(time) != length(out)) stop("time conversion is incorrect")
  return(out)
}


# statistical functions
stat_median_one_sigma <- function(y) {
  y = median(y)
  ymin = ifelse(median(y)-sd(y)/2 >= 0,median(y)-sd(y)/2,0) # show only positive intervals
  ymax = median(y)+sd(y)/2
  return(data.frame(y = y, ymin = ymin, ymax = ymax))
}
stat_one_sigma <- function(y) {
  ymin = ifelse(median(y)-sd(y)/2 >= 0,median(y)-sd(y)/2,0) # show only positive intervals
  ymax = median(y)+sd(y)/2
  return(data.frame(ymin = ymin, ymax = ymax))
}
stat_median <- function(y) {
  return(data.frame(y=median(y)))
}

stat_activitydecay <- function(activity,time,activity_padding,time_padding){
  data <- bind_cols(tibble(activity),tibble(time))
  interrange <- quantile(data$activity,c(0 + activity_padding,1 - activity_padding),na.rm = TRUE)
  halfactivity <- mean(interrange,na.rm = TRUE)
  data_interrange <- data %>% 
    filter(between(data$activity, halfactivity - time_padding,halfactivity + time_padding))
  halftime <- median(data_interrange$time,na.rm = TRUE)
  return(list(activity = halfactivity, time = halftime))
}

stat_wilcox <- function(label, value){
  # Wilcox test
  # a seriuous limitation is that since we are looking at time series data that there is likely serial dependence in the differences.
  # The Wilcoxon signed rank test is an alternative to the t-test that does not assume the sample to follow a normal distribution (non-parametric).
  # Ref 1: David F. Bauer (1972). Constructing confidence sets using rank statistics. Journal of the American Statistical Association 67, 687–690. doi: 10.1080/01621459.1972.10481279.
  # Ref 2: Myles Hollander and Douglas A. Wolfe (1973). Nonparametric Statistical Methods. New York: John Wiley & Sons. Pages 27–33 (one-sample), 68–75 (two-sample).
  # Or second edition (1999).
  label <- label %>% unlist() %>% droplevels()
  if(length(unique(label)) != 2){
    print("Wilcox test can only compare two distributions (not more or less)")
    return(list(pval = NA, conf_low = NA, conf_high = NA))
  }
  if(!is.factor(label)) stop("label needs to be a factor.")
  a <- value[label == levels(label)[1]] %>% unlist()
  b <- value[label == levels(label)[2]] %>% unlist()
  wilcox <- wilcox.test(a,b,conf.int = TRUE)
  pval <- wilcox$p.value
  conf <- wilcox$conf.int
  attributes(conf) <- NULL
  return(list(pval = pval, conf_low = conf[1], conf_high = conf[2]))
}

stat_student_t_pval <- function(df,variable_name,experiment_vals){
  control <- df[df$CONTROL,variable_name] %>% unlist() %>% unname()
  ttest <- t.test(x = control,y = experiment_vals)
  return(ttest$p.value)
}

stat_ste <- function(vector){
  ste <- sd(vector,na.rm = TRUE)/sqrt(length(vector[!is.na(vector)]))
  return(ste)
}



stat_change_two_levels <- function(label, value){
  label <- label %>% unlist() %>% droplevels()
  if(length(unique(label)) != 2){
    print("This test can only compare two distributions (not more or less)")
    return(list(pval = NA, conf_low = NA, conf_high = NA))
  }
  if(!is.factor(label)) stop("label needs to be a factor.")
  mean_first <- value[label == levels(label)[1]] %>% unlist() %>% mean(.,na.rm = T)
  mean_second <- value[label == levels(label)[2]] %>% unlist() %>% mean(.,na.rm = T)
  chg_abs <- mean_second - mean_first
  chg_rel <- mean_second / mean_first
  return(list(chg_abs = chg_abs, chg_rel = chg_rel))
}
# stat_boltzmann_sigmoid <- function(x,y,activity_padding){
#   # initial values
#   top = quantile((unlist(y)),probs = 1 - activity_padding/2) %>% unname()
#   bottom = 0
#   halftime = stat_activitydecay(activity = unlist(y),time = unlist(x),activity_padding = activity_padding,time_padding = time_padding)
#   V50 = halftime$time
#   slope = 50
#   # optimize parameters if possible
#   params = tryCatch({
#     ZABfit <- nls(y ~ bolz_sig(x,top,bottom,V50,slope),start=list(top = top,bottom = bottom,V50 = V50,slope = slope))
#     coefs <- summary(ZABfit)$coefficients
#     list(top = coefs["top","Estimate"], bottom = coefs["bottom","Estimate"], V50 = coefs["V50","Estimate"], slope = coefs["slope","Estimate"])
#   },
#   error = function(c) list(top = top, bottom = bottom, V50 = V50, slope = slope),
#   warning = function(c) "warning",
#   message = function(c) "message")
#   # return params
#   y_pred <- bolz_sig(x,top = params$top,bottom = params$bottom,V50 = params$V50,slope = params$slope)
#   return(y_pred)
# }


stat_bolz_sig <- function(time,top,bottom,V50,slope){
  pred = top + ((bottom-top) / (1 + exp((V50 - time)/slope)))
  return(pred)
}

wrangle_bolz_sig <- function(activity,
                             time,
                             activity_padding,
                             time_padding,
                             top_free = FALSE){
  # best estimates
  print("start")
  top_free <- as.logical(top_free)
  top_est = quantile(activity,1- (activity_padding/2))
  bottom_est = ifelse(quantile(activity,(activity_padding/2)) > 0,quantile(activity,(activity_padding/2)),0)
  V50_est = stat_activitydecay(activity = activity,time = time, activity_padding = activity_padding, time_padding = time_padding)$time
  slope_est = V50_est/top_est
  
  # skip NA estimates
  print(paste("top",top_est,"bottom",bottom_est,"V50",V50_est,"slope_est",slope_est))
  is_na <- sapply(c(top_est,bottom_est,V50_est,slope_est), is.na)
  no_na <- !any(is_na)
  
  # can coefs be optimized?
  if (!top_free) feasible = tryCatch({
    nls(activity ~ stat_bolz_sig(time,top,bottom,V50,slope),start=list(top = top_est,bottom = bottom_est,V50 = V50_est,slope = slope_est))
    "converges"
    },error = function(c) "error")
  if (top_free) feasible = tryCatch({
    nls(activity ~ stat_bolz_sig(time,top_est,bottom,V50,slope),start=list(bottom = bottom_est,V50 = V50_est,slope = slope_est))
    "converges"
  },error = function(c) print("error"))
  
  # optimize - fully
  if(feasible == "converges" & no_na & !top_free){
    print("top not fixed")
    model <- nls(activity ~ stat_bolz_sig(time,top,bottom,V50,slope),start=list(top = top_est,bottom = bottom_est,V50 = V50_est,slope = slope_est))
    coefs <- summary(model)$coefficients
    comp <- list(top = coefs["top","Estimate"],bottom = coefs["bottom","Estimate"],V50 = coefs["V50","Estimate"],slope = coefs["slope","Estimate"])
  }
  if(feasible == "converges" & no_na & top_free){
    print("top fixed")
    model <- nls(activity ~ stat_bolz_sig(time,top_est,bottom,V50,slope),start=list(bottom = bottom_est, V50 = V50_est, slope = slope_est))
    coefs <- summary(model)$coefficients
    comp <- list(top = top_est, bottom = coefs["bottom","Estimate"],V50 = coefs["V50","Estimate"],slope = coefs["slope","Estimate"])
  }
  # don't optimize
  if(feasible == "error" & no_na){
    print("entered error branch")
    # model <- nls(activity ~ stat_bolz_sig(time,top_est,bottom_est,V50_est,slope),start=list(slope = slope_est))
    # coefs <- summary(model)$coefficients
    # comp <- list(top = top_est, bottom = bottom_est, V50 = V50_est, slope = slope_est)
    # comp <- list(top = 100, bottom = 0, V50 = V50_est, slope = slope_est)
  } 
  # apply model if possible
  if (!no_na | feasible == "error"){
    print("add NA")
    y_pred <- rep_len(NA,length.out = length(activity))
  } else {
    print("add model info")
    y_pred <- stat_bolz_sig(time,top = comp$top, bottom = comp$bottom, V50 = comp$V50,slope = comp$slope)
  }
  
  return(y_pred)
}


# plot functions
plot_modif <- function(ggobj,p_str,size_scale = 1,style_legend){
  adj_size <- map_dbl(p_str$size,.f = ~ .x * size_scale)  
  if(p_str$modify_axis) ggobj <- ggobj + coord_cartesian(ylim = p_str$y_coord_range, xlim = p_str$x_coord_range)
  if(p_str$facet_info != ". ~ .") ggobj <- ggobj + facet_grid(p_str$facet_info)
  ggobj <- ggobj + scale_size_manual(values = adj_size,guide = "none")
  ggobj <- ggobj + scale_alpha_manual(values = c(0.1,1),guide = "none")
  ggobj <- ggobj + scale_color_manual(values = p_str$pal)
  ggobj <- ggobj + scale_linetype_manual(values = p_str$lty)
  ggobj <- ggobj + style_legend
  ggobj
}

plot_shiny_theme <- theme(
  panel.background = element_blank(),
  plot.background = element_blank(),
  legend.key = element_blank(),
  legend.background = element_blank(),
  legend.box.background = element_blank()
)

# export functions
export_plot <- function(plot,p_str,appearance_selected){
  filename <- paste0(names(plot),".",p_str$plot.format)
  if(p_str$filename_timestamp) filename <- paste0(Sys.Date(),"_",filename)
  # bg <- ifelse(appearance_selected == "appearance_pres_trans" & p_str$plot.format %in% c("png","svg"), "transparent",NA)
  bg <- "transparent"
  plot <- plot[[1]]
  ggsave(filename = filename,
         plot = plot,
         device = p_str$plot.format,
         dpi = p_str$plot.dpi, 
         width = p_str$plot.width,
         height = p_str$plot.height,
         bg = bg,
         units = "cm")
  return(filename)
}


export_tables <- function(ctrvsexp,processed_movement,raw_meta,file,filename_timestamp){
  wb <- createWorkbook()
  addWorksheet(wb, "Arsenite for each linetype")
  addWorksheet(wb, "Arsenite for each color")
  addWorksheet(wb, "M9 for each linetype" )
  addWorksheet(wb, "M9 for each color" )
  addWorksheet(wb, "All for each linetype")
  addWorksheet(wb, "All for each color")
  addWorksheet(wb, "Processed movement")
  addWorksheet(wb, "Meta information")
  
  addWorksheet(wb, "pval_values")
  addWorksheet(wb, "pval_symbol")
  addWorksheet(wb, "pval_stressor_values")
  addWorksheet(wb, "pval_stressor_symbol")

    
  # Define general styling
  italics <- createStyle(textDecoration = "italic")
  bold <- createStyle(textDecoration = "bold")
  headerStyle <- createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#4F81BD", border="TopBottom", borderColour = "#4F81BD")
  
  # Fill divisions & format headers
  ncol_style <- 8
  addStyle(wb, sheet = "Arsenite for each linetype", createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#9c0505", border="TopBottom", borderColour = "#4F81BD"), rows = 1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "Arsenite for each color", createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#9c0505", border="TopBottom", borderColour = "#4F81BD"), rows = 1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "M9 for each linetype" , createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#9c0505", border="TopBottom", borderColour = "#4F81BD"), rows = 1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "M9 for each color" , createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#9c0505", border="TopBottom", borderColour = "#4F81BD"), rows = 1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "All for each linetype", createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#9c0505", border="TopBottom", borderColour = "#4F81BD"), rows = 1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "All for each color", createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#9c0505", border="TopBottom", borderColour = "#4F81BD"), rows = 1, cols = 1:ncol_style, gridExpand = TRUE)
 
  addStyle(wb, sheet = "pval_values", createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center", valign = "center", wrapText = TRUE, fgFill = "#9c0505", border="TopBottom", borderColour = "#4F81BD"), rows = 1, cols = 1:ncol(ctrvsexp$pval_values), gridExpand = TRUE)
  addStyle(wb, sheet = "pval_values", createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center", valign = "center", fgFill = "#f7f2f2", border="TopBottom", borderColour = "#4F81BD"), rows = 1:nrow(ctrvsexp$pval_values) + 1, cols = 1, gridExpand = TRUE)
  addStyle(wb, sheet = "pval_values", createStyle(halign = "center", valign = "center"), rows = 1:nrow(ctrvsexp$pval_values) + 1, cols = 1:ncol(ctrvsexp$pval_values), gridExpand = TRUE)
  
   # addStyle(wb, sheet = "pval_symbol", createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#9c0505", border="TopBottom", borderColour = "#4F81BD"), rows = 1, cols = 1:ncol_style, gridExpand = TRUE)
  # addStyle(wb, sheet = "pval_stressor_values", createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#9c0505", border="TopBottom", borderColour = "#4F81BD"), rows = 1, cols = 1:ncol_style, gridExpand = TRUE)
  # addStyle(wb, sheet = "pval_stressor_symbol", createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#9c0505", border="TopBottom", borderColour = "#4F81BD"), rows = 1, cols = 1:ncol_style, gridExpand = TRUE)
  
  writeData(wb, "Arsenite for each linetype", ctrvsexp$tab_summary_lty_Ars)
  writeData(wb, "Arsenite for each color", ctrvsexp$tab_summary_col_Ars)
  writeData(wb, "M9 for each linetype" , ctrvsexp$tab_summary_lty_M9)
  writeData(wb, "M9 for each color" , ctrvsexp$tab_summary_col_M9)
  writeData(wb, "All for each linetype", ctrvsexp$tab_summary_lty_all)
  writeData(wb, "All for each color", ctrvsexp$tab_summary_col_all)
  addStyle(wb, sheet = "Processed movement", createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#8c95a3", border="TopBottom", borderColour = "#4F81BD"), rows = 1, cols = 1:ncol(processed_movement), gridExpand = TRUE)
  writeData(wb, "Processed movement", processed_movement)
  # addStyle(wb, sheet = "Meta information", createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",fgFill = "#7f4646", border="TopBottom", borderColour = "#4F81BD"), rows = 1, cols = 1:ncol(raw_meta), gridExpand = TRUE)
  writeData(wb, "Meta information", raw_meta)
  writeData(wb, "pval_values", ctrvsexp$pval_values)
  writeData(wb, "pval_symbol", ctrvsexp$pval_symbol)
  writeData(wb, "pval_stressor_values", ctrvsexp$pval_stressor_values)
  writeData(wb, "pval_stressor_symbol", ctrvsexp$pval_stressor_symbol)
  
  
  # set multi-sheet formatting:
  col_width <- 16
  setColWidths(wb, "Arsenite for each linetype", cols = 1:ncol_style, widths = c(rep(col_width,ncol_style-1),100))
  setColWidths(wb, "Arsenite for each color", cols = 1:ncol_style, widths = c(rep(col_width,ncol_style-1),100))
  setColWidths(wb, "M9 for each linetype" , cols = 1:ncol_style, widths = c(rep(col_width,ncol_style-1),100))
  setColWidths(wb, "M9 for each color" , cols = 1:ncol_style, widths = c(rep(col_width,ncol_style-1),100))
  setColWidths(wb, "All for each linetype", cols = 1:ncol_style, widths = c(rep(col_width,ncol_style-1),100))
  setColWidths(wb, "All for each color", cols = 1:ncol_style, widths = c(rep(col_width,ncol_style-1),100))
  setColWidths(wb, "Processed movement", cols = 1:ncol(processed_movement), widths = c(14,14,14,14,14,rep(4,ncol(processed_movement)-5))) # column widths
  addStyle(wb, "Processed movement", italics, rows = 1:nrow(processed_movement)+1, cols = c(2,4), gridExpand = TRUE)
  
  setRowHeights(wb, "pval_values", rows = 0:nrow(ctrvsexp$pval_values)+1, heights = c(30,rep(60,nrow(ctrvsexp$pval_values)))) # column widths
  setColWidths(wb, "pval_values", cols = 1:ncol(ctrvsexp$pval_values), widths = c(25,rep(12,ncol(ctrvsexp$pval_values)-1))) # column widths

  # Coloring
  colors <- c("#f7f2f2", "#bdaaaa")
  wrangle_color_grouped_tables_by_firstfactors <- function(df, exclude_n_last_factors = 1){
    fcts <- df %>% map_lgl(.,is.factor) %>% which()
    col_cols <- fcts[1:(length(fcts)-exclude_n_last_factors)]
    idx_col <- df %>% select(col_cols) %>% unite(col = idx) 
    x <- !duplicated(idx_col)
    idx <- x %>% cumsum() %% 2 == 1
    return(idx)
  }
  addStyle(wb, sheet = "Arsenite for each linetype", createStyle(fgFill = colors[[1]]), rows = which(wrangle_color_grouped_tables_by_firstfactors(ctrvsexp$tab_summary_lty_Ars))+1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "Arsenite for each linetype", createStyle(fgFill = colors[[2]]), rows = which(!wrangle_color_grouped_tables_by_firstfactors(ctrvsexp$tab_summary_lty_Ars))+1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "Arsenite for each color", createStyle(fgFill = colors[[1]]), rows = which(wrangle_color_grouped_tables_by_firstfactors(ctrvsexp$tab_summary_col_Ars))+1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "Arsenite for each color", createStyle(fgFill = colors[[2]]), rows = which(!wrangle_color_grouped_tables_by_firstfactors(ctrvsexp$tab_summary_col_Ars))+1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "All for each linetype", createStyle(fgFill = colors[[1]]), rows = which(wrangle_color_grouped_tables_by_firstfactors(ctrvsexp$tab_summary_lty_all))+1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "All for each linetype", createStyle(fgFill = colors[[2]]), rows = which(!wrangle_color_grouped_tables_by_firstfactors(ctrvsexp$tab_summary_lty_all))+1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "All for each color", createStyle(fgFill = colors[[1]]), rows = which(wrangle_color_grouped_tables_by_firstfactors(ctrvsexp$tab_summary_col_all))+1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "All for each color", createStyle(fgFill = colors[[2]]), rows = which(!wrangle_color_grouped_tables_by_firstfactors(ctrvsexp$tab_summary_col_all))+1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "M9 for each linetype", createStyle(fgFill = colors[[1]]), rows = which(wrangle_color_grouped_tables_by_firstfactors(ctrvsexp$tab_summary_lty_M9))+1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "M9 for each linetype", createStyle(fgFill = colors[[2]]), rows = which(!wrangle_color_grouped_tables_by_firstfactors(ctrvsexp$tab_summary_lty_M9))+1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "M9 for each color", createStyle(fgFill = colors[[1]]), rows = which(wrangle_color_grouped_tables_by_firstfactors(ctrvsexp$tab_summary_col_M9))+1, cols = 1:ncol_style, gridExpand = TRUE)
  addStyle(wb, sheet = "M9 for each color", createStyle(fgFill = colors[[2]]), rows = which(!wrangle_color_grouped_tables_by_firstfactors(ctrvsexp$tab_summary_col_M9))+1, cols = 1:ncol_style, gridExpand = TRUE)
  
  addStyle(wb, sheet = "Processed movement", createStyle(fgFill = colors[[2]]), rows = 2:(nrow(processed_movement)+1), cols = 1:5, gridExpand = TRUE)
  addStyle(wb, sheet = "Processed movement", createStyle(fgFill = colors[[1]]), rows = 2:(nrow(processed_movement)+1), cols = 6:ncol(processed_movement), gridExpand = TRUE)
  
  saveWorkbook(wb, file = file, TRUE)
  }
  
import_file <- function(path){
  if (is.null(path)) return(NULL)
  if(grepl(pattern = ".csv",x = path$datapath)) table <- read.csv(file = path$datapath,header = FALSE,stringsAsFactors = FALSE)
  if(grepl(pattern = ".xlsx",x = path$datapath)) table <- xlsx::read.xlsx(file = path$datapath,sheetIndex = 1,header = FALSE, stringsAsFactors = FALSE)
  return(table)
}

wrangle_legendtitles <- function(meta){
  labels <- meta %>% dplyr::pull(1)
  color_lab <- labels[grep(pattern = "Color",x = labels)]
  color_lab <- str_split(color_lab,pattern = ":") %>% unlist() %>% .[2] %>% str_trim()
  linetype_lab <- labels[grep(pattern = "Linetype",x = labels)]
  linetype_lab <- str_split(linetype_lab,pattern = ":") %>% unlist() %>% .[2] %>% str_trim()
  return(list(color = color_lab, linetype = linetype_lab))
}

# Bibliography & References
wrangle_full_bib <- function(packages){
  bib <- map(packages,citation)
  write_bib(x = packages,file = "ref.bib")
  df <- bib2df("ref.bib", separate_names = FALSE)
  df <- df %>% dplyr::select(AUTHOR, TITLE, YEAR, URL, NOTE)
  Author = df %>% dplyr::select(AUTHOR) %>% pmap(.f = function(AUTHOR) {
    paste(AUTHOR, collapse = ", ")
  }) %>% unlist() %>% enframe("idx","Author") %>% dplyr::select(Author)
  df <- cbind(Author, df) %>% dplyr::select(-AUTHOR) %>% 
    rename(Title = TITLE, Year = YEAR, Weblink = URL, `Package version` = NOTE) %>% arrange(Author)
  return(df)
}

### Theme publication
theme_publ <- function(base_size=14, base_family="Helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.height= unit(0.57, "cm"), #spacing
            legend.key.width = unit(0.8,"cm"),
            legend.margin = unit(0.8, "cm"),
            legend.title = element_text(face="bold"),
            legend.text = element_text(size = rel(1),face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}


theme_pres_white <- function(base_size=17, base_family="Helvetica") {
  fill_bg <- "white"
  col_label <- "black"
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",size = rel(1.5), hjust = 0.5,colour = col_label),
           text = element_text(),
           panel.background = element_rect(colour = NA,fill = fill_bg),
           plot.background = element_rect(colour = NA,fill = fill_bg),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1.3),colour = col_label),
           axis.title.y = element_text(angle=90,vjust =2,colour = col_label),
           axis.title.x = element_text(vjust = -0.2,colour = col_label),
           axis.text = element_text(size = rel(1.2),colour = col_label), 
           axis.text.x = element_text(angle = 45, hjust = 1,colour = col_label),
           axis.line = element_line(colour = col_label),
           axis.ticks = element_line(colour = col_label),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA, fill = fill_bg),
           legend.position = "right",
           legend.direction = "vertical",
           legend.key.height= unit(0.6, "cm"), #spacing
           legend.key.width = unit(1.4,"cm"),
           legend.margin = unit(0.8, "cm"),
           legend.background = element_rect(fill=fill_bg, colour = NA),
           legend.title = element_text(face="bold",size = rel(1.2),colour = col_label),
           legend.text = element_text(size = rel(1.1),face="italic",colour = col_label),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold",colour = col_label)
   ))
}


theme_pres_black <- function(base_size=17, base_family="Helvetica") {
  fill_bg <- "black"
  col_label <- "white"
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",size = rel(1.5), hjust = 0.5,colour = col_label),
           text = element_text(),
           panel.background = element_rect(colour = NA,fill = fill_bg),
           plot.background = element_rect(colour = NA,fill = fill_bg),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1.3),colour = col_label),
           axis.title.y = element_text(angle=90,vjust =2,colour = col_label),
           axis.title.x = element_text(vjust = -0.2,colour = col_label),
           axis.text = element_text(size = rel(1.2),colour = col_label), 
           axis.text.x = element_text(angle = 45, hjust = 1,colour = col_label),
           axis.line = element_line(colour = col_label),
           axis.ticks = element_line(colour = col_label),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA, fill = fill_bg),
           legend.position = "right",
           legend.direction = "vertical",
           legend.key.height= unit(0.6, "cm"), #spacing
           legend.key.width = unit(1.4,"cm"),
           legend.margin = unit(0.8, "cm"),
           legend.background = element_rect(fill=fill_bg, colour = NA),
           legend.title = element_text(face="bold",size = rel(1.2),colour = col_label),
           legend.text = element_text(size = rel(1.1),face="italic",colour = col_label),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold",colour = col_label)
   ))
}



theme_pres_trans <- function(base_size=17, base_family="Helvetica") {
  fill_bg <- "transparent"
  col_label <- "black"
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",size = rel(1.5), hjust = 0.5,colour = col_label),
           text = element_text(),
           panel.background = element_rect(colour = NA,fill = fill_bg),
           plot.background = element_rect(colour = NA,fill = fill_bg),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1.3),colour = col_label),
           axis.title.y = element_text(angle=90,vjust =2,colour = col_label),
           axis.title.x = element_text(vjust = -0.2,colour = col_label),
           axis.text = element_text(size = rel(1.2),colour = col_label), 
           axis.text.x = element_text(angle = 45, hjust = 1,colour = col_label),
           axis.line = element_line(colour = col_label),
           axis.ticks = element_line(colour = col_label),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA, fill = fill_bg),
           legend.position = "right",
           legend.direction = "vertical",
           legend.key.height= unit(0.6, "cm"), #spacing
           legend.key.width = unit(1.4,"cm"),
           legend.margin = unit(0.8, "cm"),
           legend.background = element_rect(fill=fill_bg, colour = NA),
           legend.title = element_text(face="bold",size = rel(1.2),colour = col_label),
           legend.text = element_text(size = rel(1.1),face="italic",colour = col_label),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold",colour = col_label)
   ))
}

# theme_pres <- function(base_size=17, base_family="Helvetica") {
#   (theme_foundation(base_size=base_size, base_family=base_family)
#    + theme(plot.title = element_text(face = "bold",size = rel(1.5), hjust = 0.5),
#            text = element_text(),
#            panel.background = element_rect(colour = NA),
#            plot.background = element_rect(colour = NA),
#            panel.border = element_rect(colour = NA),
#            axis.title = element_text(face = "bold",size = rel(1.3)),
#            axis.title.y = element_text(angle=90,vjust =2),
#            axis.title.x = element_text(vjust = -0.2),
#            axis.text = element_text(size = rel(1.2)), 
#            axis.text.x = element_text(angle = 45, hjust = 1),
#            axis.line = element_line(colour="black"),
#            axis.ticks = element_line(),
#            panel.grid.major = element_line(colour="#f0f0f0"),
#            panel.grid.minor = element_blank(),
#            legend.key = element_rect(colour = NA),
#            legend.position = "right",
#            legend.direction = "vertical",
#            legend.key.height= unit(0.6, "cm"), #spacing
#            legend.key.width = unit(1.4,"cm"),
#            legend.margin = unit(0.8, "cm"),
#            legend.title = element_text(face="bold",size = rel(1.2)),
#            legend.text = element_text(size = rel(1.1),face="italic"),
#            plot.margin=unit(c(10,5,5,5),"mm"),
#            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#            strip.text = element_text(face="bold")
#    ))
# }


themes <- list("appearance_publ" = theme_publ,
               "appearance_pres_black" = theme_pres_black,
               "appearance_pres_white" = theme_pres_white,
               "appearance_pres_trans" = theme_pres_trans
               )

set_legend_style <- function(keywidth, legend_titles){
  legend_publ <- guides(colour = guide_legend(title = legend_titles[["color"]],override.aes = list(size=1.7, fill = NA), keywidth = keywidth, reverse = TRUE,order = 1), #fill = NA to remove grey legend bg from geom_smooth confidence intervals
                        linetype = guide_legend(title = legend_titles[["linetype"]], override.aes = list(size=1.7, fill = NA, color = "black"), keywidth = keywidth, reverse = TRUE,order = 0)  #fill = NA to remove grey legend bg from geom_smooth confidence intervals
  )
  legend_pres_black <- guides(colour = guide_legend(title = legend_titles[["color"]],override.aes = list(size=1.7, fill = NA), keywidth = keywidth, reverse = TRUE,order = 1), #fill = NA to remove grey legend bg from geom_smooth confidence intervals
                              linetype = guide_legend(title = legend_titles[["linetype"]], override.aes = list(size=1.7, fill = NA, color = "white"), keywidth = keywidth, reverse = TRUE,order = 0)  #fill = NA to remove grey legend bg from geom_smooth confidence intervals
  )
  legend_pres_white <- guides(colour = guide_legend(title = legend_titles[["color"]],override.aes = list(size=1.7, fill = NA), keywidth = keywidth, reverse = TRUE,order = 1), #fill = NA to remove grey legend bg from geom_smooth confidence intervals
                              linetype = guide_legend(title = legend_titles[["linetype"]], override.aes = list(size=1.7, fill = NA, color = "black"), keywidth = keywidth, reverse = TRUE,order = 0)  #fill = NA to remove grey legend bg from geom_smooth confidence intervals
  )
  legend_pres_trans <- guides(colour = guide_legend(title = legend_titles[["color"]], override.aes = list(size=1.7, fill = NA), keywidth = keywidth, reverse = TRUE,order = 1), #fill = NA to remove grey legend bg from geom_smooth confidence intervals
                              linetype = guide_legend(title = legend_titles[["linetype"]], override.aes = list(size=1.7, fill = NA, color = "black"), keywidth = keywidth, reverse = TRUE,order = 0)  #fill = NA to remove grey legend bg from geom_smooth confidence intervals
  )
  
  style_legends <- list("appearance_publ" = legend_publ,
                        "appearance_pres_white" = legend_pres_white,
                        "appearance_pres_black" = legend_pres_black,
                        "appearance_pres_trans" = legend_pres_trans
  )
  return(style_legends)
}





##### Styling button (from Dean Attali)
# Set up a button to have an animated loading indicator and a checkmark
# for better user experience
# Need to use with the corresponding `withBusyIndicator` server function
withBusyIndicatorUI <- function(button) {
  id <- button[['attribs']][['id']]
  div(
    `data-for-btn` = id,
    button,
    div(class = "btn-loading-container",
      hidden(
        # img(src = "ajax-loader-bar.gif", class = "btn-loading-indicator"),
        img(src = "circle_small.gif", class = "btn-loading-indicator"),
        icon("check", class = "btn-done-indicator")
      )
    ),
    hidden(
      div(class = "btn-err",
          div(icon("exclamation-circle"),
              tags$b("Error: "),
              span(class = "btn-err-msg")
          )
      )
    )
  )
}

# Call this function from the server with the button id that is clicked and the
# expression to run when the button is clicked
withBusyIndicatorServer <- function(buttonId, expr) {
  # UX stuff: show the "busy" message, hide the other messages, disable the button
  loadingEl <- sprintf("[data-for-btn=%s] .btn-loading-indicator", buttonId)
  doneEl <- sprintf("[data-for-btn=%s] .btn-done-indicator", buttonId)
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  shinyjs::disable(buttonId)
  shinyjs::show(selector = loadingEl)
  shinyjs::hide(selector = doneEl)
  shinyjs::hide(selector = errEl)
  on.exit({
    shinyjs::enable(buttonId)
    shinyjs::hide(selector = loadingEl)
  })
  
  # Try to run the code when the button is clicked and show an error message if
  # an error occurs or a success message if it completes
  tryCatch({
    value <- expr
    shinyjs::show(selector = doneEl)
    shinyjs::delay(2000, shinyjs::hide(selector = doneEl, anim = TRUE, animType = "fade",
                                       time = 0.5))
    value
  }, error = function(err) { errorFunc(err, buttonId) })
}

# When an error happens after a button click, show the error
errorFunc <- function(err, buttonId) {
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  errElMsg <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)
  errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)
  shinyjs::html(html = errMessage, selector = errElMsg)
  shinyjs::show(selector = errEl, anim = TRUE, animType = "fade")
}

appCSS <- "
.btn-loading-container {
margin-left: 18px;
color: green;
font-size: 1.2em;
}
.btn-done-indicator {
color: green;
}
.btn-err {
margin-top: 10px;
color: red;
}
"

AppJS <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
"

