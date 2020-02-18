# Functions to obtain ROCs and AUCs from gene scores and gene types
# Both at cohort-level and sample-level


# Top-level function
get_gene_rocs = function(scores_df, level = c("cohort", "sample"), TP_FP_comparisons, verbose = F){
  
  # Required packages
  require(tidyr)
  require(dplyr)
  
  
  # Prepare input data frame with scores and ranks
  if (verbose) cat("Preparing input\n")
  scores_df = .prepare_input(scores_df, TP_FP_comparisons, verbose = verbose)
  
  
  # Either run at cohort or sample level
  level = match.arg(level)
  if (level == "cohort"){
    if (verbose) cat("Getting ROCs at cohort-level\n")
    roc = .get_cohortLevel_rocs(scores_df, TP_FP_comparisons, verbose = verbose)
  } else if (level == "sample"){
    if (verbose) cat("Getting ROCs at sample-level\n")
    roc = .get_sampleLevel_rocs(scores_df, TP_FP_comparisons, verbose = verbose)
  }
  
  
  # Output
  if (verbose) cat("Finished\n")
  return(roc)
}



# Prepare input
# Load from file name if required
# Check/standardise columns
# Check TP_FP_comparisons are appropriate
.prepare_input = function(scores_df, TP_FP_comparisons, verbose = F){

  # If input is a string, assume file name and try to load
  if (is.character(scores_df)){
    
    # Get file extension
    file_extension = strsplit(basename(scores_df), split = ".", fixed = T)[[1]]
    file_extension = file_extension[length(file_extension)]
    if (verbose) cat("Loading from .", file_extension, " file\n", sep = "")
    
    
    # Load appropriately
    if (file_extension == "rds"){
      scores_df = readRDS(scores_df)
    } else if (file_extension %in% c("txt", "tsv")){
      require(readr)
      scores_df = read_tsv(scores_df, col_types = cols())
    }
  }
  
  
  # Input must actually be a data frame (or tibble)
  if (!("data.frame" %in% class(scores_df))) stop("Please provide either a data.frame or a file name of a data.frame")
  
  
  # Standardise columns
  # Sample
  if (!("sample" %in% colnames(scores_df))) stop("Input must have a sample column")
  # Score
  if (!("score" %in% colnames(scores_df))) stop("Input must have a score column")
  # Gene ID
  if (!("gene" %in% colnames(scores_df))){
    if ("entrez" %in% colnames(scores_df)){
      scores_df = scores_df %>% rename(gene = entrez)
      if (verbose) cat("Using 'entrez' as gene ID column\n")
    } else if ("symbol" %in% colnames(scores_df)){
      scores_df = scores_df %>% rename(gene = symbol)
      if (verbose) cat("Using 'symbol' as gene ID column\n")
    } else {
      stop("Input must have a gene/entrez/symbol column")
    }
  }
  # Gene type
  if (!("gene_type" %in% colnames(scores_df))){
    if ("type" %in% colnames(scores_df)){
      scores_df = scores_df %>% rename(gene_type = type)
      if (verbose) cat("Using 'type' as gene type column\n")
    } else {
      stop("Input must have a type/gene_type column")
    }
  }
  
  
  # Get minimal input
  scores_df = scores_df %>% select(sample, gene, gene_type, score)
  
  
  # Check NAs
  if (any(is.na(scores_df))){
    scores_df$NA_present = apply(scores_df, 1, function(x) any(is.na(x)))
    warning(paste0(sum(scores_df$NA_present), " rows of scores data.frame have NA values, removing"))
    scores_df = scores_df %>% filter(!NA_present) %>% select(-NA_present)
  }
  
  
  # Get sample ranks
  scores_df = scores_df %>% 
    arrange(sample, desc(score)) %>%
    group_by(sample) %>%
    mutate(sample_rank = rank(-score)) %>%
    ungroup
  
  
  # Make sure that all TP and FP types are represented correctly
  if (!is.list(TP_FP_comparisons) | 
      !identical(2L, unique(lengths(TP_FP_comparisons))) | 
      !identical("character", unique(unlist(lapply(TP_FP_comparisons, class))))){
    stop("TP_FP_comparisons must be a list of length 2 character vectors, corresponding to different gene types")
  }
  if (any(!(unique(unlist(TP_FP_comparisons)) %in% scores_df$gene_type))){
    missing_gene_types = setdiff(unique(unlist(TP_FP_comparisons)), scores_df$gene_type) %>% paste(collapse = ", ")
    if (verbose) warning(paste0("TP/FP gene type(s) ", missing_gene_types, " are not present in input data frame"))
  }
  
  
  # Return input
  return(scores_df)
}



# Backfill a vector which has NAs
.backfill = function(x){
  if (is.na(x[1])) x[1] = 0
  if (length(x) > 1){
    for (i in 2:length(x)){
      if (is.na(x[i])) x[i] = x[i-1]
    }
  }
  return(x)
}



# Cohort-level ROCs
.get_cohortLevel_rocs = function(scores_df, TP_FP_comparisons, verbose = F){
  
  
  # Double check that rows of scores_df are in the right order
  # This was done in .prepare_data, but better safe than sorry
  scores_df = scores_df %>% arrange(sample, desc(score))
  
  
  # Output is a list
  cohortLevel_rocs = list()
  
  
  # Smallest rank for each gene
  min_ranks = scores_df %>%
    group_by(gene, gene_type) %>%
    summarise(min_rank = min(sample_rank)) %>%
    ungroup
  
  
  # Numbers of altered genes of each type
  type_counts = scores_df %>%
    group_by(gene_type) %>%
    summarise(n_altered_tot = n_distinct(gene))
  
  
  # Bin genes by type and min_rank, calculate recall for each gene type
  recall = min_ranks %>%
    count(top_n = min_rank, gene_type) %>%
    group_by(gene_type) %>%
    mutate(n = cumsum(n)) %>%
    ungroup %>%
    left_join(type_counts, by = "gene_type") %>%
    mutate(recall = n / n_altered_tot)
  recall = rbind(type_counts %>% mutate(top_n = 0, n = 0, recall = 0), recall)
  recall = rbind(recall, type_counts %>% mutate(top_n = max(scores_df$sample_rank), n = 0, recall = 1))
  
  
  # Loop through TP_FP_comparisons
  for (TP_FP_types in TP_FP_comparisons){

    # Extract TP and FP gene types for comparison
    TP_type = TP_FP_types[1]
    FP_type = TP_FP_types[2]
    if (verbose) cat("Comparing", TP_type, "to", FP_type, "\n")
    this_scores = scores_df %>% subset(gene_type %in% TP_FP_types)
    comparison_name = paste(TP_FP_types, collapse = "_to_")
    
    # If we don't have both gene types, return dummy result
    if (!all(TP_FP_types %in% this_scores$gene_type)){
      this_roc = data.frame(top_n = NA, TP_recall = NA, FP_recall = NA, auc = NA)
    } else {
      # Get ROC
      this_roc = recall %>%
        subset(gene_type %in% c(TP_type, FP_type)) %>%
        arrange(top_n) %>%
        mutate(top_n = factor(top_n)) %>%
        select(-n_altered_tot, -n) %>%
        unique %>%
        spread(gene_type, recall)
      
      
      # Backfill
      colnames(this_roc)[colnames(this_roc) == TP_type] = "TP_recall"
      colnames(this_roc)[colnames(this_roc) == FP_type] = "FP_recall"
      this_roc$TP_recall = .backfill(this_roc$TP_recall)
      this_roc$FP_recall = .backfill(this_roc$FP_recall)
      
      
      # AUC - average of upper and lower integrals
      this_auc = this_roc %>% 
        mutate(d_x = FP_recall - lag(FP_recall), ave_y = (TP_recall + lag(TP_recall)) / 2) %>%
        mutate(auc_segment = d_x * ave_y) %>%
        .$auc_segment %>%
        sum(na.rm = T)
      this_roc$auc = this_auc
      
    }
    # Join to output
    cohortLevel_rocs[[comparison_name]] = this_roc
  }
  
  return(cohortLevel_rocs)
}



# Sample-level ROCs
.get_sampleLevel_rocs = function(scores_df, TP_FP_comparisons, verbose = F){
  
  
  # Double check that rows of scores_df are in the right order
  # This was done in .prepare_data, but better safe than sorry
  scores_df = scores_df %>% arrange(sample, desc(score))
  
  
  # Output is a list
  sampleLevel_rocs = list()
  
  
  # Loop through TP_FP_comparisons
  for (TP_FP_types in TP_FP_comparisons){
    
    # Extract TP and FP gene types for comparison
    TP_type = TP_FP_types[1]
    FP_type = TP_FP_types[2]
    if (verbose) cat("Comparing", TP_type, "to", FP_type, "\n")
    this_scores = scores_df %>% subset(gene_type %in% TP_FP_types)
    comparison_name = paste(TP_FP_types, collapse = "_to_")
  
    
    # If we don't have both gene types in any samples, return dummy result
    if (!all(TP_FP_types %in% this_scores$gene_type)){
      this_perSample_rocs = data.frame(sample = unique(scores_df$sample), TP_recall = NA, FP_recall = NA, auc = NA)
      this_median_roc = data.frame(TP_recall = NA, FP_recall = NA, auc = NA)
      no_TP_FP_samples = unique(scores_df$sample)
    } else {
      # Get recall of both types per sample
      this_recall = this_scores %>% 
        select(sample, gene_type) %>%
        group_by(sample) %>% 
        mutate(n_TP_tot = sum(gene_type == TP_type), n_FP_tot = sum(gene_type == FP_type),
               n_TP_cum = cumsum(gene_type == TP_type), n_FP_cum = cumsum(gene_type == FP_type)) %>%
        mutate(TP_recall = n_TP_cum / n_TP_tot, FP_recall = n_FP_cum / n_FP_tot)
      
      
      # Remove samples with no TP or no FP genes
      no_TP_FP_samples = this_recall %>%
        subset(n_TP_tot == 0 | n_FP_tot == 0) %>%
        .$sample %>%
        unique
      this_recall = this_recall %>% subset(n_TP_tot > 0 & n_FP_tot > 0)
      
      
      # Get ROC curves for each sample
      this_perSample_rocs = this_recall %>%
        mutate(TP_change = TP_recall == lead(TP_recall), FP_change = FP_recall == lead(FP_recall)) %>%
        mutate(last_TP_change = TP_change != lag(TP_change), last_FP_change = FP_change != lag(FP_change)) %>%
        subset(TP_change != lag(TP_change) | is.na(TP_change) | is.na(last_TP_change)) %>%
        mutate(d_x = FP_recall - lag(FP_recall)) %>%
        mutate(auc_segment = d_x * TP_recall) %>%
        mutate(auc = sum(auc_segment, na.rm = T)) %>%
        select(sample, TP_recall, FP_recall, auc) %>%
        ungroup
      # Add FP=0 and TP=0 row for each sample (helps plotting)
      this_perSample_rocs = rbind(
        this_perSample_rocs,
        this_perSample_rocs %>% subset(!duplicated(sample)) %>% mutate(TP_recall = 0, FP_recall = 0)
      ) %>%
        arrange(sample, TP_recall, FP_recall) %>%
        unique # Just in case we already had this row; otherwise get error with spread
      
      
      # Get the median ROC across samples
      # Expand and average the TP_recall over samples for each observed value of FP_recall
      this_median_roc = this_perSample_rocs %>%
        select(-auc) %>%
        subset(!(FP_recall == lead(FP_recall) & sample == lead(sample))) %>%
        spread(sample, TP_recall) %>%
        arrange(FP_recall)
      sample_cols = colnames(this_median_roc)[2:ncol(this_median_roc)]
      for (sample_col in sample_cols) this_median_roc[[sample_col]] = .backfill(this_median_roc[[sample_col]])
      this_median_roc$TP_recall = apply(this_median_roc %>% select(sample_cols), 1, median)
      this_median_roc = this_median_roc %>%
        select(FP_recall, TP_recall) %>%
        subset(TP_recall != lag(lag(TP_recall)) & FP_recall != lag(lag(FP_recall)) | 
               is.na(lag(lag(TP_recall))) | is.na(lag(lag(FP_recall))) | is.na(lead(FP_recall)))
      
      
      # Add FP=0 and TP=0 row for each sample (helps plotting)
      this_median_roc = rbind(
        this_median_roc,
        this_median_roc %>% head(1) %>% mutate(TP_recall = 0, FP_recall = 0)
      ) %>%
        arrange(TP_recall, FP_recall) %>%
        unique # Just in case we already had this row
      
      
      # Median of per-sample AUCs
      this_median_roc$auc = this_perSample_rocs %>%
        select(sample, auc) %>%
        unique %>%
        .$auc %>%
        median
    }
    
    # Join to output
    sampleLevel_rocs[[comparison_name]] = list(
      perSample_rocs = this_perSample_rocs,
      median_roc = this_median_roc,
      no_TP_FP_samples = no_TP_FP_samples
    )
  }
  
  return(sampleLevel_rocs)
}




# # Sample-level ROCs when not all genes are scored
# scores_df = oncoIMPACT; TP_FP_comparisons = list(c("CG", "FP"), c("CG", "Non-CG")); verbose = T
# .get_partial_sampleLevel_rocs = function(scores_df, n_damaged, TP_FP_comparisons, verbose = F){
#   
#   
#   # Double check that rows of scores_df are in the right order
#   # This was done in .prepare_data, but better safe than sorry
#   scores_df = scores_df %>% arrange(sample, desc(score))
#   
#   
#   # Output is a list
#   sampleLevel_rocs = list()
#   
#   
#   # Loop through TP_FP_comparisons
#   for (TP_FP_types in TP_FP_comparisons){
#     
#     # Extract TP and FP gene types for comparison
#     TP_type = TP_FP_types[1]
#     FP_type = TP_FP_types[2]
#     if (verbose) cat("Comparing", TP_type, "to", FP_type, "\n")
#     this_scores = scores_df %>% subset(gene_type %in% TP_FP_types)
#     comparison_name = paste(TP_FP_types, collapse = "_to_")
#     
#     
#     # Get recall of both types per sample
#     this_n_damaged = n_damaged %>%
#       subset(gene_type %in% c(TP_type, FP_type)) %>%
#       spread(gene_type, n, fill = 0)
#     colnames(this_n_damaged)[colnames(this_n_damaged) == TP_type] = "n_TP_tot"
#     colnames(this_n_damaged)[colnames(this_n_damaged) == FP_type] = "n_FP_tot"
#     this_recall = this_scores %>% 
#       select(sample, gene_type) %>%
#       left_join(this_n_damaged, by = "sample") %>%
#       group_by(sample) %>% 
#       mutate(n_TP_cum = cumsum(gene_type == TP_type), n_FP_cum = cumsum(gene_type == FP_type)) %>%
#       mutate(TP_recall = n_TP_cum / n_TP_tot, FP_recall = n_FP_cum / n_FP_tot)
#     
#     
#     # Remove samples with no TP or no FP genes
#     no_TP_FP_samples = n_damaged %>%
#       subset(gene_type %in% c(TP_type, FP_type) & n == 0) %>%
#       .$sample %>%
#       unique
#     this_recall = this_recall %>% subset(!(sample %in% no_TP_FP_samples))
#     
#     
#     # Get ROC curves for each sample
#     this_perSample_rocs = this_recall %>%
#       mutate(TP_change = TP_recall == lead(TP_recall), FP_change = FP_recall == lead(FP_recall)) %>%
#       mutate(last_TP_change = TP_change != lag(TP_change), last_FP_change = FP_change != lag(FP_change)) %>%
#       subset(TP_change != lag(TP_change) | is.na(TP_change) | is.na(last_TP_change)) %>%
#       mutate(d_x = FP_recall - lag(FP_recall)) %>%
#       mutate(auc_segment = d_x * TP_recall) %>%
#       mutate(auc = sum(auc_segment, na.rm = T)) %>%
#       select(sample, TP_recall, FP_recall, auc) %>%
#       ungroup
#     this_perSample_rocs[is.nan(this_perSample_rocs$TP_recall), "TP_recall"] = 0
#     this_perSample_rocs[is.nan(this_perSample_rocs$FP_recall), "FP_recall"] = 0
#     
#     
#     plot_samples = base::sample(unique(this_perSample_rocs$sample), size = 6)
#     this_perSample_rocs %>%
#       subset(sample %in% plot_samples) %>%
#       ggplot(aes(FP_recall, TP_recall)) +
#       geom_step() +
#       facet_wrap(~sample, nrow = 2, scales = "free")
#     
#     
#     
#     # Add FP=0 and TP=0 row for each sample (helps plotting)
#     this_perSample_rocs = rbind(
#       this_perSample_rocs,
#       this_perSample_rocs %>% subset(!duplicated(sample)) %>% mutate(TP_recall = 0, FP_recall = 0)
#     ) %>%
#       arrange(sample, TP_recall, FP_recall) %>%
#       unique # Just in case we already had this row; otherwise get error with spread
#     
#     
#     # Get the median ROC across samples
#     # Expand and average the TP_recall over samples for each observed value of FP_recall
#     this_median_roc = this_perSample_rocs %>%
#       select(-auc) %>%
#       subset(!(FP_recall == lead(FP_recall) & sample == lead(sample))) %>%
#       spread(sample, TP_recall) %>%
#       arrange(FP_recall)
#     sample_cols = colnames(this_median_roc)[2:ncol(this_median_roc)]
#     for (sample_col in sample_cols) this_median_roc[[sample_col]] = .backfill(this_median_roc[[sample_col]])
#     this_median_roc$TP_recall = apply(this_median_roc %>% select(sample_cols), 1, median)
#     this_median_roc = this_median_roc %>%
#       select(FP_recall, TP_recall) %>%
#       subset(TP_recall != lag(lag(TP_recall)) & FP_recall != lag(lag(FP_recall)) | 
#                is.na(lag(lag(TP_recall))) | is.na(lag(lag(FP_recall))) | is.na(lead(FP_recall)))
#     
#     
#     # Add FP=0 and TP=0 row for each sample (helps plotting)
#     this_median_roc = rbind(
#       this_median_roc,
#       this_median_roc %>% head(1) %>% mutate(TP_recall = 0, FP_recall = 0)
#     ) %>%
#       arrange(TP_recall, FP_recall) %>%
#       unique # Just in case we already had this row
#     
#     
#     # Median of per-sample AUCs
#     this_median_roc$auc = this_perSample_rocs %>%
#       select(sample, auc) %>%
#       unique %>%
#       .$auc %>%
#       median
#     
#     
#     # Join to output
#     sampleLevel_rocs[[comparison_name]] = list(
#       perSample_rocs = this_perSample_rocs,
#       median_roc = this_median_roc,
#       no_TP_FP_samples = no_TP_FP_samples
#     )
#   }
#   
#   return(sampleLevel_rocs)
# }
