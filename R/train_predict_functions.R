# These functions are for the 'clean' script, accompanying the methods paper


# Mark training/prediction sets - provide cancer genes and their corresponding cancer-promoting alteration types
mark_truePositives = function(df, truePositive_drivers){
  
  # The altCond column in truePositive_drivers is a string which can be evaluated to a Boolean on each row of df
  # But we have potentially different altConds for different driver genes; loop through each type of altCond 
  for (driverType in unique(truePositive_drivers$altCond)){
    # The genes with this kind of driver alteration
    driverGenes = truePositive_drivers %>% filter(altCond == driverType) %>% .$entrez
    
    # Which rows of df match this for these genes
    driverRows = df %>% mutate(TP = eval(parse(text = driverType)) & entrez %in% driverGenes) %>% .$TP
    
    # Mark these rows
    df[driverRows, "type"] = "TP"
  }
  
  
  # We end up with a new column, a factor with the levels "TP" and "prediction"
  df[is.na(df$type), "type"] = "prediction"
  df$type = factor(df$type)
  
  
  return(df)
}





# Scale data
get_scaledTable = function(df, scaling_factors = NULL, 
                           mean_center = c("none", "all", "training", "prediction"),
                           sd_scale = c("all", "none", "training", "prediction")){

  # If scaling_factors are not provided (i.e. we are training rather than predicting), we need to determine them first
  if (is.null(scaling_factors)){
    
    
    # Only scale numeric features, but obviously not the sample/entrez columns if these are numeric!
    cols_toScale = setdiff(names(which(sapply(df, function(x) class(x) %in% c("numeric", "integer")))), c("sample", "entrez"))
    
    
    # Calculate mean centering factors
    mean_center = match.arg(mean_center)
    if (mean_center == "all"){
      m = apply(df %>% select(cols_toScale), 2, mean)
    } else if (mean_center == "none"){
      m = rep(0, length(cols_toScale))
      names(m) = cols_toScale
    } else if (mean_center == "training"){
      m = apply(df %>% subset(type == "TP") %>% select(cols_toScale), 2, mean)
    } else if (mean_center == "prediction"){
      m = apply(df %>% subset(type == "prediction") %>% select(cols_toScale), 2, mean)
    }
    
    
    # Calculate SD scaling factors
    sd_scale = match.arg(sd_scale)
    if (sd_scale == "all"){
      s = apply(df %>% select(cols_toScale), 2, sd)
    } else if (sd_scale == "none"){
      s = rep(1, length(cols_toScale))
      names(s) = cols_toScale
    } else if (sd_scale == "training"){
      s = apply(df %>% subset(type == "TP") %>% select(cols_toScale), 2, sd)
    } else if (sd_scale == "prediction"){
      s = apply(df %>% subset(type == "prediction") %>% select(cols_toScale), 2, sd)
    }
    
    
    # Compile into single list
    scaling_factors = list()
    for (col in setdiff(colnames(df), c("sample", "entrez", "type"))){
      scaling_factors[[col]] = list(
        type = class(df[[col]]),
        scale = col %in% cols_toScale
      )
      if (col %in% cols_toScale){
        scaling_factors[[col]]$factors = c(mean_translation = unname(m[col]), sd_scaling = unname(s[col]))
      }
    }
  }
  
  
  # Scale the data
  df_scaled = df
  for (col in names(scaling_factors)){
    if (scaling_factors[[col]]$scale){
      df_scaled[[col]] = (df_scaled[[col]] - scaling_factors[[col]]$factors["mean_translation"]) / scaling_factors[[col]]$factors["sd_scaling"]
    }
  }
  return(list(df_scaled = df_scaled, scaling_factors = scaling_factors))
}





# Prepare training and prediction sets
prepare_trainingPrediction = function(# Data
  cohort_data,
  truePositive_drivers,
  output_dir = NULL,
  # Exclude observations/features
  exclude_features = c("no_NSI_muts", "rnaseq_ubiquitous", "protein_low", "ppin_central", "rnaseq_selective", "length_RefSeq", "domains_InterPro", "multiDomain"),
  damaging_subsetCond = "no_TRUNC_muts + no_NTDam_muts + no_GOF_muts > 0 | CNVGain == 1 | Copy_number == 0",
  normalTissue_expressionData = NULL,
  normalTissue_name = NULL,
  # Scaling
  factorise_twoLevel = T,
  nonFactor_features = c("sample", "entrez", "no_ALL_muts", "no_NSI_muts", "no_NTDam_muts", "no_TRUNC_muts", "no_GOF_muts", "BND", "INV", "INS"),
  mean_center = c("none", "all", "training", "prediction"),
  sd_scale = c("all", "none", "training", "prediction"),
  # Housekeeping
  sample_gene_sep = "__"
){
  
  require(tibble)
  require(tidyr)
  require(dplyr)
  cat("Preparing training and prediction sets\n")
  
  
  # Load cohort_data if provided as a file name
  if (is.character(cohort_data)) cohort_data = readRDS(cohort_data)
  df = cohort_data
  
  
  # Remove non-expressed genes if normal tissue expression data (e.g. GTEx) are provided
  if (!is.null(normalTissue_expressionData)){
    nonExpressed_genes = normalTissue_expressionData %>% 
      subset(tissue %in% normalTissue_name) %>%
      group_by(entrez) %>%
      summarise(distinct_exp_level = paste(unique(exp_level), collapse = ";")) %>%
      subset(distinct_exp_level == "Not expressed") %>%
      .$entrez
    
    cat("Removing", length(intersect(df$entrez, nonExpressed_genes)), "genes not expressed in normal", paste(normalTissue_name, collapse = ", "), "\n")
    df = df %>% filter(!(entrez %in% nonExpressed_genes))
  }
  
  
  # Remove any features that we want to exclude
  if (!is.null(exclude_features)) cat("Excluding features:", paste(exclude_features, collapse = ", "), "\n")
  features = setdiff(colnames(df), exclude_features)
  df = df %>% select(one_of(c("sample", "entrez", features)))
  
  
  # Remove all alterations that are not determined to be damaging, unless damaging_subsetCond is null
  if (!is.null(damaging_subsetCond)){
    n_nonDamaging_rows = df %>% filter(!eval(parse(text = damaging_subsetCond))) %>% nrow
    cat("Removing", n_nonDamaging_rows, "alterations which are non-damaging\n")
    df = df %>% filter(eval(parse(text = damaging_subsetCond)))
  }
  
  
  # Mark which alterations are "true positive" drivers, and separate into training and prediction sets
  # Make sure that any gene included in the training set (no matter what alteration types) is not in the prediction set
  cat("Separating training and prediction sets\n")
  df = mark_truePositives(df = df, truePositive_drivers = truePositive_drivers)
  
  
  # Convert features with only two unique values to factors, unless otherwise specified
  # These features will not undergo any scaling
  if (factorise_twoLevel){
    twoLevel_cols = names(which(sapply(df, function(x) n_distinct(x) == 2)))
    factorCols = setdiff(twoLevel_cols, nonFactor_features)
    cat("Converting binary features to factors:", paste(factorCols, collapse = ", "), "\n")
    for (col in factorCols) df[[col]] = factor(df[[col]])
  }
  
  
  # Scale data
  scaling_res = get_scaledTable(df = df, mean_center = mean_center, sd_scale = sd_scale)
  df_scaled = scaling_res$df_scaled
  scaling_factors = scaling_res$scaling_factors
  
  
  # Remove sample and entrez columns, instead use sample__entrez as row name
  # Function to remove row names from a data frame; sometimes column_to_rownames doesn't work otherwise
  unrowname = function(x){
    rownames(x) = c()
    return(x)
  }
  # Unscaled data
  df = df %>% 
    subset(!(type == "prediction" & entrez %in% truePositive_drivers$entrez)) %>%
    unite("key", sample, entrez, sep = sample_gene_sep) %>%
    unrowname %>%
    column_to_rownames("key")
  training_ns = df %>% subset(type == "TP")
  prediction_ns = df %>% subset(type == "prediction") %>% select(-type)
  # Scaled data
  df_scaled = df_scaled %>% 
    subset(!(type == "prediction" & entrez %in% truePositive_drivers$entrez)) %>%
    unite("key", sample, entrez, sep = sample_gene_sep) %>%
    unrowname %>%
    column_to_rownames(var = "key")
  training = df_scaled %>% subset(type == "TP")
  prediction = df_scaled %>% subset(type == "prediction") %>% select(-type)
  
  
  # Save scaling factors, training and prediction sets
  if (!is.null(output_dir)){
    if (!dir.exists(output_dir)) dir.create(output_dir)
    saveRDS(scaling_factors, file = paste(output_dir, "scaling_factors.rds", sep = "/"))
    saveRDS(training, file = paste(output_dir, "training_set.rds", sep = "/"))
    saveRDS(training_ns, file = paste(output_dir, "training_set_noScale.rds", sep = "/"))
    saveRDS(prediction, file = paste(output_dir, "prediction_set.rds", sep = "/"))
    saveRDS(prediction_ns, file = paste(output_dir, "prediction_set_noScale.rds", sep = "/"))
  }
  
  
  # Output
  return(list(scaling_factors = scaling_factors, 
              training_set = training, training_set_ns = training_ns,
              prediction_set = prediction, prediction_set_ns = prediction_ns))
  cat("Dataset divided into", nrow(training), "training observations and", nrow(prediction), "observations for prediction\n")
}





# Carry out cross validation iterations, recording the sensitivity
# Runs in a parallel environment
run_crossValidation_par = function(inPath,
                                   outPath = NULL,
                                   folds = 3, 
                                   iters = 100,
                                   kernels_paramGrids = list(
                                     linear = list(nu = seq(0.05, 0.35, 0.05)),
                                     polynomial = list(nu = seq(0.05, 0.35, 0.05), gamma = 1, coef0 = 0, degree = c(2, 3, 4)),
                                     radial = list(nu = seq(0.05, 0.35, 0.05), gamma = 2^seq(-7, 4)),
                                     sigmoid = list(nu = seq(0.05, 0.35, 0.05), gamma = 2^seq(-7, 4), coef0 = 1)
                                   ),
                                   cores = 2, 
                                   sample_gene_sep = "__",
                                   verbose = F,
                                   parallelLib = "snow"
){
  
  require(parallel)
  require(snow)
  require(doSNOW)
  require(tibble)
  require(dplyr)
  require(e1071)
  cat("Beginning cross-validation for hyperparameter tuning\n")
  
  
  ##---- Sort out directory structure for verbose mode ----
  
  # In verbose mode, we will save the results of every iteration for every kernel-parameter combination
  # Create the directories required for this
  if (verbose){
    # Containing directory, called CV
    cv_dir = paste(outPath, "CV", sep = "/")
    system(paste("mkdir", cv_dir))
    
    
    # Just below this, one directory per kernel
    kernel_dirs = paste(cv_dir, names(kernels_paramGrids), sep = "/") %>% paste(collapse = " ")
    system(paste("mkdir", kernel_dirs))
    
    
    # Keep track of what directory is for what kernel/parameters/iteration
    kernelParamIter_to_dir = list()
    
    
    for (k in names(kernels_paramGrids)){
      # Make the directory names, keeping track of what setting they belong to
      params_iterations = kernels_paramGrids[[k]]
      params_iterations$iteration = 1:iters
      paramIter_to_dir = expand.grid(params_iterations, KEEP.OUT.ATTRS = F)
      paramIter_to_dir$dir = apply(paramIter_to_dir, 1, function(x){
        paramDir = x[names(x) != "iteration"]
        paramDir = paste(names(paramDir), paramDir, sep = ".") %>% paste(collapse = "_")
        paramDir = paste(cv_dir, k, paramDir, sep = "/")
        paste0(paramDir, "/iteration.", x["iteration"])
      })
      kernelParamIter_to_dir[[k]] = paramIter_to_dir
      
      
      # Make the directories
      paramIter_dirString = unique(dirname(paramIter_to_dir$dir)) %>% paste(collapse = " ")
      system(paste("mkdir", paramIter_dirString))
      system(paste("mkdir", paramIter_to_dir$dir %>% paste(collapse = " ")))
    }
  }
  
  
  ##---- Prepare grid search ----
  
  # Get the training data from inPath
  # To avoid confusion between this and the subsets that will actually be used to train in each CV iteration, call this truePositives
  truePositives = readRDS(paste(inPath, "training_set.rds", sep = "/"))
  # Unique genes in the true postitive set
  truePositive_genes = truePositives %>% 
    tibble::rownames_to_column("key") %>% 
    separate(key, into = c("sample", "entrez"), sep = sample_gene_sep) %>% 
    .$entrez %>% 
    unique
  
  
  # Record of sensitivity for all kernel-parameter-iteration combinations
  cv_stats_full = NULL
  
  
  # Begin parameter grid search, separately for each kernel
  kernels = names(kernels_paramGrids)
  cat("Parameter grid search for:", "\n")
  for(k in kernels){
    cat(k, "\n")
    
    # Table of parameter combinations and iterations to work through for this kernel
    # Parallisation works through the rows of this table
    kernels_paramGrids[[k]]$iteration = 1:iters
    params_iterations = expand.grid(kernels_paramGrids[[k]], KEEP.OUT.ATTRS = F)
    
    
    # For parameters that have not been specified, fix to default values
    # Some, e.g. gamma in the linear kernel, will just be ignored by the svm function
    # Give a warning to the user if they haven't specified one that might be relevant
    # Nu, used for all kernels
    nuDefault = 0.05
    if (!"nu" %in% colnames(params_iterations)){
      warning("You have not specified a grid range for nu. Fixing at ", nuDefault)
      params_iterations$nu = nuDefault
    }
    # Gamma, used for polynomial and sigmoid kernels
    gammaDefault = 1
    if (!"gamma" %in% colnames(params_iterations)){
      if (k %in% c("polynomial", "sigmoid")) warning("You have not specified a grid range for gamma. Fixing at ", gammaDefault)
      params_iterations$gamma = gammaDefault
    }
    # Coef0, used for polynomial and sigmoid kernels
    coef0Default = 0
    if (!"coef0" %in% colnames(params_iterations)){
      if (k %in% c("polynomial", "sigmoid")) warning("You have not specified a grid range for coef0. Fixing at ", coef0Default)
      params_iterations$coef0 = coef0Default
    }
    # Degree, used for polynomial kernel
    degreeDefault = 2
    if (!"degree" %in% colnames(params_iterations)){
      if (k == "polynomial") warning("You have not specified a grid range for degree Fixing at ", degreeDefault)
      params_iterations$degree = degreeDefault
    }
    
    
    # Set up parallel environment, including progress bar
    n_tasks = nrow(params_iterations)
    if (parallelLib == "snow"){
      cl = snow::makeCluster(cores, type = "SOCK") # For use on an HPC
    } else if (parallelLib == "parallel"){
      cl = parallel::makeCluster(cores, type = "SOCK") # For use on a local machine
    }
    registerDoSNOW(cl)
    pb = txtProgressBar(max = n_tasks, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress = progress)
    
    
    ##---- Parallel grid search ----
    
    # Work through parameter combinations and iterations in a parallel manner
    # Record the sensitivity on the training set in each case
    cv_stats_thisKernel = foreach(i = 1:n_tasks, 
                                  .combine = rbind, 
                                  .packages = c("e1071", "dplyr"),
                                  .options.snow = opts) %dopar%
      {
        ##---- Setup ----
        
        # Because parallel
        Sys.sleep(0.01)
        
        
        # Extract kernel parameters and iteration number
        mynu = params_iterations[i, "nu"]
        mygamma = params_iterations[i, "gamma"]
        mydegree = params_iterations[i, "degree"]
        mycoef0 = params_iterations[i, "coef0"]
        myiteration = params_iterations[i, "iteration"]
        
        
        # Identify directory for this iteration (verbose mode only)
        if (verbose){
          cols_to_match = intersect(colnames(params_iterations), colnames(kernelParamIter_to_dir[[k]]))
          j = sapply(cols_to_match, function(x) kernelParamIter_to_dir[[k]][[x]] == params_iterations[i, x])
          j = apply(j, 1, all)
          iteration_dir = kernelParamIter_to_dir[[k]][which(j), "dir"]
        }
        
        
        ##---- Separate training and test sets ----
        
        # Choose which true positive genes to use for training, and which to use to test the sensitivity
        # This is leave-one-out cross-validation based on the unique genes rather than the rows of the full training set
        traingenes = base::sample(truePositive_genes, size = ceiling( (folds-1)/folds * length(truePositive_genes) ))
        testgenes = setdiff(truePositive_genes, traingenes)
        
        
        # Extract trainset and testset
        trainset = truePositives %>% filter(sapply(rownames(.), function(key) strsplit(key, split = sample_gene_sep)[[1]][2]) %in% traingenes)
        testset = truePositives %>% filter(!sapply(rownames(.), function(key) strsplit(key, split = sample_gene_sep)[[1]][2]) %in% traingenes)
        
        
        # Save the trainset and testset (verbose mode only)
        if (verbose){
          saveRDS(trainset, file = paste(iteration_dir, "trainset.rds", sep = "/"))
          saveRDS(testset, file = paste(iteration_dir, "testset.rds", sep = "/"))
        }
        
        
        ##---- SVM ----
        
        # Train the one-class SVM
        svm.model = svm(formula = type ~ ., data = trainset, scale = FALSE,
                        type = "one-classification", kernel = k, 
                        nu = mynu, gamma = mygamma, degree = mydegree, coef0 = mycoef0)
        
        
        # Predict on the remaining fold of true positives
        prediction = predict(svm.model, newdata = testset %>% select(-type), decision.values = TRUE, probability = TRUE)
        
        
        # Save SVM model and predictions (verbose mode only)
        if (verbose){
          saveRDS(svm.model, file = paste(iteration_dir, "svmModel.rds", sep="/"))
          saveRDS(prediction, file = paste(iteration_dir, "prediction.rds", sep="/"))
        }
        
        
        ##---- Assess accuracy ----
        
        # Record the sensitivity on the test set
        cv_stats = data.frame(
          kernel = k,
          nu = mynu,
          gamma = mygamma,
          coef0 = mycoef0,
          degree = mydegree,
          iteration = myiteration,
          trainSize = nrow(trainset),
          sensitivity = sum(prediction) / length(prediction),
          stringsAsFactors = F
        )
        return(cv_stats)
      }
    
    
    # Close the parallel environment and progress bar
    close(pb)
    stopCluster(cl)
    
    
    # Add sensitivity to the records for all kernels
    cv_stats_full = rbind(cv_stats_full, cv_stats_thisKernel)
    rm(cv_stats_thisKernel)
  }
  
  
  # Write cv_stats to file
  if (!is.null(outPath)){
    cv_stats_fn = paste(outPath, "cv_stats.tsv", sep="/")
    write_tsv(cv_stats_full, path = cv_stats_fn)
    cat("CV results written to", cv_stats_fn, "\n")
  }
  
  return(cv_stats_full)
}





# Use a CV stats file to select best models according to sensitivity
selectParams_from_CVstats = function(cv_stats, output_dir = NULL){ #step = NULL, 
  
  require(readr)
  require(dplyr)
  
  
  # If file name was provided instead of the actual data frame, load into R
  if (is.character(cv_stats)) cv_stats = read_tsv(cv_stats, col_types = cols())
  
  
  # # Determine at what steps (of # of iterations) to do scoring
  # iters = max(cv_stats$iteration)
  # if (is.null(step)){
  #   steps = iters
  # } else {
  #   steps = seq(step, iters, step)
  #   if (steps[length(steps)] != iters) steps = c(steps, iters)
  # }
  # 
  # 
  # # Calculate mean and SD of sensitivities at cumulative intervals
  # cv_summary = data.frame()
  # for (st in steps){
  #   df = cv_stats %>%
  #     subset(iteration <= st) %>%
  #     group_by(kernel, nu, gamma, degree, coef0) %>%
  #     summarise(mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity)) %>%
  #     ungroup %>%
  #     mutate(iterations_cumulative = st)
  #   
  #   cv_summary = rbind(cv_summary, df)
  # }
  # rm(df)
  # Calculate mean and SD of sensitivities cumulatively at each iteration (to assess convergence)
  cv_summary = cv_stats %>%
    rename(iterations_cumulative = iteration) %>%
    arrange(iterations_cumulative) %>%
    group_by(kernel, nu, gamma, degree, coef0) %>%
    mutate(sens_cumul = cumsum(sensitivity), sens2_cumul = cumsum(sensitivity^2)) %>%
    ungroup %>%
    mutate(mean_sensitivity = sens_cumul / iterations_cumulative) %>%
    mutate(sd_sensitivity = sqrt((sens2_cumul - iterations_cumulative*mean_sensitivity^2) / (iterations_cumulative - 1)))
  
  
  # Function to convert a vector to a z-score
  z_score = function(x) (x-mean(x)) / sd(x)
  
  
  # At each iteration step for each kernel, select the parameter combination with the highest diff_z
  # diff_z measures the relative performance based on high mean and low variance of sensitivity
  best_models_cumulative = cv_summary %>%
    group_by(kernel, iterations_cumulative) %>%
    mutate(mean_z = z_score(mean_sensitivity), sd_z = z_score(sd_sensitivity)) %>%
    mutate(diff_z = mean_z - sd_z) %>%
    top_n(1, diff_z) %>%
    ungroup
  
  
  # For clarity, also get the final (i.e. after all CV iterations) best model parameters
  best_model_final = best_models_cumulative %>% 
    subset(iterations_cumulative == max(iterations_cumulative)) %>%
    rename(iterations_total = iterations_cumulative)
  
  
  # Save
  if (!is.null(output_dir)){
    write_tsv(best_models_cumulative, path = paste(output_dir, "best_models_cumulative.tsv", sep = "/"))
    write_tsv(best_model_final, path = paste(output_dir, "best_model_final.tsv", sep = "/"))
  }
  
  
  # Output
  return(list(best_model_final = best_model_final, best_models_cumulative = best_models_cumulative))
}





# Given kernel parameters (i.e. best models), train on true positives
# Default values from Supplementary Figure 3A in paper
train_sysSVM2 = function(model_parameters = list(linear = list(nu = 0.05),
                                                 polynomial = list(nu = 0.05, degree = 2),
                                                 radial = list(nu = 0.05, gamma = 2^-7),
                                                 sigmoid = list(nu = 0.05, gamma = 2)),
                         training_set, scaling_factors, output_dir = NULL){
  
  require(readr)
  require(dplyr)
  require(e1071)
  
  
  # Read model parameters and training data if they weren't provided directly
  if (is.character(model_parameters)) model_parameters = read_tsv(model_parameters)
  if (is.character(training_set)) training_set = readRDS(training_set)
  
  
  # If model parameters are provided as a list, convert to the legacy data frame format
  if (is.list(model_parameters)){
    model_parameters_list = model_parameters
    model_parameters = data.frame()
    for (k in names(model_parameters_list)){
      model_parameters = rbind(
        model_parameters,
        tibble(
          kernel = k,
          nu = ifelse("nu" %in% names(model_parameters_list[[k]]), model_parameters_list[[k]]$nu, 0.05),
          gamma = ifelse("gamma" %in% names(model_parameters_list[[k]]), model_parameters_list[[k]]$gamma, 1),
          coef0 = ifelse("coef0" %in% names(model_parameters_list[[k]]), model_parameters_list[[k]]$coef0, 0),
          degree = ifelse("degree" %in% names(model_parameters_list[[k]]), model_parameters_list[[k]]$degree, 2)
        )
      )
    }
  }
  
  
  # Train a ocSVM for each kernel in turn
  trained_sysSVM = list()
  for (i in 1:nrow(model_parameters)){
    # Extract parameters
    k = as.character(model_parameters$kernel[i])
    mynu = model_parameters$nu[i]
    mygamma = model_parameters$gamma[i]
    mycoef0 = model_parameters$coef0[i]
    mydegree = model_parameters$degree[i]
    
    
    # Name of model, for saving
    if (k == "linear"){
      modelName = paste0(k, "__nu_", mynu)
    } else if (k == "polynomial"){
      modelName = paste0(k, "__nu_", mynu, "__gamma_", mygamma, "__coef0_", mycoef0, "__degree_", mydegree)
    } else if (k %in% c("radial", "sigmoid")){
      modelName = paste0(k, "__nu_", mynu, "__gamma_", mygamma, "__coef0_", mycoef0)
    }
    
    
    # Train on full training set
    svm_model = svm(type ~ ., data = training_set, scale = FALSE, 
                    kernel = k, type = "one-classification", 
                    nu = mynu, gamma = mygamma, coef0 = mycoef0, degree = mydegree)
    
    
    # Get sensitivity mean and variance - these are used in the final score
    sensitivity_mean = model_parameters$mean_sensitivity[i]
    sensitivity_var = model_parameters$sd_sensitivity[i]^2
    
    
    # Add everything to list
    trained_sysSVM[[k]] = list(
      name = modelName,
      svm_model = svm_model,
      sensitivity_mean = sensitivity_mean,
      sensitivity_var = sensitivity_var
    )
  }
  
  
  # Add scaling factors to saved object
  trained_sysSVM = list(models = trained_sysSVM, scaling_factors = scaling_factors)
  
  
  # Save
  if (!is.null(output_dir)) saveRDS(trained_sysSVM, paste(output_dir, "trained_sysSVM.rds", sep = "/"))
  
  
  # Output
  return(trained_sysSVM)
}





# Use a trained sysSVM model to score alterations
predict_sysSVM2 = function(trained_sysSVM, 
                           molecular_data = NULL, systemsLevel_data = "example_data/systemsLevel_features_allGenes.tsv", # For use with pre-trained models
                           prediction_set = NULL, prediction_set_ns = NULL, # For use with de novo models
                           sample_gene_sep = "__", output_dir = NULL){
  
  require(tidyr)
  require(tibble)
  require(dplyr)
  require(e1071)
  
  
  # Make sure either molecular_data/systemsLevel_data or prediction_set/prediction_set_ns are provided
  pretrained = !is.null(molecular_data) & !is.null(systemsLevel_data)
  denovo = !is.null(prediction_set) & !is.null(prediction_set_ns)
  if (!pretrained & !denovo) stop("Please provide either molecular_data and systemsLevel_data (for pre-trained models), or prediction_set and prediction_set_ns (de novo models)")
  
  
  # Load trained model if not provided directly
  if (is.character(trained_sysSVM)) trained_sysSVM = readRDS(trained_sysSVM)
  
  
  # If molecular_data is provided, scale it
  if (pretrained){
    message("Applying pre-trained model")
    scaling_factors = trained_sysSVM$scaling_factors
    
    # Load data if needed
    if (is.character(molecular_data)) molecular_data = read_tsv(molecular_data, col_types = cols())
    if (is.character(systemsLevel_data)) systemsLevel_data = read_tsv(systemsLevel_data, col_types = cols())
    
    # Join molecular and systems-level data
    prediction_set_ns = inner_join(molecular_data, systemsLevel_data, by = "entrez")
    
    # Select appropriate columns
    if (length(setdiff(names(scaling_factors), colnames(prediction_set_ns))) > 0){
      stop(paste0("Missing features: ", paste(setdiff(names(scaling_factors), colnames(prediction_set_ns)), collapse = ", ")))
    }
    prediction_set_ns = prediction_set_ns[, c("sample", "entrez", names(scaling_factors))]
    prediction_set_ns = prediction_set_ns %>% 
      unite("key", sample, entrez, sep = sample_gene_sep) %>%
      column_to_rownames("key")
    
    # Convert columns to factors where necessary
    for (col in colnames(prediction_set_ns)){
      if (scaling_factors[[col]]$type == "factor") prediction_set_ns[[col]] = factor(prediction_set_ns[[col]], levels = c("0", "1"))
    }
    
    # Scale data
    prediction_set = get_scaledTable(df = prediction_set_ns, scaling_factors = trained_sysSVM$scaling_factors)$df_scaled
    
  } else {
    if (is.character(prediction_set)) prediction_set = readRDS(prediction_set)
    if (is.character(prediction_set_ns)) prediction_set_ns = readRDS(prediction_set_ns)
  }
  
  
  # Get decision values and score contributions for each kernel
  scored_data = prediction_set_ns # Record the unscaled features
  for (k in names(trained_sysSVM$models)){
    # Apply the ocSVM
    dv = predict(trained_sysSVM$models[[k]]$svm_model, newdata = prediction_set, decision.values = TRUE) # But predict on scaled features
    
    
    # Extract decision values
    scored_data[[paste(k, "decision_value", sep = "_")]] = as.numeric(attr(dv, "decision.values"))
    
    
    # Sensitivity weighting
    BMS_i = trained_sysSVM$models[[k]]$sensitivity_mean #/ trained_sysSVM$models[[k]]$sensitivity_var
    
    
    # Score contribution
    scored_data[[paste(k, "score", sep = "_")]] = scored_data %>% 
      rownames_to_column("key") %>%
      separate("key", into = c("sample", "entrez"), sep = sample_gene_sep) %>%
      select(sample, decision_value = paste(k, "decision_value", sep = "_")) %>% 
      group_by(sample) %>%
      mutate(r = rank(decision_value), N_s = n()) %>% # base::rank assigns 1 to the lowest value, n to the highest value
      mutate(R_igs = N_s - r + 1) %>%
      ungroup %>%
      mutate(score = -log10(R_igs / N_s) * BMS_i / log10(N_s)) %>%
      .$score
  }
  
  
  # Recover sample and entrez columns
  scored_data = scored_data %>%
    rownames_to_column("key") %>%
    separate("key", into = c("sample", "entrez"), sep = sample_gene_sep) %>%
    mutate(entrez = as.numeric(entrez))
  
  
  # Average scores from the different kernels
  scored_data$score = apply(scored_data %>% select(paste(names(trained_sysSVM$models), "score", sep = "_")), 1, mean)
  
  
  # Patient-specific ranks
  scored_data = scored_data %>%
    arrange(desc(score)) %>%
    group_by(sample) %>%
    mutate(sample_rank = 1:n()) %>%
    ungroup
  
  
  # Save
  if (!is.null(output_dir)) saveRDS(scored_data, paste(output_dir, "scores.rds", sep = "/"))
  
  
  # Output
  return(scored_data)
}





# Get lists of predicted drivers for each sample, by topping up canonical drivers with high-scored predictions
topUp_drivers = function(
  gene_scores, 
  canonical_drivers = "example_data/canonical_drivers.rds", 
  entrez_geneSymbol_mapping = "annotation_reference_files/gene_aliases_entrez.tsv",
  n_drivers_per_sample = 5, 
  output_dir = NULL, 
  sample_gene_sep = "__", 
  all_genes = NULL){
  
  require(tidyr)
  require(tibble)
  require(dplyr)
  require(tools)
  
  
  # Extract list of all damaged genes from list object if training denovo (need training and prediction genes altogether)
  if (!is.null(all_genes)){
    all_genes = rbind(all_genes$training_set_ns %>% select(-type), all_genes$prediction_set_ns) %>%
      rownames_to_column("id") %>%
      separate(id, into = c("sample", "entrez"), sep = sample_gene_sep) %>%
      mutate(entrez = as.numeric(entrez)) %>%
      select(sample, entrez)
  } else {
    # If we're using a pre-trained model then all genes will be in the gene_scores data frame
    all_genes = gene_scores %>%
      select(sample, entrez)
  }
  
  
  # For gene_scores, require sample, entrez, and score columns
  gene_scores = gene_scores %>% select(sample, entrez, score)
  
  
  # Join together - some scores may be missing, e.g. for genes in the training set
  all_scores = left_join(all_genes, gene_scores, by = c("sample", "entrez"))
  
  
  # For canonical drivers, can use either
  #  1. A vector of Entrez IDs
  #  2. A data frame with an entrez column
  #  3. A list of Entrez IDs
  # Convert whatever is provided into a vector
  if (is.character(canonical_drivers) & length(canonical_drivers) == 1){
    extension = tools::file_ext(canonical_drivers)
    if (extension == "rds"){
      canonical_drivers = readRDS(canonical_drivers)
    } else {
      canonical_drivers = read_tsv(canonical_drivers)
    }
  }
  if (is.data.frame(canonical_drivers)){
    canonical_drivers = canonical_drivers %>% pull(entrez)
  } else if (is.list(canonical_drivers)){
    canonical_drivers = unlist(canonical_drivers)
  }
  if (!is.numeric(canonical_drivers)) stop("Provided canonical_drivers must be numeric entrez IDs")
  
  
  # Annotate all_scores with canonical driver statuses
  all_scores = all_scores %>% mutate(canonical_driver = entrez %in% canonical_drivers) %>% arrange(sample, desc(score))
  
  
  # Top-up procedure
  topUp = all_scores %>%
    subset(canonical_driver | !is.na(score)) %>%
    group_by(sample) %>%
    mutate(n_canonical_thisSample = sum(canonical_driver),
           less5_flag = n() <= 5,
           n_topUp_thisSample = max(0, n_drivers_per_sample - n_canonical_thisSample),
           rank = rank(-score, na.last = "keep")) %>%
    # If using a pre-trained model, scores will not be NA so handle canonical drivers explicitly
    mutate(
      rank_nonCanonicalOnly = case_when(canonical_driver ~ NA_real_, T ~ rank),
      rank_nonCanonicalOnly = rank(rank_nonCanonicalOnly, na.last = "keep")
      ) %>%
    ungroup %>%
    subset(canonical_driver | rank_nonCanonicalOnly <= n_topUp_thisSample)
  
  
  # Annotate gene symbols
  if (is.character(entrez_geneSymbol_mapping)) entrez_geneSymbol_mapping = read_tsv(entrez_geneSymbol_mapping, col_types = cols())
  entrez_geneSymbol_mapping = entrez_geneSymbol_mapping %>% select(entrez, symbol) %>% unique
  topUp = topUp %>% left_join(entrez_geneSymbol_mapping, by = "entrez")
  
  
  # Save
  if (!is.null(output_dir)) saveRDS(topUp, paste(output_dir, "drivers_toppedUp.rds", sep = "/"))
  
  
  # Output
  return(topUp)
}





##---- DEPRECATED ----
# 
# # Calculate scores based on decision values and mean/variance of sensitivity
# get_scores = function(preds, best_models){
#   
#   scores = data.frame()
#   for (i in 1:nrow(best_models)){
#     # Sensitivity for this kernel
#     k = as.character(best_models$kernel[i])
#     sensitivity_mean = best_models$mean[i]
#     sensitivity_var = best_models$var[i]
#     BMS_i = sensitivity_mean / sensitivity_var
#     
#     
#     # Calculate scores
#     scores_thisKernel = preds %>% 
#       mutate(kernel = k) %>%
#       select(sample, entrez, kernel, decisionValue = paste("decisionValue", k, sep = "_")) %>% 
#       group_by(sample) %>%
#       arrange(desc(decisionValue)) %>%
#       mutate(R_igs = 1:n(), N_s = n()) %>%
#       ungroup %>%
#       mutate(score = -log10(R_igs / N_s) * BMS_i / log10(N_s))
#     scores = rbind(scores, scores_thisKernel)
#   }
#   
#   
#   # Average over kernels
#   scores = scores %>% 
#     group_by(sample, entrez) %>%
#     summarise(score = mean(score), kernels_predicted = paste(kernel[decisionValue > 0], collapse = ";"), kernels_predicted_no = sum(decisionValue > 0)) %>%
#     ungroup
#   return(scores)
# }
# 
# 
# 
# 
# 
# # Rank prediction genes, join to training genes, remove non-expressed (optional)
# get_sysCans = function(scores, preds, normalTissue_expressionData = NULL, normalTissue_name = "Esophagus"){
#   
#   # Join all information
#   syscan = left_join(preds, scores, by = c("sample", "entrez"))
#   
#   
#   # Remove non-expressed genes (optional)
#   if (!is.null(normalTissue_expressionData)){
#     nonExpressed_genes = normalTissue_expressionData %>% 
#       filter(exp_level=="Not expressed" & tissue == normalTissue_name) %>%
#       .$entrez %>%
#       unique
#     syscan = syscan %>% filter(!(entrez %in% nonExpressed_genes))
#   }
#   
#   
#   # Function to rank an ordered vector of Booleans
#   # Increment by 1 for each TRUE, put NA for each FALSE
#   rank_cond = function(x){
#     y = cumsum(x)
#     y[!x] = NA
#     return(y)
#   }
#   
#   
#   # Add patient-specific ranks for non-true positives only
#   syscan = syscan %>%
#     arrange(sample, desc(score)) %>%
#     group_by(sample) %>%
#     mutate(patient_rank = rank_cond(type != "TP")) %>%
#     ungroup
#   
#   
#   return(syscan)
# }
# 
# 
# 
# 
# 
# # Predict and score, cumulatively every step iterations
# score_cumulatively = function(path,
#                               cv_stats = NULL, 
#                               output_dir = NULL,
#                               step = NULL,
#                               normalTissue_expressionData = NULL, 
#                               normalTissue_name = "Esophagus"){
#   
#   require(readr)
#   require(tibble)
#   require(dplyr)
#   require(e1071)
#   cat("Summarising cross-validation\n")
#   
#   
#   # Get cv_stats if it wasn't passed to this function
#   if (is.null(cv_stats)){
#     # Typically quite a large file, use readr's read_tsv function
#     cv_stats = read_tsv(paste(path, "cv_stats.tsv", sep = "/"), col_types = cols(), progress = F) 
#   }
#   
#   
#   # Make output_dir if necessary
#   if (is.null(output_dir)) output_dir = paste(path, "Results", sep = "/")
#   if (!dir.exists(output_dir)) dir.create(output_dir)
#     
#   
#   # Load training and prediction sets
#   training = readRDS(paste(path, "training_set.rds", sep="/"))
#   training_ns = readRDS(paste(path, "training_set_noScale.rds", sep="/"))
#   prediction = readRDS(paste(path, "prediction_set.rds", sep="/"))
#   prediction_ns = readRDS(paste(path, "prediction_set_noScale.rds", sep="/"))
#   featureCols = colnames(prediction_ns)
#   prediction_ns = prediction_ns %>% 
#     tibble::rownames_to_column() %>%
#     separate(rowname, into = c("sample", "entrez"), sep="\\.") %>%
#     mutate(entrez = as.numeric(entrez)) %>%
#     select(sample, entrez, featureCols)
#   
#   
#   # Determine at what steps (of # of iterations) to do scoring
#   iters = max(cv_stats$iteration)
#   if (is.null(step)){
#     steps = iters
#   } else {
#     steps = seq(step, iters, step)
#     if (steps[length(steps)] != iters) steps = c(steps, iters)
#   }
#   
#  
#   # Do scoring at each step
#   for(st in steps){
#     cat("Scoring at", st, "iterations\n")
#     
#     
#     # Write results in a separate directory
#     step_dir = paste0(output_dir, "/", st, "iterations")
#     dir.create(step_dir)
#     
#     
#     # Summarise sensitivity distribution for CV iterations up to st
#     cv_stats_summary = cv_stats %>% 
#       subset(iteration <= st) %>%
#       group_by(kernel, nu, gamma, coef0, degree) %>% 
#       summarize(iterations = n(),
#                 min = min(sensitivity),
#                 q1 = quantile(sensitivity)[2],
#                 median = median(sensitivity),
#                 mean = mean(sensitivity),
#                 q3 = quantile(sensitivity)[4],
#                 max = max(sensitivity),
#                 var = stats::var(sensitivity)) %>% 
#       ungroup
#     write_tsv(cv_stats_summary, path = paste(step_dir, "cv_stats_summary.tsv", sep = "/"))
#     
#     
#     # Get the best model for each kernel
#     # Choose the one with the lowest variance of sensitivity, from among the top 5 by mean sensitivity
#     n_highestMean = min(5, round(iters / 2))
#     best_models = cv_stats_summary %>%
#       group_by(kernel) %>%
#       top_n(n_highestMean, mean) %>%
#       top_n(1, -var) %>%
#       ungroup
#     write_tsv(best_models, path = paste(step_dir, "best_models.tsv", sep = "/"))
#     
#     
#     # Train and predict with each best model
#     preds = rbind(
#       training_ns %>% rownames_to_column() %>% separate(rowname, into = c("sample", "entrez"), sep = "\\.") %>% mutate(entrez = as.numeric(entrez)),
#       prediction_ns %>% mutate(type = "prediction")
#     )
#     for (i in 1:nrow(best_models)){
#       # Extract parameters
#       k = as.character(best_models$kernel[i])
#       mynu = best_models$nu[i]
#       mygamma = best_models$gamma[i]
#       mycoef0 = best_models$coef0[i]
#       mydegree = best_models$degree[i]
#       
#       
#       # Name of model, for saving
#       if (k == "linear"){
#         modelName = paste0(k, "_nu.", mynu)
#       } else if (k == "polynomial"){
#         modelName = paste0(k, "_nu.", mynu, "_gamma.", mygamma, "_coef0.", mycoef0, "_degree.", mydegree)
#       } else if (k %in% c("radial", "sigmoid")){
#         modelName = paste0(k, "_nu.", mynu, "_gamma.", mygamma, "_coef0.", mycoef0)
#       }
#       
#       
#       # Train on full training set
#       svm.model = svm(type ~ ., data = training, scale = FALSE, 
#                       kernel = k, type = "one-classification", 
#                       nu = mynu, gamma = mygamma, coef0 = mycoef0, degree = mydegree)
#       # Save in text format
#       write.svm(svm.model, svm.file = paste0(step_dir, "/model_", modelName, ".svm", sep = ""), scale.file = paste(step_dir, "scale_temp.svm", sep = "/")) 
#       system(command = paste0("rm ", paste0(step_dir,"/scale_temp.svm")))
#       
#       
#       # Apply to everything, including the training set itself
#       # E.g. some training genes might not actually be predicted by some kernels
#       df = rbind(training %>% select(-type), prediction)
#       preds_thisKernel = predict(svm.model, df, decision.values = TRUE)
#       preds_thisKernel = data.frame(key = names(preds_thisKernel), 
#                                     prediction = as.vector(preds_thisKernel), 
#                                     decisionValue = as.vector(attr(preds_thisKernel, "decision.values")), 
#                                     stringsAsFactors = F) %>%
#         separate(key, into=c("sample", "entrez"), sep="\\.") %>%
#         mutate(entrez = as.numeric(entrez))
#       
#       
#       # Join to preds from other kernels
#       colnames(preds_thisKernel)[colnames(preds_thisKernel) == "prediction"] = paste("prediction", k, sep = "_")
#       colnames(preds_thisKernel)[colnames(preds_thisKernel) == "decisionValue"] = paste("decisionValue", k, sep = "_")
#       preds = left_join(preds, preds_thisKernel, by = c("sample", "entrez"))
#     }
#     saveRDS(preds, file = paste(step_dir, "predictions.rds", sep = "/"))
#     
#     
#     # Get the scores for everything
#     scores_all = get_scores(preds, best_models)
#     # And just for the prediction set
#     scores_predictionSet = get_scores(preds %>% subset(type == "prediction"), best_models)
#     # Combine
#     scores = left_join(scores_all %>% rename(score_inclTraining = score), scores_predictionSet, by = c("sample", "entrez", "kernels_predicted", "kernels_predicted_no"))
# 
#     
#     # Remove non-expressed genes (optional), rank non-true positive predictions
#     syscan = get_sysCans(scores, preds, normalTissue_expressionData, normalTissue_name)
#     saveRDS(syscan, file = paste(step_dir, "syscan.rds", sep = "/"))
#   }
#   
#   cat("Finished scoring\n")
# }
# 
# 
# 
# 
