# These functions are used to train the one-class SVM, and make predictions using a set of trained models


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
get_scaledTable = function(df, scaling_factors = NULL, scaleMean = F){
  
  # If scaling_factors are not provided (i.e. we are training rather than predicting), we need to determine them first
  if (is.null(scaling_factors)){
    # Only scale numeric features
    factorCols = names(which(sapply(df, function(x) class(x)=="factor")))
    cols_toScale = setdiff(colnames(df), c("sample", "entrez", factorCols))
    
    # Calculate mean and standard deviation
    scaling_factors = list()
    for (col in cols_toScale){
      scaling_factors[[col]] = c(
        mean_translation = ifelse(scaleMean, mean(df[[col]]), 0), # Will only translate for zero mean if scaleMean = TRUE
        sd_scaling = sd(df[[col]])
      )
    }
  }
  
  
  # Scale the data
  df_scaled = df
  for (col in names(scaling_factors)){
    df_scaled[[col]] = (df_scaled[[col]] - scaling_factors[[col]]["mean_translation"]) / scaling_factors[[col]]["sd_scaling"]
  }
  return(list(df_scaled = df_scaled, scaling_factors = scaling_factors))
}





# Prepare training and prediction sets
prepare_trainingPrediction = function(path = NULL,
                                      create_newDirectory = T,
                                      cohort_data,
                                      normalTissue_expressionData = NULL,
                                      normalTissue_name = "Esophagus",
                                      exclude_features = NULL,
                                      nonFactor_features = c("sample", "entrez", "no_ALL_muts", "no_NSI_muts", "no_NTDam_muts", "no_TRUNC_muts", "no_GOF_muts", "BND", "INV", "INS"),
                                      truePositive_drivers,
                                      damaging_subsetCond = "no_TRUNC_muts + no_NTDam_muts + no_GOF_muts > 0 | CNVGain == 1 | Copy_number == 0",
                                      scaleMean = F){
  
  require(tibble)
  require(tidyr)
  require(dplyr)
  cat("Preparing training and prediction sets\n")
  df = cohort_data
  
  
  # Remove non-expressed genes if normal tissue expression data (e.g. GTEx) are provided
  if (!is.null(normalTissue_expressionData)){
    nonExpressed_genes = normalTissue_expressionData %>% 
      filter(exp_level=="Not expressed" & tissue == normalTissue_name) %>%
      .$entrez %>%
      unique
    
    cat("Removing", length(intersect(df$entrez, nonExpressed_genes)), "genes not expressed in normal", normalTissue_name, "\n")
    df = df %>% filter(!(entrez %in% nonExpressed_genes))
  }
  
  
  # Remove any features that we want to exclude
  if (!is.null(exclude_features)) cat("Excluding features:", paste(exclude_features, collapse = ", "), "\n")
  features = setdiff(colnames(df), exclude_features)
  df = df %>% select(one_of(c("sample", "entrez", features)))

  
  # Convert features with only two unique values to factors, unless otherwise specified
  # These features will not undergo any scaling
  twoLevel_cols = names(which(sapply(df, function(x) n_distinct(x) == 2)))
  factorCols = setdiff(twoLevel_cols, nonFactor_features)
  cat("Converting binary features to factors:", paste(factorCols, collapse = ", "), "\n")
  for (col in factorCols) df[[col]] = factor(df[[col]])
  
  
  # Remove all alterations that are not determined to be damaging, unless damaging_subsetCond is null
  if (!is.null(damaging_subsetCond)){
    n_nonDamaging_rows = df %>% filter(!eval(parse(text = damaging_subsetCond))) %>% nrow
    cat("Removing", n_nonDamaging_rows, "alterations which are non-damaging\n")
    df = df %>% filter(eval(parse(text = damaging_subsetCond)))
  }
  
  
  # Scale data
  if (scaleMean){
    cat("Scaling data; mean and standard deviation\n")
  } else {
    cat("Scaling data; standard deviation only\n")
  }
  scaled = get_scaledTable(df = df, scaleMean = scaleMean)
  df_scaled = scaled$df_scaled
  scaling_factors = scaled$scaling_factors
  rm(scaled)
  
  
  # Mark which alterations are "true positive" drivers, and separate into training and prediction sets
  cat("Separating training and prediction sets\n")
  df = mark_truePositives(df = df, truePositive_drivers = truePositive_drivers)
  df_scaled$type = df$type
  # Unscaled data
  df = df %>%
    subset(!(type == "prediction" & entrez %in% truePositive_drivers$entrez)) %>%
    mutate(key = paste(sample, entrez, sep = ".")) %>%
    column_to_rownames(var = "key") %>%
    select(-sample, -entrez)
  training_ns = df %>% subset(type == "TP")
  prediction_ns = df %>% subset(type == "prediction") %>% select(-type)
  # Scaled data
  df_scaled = df_scaled %>%
    subset(!(type == "prediction" & entrez %in% truePositive_drivers$entrez)) %>%
    mutate(key = paste(sample, entrez, sep = ".")) %>%
    column_to_rownames(var = "key") %>%
    select(-sample, -entrez)
  training = df_scaled %>% subset(type == "TP")
  prediction = df_scaled %>% subset(type == "prediction") %>% select(-type)
  
  
  # Save scaling factors, training and prediction sets
  if (create_newDirectory) dir.create(path)
  saveRDS(scaling_factors, file = paste(path, "scaling_factors.rds", sep = "/"))
  saveRDS(training, file = paste(path, "training_set.rds", sep = "/"))
  saveRDS(training_ns, file = paste(path, "training_set_noScale.rds", sep = "/"))
  saveRDS(prediction, file = paste(path, "prediction_set.rds", sep = "/"))
  saveRDS(prediction_ns, file = paste(path, "prediction_set_noScale.rds", sep = "/"))
  cat("Dataset divided into", nrow(training), "training observations and", nrow(prediction), "observations for prediction\n")
}





# Carry out cross validation iterations, recording the sensitivity
# Runs in a parallel environment
run_crossValidation_par = function(path,
                                   folds = 3, 
                                   iters = 100,
                                   kernels_paramGrids = list(
                                     linear = list(nu = seq(0.05, 0.35, 0.05)),
                                     polynomial = list(nu = seq(0.05, 0.35, 0.05), gamma = 1, coef0 = 0, degree = c(3, 4)),
                                     radial = list(nu = seq(0.05, 0.35, 0.05), gamma = 2^seq(-7, 4)),
                                     sigmoid = list(nu = seq(0.05, 0.35, 0.05), gamma = 2^seq(-7, 4), coef0 = 1)
                                   ),
                                   cores = 2, 
                                   verbose = F,
                                   parallelLib = "snow"
                                   ){
  
  require(parallel)
  require(snow)
  require(doSNOW)
  require(tibble)
  require(dplyr)
  require(e1071)
  
  
  ##---- Sort out directory structure for verbose mode ----
  
  # In verbose mode, we will save the results of every iteration for every kernel-parameter combination
  # Create the directories required for this
  if (verbose){
    # Containing directory, called CV
    cv_dir = paste(path, "CV", sep = "/")
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
  
  # Get the training data from path
  # To avoid confusion between this and the subsets that will actually be used to train in each CV iteration, call this truePositives
  truePositives = readRDS(paste(path, "training_set.rds", sep = "/"))
  # Unique genes in the true postitive set
  truePositive_genes = truePositives %>% 
    tibble::rownames_to_column() %>% 
    separate(rowname, into=c("sample", "entrez"), sep="\\.") %>% 
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
        trainset = truePositives %>% filter(sapply(rownames(.), function(key) strsplit(key, split = "\\.")[[1]][2]) %in% traingenes)
        testset = truePositives %>% filter(!sapply(rownames(.), function(key) strsplit(key, split = "\\.")[[1]][2]) %in% traingenes)
        

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
  cv_stats_fn = paste(path, "cv_stats.tsv", sep="/")
  write_tsv(cv_stats_full, path = cv_stats_fn)
  cat("CV results written to", cv_stats_fn, "\n")
}





# Calculate scores based on decision values and mean/variance of sensitivity
get_scores = function(preds, best_models){
  
  scores = data.frame()
  for (i in 1:nrow(best_models)){
    # Sensitivity for this kernel
    k = as.character(best_models$kernel[i])
    sensitivity_mean = best_models$mean[i]
    sensitivity_var = best_models$var[i]
    BMS_i = sensitivity_mean / sensitivity_var
    
    
    # Calculate scores
    scores_thisKernel = preds %>% 
      mutate(kernel = k) %>%
      select(sample, entrez, kernel, decisionValue = paste("decisionValue", k, sep = "_")) %>% 
      group_by(sample) %>%
      arrange(desc(decisionValue)) %>%
      mutate(R_igs = 1:n(), N_s = n()) %>%
      ungroup %>%
      mutate(score = -log10(R_igs / N_s) * BMS_i / log10(N_s))
    scores = rbind(scores, scores_thisKernel)
  }
  
  
  # Average over kernels
  scores = scores %>% 
    group_by(sample, entrez) %>%
    summarise(score = mean(score), kernels_predicted = paste(kernel[decisionValue > 0], collapse = ";"), kernels_predicted_no = sum(decisionValue > 0)) %>%
    ungroup
  return(scores)
}





# Rank prediction genes, join to training genes, remove non-expressed (optional)
get_sysCans = function(scores, preds, normalTissue_expressionData = NULL, normalTissue_name = "Esophagus"){
  
  # Join all information
  syscan = left_join(preds, scores, by = c("sample", "entrez"))
  
  
  # Remove non-expressed genes (optional)
  if (!is.null(normalTissue_expressionData)){
    nonExpressed_genes = normalTissue_expressionData %>% 
      filter(exp_level=="Not expressed" & tissue == normalTissue_name) %>%
      .$entrez %>%
      unique
    syscan = syscan %>% filter(!(entrez %in% nonExpressed_genes))
  }
  
  
  # Function to rank an ordered vector of Booleans
  # Increment by 1 for each TRUE, put NA for each FALSE
  rank_cond = function(x){
    y = cumsum(x)
    y[!x] = NA
    return(y)
  }
  
  
  # Add patient-specific ranks for non-true positives only
  syscan = syscan %>%
    arrange(sample, desc(score)) %>%
    group_by(sample) %>%
    mutate(patient.rank = rank_cond(type != "TP")) %>%
    ungroup
  
  
  return(syscan)
}





# Predict and score, cumulatively every step iterations
score_cumulatively = function(path,
                              cv_stats = NULL, 
                              output_dir = NULL,
                              step = NULL,
                              normalTissue_expressionData = NULL, 
                              normalTissue_name = "Esophagus"){
  
  require(readr)
  require(tibble)
  require(dplyr)
  require(e1071)
  cat("Summarising cross-validation\n")
  
  
  # Get cv_stats if it wasn't passed to this function
  if (is.null(cv_stats)){
    # Typically quite a large file, use readr's read_tsv function
    cv_stats = read_tsv(paste(path, "cv_stats.tsv", sep = "/"), col_types = cols(), progress = F) 
  }
  
  
  # Make output_dir if necessary
  if (is.null(output_dir)) output_dir = paste(path, "Results", sep = "/")
  if (!dir.exists(output_dir)) dir.create(output_dir)
    
  
  # Load training and prediction sets
  training = readRDS(paste(path, "training_set.rds", sep="/"))
  training_ns = readRDS(paste(path, "training_set_noScale.rds", sep="/"))
  prediction = readRDS(paste(path, "prediction_set.rds", sep="/"))
  prediction_ns = readRDS(paste(path, "prediction_set_noScale.rds", sep="/"))
  featureCols = colnames(prediction_ns)
  prediction_ns = prediction_ns %>% 
    tibble::rownames_to_column() %>%
    separate(rowname, into = c("sample", "entrez"), sep="\\.") %>%
    mutate(entrez = as.numeric(entrez)) %>%
    select(sample, entrez, featureCols)
  
  
  # Determine at what steps (of # of iterations) to do scoring
  iters = max(cv_stats$iteration)
  if (is.null(step)){
    steps = iters
  } else {
    steps = seq(step, iters, step)
    if (steps[length(steps)] != iters) steps = c(steps, iters)
  }
  
 
  # Do scoring at each step
  for(st in steps){
    cat("Scoring at", st, "iterations\n")
    
    
    # Write results in a separate directory
    step_dir = paste0(output_dir, "/", st, "iterations")
    dir.create(step_dir)
    
    
    # Summarise sensitivity distribution for CV iterations up to st
    cv_stats_summary = cv_stats %>% 
      subset(iteration <= st) %>%
      group_by(kernel, nu, gamma, coef0, degree) %>% 
      summarize(iterations = n(),
                min = min(sensitivity),
                q1 = quantile(sensitivity)[2],
                median = median(sensitivity),
                mean = mean(sensitivity),
                q3 = quantile(sensitivity)[4],
                max = max(sensitivity),
                var = stats::var(sensitivity)) %>% 
      ungroup
    write_tsv(cv_stats_summary, path = paste(step_dir, "cv_stats_summary.tsv", sep = "/"))
    
    
    # Get the best model for each kernel
    # Choose the one with the lowest variance of sensitivity, from among the top 5 by mean sensitivity
    n_highestMean = min(5, round(iters / 2))
    best_models = cv_stats_summary %>%
      group_by(kernel) %>%
      top_n(n_highestMean, mean) %>%
      top_n(1, -var) %>%
      ungroup
    write_tsv(best_models, path = paste(step_dir, "best_models.tsv", sep = "/"))
    
    
    # Train and predict with each best model
    preds = rbind(
      training_ns %>% rownames_to_column() %>% separate(rowname, into = c("sample", "entrez"), sep = "\\.") %>% mutate(entrez = as.numeric(entrez)),
      prediction_ns %>% mutate(type = "prediction")
    )
    for (i in 1:nrow(best_models)){
      # Extract parameters
      k = as.character(best_models$kernel[i])
      mynu = best_models$nu[i]
      mygamma = best_models$gamma[i]
      mycoef0 = best_models$coef0[i]
      mydegree = best_models$degree[i]
      
      
      # Name of model, for saving
      if (k == "linear"){
        modelName = paste0(k, "_nu.", mynu)
      } else if (k == "polynomial"){
        modelName = paste0(k, "_nu.", mynu, "_gamma.", mygamma, "_coef0.", mycoef0, "_degree.", mydegree)
      } else if (k %in% c("radial", "sigmoid")){
        modelName = paste0(k, "_nu.", mynu, "_gamma.", mygamma, "_coef0.", mycoef0)
      }
      
      
      # Train on full training set
      svm.model = svm(type ~ ., data = training, scale = FALSE, 
                      kernel = k, type = "one-classification", 
                      nu = mynu, gamma = mygamma, coef0 = mycoef0, degree = mydegree)
      # Save in text format
      write.svm(svm.model, svm.file = paste0(step_dir, "/model_", modelName, ".svm", sep = ""), scale.file = paste(step_dir, "scale_temp.svm", sep = "/")) 
      system(command = paste0("rm ", paste0(step_dir,"/scale_temp.svm")))
      
      
      # Apply to everything, including the training set itself
      # E.g. some training genes might not actually be predicted by some kernels
      df = rbind(training %>% select(-type), prediction)
      preds_thisKernel = predict(svm.model, df, decision.values = TRUE)
      preds_thisKernel = data.frame(key = names(preds_thisKernel), 
                                    prediction = as.vector(preds_thisKernel), 
                                    decisionValue = as.vector(attr(preds_thisKernel, "decision.values")), 
                                    stringsAsFactors = F) %>%
        separate(key, into=c("sample", "entrez"), sep="\\.") %>%
        mutate(entrez = as.numeric(entrez))
      
      
      # Join to preds from other kernels
      colnames(preds_thisKernel)[colnames(preds_thisKernel) == "prediction"] = paste("prediction", k, sep = "_")
      colnames(preds_thisKernel)[colnames(preds_thisKernel) == "decisionValue"] = paste("decisionValue", k, sep = "_")
      preds = left_join(preds, preds_thisKernel, by = c("sample", "entrez"))
    }
    saveRDS(preds, file = paste(step_dir, "predictions.rds", sep = "/"))
    
    
    # Get the scores for the prediction set
    scores = get_scores(preds, best_models)
    # saveRDS(scores, file = paste(step_dir, "scores.rds", sep = "/"))  # syscan now has all the score information
    
    
    # Remove non-expressed genes (optional), rank non-true positive predictions
    syscan = get_sysCans(scores, preds, normalTissue_expressionData, normalTissue_name)
    saveRDS(syscan, file = paste(step_dir, "syscan.rds", sep = "/"))
  }
  
  cat("Finished scoring\n")
}





# Predict on new data, using pre-trained sysSVM models
predict_sysSVM = function(output_dir,
                          new_data,
                          normalTissue_expressionData = NULL,
                          normalTissue_name = "Esophagus",
                          exclude_features = NULL,
                          nonFactor_features = c("sample", "entrez", "no_ALL_muts", "no_NSI_muts", "no_NTDam_muts", "no_TRUNC_muts", "no_GOF_muts", "BND", "INV", "INS"),
                          scaling_factors,
                          trained_models){
  
  require(tidyr)
  require(dplyr)
  cat("Preparing data for sysSVM\n")
  df = new_data
  
  
  # Remove non-expressed genes if normal tissue expression data (e.g. GTEx) are provided
  if (!is.null(normalTissue_expressionData)){
    nonExpressed_genes = normalTissue_expressionData %>% 
      filter(exp_level=="Not Expressed" & tissue == normalTissue_name) %>%
      .$entrez %>%
      unique
    
    cat("Removing", length(nonExpressed_genes), "genes not expressed in normal", normalTissue_name, "\n")
    df = df %>% filter(!(entrez %in% nonExpressed_genes))
  }
  
  
  # Remove any features that we want to exclude
  if (!is.null(exclude_features)) cat("Excluding features:", paste(exclude_features, collapse = ", "))
  features = setdiff(colnames(df), exclude_features)
  df = df %>% select(one_of(c("sample", "entrez", features)))
  
  
  # Convert features with only two unique values to factors, unless otherwise specified
  # These features will not undergo any scaling
  twoLevel_cols = names(which(sapply(df, function(x) n_distinct(x) == 2)))
  factorCols = setdiff(twoLevel_cols, nonFactor_features)
  cat("Converting binary features to factors:", paste(factorCols, collapse = ", "))
  for (col in factorCols) df[[col]] = factor(df[[col]])
  
  
  # Scale the new data
  cat("Scaling data\n")
  scaled_df = getScaledTable(df = df, scaling.factors = scaling_factors)
  
  
  # Save the scaled data
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }else{
    stop("Output directory exists, can't overwrite")
  }
  save(df, file = paste(output_dir,"prediction_set_noScale.Rdata", sep="/"))
  save(scaled_df[["df_scaled"]], file = paste(output_dir, "prediction_set.Rdata", sep="/"))
  
  
  # Make predictions
  
}




