##---- Setup ----

library(readr)
library(tibble)
library(tidyr)
library(ggplot2)
library(plotly)
library(corrplot)
library(e1071)
library(dplyr)

setwd("/Volumes/lab-ciccarellif/working/Joel/OAC_sysSVM_2.0/methods_paper/sysSVM_implementation")


##---- Make gene types table ----

# Load files
tsg_og = readRDS("cancerAltTypes_TSG_OG.rds")
geneInfo = read_tsv("../ncg6_allgeneinfo.txt", col_types = cols())


# Combine
geneTypes = geneInfo %>%
  select(entrez, symbol = Symbol, cancer_type) %>%
  left_join(tsg_og %>% select(entrez, tsg_og = geneType), by = "entrez") %>%
  mutate(type = coalesce(tsg_og, cancer_type)) %>%
  select(entrez, symbol, type)
geneTypes$type = gsub("can", "Candidate", geneTypes$type)
geneTypes$type = gsub("cgc", "CGC", geneTypes$type)
geneTypes$type = gsub("oncogene", "Oncogene", geneTypes$type)
geneTypes$type = gsub("rst", "Non-cancer", geneTypes$type)
geneTypes$type = gsub("tumourSuppressor", "TSG", geneTypes$type)
geneTypes$type = factor(geneTypes$type, levels = c("Oncogene", "TSG", "TP53", "CGC", "Candidate", "Non-cancer"))


##---- Plotting functions ----

# Sensitivity with parameters
plot_sensitivityVsParams = function(cv_stats, thisKernel, xParam, yParam = NULL, logx = F, logy = F, size = 300){
  
  # Function to convert a vector to a z-score
  z_score = function(x) (x-mean(x)) / sd(x)
  
  
  if (is.null(yParam)){ # Make a 2D plot - sensitivity vs a single parameter, e.g. linear sensitivity vs nu
    
    # Get mean, sd and diff_z from cv_stats
    df = cv_stats %>%
      subset(kernel == thisKernel) %>%
      select(iteration, xParam = xParam, sensitivity) %>%
      group_by(xParam) %>% 
      summarise(mean = mean(sensitivity), sd = sd(sensitivity)) %>% 
      ungroup %>%
      mutate(diff_z = z_score(mean) - z_score(sd))
    # Label the best model
    df[df$diff_z == max(df$diff_z), "type"] = "Best model"
    df[df$diff_z != max(df$diff_z), "type"] = "Other"
    df$type = factor(df$type, levels = c("Best model", "Other"))

    
    # Log2 transform variables if requested - useful for gamma in particular
    if (logx){
      df$xParam = log2(df$xParam)
      bm$xParam = log2(bm$xParam)
      xParam = paste(xParam, "(log2)")
    }
    
    
    # Make the string distinct from the column in df
    xParam_name = xParam
    rm(xParam)
    
    
    # Capitalise first letter of kernel name
    substr(thisKernel, 1, 1) = toupper(substr(thisKernel, 1, 1))
    
    
    # Y-axis minimum
    sens_min = floor(min(df$mean - df$sd)*10) / 10
    
    
    # Make plot
    p = plot_ly(data = df, 
                x = ~xParam, 
                y = ~mean, 
                type = 'scatter', 
                mode = 'markers',
                marker = list(size = 15),
                error_y = ~list(array = sd, color = '#000000'),
                color = ~diff_z,
                colors = c('#B92B27', '#1565C0'),
                symbol = ~type,
                symbols = c('x', 'circle'),
                hoverinfo = 'text',
                text = ~paste0("</br>", xParam_name, ": ", xParam,
                               "</br>Mean: ", signif(mean, 3),
                               "</br>SD: ", signif(sd, 3),
                               "</br>diff_z: ", signif(diff_z, 3))) %>%
      layout(title = thisKernel,
             xaxis = list(title = xParam_name),
             yaxis = list(title = "Sensitivity", range = c(sens_min, 1)))
    
    print(p)
    
    
  } else { # Make a 3D plot
    
    # Get mean, sd and diff_z from cv_stats
    df = cv_stats %>%
      subset(kernel == thisKernel) %>%
      select(iteration, xParam = xParam, yParam = yParam, sensitivity) %>%
      group_by(xParam, yParam) %>% 
      summarise(mean = mean(sensitivity), sd = sd(sensitivity)) %>% 
      ungroup %>%
      mutate(diff_z = z_score(mean) - z_score(sd))
    # Label the best model
    df[df$diff_z == max(df$diff_z), "type"] = "Best model"
    df[df$diff_z != max(df$diff_z), "type"] = "Other"
    df$type = factor(df$type, levels = c("Best model", "Other"))
    
    
    # Log2 transform variables if requested - useful for gamma in particular
    if (logx){
      df$xParam = log2(df$xParam)
      xParam = paste(xParam, "(log2)")
    }
    if (logy){
      df$yParam = log2(df$yParam)
      yParam = paste(yParam, "(log2)")
    }
    
    
    # Make the strings distinct from the columns in df
    xParam_name = xParam
    yParam_name = yParam
    rm(xParam, yParam)
    
    
    # Capitalise first letter of kernel name
    substr(thisKernel, 1, 1) = toupper(substr(thisKernel, 1, 1))
    
    
    # Z-axis minimum
    sens_min = floor(min(df$mean - df$sd)*10) / 10
    
    
    # Make plot
    p = plot_ly(data = df,
                x = ~xParam, 
                y = ~yParam, 
                z = ~mean, 
                type = 'scatter3d', 
                mode = 'markers',
                size = I(size),
                color = ~diff_z,
                colors = c('#B92B27', '#1565C0'),
                symbol = ~type,
                symbols = c('x', 'circle'),
                error_z = ~list(array = sd, color = '#000000'),
                hoverinfo = 'text',
                text = ~paste0("</br>", xParam_name, ": ", xParam,
                               "</br>", yParam_name, ": ", yParam,
                               "</br>Mean: ", signif(mean, 3),
                               "</br>SD: ", signif(sd, 3),
                               "</br>diff_z: ", signif(diff_z, 3))) %>%
      layout(title = thisKernel,
             scene = list(xaxis = list(title = xParam_name),
                          yaxis = list(title = yParam_name, dtick = 1),
                          zaxis = list(title = "Sensitivity", range = c(sens_min, 1))))
    
    print(p)
  }
  return(df)
}


# Convergence of parameters
plot_paramConvergence = function(cv_stats, step = 1, logGamma = T, title = NULL){
  
  # Function to convert a vector to a z-score
  z_score = function(x) (x-mean(x)) / sd(x)
  
  
  # Cumulatively assess CV stats
  cumul = cv_stats %>%
    arrange(iteration) %>%
    group_by(kernel, nu, gamma, degree, coef0) %>%
    mutate(sens_cumul = cumsum(sensitivity), sens2_cumul = cumsum(sensitivity^2)) %>%
    ungroup %>%
    mutate(mean = sens_cumul / iteration) %>%
    mutate(sd = sqrt((sens2_cumul - iteration*mean^2) / (iteration - 1))) %>%
    group_by(kernel, iteration) %>%
    mutate(mean_z = z_score(mean), sd_z = z_score(sd)) %>%
    mutate(diff_z = mean_z - sd_z) %>%
    mutate(bestModel = diff_z == max(diff_z)) %>%
    ungroup %>%
    subset(iteration %% step == 0)

  
  # Identify which parameters were tuned for which kernels
  params_tuned = list()
  for (k in unique(cv_stats$kernel)){
    params_tuned[[k]] = names(which(apply(cv_stats %>% subset(kernel == k) %>% select(nu, gamma, degree, coef0), 2, n_distinct) > 1))
  }
  
  
  # Reshape for parameter convergence
  df = cumul %>%
    subset(bestModel) %>%
    select(kernel, nu, gamma, degree, coef0, iteration) %>%
    gather("parameter", "value", c(nu, gamma, degree, coef0))
  for (k in unique(cv_stats$kernel)) df = df %>% subset(kernel != k | parameter %in% params_tuned[[k]])
  
  
  # Use a log scale for gamma
  if (logGamma){
    df[df$parameter == "gamma", "value"] = log2(df[df$parameter == "gamma", "value"])
    df[df$parameter == "gamma", "parameter"] = "gamma (log2)"
  }
  
  
  # Housekeeping  
  capFirst = function(x) {
    firstLetters = substring(x, 1, 1)
    firstLetters = toupper(firstLetters)
    otherLetters = substring(x, 2)
    paste0(firstLetters, otherLetters)
  }
  df$kernel = capFirst(df$kernel)
  df$kernel = factor(df$kernel, levels = c("Linear", "Polynomial", "Radial", "Sigmoid"))
  df$parameter = capFirst(df$parameter)
  if (!logGamma){
    df$parameter = factor(df$parameter, levels = c("Nu", "Gamma", "Degree", "Coef0"))
  } else {
    df$parameter = factor(df$parameter, levels = c("Nu", "Gamma (log2)", "Degree", "Coef0"))
  }
  
  
  # Plot
  p = ggplot(df %>% subset(iteration %% step == 0), aes(x = iteration, y = value)) +
    geom_line() +
    facet_grid(parameter ~ kernel, scales = "free_y") +
    labs(x = "Iterations (n)", y = "Parameter value", title = title) +
    scale_x_continuous(breaks = seq(0, 10000, 2500)) +
    theme_light(base_size = 14) +
    theme(legend.position = "top", panel.spacing.x=unit(1.5, "lines"), plot.title = element_text(hjust = 0.5))
  
  print(p)
  return(df)
}


# Convergence of underlying metrics (one plot per kernel)
plot_metricConvergence = function(cv_stats, step = 100, top_n = 4){
  
  # Function to convert a vector to a z-score
  z_score = function(x) (x-mean(x)) / sd(x)
  
  
  # Use a log scale for gamma
  cv_stats = cv_stats %>% mutate(logGamma = log2(gamma))
  
  
  # Cumulatively assess CV stats
  cumul = cv_stats %>%
    arrange(iteration) %>%
    group_by(kernel, nu, logGamma, degree, coef0) %>%
    mutate(sens_cumul = cumsum(sensitivity), sens2_cumul = cumsum(sensitivity^2)) %>%
    ungroup %>%
    mutate(Mean = sens_cumul / iteration) %>%
    mutate(SD = sqrt((sens2_cumul - iteration*Mean^2) / (iteration - 1))) %>%
    group_by(kernel, iteration) %>%
    mutate(mean_z = z_score(Mean), sd_z = z_score(SD)) %>%
    mutate(diff_z = mean_z - sd_z) %>%
    ungroup %>%
    subset(iteration %% step == 0)
  
  
  # Identify which parameters were tuned for which kernels
  params_tuned = list()
  for (k in unique(cv_stats$kernel)){
    params_tuned[[k]] = names(which(apply(cv_stats %>% subset(kernel == k) %>% select(nu, logGamma, degree, coef0), 2, n_distinct) > 1))
  }
  
  
  # Reshape for metric convergence
  df = cumul
  for (k in unique(cv_stats$kernel)){
    if (length(params_tuned[[k]]) == 1){
      df[df$kernel == k, "model"] = paste0(params_tuned[[k]][1], "_", df %>% subset(kernel == k) %>% .[[params_tuned[[k]][1]]])
    } else if (length(params_tuned[[k]]) == 2){
      df[df$kernel == k, "model"] = paste0(params_tuned[[k]][1], "_", df %>% subset(kernel == k) %>% .[[params_tuned[[k]][1]]], "__",
                                           params_tuned[[k]][2], "_", df %>% subset(kernel == k) %>% .[[params_tuned[[k]][2]]])
    } else {
      df[df$kernel == k, "model"] = paste0(params_tuned[[k]][1], "_", df %>% subset(kernel == k) %>% .[[params_tuned[[k]][1]]], "__",
                                           params_tuned[[k]][2], "_", df %>% subset(kernel == k) %>% .[[params_tuned[[k]][2]]], "__",
                                           params_tuned[[k]][3], "_", df %>% subset(kernel == k) %>% .[[params_tuned[[k]][3]]])
    }
  }
  df = df %>%
    select(kernel, Model = model, iteration, Mean, SD, diff_z) %>%
    gather("metric", "value", c(Mean, SD, diff_z))

  
  # Housekeeping  
  capFirst = function(x) {
    firstLetters = substring(x, 1, 1)
    firstLetters = toupper(firstLetters)
    otherLetters = substring(x, 2)
    paste0(firstLetters, otherLetters)
  }
  df$kernel = capFirst(df$kernel)
  df$kernel = factor(df$kernel, levels = c("Linear", "Polynomial", "Radial", "Sigmoid"))
  df$metric = factor(df$metric, levels = c("Mean", "SD", "diff_z"))
  
  
  # Annotate the top models (at the end of all CV iterations)
  top_models = df %>%
    subset(iteration == max(iteration) & metric == "diff_z") %>%
    group_by(kernel) %>%
    mutate(rank = n() - rank(value) + 1) %>%
    mutate(top_models = rank <= top_n) %>%
    ungroup %>%
    select(kernel, Model, rank, top_models)
  df = df %>% left_join(top_models, by = c("kernel", "Model"))
  
  
  # One plot per kernel
  for (k in unique(df$kernel)){
    p = ggplot(df %>% subset(kernel == k & !top_models), aes(x = iteration, y = value, group = Model)) +
      geom_line(color = "gray80") +
      geom_line(data = df %>% 
                  subset(kernel == k & top_models) %>%
                  arrange(rank) %>%
                  mutate(Model = factor(Model, levels = unique(Model))), 
                aes(col = Model)) +
      facet_wrap(~metric, scales = "free_y", nrow = 3) +
      labs(x = "Iterations (n)", y = "Value", title = k) +
      scale_x_continuous(breaks = seq(0, 10000, 2500)) +
      theme_light(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_brewer(type = "qual", palette = 2)
    
    print(p)
  }
  return(df)
}


# Correlation of kernel decision values
plot_interKernel_corr = function(scores){
  
  # Reshape
  decVal_corr = scores %>%
    select(sample, entrez, Linear=linear_decision_value, Polynomial=polynomial_decision_value, Radial=radial_decision_value, Sigmoid=sigmoid_decision_value) %>% 
    unite("id", sample, entrez) %>%
    column_to_rownames("id")
  

  # Get correlation values and p-values
  correlation_matrix = cor(decVal_corr)
  correlation_matrix_p_values = cor.mtest(decVal_corr, conf.level = .99)
  
  
  # Plot
  colours = colorRampPalette(c("red", "white", "white", "blue"))
  corrplot(correlation_matrix, 
           method="circle", 
           type = "upper", 
           p.mat = correlation_matrix_p_values$p,
           sig.level = .05, 
           insig = "blank",
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 1.2,
           na.label = "square", 
           na.label.col = "lightgrey",
           col = colours(20),
           diag=FALSE,
           addCoef.col= "black",
           mar = c(0.2, 0, 0.2, 0) + 0.1)
  
  return(correlation_matrix)
}


# Feature importance (RFE)
# Function to do RFE (from Thanos)
plot_rfe = function(training_set, kernel, best_models = NULL, mynu = 0.05, mygamma = 1, mydegree = 2, mycoef0 = 0, title = NULL){

  # Prepare data
  x = training_set %>% select(-type) %>% data.matrix()
  x = x[,colSums(!is.na(x)) > 0]
  y = training_set %>% .$type
  n = ncol(x)
  
  
  # Initialisations
  survivingFeaturesIndexes = seq(1:n)
  featureRankedList = vector(length=n)
  rankedFeatureIndex = n
  weight_list = list()
  count = 1
  
  
  # If best_models is provided, get parameter values from it
  thisKernel = kernel
  if (!is.null(best_models)){
    mynu = best_models %>% subset(kernel == thisKernel) %>% .$nu
    mygamma = best_models %>% subset(kernel == thisKernel) %>% .$gamma
    mydegree = best_models %>% subset(kernel == thisKernel) %>% .$degree
    mycoef0 = best_models %>% subset(kernel == thisKernel) %>% .$coef0
  }
  
  
  # Recursive feature elimination
  while(length(survivingFeaturesIndexes) > 0){
    if (length(survivingFeaturesIndexes) %% 5 == 0) message(length(survivingFeaturesIndexes), " features remaining")
    
    # Train the support vector machine
    svmModel = svm(x[, survivingFeaturesIndexes], y, kernel = kernel, type = "one-classification", scale = F, 
                   nu = mynu, gamma = mygamma, degree = mydegree, coef0 = mycoef0)

    # Compute the weight vector
    w = t(svmModel$coefs) %*% svmModel$SV
    
    # Compute ranking criteria
    rankingCriteria = w * w
    weight_list[[as.character(count)]] = rankingCriteria
    # Rank the features
    ranking = sort(rankingCriteria, index.return = T)$ix
    
    # Update feature ranked list
    featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
    rankedFeatureIndex = rankedFeatureIndex - 1
    
    # Eliminate the feature with smallest ranking criterion
    survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]]
    count = count + 1
  }
  
  
  # Make data frame for plotting
  capFirst = function(x) {
    firstLetters = substring(x, 1, 1)
    firstLetters = toupper(firstLetters)
    otherLetters = substring(x, 2)
    paste0(firstLetters, otherLetters)
  }
  df = data.frame()
  for (i in 1:(length(weight_list)-1)){
    df = rbind(df, 
               t(weight_list[[i]]) %>% as.data.frame %>% rownames_to_column("feature") %>% rename(w2 = V1) %>% mutate(elim_n = i))
  }
  df$feature = factor(df$feature, levels = rev(colnames(weight_list[[1]])[featureRankedList]))
  
  
  # Plot
  p = ggplot(df, aes(x = feature, y = w2)) + 
    geom_boxplot() +
    labs(x = "Feature", y = "Squared weight", title = capFirst(kernel)) +
    theme_light(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = element_text(hjust = 0.5)) +
    coord_flip()
  print(p)
  
  
  # Output data
  return(list(featureRankedList = featureRankedList, weights = weight_list))
}


# Score distributions for different gene types
plot_scoreByType = function(scores, geneTypes, score_column = "score_all"){
  
  # Join information
  df = scores %>%
    select(sample, entrez, score = score_column) %>%
    subset(!is.na(score)) %>%
    left_join(geneTypes, by = "entrez")
  
  
  # Plot
  p = ggplot(df, aes(type, score, fill = type)) +
    geom_boxplot() +
    stat_summary(aes(label=round(..y.., 1)), fun.y=min, geom="text", size=4, position = position_nudge(y = -20)) +
    stat_summary(aes(label=round(..y.., 1)), fun.y=max, geom="text", size=4, position = position_nudge(y = 20)) +
    stat_summary(aes(label=round(..y.., 1)), fun.y=median, geom="label", size=4, fill="white") +
    labs(x = NULL, y = "sysSVM score") +
    scale_fill_brewer(type = "qual", palette = 4) +
    theme_light(base_size = 14) +
    theme(legend.position = "none", panel.grid = element_blank())
  
  print(p)
  
  return(df)
}


# Patient-specific rank distributions for different gene types
plot_rankByType = function(scores, geneTypes, score_column = "score_all"){
  
  # Join information
  df = scores %>%
    select(sample, entrez, score = score_column) %>%
    subset(!is.na(score)) %>%
    left_join(geneTypes, by = "entrez")
  
  
  # Calculate patient-specific ranks
  df = df %>%
    arrange(desc(score)) %>%
    group_by(sample) %>%
    mutate(patient_rank = 1:n()) %>%
    ungroup

 
  # Plot
  p = ggplot(df, aes(type, patient_rank, fill = type)) +
    geom_boxplot() +
    stat_summary(aes(label=round(..y.., 1)), fun.y=min, geom="text", size=4, position = position_nudge(y = -20)) +
    stat_summary(aes(label=round(..y.., 1)), fun.y=max, geom="text", size=4, position = position_nudge(y = 20)) +
    stat_summary(aes(label=round(..y.., 1)), fun.y=median, geom="label", size=4, fill="white") +
    labs(x = NULL, y = "Patient-specific rank") +
    scale_fill_brewer(type = "qual", palette = 4) +
    theme_light(base_size = 14) +
    theme(legend.position = "none", panel.grid = element_blank())
  
  print(p)
  
  return(df)
}


##---- Make plots ----

# Directory of a sysSVM implementation
wdir = "noMeanScale/10000.iterations_in_20.batches"
wdir = "originalSettings/10000.iterations_in_20.batches"
wdir = "rmCorrFeatures_noMeanScale/10000.iterations_in_25.batches"
wdir = "noBinaryFeatures_noMeanScale/10000.iterations_in_25.batches"


# Load the necessary files
training_set = readRDS(paste(wdir, "training_set.rds", sep = "/"))
cv_stats = read_tsv(paste(wdir, "cv_stats.tsv", sep = "/"), col_types = cols())
best_models = read_tsv(paste(wdir, "best_model_final.tsv", sep = "/"), col_types = cols())
scores = readRDS(paste(wdir, "scores.rds", sep = "/"))


# Sensitivity distribution after all iterations
sensitivityVsParams_linear = plot_sensitivityVsParams(cv_stats = cv_stats, thisKernel = "linear", xParam = "nu")
sensitivityVsParams_polynomial = plot_sensitivityVsParams(cv_stats = cv_stats, thisKernel = "polynomial", xParam = "nu", yParam = "degree")
sensitivityVsParams_radial = plot_sensitivityVsParams(cv_stats = cv_stats, thisKernel = "radial", xParam = "nu", yParam = "gamma", logy = T)
sensitivityVsParams_sigmoid = plot_sensitivityVsParams(cv_stats = cv_stats, thisKernel = "sigmoid", xParam = "nu", yParam = "gamma", logy = T)


# Convergence of parameters
parameterConvergence = plot_paramConvergence(cv_stats, step = 1, logGamma = T, title = NULL)


# Convergence of underlying metrics
metricConvergence = plot_metricConvergence(cv_stats, step = 100, top_n = 4)


# Correlation of kernel decision values
interKernel_corr = plot_interKernel_corr(scores)


# Feature importance (RFE)
rfe_linear = plot_rfe(training_set = training_set, kernel = "linear", best_models = best_models)
rfe_polynomial = plot_rfe(training_set = training_set, kernel = "polynomial", best_models = best_models)
rfe_radial = plot_rfe(training_set = training_set, kernel = "radial", best_models = best_models)
rfe_sigmoid = plot_rfe(training_set = training_set, kernel = "sigmoid", best_models = best_models)


# Score distributions for different gene types
scores_byType_all = plot_scoreByType(scores, geneTypes, score_column = "score_all")
scores_byType_predOnly = plot_scoreByType(scores, geneTypes, score_column = "score_predOnly")


# Rank distributions for different gene types
ranks_byType_all = plot_rankByType(scores, geneTypes, score_column = "score_all")
ranks_byType_predOnly = plot_rankByType(scores, geneTypes, score_column = "score_predOnly")
