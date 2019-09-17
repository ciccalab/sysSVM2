# Packages
library(readr)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)



# Calculate Jaccard index between two vectors
.JaccInd = function(v1, v2)  length(intersect(v1, v2)) / length(union(v1, v2))



# Plot convergence of the unique top N predictions over iterations
plot_predictionSet_convergence = function(summaryTable, distinct_sysCans, nrow = 2, topN = 10, title = NULL, ylim = NULL, JI_breaks=seq(0, 1, 0.05), nGenes_breaks=seq(0, 2000, 100)){
  
  # Jaccard indexes to compare all predicted gene sets to each other
  pairwise_JI_allGeneSets = expand.grid(geneSet1 = distinct_sysCans$gene_set_name, geneSet2 = distinct_sysCans$gene_set_name, KEEP.OUT.ATTRS = F, stringsAsFactors = F) %>%
    left_join(distinct_sysCans %>% select(geneSet1 = gene_set_name, entrez1 = syscans), by = "geneSet1") %>%
    left_join(distinct_sysCans %>% select(geneSet2 = gene_set_name, entrez2 = syscans), by = "geneSet2")
  pairwise_JI_allGeneSets$index = apply(pairwise_JI_allGeneSets %>% select(entrez1, entrez2), 1, function(x){
    .JaccInd(x[["entrez1"]], x[["entrez2"]])
  })
  pairwise_JI_allGeneSets = pairwise_JI_allGeneSets %>% select(geneSet1, geneSet2, index)

  
  # Compile information into table
  df = summaryTable %>%
    select(reordering, iterations, gene_set_name) %>%
    mutate(n_genes = as.numeric(gsub("[^0-9]", "", .$gene_set_name)))
  for (order_num in unique(df$reordering)){
    vec = df$gene_set_name[df$reordering == order_num]
    df$previous_gene_set[df$reordering == order_num] = c(NA, vec[1:(length(vec) - 1)])
  }
  
  df = df %>%
    left_join(pairwise_JI_allGeneSets %>% dplyr::rename(index_previous = index), 
              by = c("gene_set_name"="geneSet1", "previous_gene_set"="geneSet2"))
  
  
  # Reshape for plotting
  df = df %>%
    mutate(n_genes = n_genes / 2000) %>%
    select(reordering, iterations, Num_genes=n_genes, Jaccard_index_previous_iteration=index_previous) %>%
    gather("quantity", "value", 3:4)
  
  
  # Plot
  p = ggplot(df, aes(x = iterations, y = value, col = quantity)) +
    geom_line() +
    scale_y_continuous("Jaccard index", 
                       breaks = JI_breaks, 
                       sec.axis = sec_axis(~ . * 2000, name = "Number of genes", breaks = nGenes_breaks)) +
    scale_x_continuous(breaks = seq(0, 10000, 2500)) +
    geom_hline(yintercept = 1, linetype = 2) +
    coord_cartesian(ylim = ylim) + 
    labs(x = "Iterations (n)", col = "") +
    facet_wrap(~reordering, nrow = nrow) +
    theme_light(base_size = 16) +
    theme(legend.position = "top")
  
  
  if (!is.null(title)) p = p + labs(title = title) + theme(plot.title = element_text(hjust = 0.5))
  
  print(p)
}



# Plot parameter convergence over iterations
plot_param_convergence = function(summaryTable, logGamma = T, title = NULL, exclude = c("sigmoid_coef0", "polynomial_coef0", "polynomial_gamma")){
  
  # Reshape appropriately
  all_paramColumns = c("linear_nu", "radial_nu", "radial_gamma", "sigmoid_nu", "sigmoid_gamma", "sigmoid_coef0", "polynomial_nu", "polynomial_gamma", "polynomial_coef0", "polynomial_degree")
  my_paramColumns = setdiff(all_paramColumns, exclude)
  df = summaryTable %>% 
    select(Reordering = reordering, iterations, my_paramColumns) %>%
    gather("kernel_parameter", "value", -Reordering, -iterations) %>%
    mutate(kernel = sapply(kernel_parameter, function(x) strsplit(x, split = "_")[[1]][1], USE.NAMES = F),
           parameter = sapply(kernel_parameter, function(x) strsplit(x, split = "_")[[1]][2], USE.NAMES = F)) %>%
    select(Reordering, iterations, kernel, parameter, value)

  
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
  df$Reordering = factor(df$Reordering)
  
  
  # Plot
  p = ggplot(df, aes(x = iterations, y = value, col = Reordering)) +
    geom_line() +
    facet_grid(parameter ~ kernel, scales = "free_y") +
    labs(x = "Iterations (n)", y = "Parameter value") +
    scale_x_continuous(breaks = seq(0, 10000, 2500)) +
    theme_light(base_size = 14) +
    theme(legend.position = "top", panel.spacing.x=unit(1.5, "lines"))
  
  if (!is.null(title)) p = p + labs(title = title) + theme(plot.title = element_text(hjust = 0.5))
  
  print(p)
}



# Sensitivity with parameters
plot_sensitivityVsParams = function(cv_stats, bestModel, title = NULL, logGamma = T, exclude = '((parameter=="coef0") | (kernel=="polynomial" & parameter=="gamma"))'){
  
  # Reshape
  df = cv_stats %>%
    select(iteration, sensitivity, kernel, nu, gamma, coef0, degree) %>%
    gather("parameter", "value", c("nu", "gamma", "coef0", "degree")) %>%
    subset((kernel=="linear" & parameter=="nu") |
             (kernel=="polynomial") |
             (kernel=="radial" & parameter%in%c("nu", "gamma")) |
             (kernel=="sigmoid" & parameter%in%c("nu", "gamma", "coef0"))) %>%
    subset(!eval(parse(text = exclude)))
  
  
  # Best model sensitivity
  bm = bestModel %>%
    select(-iterations, -reordering, -gene_set_name, -frequency, -consecutive) %>%
    gather("kernel_parameter", "value") %>%
    separate(kernel_parameter, into = c("kernel", "parameter"), sep = "_")
  bm = left_join(
    bm %>% subset(!grepl("Sensitivity", parameter)) %>% unique,
    bm %>% subset(parameter == "medSensitivity") %>% select(kernel, sensitivity = value),
    by = "kernel"
  )
  bm = bm %>%
    subset((kernel=="linear" & parameter=="nu") |
             (kernel=="polynomial" & parameter%in%c("nu", "gamma", "coef0", "degree")) |
             (kernel=="radial" & parameter%in%c("nu", "gamma")) |
             (kernel=="sigmoid" & parameter%in%c("nu", "gamma", "coef0"))) %>%
    subset(!eval(parse(text = exclude)))
  
  
  # Use a log scale for gamma
  if (logGamma){
    df[df$parameter == "gamma", "value"] = log2(df[df$parameter == "gamma", "value"])
    df[df$parameter == "gamma", "parameter"] = "gamma (log2)"
    bm[bm$parameter == "gamma", "value"] = log2(bm[bm$parameter == "gamma", "value"])
    bm[bm$parameter == "gamma", "parameter"] = "gamma (log2)"
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
  bm$kernel = capFirst(bm$kernel)
  bm$kernel = factor(bm$kernel, levels = c("Linear", "Polynomial", "Radial", "Sigmoid"))
  bm$parameter = capFirst(bm$parameter)
  if (!logGamma){
    bm$parameter = factor(bm$parameter, levels = c("Nu", "Gamma", "Degree", "Coef0"))
  } else {
    bm$parameter = factor(bm$parameter, levels = c("Nu", "Gamma (log2)", "Degree", "Coef0"))
  }
  
    
  # Plot
  p = ggplot(df, aes(x = factor(value), y = sensitivity)) +
    geom_boxplot() +
    geom_point(data = bm, shape = 4, col = "red") +
    facet_grid(kernel ~ parameter, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Parameter value", y = "Sensitivity") +
    theme_light(base_size = 14) +
    theme(panel.spacing.x=unit(1.5, "lines"))
  
  if (!is.null(title)) p = p + labs(title = title) + theme(plot.title = element_text(hjust = 0.5))
  
  print(p)
}



# Median and variance of sensitivity of best models over iterations
plot_sensitivityIters = function(summaryTable, title = NULL, ylim = list(c(0.8, 1), c(0, 0.005))){
  
  # Reshape
  df = summaryTable %>%
    select(-gene_set_name, -frequency, -consecutive) %>%
    gather("kernel_quantity", "value", -iterations, -reordering) %>%
    separate(kernel_quantity, into = c("kernel", "quantity"), sep = "_") %>%
    subset(quantity %in% c("medSensitivity", "varSensitivity"))
    
  
  # Housekeeping
  capFirst = function(x) {
    firstLetters = substring(x, 1, 1)
    firstLetters = toupper(firstLetters)
    otherLetters = substring(x, 2)
    paste0(firstLetters, otherLetters)
  }
  df$kernel = capFirst(df$kernel)
  df$kernel = factor(df$kernel, levels = c("Linear", "Polynomial", "Radial", "Sigmoid"))
  df = df %>% 
    rename(Reordering = reordering) %>%
    mutate(Reordering = factor(Reordering))
  df$quantity = gsub("medSensitivity", "Sensitivity median", df$quantity)
  df$quantity = gsub("varSensitivity", "Sensitivity variance", df$quantity)
  
  
  # Put in dummy points to make ylim
  dummy = df %>%
    select(kernel, quantity) %>%
    unique
  dummy = rbind(dummy, dummy) %>% 
    arrange(kernel, quantity)
  for (k in unique(dummy$kernel)){
    dummy[dummy$kernel==k & dummy$quantity=="Sensitivity median", "value"] = ylim[[1]]
    dummy[dummy$kernel==k & dummy$quantity=="Sensitivity variance", "value"] = ylim[[2]]
  }
  dummy$iterations = min(df$iterations)
  
  
  # Plot
  p = ggplot(df, aes(x = iterations, y = value, col = Reordering)) +
    geom_line() +
    geom_blank(data = dummy, inherit.aes = F, aes(x = iterations, y = value)) +
    facet_grid(quantity ~ kernel, scales = "free_y") +
    labs(x = "Iterations (n)", y = NULL) +
    theme_light(base_size = 14) +
    theme(legend.position = "top", panel.spacing.x=unit(1.5, "lines"))
  
  if (!is.null(title)) p = p + labs(title = title) + theme(plot.title = element_text(hjust = 0.5))
  
  print(p)
}


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  