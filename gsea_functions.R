##---- Packages ----

library(readr)
library(tidyr)
library(tibble)
library(dplyr)
library(gplots)
library(ggplot2)


##---- Define functions ----

# Do a gene set enrichment analysis of sysSVM predictions
# Can use the true positive genes, the top N prediction genes, or both 
do_GSEA = function(syscan, pathway_df, geneType = "prediction", topN = 10, fdr_threshold = 0.01, pathway_levelLimits = c(2, NA), exclude_pathways = NULL, pathway_sizeLimits = c(10, 500)){
  
  # Get genes to test enrichment of
  cat("Getting genes to test\n")
  # Load syscan file if file name was provided
  if (is.character(syscan)) syscan = readRDS(syscan)
  
  # Discard columns we don't need for this function (and make sure we have the ones we do need)
  syscan_gsea = syscan %>% select(sample, entrez, type, score)
  
  # Subset for the genes of interest
  if (geneType != "all") syscan_gsea = syscan_gsea %>% subset(type == geneType)
  syscan_gsea = syscan_gsea %>% arrange(sample, desc(score)) %>% group_by(sample) %>% mutate(patient_rank = 1:n()) %>% ungroup
  if (!is.null(topN)) syscan_gsea = syscan_gsea %>% subset(patient_rank <= topN)
  genes_toTest = unique(syscan_gsea$entrez)
  
  
  # Get pathway list to use
  cat("Preparing pathway list\n")
  # Load file if file name provided; make sure we have required columns
  if (is.character(pathway_df)) pathway_df = readRDS(pathway_df)
  pathway_df = pathway_df %>% select(pathwayID, pathwayName, n_genes, pathwayLevel, genes)
  
  # Level limits
  pathway_levelLimits[which(is.na(pathway_levelLimits))] = c(min(pathway_df$pathwayLevel), max(pathway_df$pathwayLevel))[which(is.na(pathway_levelLimits))]
  pathway_df = pathway_df %>% subset(between(pathwayLevel, pathway_levelLimits[1], pathway_levelLimits[2]))
  
  # Size limits
  pathway_sizeLimits[which(is.na(pathway_sizeLimits))] = c(min(pathway_df$n_genes), max(pathway_df$n_genes))[which(is.na(pathway_sizeLimits))]
  pathway_df = pathway_df %>% subset(between(n_genes, pathway_sizeLimits[1], pathway_sizeLimits[2]))
  
  # Other pathways to exclude
  if (!is.null(exclude_pathways)) pathway_df = pathway_df %>% subset(!pathwayID %in% exclude_pathways)
  
  
  # Do GSEA with Fisher (i.e. hypergeometric) test
  # NB this requires a figure for the "total" number of human genes - here we use 19549
  cat("Testing enrichment of", length(genes_toTest), "genes in", nrow(pathway_df), "pathways\n")
  gsea_res = pathway_df %>% rename(pathwayGenes_n = n_genes, pathwayGenes = genes)
  gsea_res$intersection_size = sapply(gsea_res$pathwayGenes, function(x) length(intersect(x, genes_toTest)))
  gsea_res$fisher_p.value = apply(gsea_res, 1, function(x){
    contingency.table = matrix(c(
      as.numeric(x["intersection_size"]),
      as.numeric(x["pathwayGenes_n"]) - as.numeric(x["intersection_size"]),
      length(genes_toTest) - as.numeric(x["intersection_size"]),
      19549 - as.numeric(x["pathwayGenes_n"]) - length(genes_toTest) + as.numeric(x["intersection_size"])
    ), 2, 2)
    fisher.test(contingency.table, alternative = "greater")$p.value
  })
  
  
  # Correct for FDR and threshold
  gsea_res = gsea_res %>%
    arrange(fisher_p.value) %>%
    mutate(fdr = p.adjust(fisher_p.value, method = "BH")) %>%
    mutate(enriched = fdr < fdr_threshold)
  enrichedPathways = gsea_res %>% subset(enriched) %>% .$pathwayID
  cat(length(enrichedPathways), "pathways enriched with FDR <", fdr_threshold, "\n")
  
  
  # Annotate what enriched pathways each of the sys-candidates is in
  syscan_gsea = syscan_gsea %>% 
    inner_join(pathway_df %>% select(pathwayID, genes) %>% unnest(entrez = genes), by = "entrez") %>%
    subset(pathwayID %in% enrichedPathways) %>%
    rename(enriched_pathwayID = pathwayID)
  
  
  return(list(syscan_gsea = syscan_gsea, gsea_res = gsea_res))
}



# Calculate the Jaccard index between each pair of samples, and run hierarchical clustering based on Euclidean distance
# Also return the distance object itself, as we need this for silhouette
jaccardIndex_hierClust = function(syscan_gsea){
  
  # Damaged enriched pathways in each sample, as a matrix
  sample_pathways = syscan_gsea %>%
    select(sample, enriched_pathwayID) %>%
    unique %>%
    mutate(damaged = 1) %>%
    spread(enriched_pathwayID, damaged, fill = 0) %>%
    column_to_rownames("sample") %>%
    as.matrix
  
  
  # Calculate the Jaccard indices
  cat("Calculating Jaccard indices between", nrow(sample_pathways), "samples\n")
  # Intersection of pathways between each sample
  intersection_pathways = sample_pathways %*% t(sample_pathways)
  # Union of pathways between each sample
  union_pathways = ncol(sample_pathways) - (1 - sample_pathways) %*% t(1 - sample_pathways)
  # Jaccard index
  jaccardIndex = intersection_pathways / union_pathways
  
  
  # Euclidean distance and hierarchical clustering
  cat("Calculating Euclidean distance and running hierarchical clustering\n")
  # Euclidean distance object
  distObject = dist(jaccardIndex, method = "euclidean")
  # Hierarchical clustering
  hclust_res = hclust(distObject)
  
  
  return(list(jaccardIndex = jaccardIndex, distObject = distObject, hclust_res = hclust_res))
}



# Make a silhouette plot to determine the best number of clusters
silhouette_plot = function(hclust_res, distObject, n_clust_range = 2:20){
  distMat = as.matrix(distObject)
  silhouette_res = data.frame()
  
  # Calculate silhouette of each sample for a range of n_clust
  cat("Calculating silhouette distributions for", length(n_clust_range), "different numbers of clusters")
  for (n_clust in n_clust_range){
    silhouette_thisK = data.frame(cluster = cutree(hclust_res, k = n_clust)) %>%
      rownames_to_column("sample") %>%
      group_by(cluster) %>%
      mutate(clusterSize = n()) %>%
      ungroup %>% 
      mutate(n_clusters = n_clust)
    
    
    # Calculate silhouette for each sample - adapted from https://en.wikipedia.org/wiki/Silhouette_(clustering)
    silhouette_thisK$silhouette = apply(silhouette_thisK, 1, function(x){
      thisSample = x["sample"]
      thisCluster = as.numeric(x["cluster"])
      otherSamples_byCluster = lapply(unique(silhouette_thisK$cluster), function(y) silhouette_thisK %>% subset(sample != thisSample & cluster == y) %>% .$sample)
      
      C = as.numeric(x["clusterSize"])
      if (C != 1){
        a = sum(distMat[thisSample, otherSamples_byCluster[[thisCluster]]]) / (C - 1)
        b = c()
        for (k in setdiff(unique(silhouette_thisK$cluster), thisCluster)){
          C_k = silhouette_thisK %>% subset(cluster == k) %>% .$clusterSize %>% head(1)
          b = c(b, sum(distMat[thisSample, otherSamples_byCluster[[k]]]) / C_k)
        }
        b = min(b)
        
        s = (b - a) / max(a, b)
      } else {
        s = 1
      }
      return(s)
    })
    
    silhouette_res = rbind(silhouette_res, silhouette_thisK)
  }
  
  
  # Plot median silhouette over n_clust
  median_silhouette = silhouette_res %>% group_by(n_clusters) %>% summarise(silhouette = median(silhouette))
  p = ggplot(median_silhouette, aes(x = n_clusters, y = silhouette)) +
    geom_line() +
    labs(x = "Clusters (n)", y = "Median silhouette value") +
    theme_gray(base_size = 14)
  print(p)
  
  
  return(list(silhouette_res = silhouette_res, median_silhouette = median_silhouette))
}



# Given a choice for the number of clusters, do the clustering and plot heatmap
cluster_plot = function(jaccardIndex, hclust_res, k){
  
  # Cut dendrogram to get k clusters
  clusters = cutree(hclust_res, k)
  myDendrogram = as.dendrogram(hclust_res)
  
  # Colours for each cluster
  clusterColours = rainbow(length(unique(clusters)), start=0.1, end=0.9)
  clusterColours = clusterColours[as.vector(clusters)] 
  
  # Colour scale for Jaccard index
  mycol = colorpanel(40, "darkblue", "white", "darkred")
  
  # Plot
  heatmap.2(jaccardIndex, 
            Rowv = rev(myDendrogram), 
            Colv = myDendrogram, 
            col = mycol, 
            # breaks = seq(0, 1, 0.2), # Can't get this to work yet
            density.info = "none", 
            trace = "none", 
            RowSideColors = clusterColours,
            ColSideColors = clusterColours,
            keysize = 1,
            labRow = "", # Remove this to print sample names
            labCol = "")
  
  
  # Get list of what samples are in what clusters - samples and clusters in the same order as the heatmap
  sampleOrder = hclust_res$labels[order.dendrogram(myDendrogram)]
  
  clusterAssignments = left_join(
    data.frame(sample = sampleOrder, stringsAsFactors = F),
    data.frame(cluster = clusters) %>% rownames_to_column("sample"),
    by = "sample"
  )
  
  cluster_renumbering = clusterAssignments %>%
    mutate(cluster_old = factor(cluster, levels = unique(cluster))) %>%
    group_by(cluster_old) %>%
    summarise(cluster_size = n()) %>%
    mutate(cluster_new = 1:nrow(.), cluster_old = as.numeric(as.character(cluster_old)))
  
  clusterAssignments = left_join(clusterAssignments, cluster_renumbering, by = c("cluster" = "cluster_old")) %>%
    select(sample, cluster = cluster_new, cluster_size)
  
  return(list(clusters = clusterAssignments))
}



##---- Example implementation of functions ----


gsea = do_GSEA(syscan = "/Volumes/lab-ciccarellif/working/Joel/OAC_sysSVM_2.0/methods_paper/kernel_tests/noMeanScale/10000.iterations_in_20.batches/score_order3/9900iterations/syscan.rds",
               pathway_df = "/Volumes/lab-ciccarellif/working/Joel/OAC_sysSVM_2.0/methods_paper/GSEA/reactome_v70_human.rds",
               geneType = "prediction",
               topN = 10,
               pathway_levelLimits = c(2, NA),
               pathway_sizeLimits = c(10, 500),
               fdr_threshold = 0.01)


jaccInd = jaccardIndex_hierClust(gsea$syscan_gsea)


silhouette = silhouette_plot(jaccInd$hclust_res, jaccInd$distObject, n_clust_range = 2:10)


clustering = cluster_plot(jaccInd$jaccardIndex, jaccInd$hclust_res, k = 5)
