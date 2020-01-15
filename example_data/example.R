##---- Setup ----

# Top-level input directory
setwd("/Volumes/lab-ciccarellif/working/Joel/OAC_sysSVM_2.0/methods_paper")
# Directory to work in
test_dir = "/Volumes/lab-ciccarellif/working/Joel/OAC_sysSVM_2.0/methods_paper/sysSVM_implementation/test_example"


# Functions
source("sysSVM_implementation/run_sysSVM/train_predict_functions.R")


# Load data
molecular_data = readRDS("simulate_pancan/totalTable_1000simulatedSamples_onlyFiltered_seed1.rds") %>% select(-ploidy)
systemsLevel_data = readRDS("sysProperties/ncg6_sysProp_mmImputed_expanded.rds") %>% select(-gene_type)
truePositive_drivers = readRDS("sysSVM_implementation/cancerAltTypes_TSG_OG.rds")


# Make cohort smaller
my_samples = paste0("sim", 1:100)
molecular_data = molecular_data %>% subset(sample %in% my_samples)


# Kernel parameter grids to search over
kernels_paramGrids = list(
  linear     = list(nu = c(0.05, 0.2)),
  polynomial = list(nu = c(0.05, 0.2), gamma = 1, coef0 = 0, degree = c(2, 4)),
  radial     = list(nu = c(0.05, 0.2), gamma = 2^c(-7, 0, 4)),
  sigmoid    = list(nu = c(0.05, 0.2), gamma = 2^c(-7, 0, 4), coef0 = 0)
)


##---- Run example ----

# 1. Join molecular and systems-level properties
cohort_data = inner_join(molecular_data, systemsLevel_data, by = "entrez")


# 2. Prepare data
# # marks default
sysSVM_data = prepare_trainingPrediction(
  # Data
  cohort_data = cohort_data,
  truePositive_drivers = truePositive_drivers,
  output_dir = test_dir, # 
  # Exclude observations/features
  exclude_features = c("no_NSI_muts", "young", "rnaseq_ubiquitous", "protein_low", "ppin_central", "rnaseq_selective", "length_RefSeq", "domains_InterPro", "multiDomain"),
  damaging_subsetCond = "no_TRUNC_muts + no_NTDam_muts + no_GOF_muts > 0 | CNVGain == 1 | Copy_number == 0", #
  # Scaling
  factorise_twoLevel = T, # 
  nonFactor_features = c("sample", "entrez", "no_ALL_muts", "no_NSI_muts", "no_NTDam_muts", "no_TRUNC_muts", "no_GOF_muts"), # 
  mean_center = "none", #
  sd_scale = "training", #
  # Housekeeping
  sample_gene_sep = "__" #
)


# 3. CV iterations to measure sensitivity
# Takes MUCH longer with a dense parameter grid, 10k iterations and a cohort 10x as large
cv_stats = run_crossValidation_par(inPath = test_dir,
                                   outPath = test_dir,
                                   folds = 3, 
                                   iters = 2,
                                   kernels_paramGrids = kernels_paramGrids,
                                   cores = 1, 
                                   verbose = F,
                                   parallelLib = "parallel")


# 4. Select best model
best_model = selectParams_from_CVstats(cv_stats, output_dir = test_dir)


# 5. Train best model on full training set
trained_sysSVM = train_sysSVM(model_parameters = best_model$best_model_final, training_set = sysSVM_data$training_set, output_dir = test_dir)


# 6. Predict
predictions = predict_sysSVM(trained_sysSVM, prediction_set = sysSVM_data$prediction_set, prediction_set_ns = sysSVM_data$prediction_set_ns, output_dir = test_dir)


# Explore
predictions %>%
  arrange(sample, desc(score)) %>%
  select(sample, entrez, linear_decision_value:sample_rank) %>% 
  View
