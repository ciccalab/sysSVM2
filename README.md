# sysSVM2

## Description
sysSMV2 is a computational tool for patient-specific cancer driver gene prioritisation. It is based on the principle that driver genes are characterised by particular molecular properties (*e.g.* mutations, copy number variants) and systems-level properties (*e.g.* evolutionary origin, breadth of expression). It works by identifying genes with similar properties to canonical drivers, using a one-class Support Vector Machine framework<sup>1</sup>.
\
\
Models that have been trained on data from TCGA are available to download from the appropriate [directory](trained_models). For smaller TCGA cohorts (N <200 samples), we recommend using the model trained on pan-cancer data. 
\
\
Users may also train their own sysSVM2 models. 

## Download/installation
To download sysSVM2, clone this repository. sysSVM2 is implemented in R. In this guide, we assume that the user's working directory corresponds to a clone of this repository.

## Running a pre-trained sysSVM2 model on your data
If you have your own data and want to predict driver genes without training a new model, there are two steps:
1. **Data annotation**: sysSVM2 requires mutation and copy number data in a particular format
1. **Prediction**: use one of the [pre-trained models](trained_models) to identify drivers in individual samples.
[//]: # (end list)

### 1. Data annotation
sysSVM2 uses somatic mutation and CNV data summarised at the gene-level (molecular properties). These gene-level summaries can be produced from standard variant call file formats using the annotation functions provided in this repository. In R, first source these functions:
```
source("R/annotation_functions.R")
```
Annotation of small somatic mutations (SSMs, *i.e.* SNVs and indels) requires [ANNOVAR](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/) to be installed. A VCF file (aligned to  hg19 in this example) can be converted to annotated SSMs as follows:
```
ssms_annotated = annotate_ssms(
  vcf_fn, annovar_dir, 
  genome_version = "hg19", 
  gene_aliases_entrez = "annotation_reference_files/gene_aliases_entrez.tsv", 
  hotspots = "annotation_reference_files/tcga_pancancer_hotspots_oncodriveclust.tsv"
)
```
Note that, by default, mutation hotspots identified in TCGA pan-cancer data are mapped to SSMs in order to identify gain-of-function mutations.
\
\
Annotation of CNV segments requires [bedtools](https://bedtools.readthedocs.io/en/latest/) to be installed. If available, sample ploidy values can be used to filter copy number gains. To intersect CNV segments with gene coordinates and identify gene gains and losses, run the following:
```
cnvs_annotated = annotate_cnvs(
  cnv_segments, bedtools_bin_dir,
  ploidy, 
  gene_coords = "annotation_reference_files/gene_coords_hg19.tsv"
)
```
Finally, the annotated SSMs and CNVs can be combined into a unified format for sysSVM2:
```
molecular_data = make_totalTable(
  ssms_annotated, cnvs_annotated, 
  canonical_drivers = "example_data/canonical_drivers.rds"
)
```

### 2. Prediction
A pre-trained model can now be used to make predictions on your annotated data. First, load one of the provided pre-trained models. For example, the model trained on pan-cancer data:
```
trained_sysSVM = readRDS("trained_models/PANCAN_trained_sysSVM.rds")
```
To run an example without annotating variant call data, you can use the provided (annotated) small cohort:
```
molecular_data = read_tsv("example_data/molecular_features_100samples.tsv")
```
Then simply run the trained model on the annotated data:
```
predictions = predict_sysSVM2(
  trained_sysSVM, 
  molecular_data = molecular_data, 
  systemsLevel_data = "example_data/systemsLevel_features_allGenes.tsv"
)
```
The output is a ranked list of damaged genes in each patient, with high scores/ranks corresponding to putative driver genes. 


## Training a new sysSVM2 model
Training sysSVM2 requires a cohort of cancer samples. There are four broad steps to the algorithm: 
1. **Feature mapping**: identify the molecular and systems-level properties of damaged genes in the cohort; mark canonical drivers as a training set
1. **Model selection**: tune SVM parameters to optimise performance, based on the sensitivity on the training set. This step can also be skipped if you want to use the default kernel parameters.
1. **Training**: train the model with the selected parameters
1. **Prediction**: predict on new samples/genes, assigning a score to each. If necessary, extract a list of the canonical drivers in each sample, 'topped up' with the highest-scoring predictions from the algorithm.

[//]: # (end list)

In this guide, we assume that the user's working directory corresponds to a clone of this repository. We also assume that results are output to a directory called ```~/test_sysSVM2```. 


### 1. Feature mapping
After formatting the data for sysSVM2 (see **Data annotation** above), the process of training and prediction can begin. First, source the relevant R functions:
```
source("R/train_predict_functions.R")
```
An example sysSVM2 input file, for a cohort of 100 simulated pan-cancer samples, is provided: 
```
molecular_data = read_tsv("example_data/molecular_features_100samples.tsv")
```
The required ID columns are
* ```sample```: Sample identifiers
* ```entrez```: Gene Entrez IDs,

[//]: # (end list)

and the recommended molecular feature columns (produced by ```R/annotation_functions.R```) are
* ```no_ALL_muts```: Total number of mutations in a gene
* ```no_NTDam_muts```: Number of truncating mutations
* ```no_NTDam_muts```: Number of non-truncating mutations
* ```no_GOF_muts```: Number of gain-of-function/hotspot mutations
* ```Copy_number```: Total copy number
* ```CNVGain```: Binary 0/1 indicating gene amplification (recommended copy number > 2 * ploidy)
* ```CNVLoss```: Binary 0/1 indicating gene loss (recommended copy number <= 1).

[//]: # (end list)

To complete the feature mapping of the cohort, the systems-level properties of the genes are also required. A compendium of 25 of these properties, each of which distinguish cancer genes from the rest of human genes, is provided: 
```
systemsLevel_data = read_tsv("example_data/systemsLevel_features_allGenes.tsv")
``` 
Finally, join the two tables to create the sysSVM2 input file: 
```
sysSVM2_input = inner_join(molecular_data, systemsLevel_data, by = "entrez")
```
### 2. Model selection
After preparing the input file, the next step is to separate the training and prediction sets (*i.e.* canonical drivers, and the rest of genes), and perform data normalisation. A list of canonical driver genes (tumour suppressor genes/oncogenes), along with their respective driver alteration types (loss/gain of function) is provided. By default, sysSVM2 scales numeric features to have unit standard deviation.
```
canonical_drivers = readRDS("example_data/canonical_drivers.rds")
sysSVM_data = prepare_trainingPrediction(sysSVM2_input, canonical_drivers, output_dir = "~/test_sysSVM2")
```
SVM model parameters then need to be tuned. sysSVM2 does this automatically, using three-fold Cross-Validation (CV) over a pre-determined grid of parameter combinations. However, this is a computationally intensive process. If sufficient computational resources are not available, the parameters chosen on simulated pan-cancer data, provided in Supplementary Figure 2 of the sysSVM2 publication, can be used. Skip to Training below if you want to do this. Alternatively, if users want to train a model on a similar dataset to one for which they have already carried out CV iterations, kernel parameters might be 're-used' to save on computation time.  
\
To run the CV iterations as efficiently as possible, the code is designed to run in a parallel environment. For example, to perform 10 CV iterations using 4 cores, run
```
cv_stats = run_crossValidation_par(iters = 10,
                                   cores = 4, 
                                   inPath = "~/test_sysSVM2",
                                   outPath = "~/test_sysSVM2",
                                   parallelLib = "parallel")
```
The ```parallelLib``` argument should be set to either ```"parallel"``` or ```"snow"```, depending on the user's environment. In practice, we recommend using at least 1,000 CV iterations to ensure convergence of parameter selections. 
\
\
After the CV iterations have been run, their results are assessed to identify the best parameter combinations:
```
model_selection = selectParams_from_CVstats(cv_stats, output_dir = "~/test_sysSVM2")
```
### 3. Training
After model selection, the entire training set is used to train the final sysSVM2 model:
```
trained_sysSVM2 = train_sysSVM2(model_parameters = model_selection$best_model_final, 
                                training_set = sysSVM_data$training_set, 
                                scaling_factors = sysSVM_data$scaling_factors,
                                output_dir = "~/test_sysSVM2")
```
If you want to use the default SVM parameters learned from simulated pan-cancer data, simply do not provide a ```model_parameters``` argument.
### 4. Prediction
The trained model can now be used to make predictions, either on the same cohort or on new samples. To use the prediction set from the same cohort used for training:
```
predictions = predict_sysSVM2(trained_sysSVM, 
                              prediction_set = sysSVM_data$prediction_set, 
                              prediction_set_ns = sysSVM_data$prediction_set_ns, 
                              output_dir = "~/test_sysSVM2")
```
The output is a ranked list of damaged genes in each patient, with high scores/ranks corresponding to putative driver genes. Using the above code, it will also be saved in a file called ```scores.rds```.
\
It may be desirable to extract a list of predicted drivers in each sample, rather than a scored/rankded list of all damaged genes. This can be achieved using a 'top-up' procedure, in which (1) all canonical drivers are considered to be predicted drivers, and (2) the highest-scored predictions from sysSVM2 are used to top up the predictions in each sample, up to a specified minimum number of genes per sample. In this example, we can obtain lists of at least five predicted drivers per sample as follows:
```
drivers_toppedUp = topUp_drivers(all_genes = sysSVM_data,
                                 gene_scores = predictions,
                                 canonical_drivers = canonical_drivers,
                                 n_drivers_per_sample = 5,
                                 output_dir = "~/test_sysSVM2")
```
The output will be saved in a file called ```drivers_toppedUp.rds```.


## Reference
1. B. Shoelkopf, J. Platt, J. Shawe-Taylor, A. Smola and R. Williamson. Estimating the support of a high-dimensional distribution. *Neural Comput.* **13** (2001). 
2. K. Wang, M. Li, H. Hakonarson. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data. *Nucleic Acids Research* **38** (2010). [Webpage](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/)

[//]: # (end list)


## Acknowledgements
Hrvoje Misetic\
Christopher Yau\
Thanos Mourikis\
Damjan Temelkovski

