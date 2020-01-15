# sysSVM2-NN

## Description
sysSMV2-NN is a computational tool for patient-specific cancer driver gene prioritisation. It is based on the principle that driver genes are characterised by particular molecular properties (*e.g.* mutations, copy number variants) and systems-level properties (*e.g.* evolutionary origin, breadth of expression). It works by identifying genes with similar properties to canonical drivers.
\
\
The first part of the method (sysSVM2) requires a cohort of cancer samples. There are four broad steps to the algorithm: 
1. **Feature mapping**: identify the molecular and systems-level properties of damaged genes in the cohort; mark canonical drivers as a training set
1. **Model selection**: tune SVM parameters to optimise performance, based on the sensitivity on the training set 
1. **Training**: train the model with the selected parameters
1. **Prediction**: predict on new samples/genes, assigning a score to each
<a/>
The second part of the method uses a Neural Network (NN) to incorporate additional training samples to expand the initial sysSVM2 model. This can be useful in settings where new data arrive sporadically, as it saves the user from having to re-train sysSVM2 *de novo*.  

## Download/installation
TO DO


## Running sysSVM2 on an initial cohort
In R, first source the file ```train_predict_functions.R```. This contains the functions that execute sysSVM2. 
\
\
The cohort data must be formatted for sysSVM. An example sysSVM2 input file, for a cohort of 100 simulated pan-cancer samples, is provided: ```molecular_data = read_tsv("sysSVM_NN/example_data/molecular_features_100samples.tsv")```. The required ID columns are
* sample: Sample identifiers
* entrez: Gene Entrez IDs
and the recommended molecular feature columns are
* no_ALL_muts: Total number of mutations in a gene
* no_NTDam_muts: Number of truncating mutations
* no_NTDam_muts: Number of non-truncating mutations
* no_GOF_muts: Number of gain-of-function/hotspot mutations
* Copy_number: Total copy number
* CNVGain: Binary 0/1 indicating gene amplification (recommended copy number > 2 * ploidy)
* CNVLoss: Binary 0/1 indicating gene loss (recommended copy number <= 1)
<a/>
To complete the feature mapping of the cohort, the systems-level properties of the genes are also required. A compendium of 25 of these properties, each of which distinguish cancer genes from the rest of human genes, is provided: ```systemsLevel_data = read_tsv("sysSVM_NN/example_data/systemsLevel_features_allGenes.tsv")```. 
\
\
Finally, join the two tables to create the sysSVM2 input file: ```sysSVM2_input = inner_join(molecular_data, systemsLevel_data, by = "entrez")```.

After preparing the input file, further data formatting is done using ```prepare_trainingPrediction```. One of the main purposes of this is to separate the training and prediction sets (*i.e.* canonical drivers, and the rest of genes). 
\
\
The next step is model selection, in which SVM parameters are tuned to the training set. First, a grid of parameter combinations is assessed using three-fold cross-validation. This is performed using ```run_crossValidation_par```, which can run in parallel envrionments. We recommend at least 1000 iterations in practice, and this is fairly computationally intensive. After the cross-validation iterations have been run, they are assessed using ```selectParams_from_CVstats```, which identifies the best parameter combinations. 
\
\
We can now train the SVMs with selected parameters, with ```train_sysSVM```. The trained model lcan be used to make predictions with ```predict_sysSVM```.


## Expanding a cohort with sysSVM2-NN
Run ```some code.py```


## Reference
TO DO


## Acknowledgements
Thanos Mourikis\
Damjan Temelkovski\
Christopher Yau
