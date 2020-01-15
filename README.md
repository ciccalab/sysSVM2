# sysSVM2-NN

## Description
sysSMV2-NN is a computational tool for patient-specific cancer driver gene prioritisation. It is based on the principle that driver genes are characterised by particular molecular properties (*e.g.* mutations, copy number variants) and systems-level properties (*e.g.* evolutionary origin, breadth of expression). It works by identifying genes with similar properties to canonical drivers.
\
\
The first part of the method (sysSVM2) requires a cohort of cancer samples. There are four broad steps to the algorithm: 
1. Feature mapping: identify the molecular and systems-level properties of damaged genes in the cohort; mark canonical drivers as a training set
1. Model selection: tune SVM parameters to optimise performance, based on the sensitivity on the training set 
1. Training: train the model with the selected parameters
1. Prediction: predict on new samples/genes, assigning a score to each
<a/>
The second part of the method uses a Neural Network (NN) to incorporate additional training samples to expand the initial sysSVM2 model. This can be useful in settings where new data arrive sporadically, as it saves the user from having to re-train sysSVM2 *de novo*.  

## Download/installation
TO DO


## Running sysSVM2 on an initial cohort
Run ```some code.R```


## Expanding a cohort with sysSVM2-NN
Run ```some code.py```


## Reference
TO DO


## Acknowledgements
Thanos Mourikis\
Damjan Temelkovski\
Christopher Yau
