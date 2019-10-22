#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 11:52:37 2019

@author: nulsenj
"""

#%% Parse arguments

# Set up argument parser
from argparse import ArgumentParser
parser = ArgumentParser(description = 'Secondary training of AAE to incorporate new data')



# Required arguments (positional, no defaults)
parser.add_argument("original_data_file", help = "Location of original training data (systems-level/molecular features and scores)")
parser.add_argument("new_data_file", help = "Location of new training data (no scores, systems-level/molecular features only)")
parser.add_argument("scaling_factors_file", help = "Location of scaling factors (output of initial training)")
parser.add_argument("trained_model_file", help = "Location of pre-trained AAE model")
parser.add_argument("output_dir", help = "Directory to write output to")



# Optional arguments
parser.add_argument("--val_proportion", type = float, default = 0.2, help = "Proportion of observations to use for validation")
parser.add_argument("--epochs", type = int, default = 100, help = "Number of epochs to use in training")
#parser.add_argument("--batch_size", type = int, default = 32, help = "Batch size for training; set to zero not to use batches")
#parser.add_argument("--learning_rate", type = float, default = 0.0005, help = "Learning rate for Adam optimiser")



# Extract arguments
args = parser.parse_args()



#%% Packages
print('Loading packages')

# General
import numpy as np
import pandas as pd


# Tensorflow / Keras
from tensorflow.keras.models import load_model


# Custom functions, defined in AAE_functions.py (must be in the same directory as this script)
from AAE_functions import scaleData, splitData, scoreLoss, minMax, newSampleMetric


np.random.seed(1338)



#%% Prepare data (two datasets - original with scores, and new without scores)


# File locations
originalData_fn = args.original_data_file
newData_fn = args.new_data_file
scaling_factors_fn = args.scaling_factors_file



# Load data
print('Reading', originalData_fn, 'and', newData_fn)
originalData = pd.read_csv(originalData_fn, sep = '\t', index_col = 'id')
newData = pd.read_csv(newData_fn, sep = '\t', index_col = 'id')



# Scale data
print('Scaling data using scaling factors provided at', scaling_factors_fn)
scaling_factors = pd.read_csv(scaling_factors_fn, sep = '\t', index_col = 0)
originalData = scaleData(originalData, scaling_factors = scaling_factors)
newData = scaleData(newData, scaling_factors = scaling_factors, has_scores = False)



# Concatenate into a single dataset
# Make sure column orders match what the trained NN expects
newData = newData[originalData.columns.values]
fullData = pd.concat([originalData, newData], sort = False)



# Split data into training-validation, convert to numpy
val_proportion = args.val_proportion
print('Splitting training / validation data', (1-val_proportion)*100, '/', val_proportion*100)
train_IDs, train_data_in, train_data_out, train_scores, val_IDs, val_data_in, val_data_out, val_scores  = splitData(fullData, val_proportion = val_proportion)
# Boolean indicating which of the validation samples are from the new dataset
val_newData_bool = np.array([(x in newData.index.values) for x in val_IDs])



#%% Load a saved (trained) model

# scoreLoss is defined in functions.py
# Load from HDF5 file (name given as argument)
preTrained_model_fn = args.trained_model_file
aae = load_model(preTrained_model_fn, custom_objects={'scoreLoss': scoreLoss})



#%% Train
print('Training with new samples')

# Number of training epochs
# Without batches, need substantially more epochs
my_epochs = args.epochs



# Compile with new metric to track reconstruction error of new/old samples separately
aae.compile(loss = aae.loss,
            loss_weights = aae.loss_weights,
            optimizer = aae.optimizer,
            metrics = [newSampleMetric(mask = val_newData_bool)])



# Train
training_history = aae.fit(train_data_in,
                           [train_scores, train_data_out],
                           validation_data = (val_data_in, [val_scores, val_data_out]),
                           batch_size = train_data_in.shape[0], # Batches aren't compatible with newSampleMetric
#                           verbose = 2,
                           epochs = my_epochs)


# Record the losses over epochs
training_history = pd.DataFrame.from_dict(training_history.history)
training_history['epoch'] = list(range(1, my_epochs+1))



#%% Extract score estimates for test data
print('Getting score estimates')


# Run data through the model
scorePreds_train = np.squeeze(aae.predict(train_data_in)[0])
scorePreds_val = np.squeeze(aae.predict(val_data_in)[0])



# Compile results into a single dataframe
res = pd.DataFrame({'id': list(train_IDs) + list(val_IDs),
                    'true_scaledScore': list(np.squeeze(train_scores)) + list(np.squeeze(val_scores)),
                    'predicted_scaledScore': list(scorePreds_train) + list(scorePreds_val)})



# Unscale the scores
res['true_score'] = minMax(res['true_scaledScore'],
   newMin = scaling_factors['lower']['score'], newMax = scaling_factors['upper']['score'],
   l = -0.8, u = 0.8)

res['predicted_score'] = minMax(res['predicted_scaledScore'],
   newMin = scaling_factors['lower']['score'], newMax = scaling_factors['upper']['score'],
   l = -0.8, u = 0.8)




#%% Save results and trained model
print('Saving results')


# Fix output directory so it has a forward slash at the end
output_dir = args.output_dir
if output_dir[len(output_dir)-1] != '/':
    output_dir = output_dir + '/'



# Score predictions
res.to_csv(output_dir + 'aaeExtended_preds.tsv', sep = '\t', index = False)



# Training history
training_history.to_csv(output_dir + 'aaeExtended_trainingHistory.tsv', sep = '\t', index = False)



# Trained model
aae.save(output_dir + 'aaeExtended_model.h5')



print('Finished')



                                                                                                                                                                                                                                                                                                                                                                                                                                              
