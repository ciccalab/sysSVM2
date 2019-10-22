#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 13:11:08 2019

@author: nulsenj

Based on code from https://keras.io/examples/mnist_denoising_autoencoder/
"""

#%% Parse arguments

# Set up argument parser
from argparse import ArgumentParser
parser = ArgumentParser(description = 'Initial training of AAE to reproduce sysSVM scores')



# Required arguments (positional, no defaults)
parser.add_argument("input_file", help = "Location of training data (systems-level/molecular features and scores)")
parser.add_argument("output_dir", help = "Directory to write output to")



# Optional arguments
parser.add_argument("--val_proportion", type = float, default = 0.2, help = "Proportion of observations to use for validation")
parser.add_argument("--epochs", type = int, default = 50, help = "Number of epochs to use in training")
parser.add_argument("--score_lossWeight", type = float, default = 5, help = "Relative weight of score loss compared to reconstruction loss")
parser.add_argument("--batch_size", type = int, default = 32, help = "Batch size for training; set to zero not to use batches")
parser.add_argument("--learning_rate", type = float, default = 0.0005, help = "Learning rate for Adam optimiser")



# Extract arguments
args = parser.parse_args()



#%% Packages
print('Loading packages')

# General
import numpy as np
import pandas as pd
import os


# Tensorflow / Keras
import tensorflow as tf
from tensorflow.keras.layers import Dense, Input
from tensorflow.keras.models import Model
from tensorflow.compat.v1.keras.initializers import RandomNormal
from tensorflow.keras.losses import mse, mae
from tensorflow.keras.optimizers import Adam
from tensorflow.compat.v1 import Session, reshape
import tensorflow.keras.backend as K


# Custom functions, defined in AAE_functions.py (must be in the same directory as this script)
from AAE_functions import scaleData, splitData, scoreLoss, minMax


np.random.seed(1337)



#%% Prepare data


# Load input data (systems-level and molecular features, as well as sysSVM scores)
sysMol_scores_fn = args.input_file
print('Reading', sysMol_scores_fn)
sysMol_scores = pd.read_csv(sysMol_scores_fn, sep = '\t', index_col = 'id')



# Scale data
print('Scaling data')
sysMol_scores, scaling_factors = scaleData(sysMol_scores, verbose = True)



# Split data into training-validation, convert to numpy
val_proportion = args.val_proportion
print('Splitting training / validation data', (1-val_proportion)*100, '/', val_proportion*100)
train_IDs, train_data_in, train_data_out, train_scores, val_IDs, val_data_in, val_data_out, val_scores = splitData(sysMol_scores, val_proportion = val_proportion)




#%% Network settings

# Dimensions of AE layers
in_dim = train_data_in.shape[1]
out_dim = train_data_out.shape[1]
intermediate_dims = [28, 20, 14] # Will be reversed for the decoder
latent_dim = 10



# Dimensions of score-producing network layers (takes latent RDR as input)
score_dims = [5, 3, 2]



# Activation function to use for internal layers
my_activationFunction = 'tanh'



# Initialisation distribution for weights
initSD = 0.5
my_init = RandomNormal(stddev = initSD)



# Loss function for scores, and loss function for sys/mol feature reconstruction
# Combine in a weighted fashion
my_lossFunctions = [scoreLoss, mse]
score_lossWeight = args.score_lossWeight # Loss function weight for score component, relative to reconstruction loss



# Optimiser
my_learning_rate = args.learning_rate
my_optimiser = Adam(lr = my_learning_rate)



# Training epochs
my_epochs = args.epochs



# Batch size
my_batch_size = args.batch_size
if my_batch_size == 0: # If set to zero, don't use batches
    my_batch_size = train_data_in.shape[0]



#%% Build network
print('Building network')


# Build encoder model
inputs = Input(shape = (in_dim,), name = 'encoder_input')
x0 = Dense(intermediate_dims[0], activation = my_activationFunction, kernel_initializer = my_init)(inputs)
x1 = Dense(intermediate_dims[1], activation = my_activationFunction, kernel_initializer = my_init)(x0)
x2 = Dense(intermediate_dims[2], activation = my_activationFunction, kernel_initializer = my_init)(x1)

# Generate the latent representation
latent = Dense(latent_dim, activation = my_activationFunction, name = 'latent_vector')(x2)

# Instantiate Encoder Model
encoder = Model(inputs, latent, name = 'encoder')
#encoder.summary()



# Build decoder model
latent_inputs = Input(shape = (latent_dim,), name = 'decoder_input')
y0 = Dense(intermediate_dims[2], activation = my_activationFunction, kernel_initializer = my_init)(latent_inputs)
y1 = Dense(intermediate_dims[1], activation = my_activationFunction, kernel_initializer = my_init)(y0)
y2 = Dense(intermediate_dims[0], activation = my_activationFunction, kernel_initializer = my_init)(y1)
outputs_reconstruction = Dense(out_dim, kernel_initializer = my_init, activation = 'tanh')(y2) # Make sure this activation function is compatible with data

# Instantiate the decoder model
decoder = Model(latent_inputs, outputs_reconstruction, name = 'decoder')
#decoder.summary()



# Build score prediction model
latent_inputs_scores = Input(shape = (latent_dim,), name = 'score_input')
z0 = Dense(score_dims[0], activation = my_activationFunction, kernel_initializer = my_init)(latent_inputs_scores)
z1 = Dense(score_dims[1], activation = my_activationFunction, kernel_initializer = my_init)(z0)
z2 = Dense(score_dims[2], activation = my_activationFunction, kernel_initializer = my_init)(z1)
outputs_scores = Dense(1, kernel_initializer = my_init, activation = 'tanh')(z2) # Make sure this activation function is compatible with data

# Instantiate the decoder model
score_nn = Model(latent_inputs_scores, outputs_scores, name = 'score_nn')
#score_nn.summary()



# Augmented autoencoder = Encoder + Decoder + Score predictor
# Instantiate the augmented autoencoder
aae = Model(inputs, [score_nn(encoder(inputs)), decoder(encoder(inputs))], name = 'aae')
aae.summary()



#%% Compile and train
print('Training augmented autoencoder')

# Compile
aae.compile(loss = my_lossFunctions,
            loss_weights=[score_lossWeight, 1],
            metrics = [mae, mse],
            optimizer = my_optimiser)



# Train the autoencoder - this will also fit parameters of the submodels
training_history = aae.fit(train_data_in,
                           [train_scores, train_data_out],
                           batch_size = my_batch_size, # Smaller batch size takes longer but is more accurate
                           validation_data = (val_data_in, [val_scores, val_data_out]),
#                           verbose = 2,
                           epochs = my_epochs)



# Record the losses over epochs
training_history = pd.DataFrame.from_dict(training_history.history)
training_history['epoch'] = list(range(1, my_epochs+1))



#%% Extract score estimates and RDR for test data
print('Getting score estimates and RDR for test data')


# Run test data through the model
rdr = encoder.predict(val_data_in)
score_estimates = np.squeeze(score_nn.predict(rdr))



# Convert RDR to dataframe
rdr_colNames = []
for i in list(range(latent_dim)):
    col = 'rdr_' + str(i+1)
    rdr_colNames = rdr_colNames + [col]

rdr_df = pd.DataFrame(data = rdr, columns = rdr_colNames)



# Compile results into a single dataframe
res = pd.DataFrame({'id': val_IDs,
                    'true_scaledScore': np.squeeze(val_scores),
                    'reproducedScore': score_estimates})
res = res.reset_index(drop = True)
res = pd.concat([res, rdr_df], axis = 1)



#%% Save results and trained model
print('Saving results')


# Fix output directory so it has a forward slash at the end
output_dir = args.output_dir
if output_dir[len(output_dir)-1] != '/':
    output_dir = output_dir + '/'



# Scaling factors
scaling_factors_df = pd.DataFrame.from_dict(scaling_factors, orient='index', columns=['lower', 'upper'])
scaling_factors_df.to_csv(output_dir + 'scalingFactors.tsv', sep = '\t', index = True)



# Score predictions, and RDR of un-noisy test data
res.to_csv(output_dir + 'aae_preds_rdr.tsv', sep = '\t', index = False)



# Training history
training_history.to_csv(output_dir + 'aae_trainingHistory.tsv', sep = '\t', index = False)



# Trained model
aae.save(output_dir + 'aae_model.h5')



print('Finished')



                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
