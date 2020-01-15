#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 16:11:42 2019

@author: nulsenj
"""

#%% Packages

import numpy as np
import tensorflow as tf
import tensorflow.keras.backend as K
import tensorflow.keras.metrics as metrics


#%% Scale data

# Function to do min-max scaling
# Must be compatible with activation function of autoencoder output
# E.g. between +/- 0.8 for tanh, or between 0.1 to 0.9 for sigmoid
# Keep track of the scaling factors used for new data, or work from provided scaling factors
def minMax(x, newMin, newMax, l = None, u = None):
    if l is None and u is None:
        l = min(x)
        u = max(x)
        y = (x - l) / (u - l) * (newMax - newMin) + newMin
        return y, [l, u]
    else:
        y = (x - l) / (u - l) * (newMax - newMin) + newMin
        return y



# Carry out scaling on a dataframe
# Log-transform required columns
# Min-max scale, either from provided scaling_factors or recording new ones
# Add dummy scores if missing
def scaleData(df,
              logCols = ['length_RefSeq', 'ppin_degree', 'ppin_betweenness', 'complexes', 'mirna'],
              verbose = False,
              minMax_limits = [-0.8, 0.8],
              scaling_factors = None,
              has_scores = True):

    # Log transform columns where appropriate
    if logCols is not None:
        if verbose: print('Log-transforming:', *logCols)

        for col in logCols:
            df[col] = np.log10(df[col] + 1)


    # First dataset, no scaling factors provided
    if scaling_factors is None:
        new_scaling_factors = {}
        if verbose: print('Min-max scaling data between', minMax_limits[0], 'and', minMax_limits[1])

        for col in df.columns.values:
            df[col], new_scaling_factors[col] = minMax(x = df[col], newMin = minMax_limits[0], newMax = minMax_limits[1])
    # Secondary/repeated dataset, use provided scaling factors
    else:
        for col in df.columns.values:
            df[col] = minMax(x = df[col], l = scaling_factors['lower'][col], u = scaling_factors['upper'][col], newMin = minMax_limits[0], newMax = minMax_limits[1])


    # Add dummy scores if necessary
    if not has_scores:
        df['score'] = -10000


    # Return output
    if scaling_factors is None:
        return df, new_scaling_factors
    else:
        return df



#%% Split data

# Training-validation
# Input-output
# Convert to numpy
def splitData(df, val_proportion = 0.2):

    # First split the sample_entrez IDs
    all_IDs = df.index.values
    n_total_obs = len(all_IDs)
    n_test_obs = np.int_(val_proportion * n_total_obs)
    val_IDs = list(np.random.choice(all_IDs, size = n_test_obs, replace = False))
    train_IDs = list(set(all_IDs) - set(val_IDs))


    # Then use these to split the actual data and separate the score column
    # Additionally, the genes_thisSample should only be used as input and should not be something to decode
    # Training
    train_data_in = df.loc[train_IDs].drop(columns=['score'])
    train_data_out = df.loc[train_IDs].drop(columns=['score', 'genes_thisSample'])
    train_scores = df.loc[train_IDs, ['score']]
    # Validation
    val_data_in = df.loc[val_IDs].drop(columns=['score'])
    val_data_out = df.loc[val_IDs].drop(columns=['score', 'genes_thisSample'])
    val_scores = df.loc[val_IDs, ['score']]


    # Finally, convert to numpy arrays
    # Training
    train_data_in = train_data_in.to_numpy()
    train_data_out = train_data_out.to_numpy()
    train_scores = train_scores.to_numpy()
    # Validation
    val_data_in = val_data_in.to_numpy()
    val_data_out = val_data_out.to_numpy()
    val_scores = val_scores.to_numpy()


    # Return output
    return train_IDs, train_data_in, train_data_out, train_scores, val_IDs, val_data_in, val_data_out, val_scores



#%% Custom loss function for scores

# We need to train a NN with data where some sysSVM scores are present and others are missing
# The general loss function is MAE + MSE
# Additionally, this function ignores missing score values, encoded by -10000
def scoreLoss(yTrue, yPred):

    # Remove values where yTrue is -10000
    yTrue_subset = tf.boolean_mask(yTrue, tf.not_equal(yTrue, -10000))
    yPred_subset = tf.boolean_mask(yPred, tf.not_equal(yTrue, -10000))

    # Add mae and mse
    return K.mean(K.square(yTrue_subset - yPred_subset) + K.abs(yTrue_subset - yPred_subset), axis=-1)



#%% Custom metric for new samples

# During secondary training we want to track the reconstruction error on the new samples separately
# General loss function is MSE
# This function subsets according to a Boolean mask
def newSampleMetric(mask):

    # Define the actual loss function
    def loss(yTrue, yPred):
        # Subset according to mask
        yTrue_subset = tf.boolean_mask(yTrue, mask)
        yPred_subset = tf.boolean_mask(yPred, mask)
        # MSE
        return metrics.mse(yTrue_subset, yPred_subset)


    # Return the function
    return loss
  # MSE
        return metrics.mse(yTrue_subset, yPred_subset)
    
    
    # Return the function
    return loss


