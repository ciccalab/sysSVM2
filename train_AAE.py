#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 13:11:08 2019

@author: nulsenj

Based on code from https://keras.io/examples/mnist_denoising_autoencoder/
"""

#################### Parse arguments ####################

# Set up argument parser
from argparse import ArgumentParser
parser = ArgumentParser(description = "Train an AAE to reproduce sysSVM scores, either de novo or from a warm start")


# Required arguments relating to I/O (positional, no defaults)
parser.add_argument("input_files", help = "Locations of training data (comma-separated; systems-level/molecular features and scores)")
parser.add_argument("output_dir", help = "Directory to write output to")


# Arguments required if warm-starting from a saved model
parser.add_argument("--trained_model_h5", default = None, help = "Location of pre-trained AAE model, saved in h5 format. Leave unspecified to start de novo.")
parser.add_argument("--scaling_factors", default = None, help = "Location of scaling factors (output of initial training). Required if starting from a saved model.")


# Optional arguments if starting de novo
parser.add_argument("--internal_minMax", nargs = 2, type = float, default = [-0.8, 0.8], help = "Internal range to scale data if starting de novo")
parser.add_argument("--linear_preLayer", type = int, default = 20, help = "Include a first layer with linear activation. Zero not to use, otherwise an integer specifying layer dimension.")
parser.add_argument("--encoder_layer_dims", nargs = "+", type = int, default = [20, 14, 10, 7], help = "Dimensions of encoder layers, including latent dimension. Intermediate dimensions are reveresed for the decoder.")
parser.add_argument("--score_layer_dims", nargs = "*", type = int, default = [5, 3, 2], help = "Dimensions of layers between latent representation and score output.")
parser.add_argument("--score_lossWeight", type = float, default = 5, help = "Relative weight of score loss compared to reconstruction loss")
parser.add_argument("--init_SD", type = float, default = 0.5, help = "Standard deviation of initial weight values")
parser.add_argument("--learning_rate", type = float, default = 0.0005, help = "Learning rate for Adam optimiser")


# General use arguments
parser.add_argument("--n_threads", type = int, default = 0, help = "Thread count for parallelisation")
parser.add_argument("--verbose", type = bool, default = False, help = "Verbose output")
parser.add_argument("--epochs", type = int, default = 100, help = "Number of epochs to use in training")
parser.add_argument("--batch_size", type = int, default = 0, help = "Batch size for training; set to zero not to use batches")
parser.add_argument("--val_proportion", type = float, default = 0.2, help = "Proportion of observations to use for validation")
parser.add_argument("--seed", type = int, default = None, help = "Seed for numpy psuedo-random number generator")
parser.add_argument("--output_name", default = "", help = "Optional name to add to all output files")


# Extract arguments
args = parser.parse_args()
warm_start = args.trained_model_h5 is not None
if not warm_start:
    print("Training AAE de novo")
else:
    print("Warm-starting pre-trained AAE")


# # For testing
# from argparse import Namespace
# wdir = "/camp/home/nulsenj/working/Joel/OAC_sysSVM_2.0/methods_paper/augmented_AE/"
# args = Namespace(input_files = [wdir + "input_data/sim_1000_plus_10.tsv"],
#                  output_dir = wdir + "test_newScripts",
#                  output_name = "testing",
#                  trained_model_h5 = None, scaling_factors = None, # de novo
#                  # trained_model_h5 = wdir + "initial_training/weight5_10kEpochs/trained_aae_model.h5", scaling_factors = wdir + "initial_training/weight5_10kEpochs/scaling_factors.tsv", # warm start
#                  internal_minMax = [-0.8, 0.8], linear_preLayer = 0, encoder_layer_dims = [14, 7, 4], score_layer_dims = [3, 2],
#                  score_lossWeight = 5, init_SD = 0.5, learning_rate = 0.0005,
#                  verbose = True, epochs = 5, val_proportion = 0.2, seed = 1337, n_threads = 1)


# Check argument consistency
if warm_start and args.scaling_factors is None:
    print("Error: you must provide scaling factors if warm-starting from a saved model")
    exit()

if not warm_start and args.internal_minMax[0] >= args.internal_minMax[1]:
    print("Error: the internal scaling min must be smaller than internal scaling max")
    exit()





#################### Packages ####################
print("Loading packages")

# General
import numpy as np
import pandas as pd
import re


# Tensorflow / Keras
import tensorflow as tf
from tensorflow.keras.layers import Dense, Input
from tensorflow.keras.models import Model, load_model
from tensorflow.compat.v1.keras.initializers import RandomNormal
from tensorflow.keras.losses import mse, mae
from tensorflow.keras.optimizers import Adam
from tensorflow.compat.v1 import Session, reshape
import tensorflow.keras.backend as K


# Custom functions, defined in AAE_functions.py (must be in the same directory as this script)
from AAE_functions import minMax, scaleData, splitData, scoreLoss#, newSampleMetric


# Configure tensorflow session to have appropriate number of threads
session_conf = tf.compat.v1.ConfigProto(intra_op_parallelism_threads = args.n_threads, \
                                        inter_op_parallelism_threads = args.n_threads, \
                                        allow_soft_placement=True, \
                                        device_count = {'CPU': args.n_threads})
sess = tf.compat.v1.Session(config = session_conf)
tf.compat.v1.keras.backend.set_session(sess)





#################### Prepare data ####################

# For reproducibility
if args.seed is not None:
    np.random.seed(args.seed)


# Load input data (systems-level and molecular features, as well as sysSVM scores)
print("Reading input data from", args.input_files)
input_files = args.input_files.split(",")
features_scores_df_list = []
for input_file in input_files:
    features_scores_df_list.append(pd.read_csv(input_file, sep = "\t", index_col = "id"))
# features_scores_df = pd.read_csv(args.input_file, sep = "\t", index_col = "id")
features_scores_df = pd.concat(features_scores_df_list)


# Scale data
if not warm_start:
    features_scores_df_scaled, scaling_factors_df = scaleData(df = features_scores_df.copy(),
                                                              new_minMax = args.internal_minMax,
                                                              verbose = True)
else:
    scaling_factors_df = pd.read_csv(args.scaling_factors, sep = "\t", index_col = "feature")
    features_scores_df = features_scores_df[scaling_factors_df.index.values] # Ensures columns are in the correct order
    features_scores_df_scaled = scaleData(df = features_scores_df.copy(),
                                          scaling_factors = scaling_factors_df.copy(),
                                          verbose = True)


# Replace missing score values with -10000
features_scores_df_scaled.fillna(value = -10000, inplace = True)


# Split data into training-validation, convert to numpy
print("Splitting training / validation data", (1-args.val_proportion)*100, "/", args.val_proportion*100)
training_IDs, training_np, training_scores, validation_IDs, validation_np, validation_scores = splitData(df = features_scores_df_scaled, val_proportion = args.val_proportion)





#################### Network settings ####################

if not warm_start:
    # Dimensions of AE layers
    feature_dim = training_np.shape[1]


    # Activation function to use for internal layers
    my_activationFunction = "tanh"


    # Initialisation distribution for weights
    my_init = RandomNormal(stddev = args.init_SD)


    # Loss function for scores, and loss function for sys/mol feature reconstruction
    # Relative weight of the score loss is set in args.score_lossWeight
    my_lossFunctions = [scoreLoss, mse]


    # Optimiser
    my_optimiser = Adam(lr = args.learning_rate)


# Batch size; CoLCC runs out of memory
if args.batch_size == 0:
    my_batch_size = training_np.shape[0]
else:
    my_batch_size = args.batch_size





#################### Build network (de novo only) ####################
if not warm_start:
    print("Building model architecture")


    # Encoder model
    # Input
    encoder_layers = [ Input(shape = (feature_dim,), name = "encoder_input") ]

    # Optional linear pre-layer
    if args.linear_preLayer != 0:
        encoder_layers.append( Dense(args.linear_preLayer, activation = "linear", kernel_initializer = my_init, name = "linear_preLayer")(encoder_layers[0]) )

    # Other layers
    for d in args.encoder_layer_dims:
        encoder_layers.append( Dense(d, activation = my_activationFunction, kernel_initializer = my_init)(encoder_layers[-1]) )

    # Instantiate
    encoder = Model(encoder_layers[0], encoder_layers[-1], name = "encoder")
    if args.verbose:
        encoder.summary()





    # Decoder model
    # Input (final layer of the encoder)
    decoder_layers = [ Input(shape = (args.encoder_layer_dims[-1],), name = "decoder_input") ]

    # Intermediate layers; reverse the dimensions of the encoder
    for d in args.encoder_layer_dims[-2::-1]:
        decoder_layers.append( Dense(d, activation = my_activationFunction, kernel_initializer = my_init)(decoder_layers[-1]) )

    # Output (feature reconstruction)
    decoder_layers.append( Dense(feature_dim, activation = "tanh", kernel_initializer = my_init, name = "feature_reconstruction")(decoder_layers[-1]) )

    # Instantiate
    decoder = Model(decoder_layers[0], decoder_layers[-1], name = "decoder")
    if args.verbose:
        decoder.summary()





    # Score prediction model
    # Input (final layer of the encoder)
    score_layers = [ Input(shape = (args.encoder_layer_dims[-1],), name = "score_input") ]

    # Intermediate layers
    for d in args.score_layer_dims:
        score_layers.append( Dense(d, activation = my_activationFunction, kernel_initializer = my_init)(score_layers[-1]) )

    # Output (predicted score)
    score_layers.append( Dense(1, activation = "tanh", kernel_initializer = my_init, name = "score_prediction")(score_layers[-1]) )

    # Instantiate
    score_nn = Model(score_layers[0], score_layers[-1], name = "score_nn")
    if args.verbose:
        score_nn.summary()





    # Augmented autoencoder = Encoder + Decoder + Score predictor
    # Instantiate the full augmented autoencoder
    aae = Model(encoder_layers[0], [score_nn(encoder(encoder_layers[0])), decoder(encoder(encoder_layers[0]))], name = "aae")
    if args.verbose:
        aae.summary()


    # Compile
    aae.compile(loss = my_lossFunctions,
                loss_weights=[args.score_lossWeight, 1],
                metrics = [mae, mse],
                optimizer = my_optimiser)





#################### Load trained model (warm start only) ####################
if warm_start:
    print("Loading pre-trained model from", args.trained_model_h5)

    # Load from HDF5 file
    # scoreLoss is defined in functions.py
    aae = load_model(args.trained_model_h5, custom_objects = {"scoreLoss": scoreLoss})
    if args.verbose:
        aae.summary()





#################### Compile and train ####################
print("Training augmented autoencoder")


# Train the autoencoder - this will also fit parameters of the submodels
training_history = aae.fit(training_np,
                           [training_scores, training_np],
                           batch_size = my_batch_size,
                           validation_data = (validation_np, [validation_scores, validation_np]),
                           verbose = int(args.verbose),
                           epochs = args.epochs)


# Record the losses and metrics over epochs
training_history = pd.DataFrame.from_dict(training_history.history)
training_history["epoch"] = list(range(1, args.epochs+1))
training_history.rename(columns = {"loss": "Total_loss", "val_loss": "val_Total_loss",
                                   "score_nn_loss": "Score_MSE+MAE", "val_score_nn_loss": "val_Score_MSE+MAE",
                                   "score_nn_mean_squared_error": "Score_MSE", "val_score_nn_mean_squared_error": "val_Score_MSE",
                                   "score_nn_mean_absolute_error": "Score_MAE", "val_score_nn_mean_absolute_error": "val_Score_MAE",
                                   "decoder_loss": "Reconstruction_MSE", "val_decoder_loss": "val_Reconstruction_MSE",
                                   "decoder_mean_absolute_error": "Reconstruction_MAE", "val_decoder_mean_absolute_error": "val_Reconstruction_MAE"},
                        inplace = True)
training_history.drop(columns = ["decoder_mean_squared_error", "val_decoder_mean_squared_error"], inplace = True)





#################### Save trained model and auxiliary files ####################
if not warm_start:
    print("Saving model, scaling factors and training history")
else:
    print("Saving model and training history")


# Fix output directory so it has a forward slash at the end
output_dir = args.output_dir
if output_dir[len(output_dir)-1] != "/":
    output_dir = output_dir + "/"


# Fix output name so it has an underscore at the end
output_name = args.output_name
if output_name != "":
    if output_name[len(output_name)-1] != "_":
        output_name = output_name + "_"


# Scaling factors
if not warm_start:
    scaling_factors_df.to_csv(output_dir + output_name + "scaling_factors.tsv", sep = "\t", index = True, index_label = "feature")


# Training history
training_history.to_csv(output_dir + output_name + "training_history.tsv", sep = "\t", index = False)


# Trained model
aae.save(output_dir + output_name + "AAE.h5")



print("Finished")

#
