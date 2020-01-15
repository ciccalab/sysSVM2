#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 11:52:37 2019

@author: nulsenj
"""

#################### Parse arguments ####################

# Set up argument parser
from argparse import ArgumentParser
parser = ArgumentParser(description = "Get predictions on data from a saved AAE model")


# Required arguments (positional, no defaults)
parser.add_argument("input_file", help = "Data to predict on (systems-level/molecular features)")
parser.add_argument("scaling_factors_file", help = "Scaling factors (output of initial training)")
parser.add_argument("trained_model_file", help = "Trained AAE model")
parser.add_argument("output_dir", help = "Directory to write output to")


# Optional arguments
parser.add_argument("--n_threads", type = int, default = 0, help = "Thread count for parallelisation")
parser.add_argument("--output_name", default = "", help = "Optional name to add to all output files")
parser.add_argument("--save_RDR_reconstruction", type = bool, default = False, help = "Save the RDR and feature reconstruction in addition to driver scores?")


# Extract arguments
args = parser.parse_args()


# # For testing
# from argparse import Namespace
# wdir = "/camp/home/nulsenj/working/Joel/OAC_sysSVM_2.0/methods_paper/augmented_AE/"
# args = Namespace(input_file = wdir + "input_data/sim_orig_sysMol_scores.tsv",
#                  scaling_factors_file = wdir + "test_newScripts/first20_scaling_factors.tsv",
#                  trained_model_file = wdir + "test_newScripts/first20_AAE.h5",
#                  output_dir = wdir + "test_newScripts",
#                  output_name = "preds",
#                  save_RDR_reconstruction = True,
#                  n_threads = 1)





#################### Packages ####################
print("Loading packages")

# General
import numpy as np
import pandas as pd


# Tensorflow / Keras
import tensorflow as tf
import tensorflow.keras.backend as K
from tensorflow.keras.models import load_model


# Custom functions, defined in AAE_functions.py (must be in the same directory as this script)
from AAE_functions import minMax, scaleData, scoreLoss


# Configure tensorflow session to have appropriate number of threads
session_conf = tf.compat.v1.ConfigProto(intra_op_parallelism_threads = args.n_threads, inter_op_parallelism_threads = args.n_threads)
sess = tf.compat.v1.Session(config = session_conf)
tf.compat.v1.keras.backend.set_session(sess)





#################### Prepare data (two datasets - 'original' with scores, and 'new' without scores) ####################

# Load data
print("Reading input from", args.input_file)
features_df = pd.read_csv(args.input_file, sep = "\t", index_col = "id")


# Scale data
print("Scaling data using scaling factors provided at", args.scaling_factors_file)
scaling_factors_df = pd.read_csv(args.scaling_factors_file, sep = "\t", index_col = "feature")
features_df = features_df[scaling_factors_df.drop(index = ["score"]).index.values] # Ensures columns are in the correct order
features_df_scaled = scaleData(df = features_df.copy(),
                               scaling_factors = scaling_factors_df.drop(index = ["score"]),
                               verbose = True)


# Convert to numpy
prediction_IDs = np.array(features_df_scaled.index.values)
prediction_np = features_df_scaled.to_numpy()





#################### Extract score estimates for test data ####################
print("Getting score estimates, RDR and feature reconstruction")

# Load pre-trained from HDF5 file
# scoreLoss is defined in AAE_functions.py
aae = load_model(args.trained_model_file, custom_objects = {"scoreLoss": scoreLoss})


# Run data through the model (and submodels)
rdr = aae.get_layer("encoder").predict(prediction_np)
score_estimates = aae.get_layer("score_nn").predict(rdr)


# Convert scores to dataframe and unscale
scores_df = pd.DataFrame(data = score_estimates, index = prediction_IDs, columns = ["score"])
scores_df["score"] = minMax(x = scores_df["score"],
			    old = [ scaling_factors_df["internal_min"]["score"], scaling_factors_df["internal_max"]["score"] ],
			    new = [ scaling_factors_df["external_min"]["score"], scaling_factors_df["external_max"]["score"] ])


# If you want to save the RDR and feature reconstruction, do similarly for these
if args.save_RDR_reconstruction:
    rdr_df = pd.DataFrame(data = rdr, index = prediction_IDs, columns = ["rdr_" + str(i) for i in range(1, rdr.shape[1]+1)])
    feature_reconstruction = aae.get_layer("decoder").predict(rdr)
    reconstruction_df = pd.DataFrame(data = feature_reconstruction, index = prediction_IDs, columns = list(features_df.columns.values))
    for col in reconstruction_df.columns.values:
        reconstruction_df[col] = minMax(x = reconstruction_df[col],
	   				old = [ scaling_factors_df["internal_min"][col], scaling_factors_df["internal_max"][col] ],
					new = [ scaling_factors_df["external_min"][col], scaling_factors_df["external_max"][col] ])





#################### Save results and trained model ####################
print("Saving results")


# Fix output directory so it has a forward slash at the end
output_dir = args.output_dir
if output_dir[len(output_dir)-1] != "/":
    output_dir = output_dir + "/"


# Fix output name so it has an underscore at the end
output_name = args.output_name
if output_name != "":
    if output_name[len(output_name)-1] != "_":
        output_name = output_name + "_"


# Scores
scores_df.to_csv(output_dir + output_name + "scores.tsv", sep = "\t", index = True, index_label = "id")


# RDR and feature reconstruction
if args.save_RDR_reconstruction:
    rdr_df.to_csv(output_dir + output_name + "rdr.tsv", sep = "\t", index = True, index_label = "id")
    reconstruction_df.to_csv(output_dir + output_name + "reconstruction_scores.tsv", sep = "\t", index = True, index_label = "id")


print("Finished")

#
