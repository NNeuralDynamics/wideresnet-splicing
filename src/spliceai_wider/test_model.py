###############################################################################
# This file contains code to test the SpliceAI model.
###############################################################################

import numpy as np
import pandas as pd
import sys
import time
import h5py
from keras.models import load_model
from sklearn.metrics import confusion_matrix, classification_report, average_precision_score
from utils import *
from constants import *
from spliceai import *

###############################################################################
# Ensure valid input for context length (CL)
###############################################################################

assert int(sys.argv[1]) in [80, 400, 2000, 10000]
CL = int(sys.argv[1])

###############################################################################
# Load model and test data
###############################################################################

BATCH_SIZE = 28

model = load_model('../OG_canonical_models/SpliceAI10000_c1.h5')

# Load test data
h5f = h5py.File('../ml_data/dataset_k562_rna_test_0.h5', 'r')

# if h5 has GC datsets make 3 , if not 2
num_idx = len(h5f.keys()) // 3  

###############################################################################
# Model testing
###############################################################################

start_time = time.time()
print("\033[1mTest set metrics:\033[0m")

Y_true_class_A, Y_true_class_D = [], []
Y_pred_class_A, Y_pred_class_D = [], []
Y_true_sse_A, Y_true_sse_D = [], []
Y_pred_prob_sse_A, Y_pred_prob_sse_D = [], []
Y_true_all, Y_pred_prob_all = [], []

# Iterate over all dataset indices
for idx in range(num_idx):
    Xc = h5f[f'X{idx}'][:]  # Load input
    Yc = h5f[f'Y{idx}'][:, :, :4]  # Load all four channels of Y (assuming SSU is in the 4th channel)

    # Predict classification probabilities
    Yp_class = model.predict(Xc, batch_size=BATCH_SIZE)
    if not isinstance(Yp_class, list):
        Yp_class = [Yp_class]

    # Collect true and predicted classification values
    Y_true_class_A.extend(Yc[:, :, 1].flatten())  # True acceptor labels
    Y_true_class_D.extend(Yc[:, :, 2].flatten())  # True donor labels
    Y_pred_class_A.extend(Yp_class[0][:, :, 1].flatten())  # Predicted acceptor probabilities
    Y_pred_class_D.extend(Yp_class[0][:, :, 2].flatten())  # Predicted donor probabilities

    # Collect SSU values for true acceptor and donor sites
    Y_true_sse_A.extend(Yc[:, :, 3][Yc[:, :, 1] == 1].flatten())  # True SSU for acceptor sites
    Y_true_sse_D.extend(Yc[:, :, 3][Yc[:, :, 2] == 1].flatten())  # True SSU for donor sites

    # Collect predicted probabilities (classification) for SSU comparison
    Y_pred_prob_sse_A.extend(Yp_class[0][:, :, 1][Yc[:, :, 1] == 1].flatten())
    Y_pred_prob_sse_D.extend(Yp_class[0][:, :, 2][Yc[:, :, 2] == 1].flatten())

    # Print shapes and last five values for debugging
    # print("Y_pred_class_A shape:", np.array(Y_pred_class_A).shape)
    # print("Last 5 values of Y_pred_class_A:", Y_pred_class_A[-5:])
    # print("Y_pred_class_D shape:", np.array(Y_pred_class_D).shape)
    # print("Last 5 values of Y_pred_class_D:", Y_pred_class_D[-5:])

    # Y_true_sse_A.extend(Yc[:, :, 3][Yc[:, :, 1] == 1].flatten())
    # Y_pred_prob_sse_A.extend(Yp_class[:, :, 1][Yc[:, :, 1] == 1].flatten())

    # Y_true_sse_D.extend(Yc[:, :, 3][Yc[:, :, 2] == 1].flatten())
    # Y_pred_prob_sse_D.extend(Yp_class[:, :, 2][Yc[:, :, 2] == 1].flatten())

    # Y_true_all.extend(Yc[:, :, 3].flatten())
    # Y_pred_prob_all.extend((Yp_class[:, :, 1] + Yp_class[:, :, 2]).flatten())

h5f.close()

# print(f"Y_true_all shape: {np.array(Y_true_all).shape}")
# print(f"Y_pred_prob_all shape: {np.array(Y_pred_prob_all).shape}")

###############################################################################
# Print Metrics and Compare SSE with Classification Probabilities
###############################################################################

# Y_true_stack = np.stack([1 - (np.array(Y_true_class_A) + np.array(Y_true_class_D)),
#                          np.array(Y_true_class_A), 
#                          np.array(Y_true_class_D)], axis=1)

# Y_pred_stack = np.stack([1 - (np.array(Y_pred_class_A) + np.array(Y_pred_class_D)),
#                          np.array(Y_pred_class_A), 
#                          np.array(Y_pred_class_D)], axis=1)

# Y_true_class = np.argmax(Y_true_stack, axis=1)
# Y_pred_class = np.argmax(Y_pred_stack, axis=1)

# print("\nConfusion Matrix (rows = Actual, columns = Predicted):")
# conf_matrix = confusion_matrix(Y_true_class, Y_pred_class, labels=[0, 1, 2])
# print(conf_matrix)

# print("\n\033[1mClassification Report:\033[0m")
# print(classification_report(Y_true_class, Y_pred_class, target_names=['Null', 'Acceptor', 'Donor']))

###############################################################################
# Call `print_topl_statistics` for All Predictions
###############################################################################

print("\n\033[1mTop-L Statistics for Acceptor Sites:\033[0m")
print_topl_statistics(
    np.array(Y_true_class_A),  # Binary true acceptor labels
    np.array(Y_pred_class_A),  # Predicted acceptor probabilities
    np.array(Y_true_sse_A)     # Continuous SSU values for acceptor sites
)

print("\n\033[1mTop-L Statistics for Donor Sites:\033[0m")
print_topl_statistics(
    np.array(Y_true_class_D),  # Binary true donor labels
    np.array(Y_pred_class_D),  # Predicted donor probabilities
    np.array(Y_true_sse_D)     # Continuous SSU values for donor sites
)

###############################################################################
# Call `print_regression_statistics` for All Predictions
###############################################################################

# Filter predicted classification probabilities to match true acceptor and donor sites
Y_pred_class_A_true_sites = np.array(Y_pred_class_A)[np.array(Y_true_class_A) == 1]
Y_pred_class_D_true_sites = np.array(Y_pred_class_D)[np.array(Y_true_class_D) == 1]

# Print Regression Metrics for Acceptor and Donor using Classification Probabilities
print("\n\033[1mAcceptor Regression (SSU) using Classification Probabilities - True Sites:\033[0m")
print_regression_statistics(np.array(Y_true_sse_A), Y_pred_class_A_true_sites)

print("\n\033[1mDonor Regression (SSU) using Classification Probabilities - True Sites:\033[0m")
print_regression_statistics(np.array(Y_true_sse_D), Y_pred_class_D_true_sites)

###############################################################################
# Final Output
###############################################################################

print(f"--- {time.time() - start_time:.2f} seconds ---")
print("--------------------------------------------------------------")
