import numpy as np
import sys
import time
import h5py
import tensorflow as tf
from keras.models import load_model
import tensorflow.keras.backend as kb
from spliceai import *
from utils import *
from multi_gpu import *
from constants import *
import multiprocessing

assert int(sys.argv[1]) in [80, 400, 2000, 10000]

print("Setting up...")
physical_devices = tf.config.list_physical_devices('GPU')
print("Number of physical GPUs:", len(physical_devices))

L = 32
N_CPUS = multiprocessing.cpu_count()
N_GPUS = len(physical_devices)
# Optionally set the number of threads
tf.config.threading.set_intra_op_parallelism_threads(N_CPUS) 
tf.config.threading.set_inter_op_parallelism_threads(N_CPUS)

# Print current settings (optional)
print(f"Intra-op parallelism threads: {tf.config.threading.get_intra_op_parallelism_threads()}")
print(f"Inter-op parallelism threads: {tf.config.threading.get_inter_op_parallelism_threads()}")

# Define window size (W), atrous rate (AR), and batch size based on sys.argv[1]
if int(sys.argv[1]) == 80:
    W = np.asarray([11, 11, 11, 11])
    AR = np.asarray([1, 1, 1, 1])
    BATCH_SIZE = 18 * 28  # Adjust batch size for CPU parallelism
elif int(sys.argv[1]) == 400:
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4])
    BATCH_SIZE = 18 * 28  # Adjust batch size for CPU parallelism
elif int(sys.argv[1]) == 2000:
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                    21, 21, 21, 21])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
                     10, 10, 10, 10])
    BATCH_SIZE = 12 * 28  # Adjust batch size for CPU parallelism
elif int(sys.argv[1]) == 10000:
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                    21, 21, 21, 21, 41, 41, 41, 41])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
                     10, 10, 10, 10, 25, 25, 25, 25])
    BATCH_SIZE = 28  # Adjust batch size for CPU parallelism

CL = 2 * np.sum(AR * (W - 1))
assert CL <= CL_max and CL == int(sys.argv[1])
print("\033[1mContext nucleotides: %d\033[0m" % (CL))
print("\033[1mSequence length (output): %d\033[0m" % (SL))

model = SpliceAI(L, W, AR)
model_m = make_parallel(model, N_GPUS)
model_m.compile(loss=categorical_crossentropy_2d, optimizer='adam')

#model_m.load_weights('./Models/SpliceAI10000_c5.h5')

#kb.set_value(model_m.optimizer.lr, 0.0005)
# Open the HDF5 file for training data

######################################
# ensure correct data
#h5f = h5py.File(data_dir + 'dataset' + '_k562_annotated_' + 'train' + '_' + 'all' + '.h5', 'r')
#h5f = h5py.File(data_dir + 'dataset' + '_k562_rna_seq_' + 'train' + '_' + 'all' + '.h5', 'r')
#h5f = h5py.File(data_dir + 'dataset' + '_k562_norm_annotated_' + 'train' + '_' + 'all' + '.h5', 'r')
h5f = h5py.File(data_dir + 'dataset' + '_norm_' + 'train' + '_' + 'all' + '.h5', 'r')

######################################

# Split indices for training and validation
num_idx = len(h5f.keys()) // 3
idx_all = np.random.permutation(num_idx)
idx_train = idx_all[:int(0.9 * num_idx)]
idx_valid = idx_all[int(0.9 * num_idx):]
EPOCH_NUM = 10

print("Learning rate: %.5f" % (kb.get_value(model_m.optimizer.lr)))

start_epoch = 0
start_time = time.time()

print(f"Training started for {EPOCH_NUM} epochs, with {len(idx_train)} training samples per epoch")

# Enable debug logging for device placement (optional)
tf.debugging.set_log_device_placement(True)

# Training loop
for epoch_num in range(start_epoch, EPOCH_NUM):
    print(f"Starting epoch {epoch_num + 1} of {EPOCH_NUM}")
    
    # Shuffle training indices at the start of each epoch
    np.random.shuffle(idx_train)
    
    # Iterate through shuffled indices
    for idx_num, idx in enumerate(idx_train, start=1):
        print(f"Processing training sample {idx_num} of {len(idx_train)} in epoch {epoch_num + 1}")
        X = h5f['X' + str(idx)][:]
        Y = h5f['Y' + str(idx)][:, :, :3]  # Extract only the first three values, ignoring SSU

                
        # Train the model batch by batch
        #print(f"Shape of Xc: {Xc.shape}, Shape of Yc: {Yc.shape}")
        model_m.fit(X, Y, batch_size=BATCH_SIZE, verbose=0)
        print("model.fit complete.....epoch:" + str(epoch_num + 1) + ", idx: " + str(idx))

    # Save the model after each epoch
    model.save('./Models/SpliceAI' + sys.argv[1] + '_c' + sys.argv[2] + '.h5')
    print("model saved")

    if (epoch_num + 1) >= 6:
        # Apply decay every epoch after the 6th epoch
        kb.set_value(model_m.optimizer.lr, 0.5 * kb.get_value(model_m.optimizer.lr))
        new_lr = kb.get_value(model_m.optimizer.lr)
        print(f"Learning rate changed to: {new_lr:.5f}")

    # Validation and learning rate decay logic remains the same
    if epoch_num == EPOCH_NUM - 1:
        # Validation set metrics
        print("--------------------------------------------------------------")
        print("\n\033[1mValidation set metrics:\033[0m")

        Y_true_A = []
        Y_true_D = []
        Y_pred_A = []
        Y_pred_D = []

        for idx in idx_valid:
            X = h5f['X' + str(idx)][:]
            Y = h5f['Y' + str(idx)][:, :, :3]  # Extract only the first three values, ignoring SSU

            Yp = model_m.predict(X, batch_size=BATCH_SIZE)

            if not isinstance(Yp, list):
                Yp = [Yp]

            Y_true_A.extend(Y[:, :, 1].flatten())
            Y_true_D.extend(Y[:, :, 2].flatten())
            Y_pred_A.extend(Yp[0][:, :, 1].flatten())
            Y_pred_D.extend(Yp[0][:, :, 2].flatten())

        print("\n\033[1mAcceptor:\033[0m")
        print_topl_statistics(np.asarray(Y_true_A), np.asarray(Y_pred_A))

        print("\n\033[1mDonor:\033[0m")
        print_topl_statistics(np.asarray(Y_true_D), np.asarray(Y_pred_D))

        print("\n\033[1mTraining set metrics:\033[0m")

        Y_true_1 = []
        Y_true_2 = []
        Y_pred_1 = []
        Y_pred_2 = []

        for idx in idx_train[:len(idx_valid)]:
            X = h5f['X' + str(idx)][:]
            Y = h5f['Y' + str(idx)][:, :, :3]  # Extract only the first three values, ignoring SSU

            #Xc, Yc = clip_datapoints(X, Y, CL, N_GPUS)
            Yp = model_m.predict(X, batch_size=BATCH_SIZE)

            if not isinstance(Yp, list):
                Yp = [Yp]

            Y_true_1.extend(Y[:, :, 1].flatten())
            Y_true_2.extend(Y[:, :, 2].flatten())
            Y_pred_1.extend(Yp[0][:, :, 1].flatten())
            Y_pred_2.extend(Yp[0][:, :, 2].flatten())

        print("\n\033[1mAcceptor:\033[0m")
        print_topl_statistics(np.asarray(Y_true_1), np.asarray(Y_pred_1))

        print("\n\033[1mDonor:\033[0m")
        print_topl_statistics(np.asarray(Y_true_2), np.asarray(Y_pred_2))

        print("Learning rate: %.5f" % (kb.get_value(model_m.optimizer.lr)))
        print("--- %s seconds ---" % (time.time() - start_time))
        start_time = time.time()

        print("--------------------------------------------------------------")

        # Save model
        model.save('./Models/SpliceAI' + sys.argv[1] + '_c' + sys.argv[2] + '.h5')

   

h5f.close()
