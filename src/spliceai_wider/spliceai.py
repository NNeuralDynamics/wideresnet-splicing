###############################################################################
# This file has the functions necessary to create the SpliceAI model.
###############################################################################
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Activation,Conv1D, Cropping1D, BatchNormalization, add
import tensorflow.keras.backend as kb
import numpy as np


def ResidualUnit(l, w, ar):
    # Residual unit proposed in "Identity mappings in Deep Residual Networks"
    # by He et al.

    def f(input_node):
        bn1 = BatchNormalization()(input_node)
        act1 = Activation('relu')(bn1)
        conv1 = Conv1D(filters=l, kernel_size=(w,), dilation_rate=(ar,), padding='same')(act1)
        bn2 = BatchNormalization()(conv1)
        act2 = Activation('relu')(bn2)
        conv2 = Conv1D(l, kernel_size=(w,), dilation_rate=(ar,), padding='same')(act2)
        output_node = add([conv2, input_node])

        return output_node

    return f


def SpliceAI(L, W, AR):
    # L: Number of convolution kernels
    # W: Convolution window size in each residual unit
    # AR: Atrous rate in each residual unit

    assert len(W) == len(AR)

    CL = 2 * np.sum(AR*(W-1))

    input0 = Input(shape=(None, 4))
    conv = Conv1D(L, 1)(input0)
    skip = Conv1D(L, 1)(conv)

    for i in range(len(W)):
        conv = ResidualUnit(L, W[i], AR[i])(conv)
        
        if (((i+1) % 4 == 0) or ((i+1) == len(W))):
            # Skip connections to the output after every 4 residual units
            dense = Conv1D(L, 1)(conv)
            skip = add([skip, dense])

    skip = Cropping1D(int(CL/2))(skip)

    output0 = [[] for t in range(1)]

    for t in range(1):
        # Add batch normalization before the final Conv1D layer
        bn = BatchNormalization()(skip)
        output0[t] = Conv1D(3, 1, activation='softmax')(bn)
    
    model = Model(inputs=input0, outputs=output0)

    return model

def categorical_crossentropy_2d(y_true, y_pred):
    y_true = tf.cast(y_true, tf.float32)
    y_pred = tf.cast(y_pred, tf.float32)
    # Standard categorical cross entropy for sequence outputs

    return - kb.mean(y_true[:, :, 0]*kb.log(y_pred[:, :, 0]+1e-10)
                   + y_true[:, :, 1]*kb.log(y_pred[:, :, 1]+1e-10)
                   + y_true[:, :, 2]*kb.log(y_pred[:, :, 2]+1e-10))
                   