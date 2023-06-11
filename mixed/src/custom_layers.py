from keras import backend as K
from keras.layers import Layer
from typing import Tuple


class RowNormalization(Layer):
    def __init__(self, **kwargs):
        super(RowNormalization, self).__init__(**kwargs)
    
    def build(self, input_shape: Tuple[int, int]) -> None:
        super(RowNormalization, self).build(input_shape)

    def call(self, x):
        row_sums = K.sum(x, axis=2, keepdims=True)
        normalized_data = x / (row_sums + K.epsilon())
        return normalized_data
    
    def compute_output_shape(self, input_shape: Tuple[int, int]) -> Tuple[int, int]:
        return input_shape
    
class ColumnNormalization(Layer):
    def __init__(self, **kwargs):
        super(ColumnNormalization, self).__init__(**kwargs)
    
    def build(self, input_shape: Tuple[int, int]) -> None:
        super(ColumnNormalization, self).build(input_shape)

    def call(self, inputs):
        # Calculate the mean and standard deviation for each column
        # mean, variance = K.nn.moments(inputs, axes=[0])
        mean = K.mean(inputs, axis=0)
        variance = K.var(inputs, axis=0)
        # Normalize the columns
        normalized = (inputs - mean) / K.sqrt(variance + 1e-8)
        return normalized
    
    def compute_output_shape(self, input_shape: Tuple[int, int]) -> Tuple[int, int]:
        return input_shape
    
__all__ = ['RowNormalization', 'ColumnNormalization']