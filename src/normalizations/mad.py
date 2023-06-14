import numpy as np
import pandas as pd

def letter_to_int(c):
    dict = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    if c in ["A", "T", "C", "G"]:
        return dict[c]
    else:
        return -1
    

def normalize_sequence(sequence):
    # Create a DataFrame with the sequence as a single column
    sequence_array = np.array([letter_to_int(letter) for letter in sequence], dtype=float)
    df = pd.DataFrame({'Sequence': sequence_array})

    # Compute the mean absolute deviation
    mad = df['Sequence'].mad()
    epsilon = 1e-8

    # Normalize the sequence by subtracting the mean and dividing by the mean absolute deviation
    normalized_sequence = (df['Sequence'] - df['Sequence'].mean()) / (mad + epsilon)

    return normalized_sequence.tolist()
