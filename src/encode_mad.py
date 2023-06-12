import numpy as np

def normalize_sequence(sequence):
    # Define a mapping for the genomic letters
    letter_to_int = {'A': 0, 'T': 1, 'C': 2, 'G': 3}

    # Convert the sequence string to a numpy array of integers
    sequence_array = np.array([letter_to_int[letter] for letter in sequence])

    # Calculate the median of the sequence
    median = np.median(sequence_array)

    # Calculate the median absolute deviation (MAD) of the sequence
    mad = np.median(np.abs(sequence_array - median))

    # Normalize the sequence by subtracting the median and dividing by the MAD
    normalized_sequence = (sequence_array - median) / mad

    return normalized_sequence

sequence = "ATCGATCG"
normalized_sequence = normalize_sequence(sequence)
print(normalized_sequence)
