from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
import numpy as np

fasta_file_path = sys.argv[1]
contig_type = sys.argv[2]
out_dir = os.path.join(os.path.dirname(fasta_file_path), "encode_mad")
os.makedirs(out_dir, exist_ok=True)

def normalize_sequence(sequence):
    # Define a mapping for the genomic letters
    letter_to_int = {'A': 0, 'T': 1, 'C': 2, 'G': 3}

    # Convert the sequence string to a numpy array of integers
    sequence_array = np.array([letter_to_int[letter] for letter in sequence if letter in ["A", "T", "C", "G"]])

    # Calculate the median of the sequence
    median = np.median(sequence_array)

    # Calculate the median absolute deviation (MAD) of the sequence
    mad = np.median(np.abs(sequence_array - median))

    # Normalize the sequence by subtracting the median and dividing by the MAD
    normalized_sequence = (sequence_array - median) / mad

    return normalized_sequence

normalized_sequences = []
normalized_sequencesbw = []

with open(fasta_file_path) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        sequencebw = Seq(sequence).reverse_complement()

        normalized_sequence = normalize_sequence(sequence)
        normalized_sequences.append(normalized_sequence)
        normalized_sequencebw = normalize_sequence(sequencebw)
        normalized_sequencesbw.append(normalized_sequencebw)
fasta_file.close()
print(np.shape(normalized_sequences[0]))
print(normalized_sequences[0])