import numpy as np
from Bio import SeqIO
import os
import sys
np.random.seed(42)

input_file = sys.argv[1]
seq_length = int(sys.argv[2])
contig_type = sys.argv[3]
direction = sys.argv[4]
out_dir = os.path.join(os.path.dirname(input_file), "encode", contig_type + str(seq_length) + direction + ".npy")
os.makedirs(os.path.join(os.path.dirname(input_file), "encode"), exist_ok=True)

def encode_sequence(sequence):
    encoding = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1]}
    encoded_sequence = [encoding.get(base, [1/4, 1/4, 1/4, 1/4]) for base in sequence]
    return np.array(encoded_sequence)

def fragment_sequence(sequence, fragment_length):
    sequence_length = len(sequence)
    start = np.random.randint(0, sequence_length - fragment_length)
    end = start + fragment_length
    fragment = sequence[start:end]
    return fragment

def encode_and_fragment_sequences(fasta_file, fragment_length, direction):
    sequences = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequence = str(record.seq) if direction == "fw" else str(record.seq.reverse_complement())
        fragment = fragment_sequence(sequence, fragment_length)
        encoded_fragment = encode_sequence(fragment)
        sequences.append(encoded_fragment)

    encoded_sequences = np.array(sequences)
    np.save(out_dir, encoded_sequences)
    print(f"Encodings saved in {out_dir}")

encode_and_fragment_sequences(input_file, seq_length, direction)
