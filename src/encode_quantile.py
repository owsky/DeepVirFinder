import numpy as np
from sklearn.preprocessing import quantile_transform
from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys

fasta_file_path = sys.argv[1]
contig_type = sys.argv[2]
out_dir = os.path.join(os.path.dirname(fasta_file_path), "encode_quantile")
os.makedirs(out_dir, exist_ok=True)

def calculate_kmer_counts(sequence, k):
    kmer_counts = {}
    n = len(sequence)
    for i in range(n - k + 1):
        kmer = sequence[i:i+k]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    return kmer_counts

def normalize_kmer_counts(sequence, k):
    kmer_counts = calculate_kmer_counts(sequence, k)
    # kmers = list(kmer_counts.keys())
    counts = np.array(list(kmer_counts.values())).reshape(-1, 1)
    
    # Apply quantile normalization using sklearn's quantile_transform
    normalized_counts = quantile_transform(counts, output_distribution='normal', random_state=0)

    # normalized_kmer_counts = dict(zip(kmers, normalized_counts.flatten()))
    return normalized_counts

normalized_sequences = []
normalized_sequencesbw = []
k = 3
with open(fasta_file_path) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        sequencebw = Seq(sequence).reverse_complement()

        normalized_sequence = normalize_kmer_counts(sequence, k)
        normalized_sequences.append(normalized_sequence)
        normalized_sequencebw = normalize_kmer_counts(sequencebw, k)
        normalized_sequencesbw.append(normalized_sequencebw)
fasta_file.close()
print(np.shape(normalized_sequences[0]))
print(normalized_sequences[0])
print(len(fasta_file))