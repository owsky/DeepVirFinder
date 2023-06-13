from collections import Counter
import math
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
import numpy as np
from kmer_gen import gen_all_kmers

fasta_file_path = sys.argv[1]
contig_type = sys.argv[2]
out_dir = os.path.join(os.path.dirname(fasta_file_path), "encode_log")
os.makedirs(out_dir, exist_ok=True)

def log_normalization(sequence, k):
    # Count the occurrences of each k-mer
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmer_counts = Counter(kmers)

    all_kmers = gen_all_kmers(k)
    
    # Get the frequencies of all possible k-mers (including those with zero counts)
    frequencies = [kmer_counts.get(kmer, 0) for kmer in all_kmers]

    # Apply log transformation to the frequencies
    log_frequencies = [math.log(freq) if freq != 0 else 0 for freq in frequencies]

    return log_frequencies

k = 4
normalized_sequences = []
normalized_sequencesbw = []

with open(fasta_file_path) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        sequencebw = Seq(sequence).reverse_complement()

        normalized_sequence = log_normalization(sequence, k)
        normalized_sequences.append(normalized_sequence)
        normalized_sequencebw = log_normalization(sequencebw, k)
        normalized_sequencesbw.append(normalized_sequencebw)
fasta_file.close()

normalized_sequences = np.array(normalized_sequences)
normalized_sequences = np.reshape(normalized_sequences, (*normalized_sequences.shape, 1))
normalized_sequencesbw = np.array(normalized_sequencesbw)
normalized_sequencesbw = np.reshape(normalized_sequencesbw, (*normalized_sequencesbw.shape, 1))
np.save(os.path.join(out_dir, contig_type + "_0k_codefw.npy"), normalized_sequences)
np.save(os.path.join(out_dir, contig_type + "_0k_codebw.npy"), normalized_sequencesbw)