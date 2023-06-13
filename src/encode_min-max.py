from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
import numpy as np
from kmer_gen import gen_all_kmers

fasta_file_path = sys.argv[1]
contig_type = sys.argv[2]
out_dir = os.path.join(os.path.dirname(fasta_file_path), "encode_min_max")
os.makedirs(out_dir, exist_ok=True)

def min_max_normalization(sequence, k):
    # Generate all possible k-mers
    all_kmers = gen_all_kmers(k)
    
    # Count the occurrences of each k-mer
    kmer_counts = Counter(all_kmers)
    
    # Get the frequencies of k-mers in the sequence
    sequence_kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmer_frequencies = [kmer_counts[kmer] for kmer in sequence_kmers]
    
    # Normalize the frequencies
    min_freq = min(kmer_frequencies)
    max_freq = max(kmer_frequencies)
    
    normalized_frequencies = [(freq - min_freq) / (max_freq - min_freq) if (max_freq - min_freq) != 0 else 0 for freq in kmer_frequencies]
    
    return normalized_frequencies


normalized_sequences = []
normalized_sequencesbw = []
k = 3
with open(fasta_file_path) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        sequencebw = Seq(sequence).reverse_complement()

        normalized_sequence = min_max_normalization(sequence, k)
        normalized_sequences.append(normalized_sequence)
        normalized_sequencebw = min_max_normalization(sequencebw, k)
        normalized_sequencesbw.append(normalized_sequencebw)
fasta_file.close()

normalized_sequences = np.array(normalized_sequences)
normalized_sequences = np.reshape(normalized_sequences, (*normalized_sequences.shape, 1))
normalized_sequencesbw = np.array(normalized_sequencesbw)
normalized_sequencesbw = np.reshape(normalized_sequencesbw, (*normalized_sequencesbw.shape, 1))
np.save(os.path.join(out_dir, contig_type + "_0k_codefw.npy"), normalized_sequences)
np.save(os.path.join(out_dir, contig_type + "_0k_codebw.npy"), normalized_sequencesbw)