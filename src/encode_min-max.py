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

def min_max_normalization(sequence, k):
    # Count the occurrences of each k-mer
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmer_counts = Counter(kmers)
    
    # Get the frequencies of k-mers
    frequencies = list(kmer_counts.values())
    
    min_freq = min(frequencies)
    max_freq = max(frequencies)
    
    normalized_frequencies = [(freq - min_freq) / (max_freq - min_freq) for freq in frequencies]
    
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
print(np.shape(normalized_sequences[0]))
print(normalized_sequences[0])