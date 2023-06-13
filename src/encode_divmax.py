from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
import numpy as np
from kmer_gen import gen_all_kmers

fasta_file_path = sys.argv[1]
contig_type = sys.argv[2]
out_dir = os.path.join(os.path.dirname(fasta_file_path), "encode_divmax")
os.makedirs(out_dir, exist_ok=True)

def normalize_kmers(sequence, k):
    # Count the occurrence of each k-mer in the sequence
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmer_counts = Counter(kmers)
    
    # Find the maximum count of k-mers in the sequence
    max_count = max(kmer_counts.values())
    
    # Normalize the k-mer counts
    normalized_kmers = np.zeros(4**k)
    for i, kmer in enumerate(gen_all_kmers(k)):
        normalized_kmers[i] = kmer_counts[kmer] / max_count
    
    return normalized_kmers

k = 6
normalized_sequences = []
normalized_sequencesbw = []
with open(fasta_file_path) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        sequencebw = Seq(sequence).reverse_complement()
        normalized_sequence = normalize_kmers(sequence, k)
        normalized_sequences.append(normalized_sequence)
        normalized_sequencebw = normalize_kmers(sequencebw, k)
        normalized_sequencesbw.append(normalized_sequencebw)
fasta_file.close()

normalized_sequences = np.array(normalized_sequences)
normalized_sequences = np.reshape(normalized_sequences, (*normalized_sequences.shape, 1))
normalized_sequencesbw = np.array(normalized_sequencesbw)
normalized_sequencesbw = np.reshape(normalized_sequencesbw, (*normalized_sequencesbw.shape, 1))
np.save(os.path.join(out_dir, contig_type + "_0k_codefw.npy"), normalized_sequences)
np.save(os.path.join(out_dir, contig_type + "_0k_codebw.npy"), normalized_sequencesbw)