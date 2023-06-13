from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
import numpy as np
import itertools

fasta_file_path = sys.argv[1]
contig_type = sys.argv[2]
out_dir = os.path.join(os.path.dirname(fasta_file_path), "encode_count")
os.makedirs(out_dir, exist_ok=True)

def normalize_kmers(sequence, k):
    # Count the occurrences of each k-mer
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmers = [kmer for kmer in kmers if set(kmer).issubset({'A', 'C', 'G', 'T'})]
    kmer_counts = Counter(kmers)
    
    # Calculate the total number of k-mers
    total_kmers = len(sequence) - k + 1
    
    # Normalize the k-mer counts
    normalized_kmers = {kmer: count/total_kmers for kmer, count in kmer_counts.items()}
    
    return normalized_kmers


def create_kmer_index(k):
    # Create a dictionary of all possible k-mers of length k
    kmers = [''.join(p) for p in itertools.product('ACGT', repeat=k)]
    kmer_index = {kmer: i for i, kmer in enumerate(kmers)}
    return kmer_index

def sequence_to_vector(sequence, kmer_index, k):
    # Calculate the normalized k-mer counts for the sequence
    normalized_kmers = normalize_kmers(sequence, k)
    
    # Create a fixed-length vector representation of the normalized k-mer counts
    vector = [0] * len(kmer_index)
    for kmer, count in normalized_kmers.items():
        index = kmer_index[kmer]
        vector[index] = count
    
    return vector

encoded_sequences = []
encoded_sequencesbw = []
k = 4
with open(fasta_file_path) as fasta_file:
    kmer_index = create_kmer_index(k)
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        sequencebw = Seq(sequence).reverse_complement()

        normalized_sequence = sequence_to_vector(sequence, kmer_index, k)
        normalized_sequencebw = sequence_to_vector(sequencebw, kmer_index, k)
        encoded_sequences.append(normalized_sequence)
        encoded_sequencesbw.append(normalized_sequencebw)

    encoded_sequences = np.array(encoded_sequences)
    encoded_sequences = np.reshape(encoded_sequences, (*encoded_sequences.shape, 1))
    encoded_sequencesbw = np.array(encoded_sequencesbw)
    encoded_sequencesbw = np.reshape(encoded_sequencesbw, (*encoded_sequencesbw.shape, 1))
    np.save(os.path.join(out_dir, contig_type + "_0k_codefw.npy"), encoded_sequences)
    np.save(os.path.join(out_dir, contig_type + "_0k_codebw.npy"), encoded_sequencesbw)
fasta_file.close()