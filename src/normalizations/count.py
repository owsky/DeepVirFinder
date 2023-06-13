from collections import Counter
import itertools


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


def normalize_sequence(sequence: str) -> list:
    k = 4
    kmer_index = create_kmer_index(k)
    return sequence_to_vector(sequence, kmer_index, k)